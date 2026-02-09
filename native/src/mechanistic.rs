use crate::mechanistic_params::{get_default_params, ParaOde, INITIAL_HORMONES};
use rand::Rng;
use rand_distr::{Distribution, Normal};
use std::f64::consts::PI;

const N_HORMONES: usize = 54;

#[derive(Clone, Debug)]
pub struct Follicle {
    pub size: f64,
    pub fsh_sensitivity: f64,
    pub destiny: i32, // -1: active/growing, -2: atretic, 1: ovulated, 3: big but anovulatory, 4: about to ovulate
    pub time_decrease: f64, // when decline started (for destiny 3 → -2 timeout & accelerated shrink)
    pub birth_time: f64,
}

pub struct MenstrualMechanisticModel {
    pub params: Vec<f64>,
    pub para_ode: ParaOde,
    pub hormones: [f64; N_HORMONES],
    pub follicles: Vec<Follicle>,
    pub t: f64,
    pub last_ovulation_t: f64,
    pub stress_input: f64,
    pub cumulative_ovulations: Vec<f64>,
    // Reusable buffers to avoid per-step allocations
    active_idx: Vec<usize>,
    fsh_sens_buf: Vec<f64>,
    y_buf: Vec<f64>,
    dy_buf: Vec<f64>,
    k1: Vec<f64>,
    k2: Vec<f64>,
    k3: Vec<f64>,
    k4: Vec<f64>,
    y_tmp: Vec<f64>,
}

impl MenstrualMechanisticModel {
    pub fn new(_seed: u64) -> Self {
        Self {
            params: get_default_params(),
            para_ode: ParaOde::default(),
            hormones: INITIAL_HORMONES,
            follicles: Vec::new(),
            t: 0.0,
            // MATLAB initializes Tovu = 14 — places the initial P4 Gaussian peak
            // at day 21, which gives a natural first-cycle profile (P4 ≈ 0 at t=0
            // since exp(-0.06 * 21^2) ≈ 0). Using -14 would put the Gaussian at -7,
            // giving a small but nonzero P4 at t=0. Match MATLAB for correctness.
            last_ovulation_t: 14.0,
            stress_input: 0.0,
            cumulative_ovulations: Vec::new(),
            active_idx: Vec::new(),
            fsh_sens_buf: Vec::new(),
            y_buf: Vec::new(),
            dy_buf: Vec::new(),
            k1: Vec::new(),
            k2: Vec::new(),
            k3: Vec::new(),
            k4: Vec::new(),
            y_tmp: Vec::new(),
        }
    }

    /// Advance the model by `dt` days. Returns true if an ovulation occurred.
    pub fn step(&mut self, dt: f64, rng: &mut impl Rng) -> bool {
        // --- 1. Recruit new follicles (Poisson process modulated by FSH) ---
        let y_fsh = self.hormones[7].max(0.0); // FSH in blood
        let fsh_imp = y_fsh.powi(3) / (y_fsh.powi(3) + 12.35f64.powi(3));
        let lambda_eff = ((20.0 / 14.0) + 4.0 * (20.0 / 14.0) * fsh_imp).max(0.0);

        if let Ok(poi) = rand_distr::Poisson::new(lambda_eff * dt) {
            let n_new = poi.sample(rng) as usize;
            if n_new > 0 {
                let normal = Normal::new(0.6, 0.55 * 0.6).unwrap();
                for _ in 0..n_new {
                    let sens: f64 = normal.sample(rng);
                    self.follicles.push(Follicle {
                        size: 4.0,
                        fsh_sensitivity: sens.max(0.1),
                        destiny: -1,
                        time_decrease: 0.0,
                        birth_time: self.t,
                    });
                }
            }
        }

        // --- 2. Collect active follicle indices and build ODE state ---
        self.active_idx.clear();
        self.fsh_sens_buf.clear();
        for (i, f) in self.follicles.iter().enumerate() {
            if f.destiny == -1 || f.destiny == 3 || f.destiny == 4 {
                self.active_idx.push(i);
                self.fsh_sens_buf.push(f.fsh_sensitivity);
            }
        }
        let nf = self.active_idx.len();
        let n = nf + N_HORMONES;

        // Build combined state: [follicle_sizes..., hormones...]
        self.y_buf.resize(n, 0.0);
        for (j, &idx) in self.active_idx.iter().enumerate() {
            self.y_buf[j] = self.follicles[idx].size;
        }
        self.y_buf[nf..].copy_from_slice(&self.hormones);

        // --- 3. RK4 integration ---
        // Use 4 substeps for stability (effective dt = dt/4 = 0.0625 days)
        let n_substeps = 4;
        let h = dt / n_substeps as f64;

        self.dy_buf.resize(n, 0.0);
        self.k1.resize(n, 0.0);
        self.k2.resize(n, 0.0);
        self.k3.resize(n, 0.0);
        self.k4.resize(n, 0.0);
        self.y_tmp.resize(n, 0.0);

        for _sub in 0..n_substeps {
            let t_sub = self.t + (_sub as f64) * h;

            // k1 = f(t, y)
            rhs(t_sub, &self.y_buf, &mut self.k1, nf, self.last_ovulation_t,
                self.stress_input, &self.para_ode, &self.params, &self.fsh_sens_buf, &self.active_idx, &self.follicles);

            // k2 = f(t + h/2, y + h/2 * k1)
            for i in 0..n { self.y_tmp[i] = self.y_buf[i] + 0.5 * h * self.k1[i]; }
            rhs(t_sub + 0.5 * h, &self.y_tmp, &mut self.k2, nf, self.last_ovulation_t,
                self.stress_input, &self.para_ode, &self.params, &self.fsh_sens_buf, &self.active_idx, &self.follicles);

            // k3 = f(t + h/2, y + h/2 * k2)
            for i in 0..n { self.y_tmp[i] = self.y_buf[i] + 0.5 * h * self.k2[i]; }
            rhs(t_sub + 0.5 * h, &self.y_tmp, &mut self.k3, nf, self.last_ovulation_t,
                self.stress_input, &self.para_ode, &self.params, &self.fsh_sens_buf, &self.active_idx, &self.follicles);

            // k4 = f(t + h, y + h * k3)
            for i in 0..n { self.y_tmp[i] = self.y_buf[i] + h * self.k3[i]; }
            rhs(t_sub + h, &self.y_tmp, &mut self.k4, nf, self.last_ovulation_t,
                self.stress_input, &self.para_ode, &self.params, &self.fsh_sens_buf, &self.active_idx, &self.follicles);

            // y = y + h/6 * (k1 + 2*k2 + 2*k3 + k4)
            for i in 0..n {
                self.y_buf[i] += h / 6.0 * (self.k1[i] + 2.0 * self.k2[i] + 2.0 * self.k3[i] + self.k4[i]);
                // Clamp negative values (concentrations can't go negative)
                if self.y_buf[i] < 0.0 { self.y_buf[i] = 0.0; }
            }
        }

        // --- 4. Write back state ---
        for (j, &idx) in self.active_idx.iter().enumerate() {
            self.follicles[idx].size = self.y_buf[j];
        }
        self.hormones.copy_from_slice(&self.y_buf[nf..]);
        self.t += dt;

        // --- 5. Follicle fate decisions (matching MATLAB testfun logic) ---
        // Evaluate growth rates for active follicles
        self.dy_buf.fill(0.0);
        rhs(self.t, &self.y_buf, &mut self.dy_buf, nf, self.last_ovulation_t,
            self.stress_input, &self.para_ode, &self.params, &self.fsh_sens_buf, &self.active_idx, &self.follicles);

        let y_lh = self.hormones[5]; // LH in blood
        let big_but_alive_timeout = 5.0; // MATLAB para(10) = 5 days
        let mut ovulated = false;

        for j in 0..self.active_idx.len() {
            let idx = self.active_idx[j];
            let f = &mut self.follicles[idx];
            let growth_rate = self.dy_buf[j];
            let age = self.t - f.birth_time;

            match f.destiny {
                -1 => {
                    // Active/growing follicle
                    if growth_rate <= 0.01 {
                        // Growth stalled → atretic
                        f.destiny = -2;
                        f.time_decrease = self.t;
                    } else if growth_rate <= 0.1 && age >= 2.0 {
                        // Slow growth after 2+ days → atretic
                        f.destiny = -2;
                        f.time_decrease = self.t;
                    } else if f.size >= self.para_ode.folmax && y_lh < 25.0 {
                        // Big enough but LH not high enough → waiting state
                        f.destiny = 3;
                        f.time_decrease = self.t;
                    } else if f.size >= self.para_ode.folmax && y_lh >= 25.0 {
                        // Ready to ovulate — enter pre-ovulation state
                        f.destiny = 4;
                        f.time_decrease = self.t;
                    }
                }
                3 => {
                    // Big but anovulatory — waiting for LH surge or timeout
                    if y_lh >= 25.0 {
                        f.destiny = 4;
                        f.time_decrease = self.t;
                    } else if self.t - f.time_decrease >= big_but_alive_timeout {
                        f.destiny = -2;
                        f.time_decrease = self.t;
                    }
                }
                4 => {
                    // About to ovulate — 0.5 day delay then ovulate
                    if self.t - f.time_decrease >= 0.5 {
                        f.destiny = 1;
                        self.last_ovulation_t = self.t;
                        self.cumulative_ovulations.push(self.t);
                        ovulated = true;
                    }
                }
                _ => {}
            }
        }

        // Apply accelerated decline for atretic follicles that are still in the ODE
        // (MATLAB: f(i) = -0.05 * y(i) * (t - TimeDecrease))
        // We do this as a post-hoc size update since atretic follicles aren't in the ODE
        for f in &mut self.follicles {
            if f.destiny == -2 && f.size > 0.1 {
                let dt_since = self.t - f.time_decrease;
                f.size *= (-0.05 * dt_since * dt).exp(); // exponential decay
                if f.size < 0.5 {
                    f.size = 0.0; // effectively dead
                }
            }
        }

        // Purge dead follicles to keep memory bounded
        if self.follicles.len() > 50 {
            self.follicles.retain(|f| f.destiny == -1 || f.destiny == 3 || f.destiny == 4);
        }

        ovulated
    }
}

/// Evaluate the RHS of the coupled follicle-hormone ODE system.
///
/// State layout: [follicle_0, ..., follicle_{nf-1}, hormone_0, ..., hormone_53]
///
/// The 54 hormone states map to MATLAB indices (1-based → 0-based):
///   hormones[0]  = baseline E2 compartment (constant in MATLAB)
///   hormones[1]  = E2 blood (algebraic in MATLAB, relaxation ODE here)
///   hormones[2]  = (unused/constant)
///   hormones[3]  = P4 blood (algebraic in MATLAB, relaxation ODE here)
///   hormones[4]  = RP_LH (LH in pituitary)
///   hormones[5]  = LH (LH in blood)
///   hormones[6]  = RP_FSH (FSH in pituitary)
///   hormones[7]  = FSH (FSH in blood)
///   hormones[8]  = FSHfoll (FSH in ovaries)
///   ...
///   hormones[35] = RLH
///   hormones[39] = RFSH (free FSH receptors)
///   hormones[40] = FSHR (bound FSH receptors)
///   hormones[41] = RFSH_des (desensitized FSH receptors)
///   ...
///   hormones[49] = GnRH
///   hormones[50] = RecGa (active GnRH receptor)
///   hormones[51] = RecGi (inactive GnRH receptor)
///   hormones[52] = GReca (active GnRH-receptor complex)
///   hormones[53] = GReci (inactive GnRH-receptor complex)
fn rhs(
    t: f64,
    y: &[f64],
    dy: &mut [f64],
    nf: usize,
    last_ovulation_t: f64,
    stress_input: f64,
    po: &ParaOde,
    p: &[f64],
    fsh_sens: &[f64],
    active_idx: &[usize],
    follicles: &[Follicle],
) {
    dy.fill(0.0);

    // Hormone indices (offset by nf)
    let ie2 = nf + 1;
    let ip4 = nf + 3;
    let irp_lh = nf + 4;
    let ilh = nf + 5;
    let irp_fsh = nf + 6;
    let ifsh = nf + 7;
    let ifshfoll = nf + 8;
    let irlh = nf + 35;
    let irfsh = nf + 39;
    let ifshr = nf + 40;
    let irfsh_des = nf + 41;
    let ignrh = nf + 49;
    let irec_ga = nf + 50;
    let irec_gi = nf + 51;
    let igrec_a = nf + 52;
    let igrec_i = nf + 53;

    let y_e2 = y[ie2].max(0.0);
    let y_p4 = y[ip4].max(0.0);
    let y_lh = y[ilh].max(0.0);
    let y_fsh = y[ifsh].max(0.0);

    // --- Follicular surface area (for E2 production) ---
    // SF = π * Σ(x^5 / (x^5 + 15^5) * x^2)  over active, non-atretic follicles
    // MATLAB zeroes out follicles with destiny == -2 or == 4 from the SF calculation
    let sf: f64 = if nf == 0 {
        0.0
    } else {
        let mut sum = 0.0;
        for j in 0..nf {
            let d = follicles[active_idx[j]].destiny;
            if d == -2 || d == 4 {
                continue; // MATLAB excludes atretic and pre-ovulation follicles
            }
            let x = y[j].max(0.0);
            let x5 = x.powi(5);
            sum += (x5 / (x5 + 15.0f64.powi(5))) * x.powi(2);
        }
        PI * sum
    };

    // --- E2 and P4 targets (algebraic constraints in MATLAB, relaxation here) ---
    let t_ov7 = last_ovulation_t + 7.0;
    let gauss = (-0.06 * (t - t_ov7).powi(2)).exp();
    let e2_target = 13.0 * (1.0 + 0.02 * sf) + 150.0 * gauss;
    let p4_target = 15.0 * gauss;

    // Fast relaxation to algebraic constraint.
    // The MATLAB treats E2/P4 as algebraic (mass matrix M[i,i]=0) so they
    // instantly track their targets. We approximate with a fast ODE.
    // Rate 50 gives τ = 0.02 days — fast enough at 0.0625-day RK4 substeps.
    dy[ie2] = (e2_target - y_e2) * 50.0;
    dy[ip4] = (p4_target - y_p4) * 50.0;

    // --- GnRH frequency and mass (MATLAB lines 39-46) ---
    let stress_factor = (1.0 - stress_input).max(0.0);
    let yg_freq = 16.0 / (1.0 + (y_p4 / p[203]).powf(p[204]))
        * (1.0 + hill_pos(y_e2, p[205], p[206]))
        * stress_factor;

    let yg_mass = hill_pos(y_e2, p[208], p[209]) + hill_neg(y_e2, p[210], p[211]);

    // --- LH pituitary production and release (MATLAB lines 52-72) ---
    let hp_e2 = hill_pos(y_e2, p[3], p[6]);
    let f_lh_prod = (p[1] + p[2] * hp_e2) / (1.0 + (y_p4 / p[4]).powf(p[7]));
    let f_lh_rel = p[16] + p[5] * hill_pos(y[igrec_a], p[8], p[9]);

    dy[irp_lh] = f_lh_prod - f_lh_rel * y[irp_lh];
    dy[ilh] = (1.0 / p[12]) * f_lh_rel * y[irp_lh]
        - p[230] * y_lh * y[irlh]
        - p[231] * y_lh;

    // --- FSH pituitary production and release (MATLAB lines 78-100) ---
    let hm_freq = 1.0 / (1.0 + (yg_freq / p[11]).powf(p[13]));
    let f_fsh_prod = (p[21] / (1.0 + (y_p4 / 5.0).powi(2))) * hm_freq;
    let f_fsh_rel = p[17] + p[28] * hill_pos(y[igrec_a], p[18], p[20]);

    dy[irp_fsh] = f_fsh_prod - f_fsh_rel * y[irp_fsh];
    dy[ifsh] = (1.0 / p[12]) * f_fsh_rel * y[irp_fsh]
        - (p[243] + p[241]) * y_fsh;

    // --- FSH in ovaries and receptor kinetics (MATLAB lines 107-126) ---
    dy[ifshfoll] = p[243] * y_fsh * p[246] - p[240] * y[ifshfoll] * y[irfsh];
    dy[irfsh] = p[242] * y[irfsh_des] - p[240] * y[ifshfoll] * y[irfsh];
    dy[ifshr] = p[240] * y[ifshfoll] * y[irfsh] - p[244] * y[ifshr];
    dy[irfsh_des] = p[244] * y[ifshr] - p[242] * y[irfsh_des];

    // --- GnRH and receptor dynamics (MATLAB lines 132-168) ---
    dy[ignrh] = p[301] * yg_mass * yg_freq
        - p[302] * y[ignrh] * y[irec_ga]
        + p[303] * y[igrec_a]
        - p[300] * y[ignrh];

    dy[irec_ga] = p[303] * y[igrec_a]
        - p[302] * y[ignrh] * y[irec_ga]
        - p[306] * y[irec_ga]
        + p[307] * y[irec_gi];

    dy[irec_gi] = p[311]
        + p[306] * y[irec_ga]
        - p[307] * y[irec_gi]
        + p[305] * y[igrec_i]
        - p[308] * y[irec_gi];

    dy[igrec_a] = p[302] * y[ignrh] * y[irec_ga]
        - p[303] * y[igrec_a]
        - p[309] * y[igrec_a]
        + p[310] * y[igrec_i];

    dy[igrec_i] = p[309] * y[igrec_a]
        - p[310] * y[igrec_i]
        - p[305] * y[igrec_i]
        - p[304] * y[igrec_i];

    // --- Follicle growth (MATLAB testfun_NormalCycle lines 36-83) ---
    let fsh_rez_comp = y[ifshr].max(0.0);
    let sum_v: f64 = (0..nf).map(|j| y[j].max(0.0).powf(po.v)).sum();

    for j in 0..nf {
        let d = follicles[active_idx[j]].destiny;
        if d == -2 {
            // Atretic: accelerated decline (handled outside ODE)
            dy[j] = 0.0;
            continue;
        }
        if d == 4 {
            // Pre-ovulation: follicle rests
            dy[j] = 0.0;
            continue;
        }

        let x = y[j].max(0.0);
        let sens = fsh_sens[j];

        // FSH responsiveness: Hill function with individual sensitivity
        let f_fsh = fsh_rez_comp / (fsh_rez_comp + sens);

        // Growth rate modulated by P4 and FSH receptor
        let gamma_p4 = po.gamma
            * ((1.0 / (1.0 + (y_p4 / 3.0).powi(5)))
                + hill_pos_raw(fsh_rez_comp, 0.59, 10.0));

        // Competition: negative Hill function for FSH
        let kappa_fsh = po.k * hill_neg_raw(fsh_rez_comp, 0.85, 25.0);

        dy[j] = f_fsh * (po.xi - x) * x * (gamma_p4 - kappa_fsh * (sum_v - po.mu * x.powf(po.v)));
    }
}

/// Hill positive: s^n / (k^n + s^n)
#[inline]
fn hill_pos(s: f64, k: f64, n: f64) -> f64 {
    if s <= 0.0 || k <= 0.0 {
        return 0.0;
    }
    let sn = s.powf(n);
    let kn = k.powf(n);
    sn / (kn + sn)
}

/// Hill negative: k^n / (k^n + s^n)
#[inline]
fn hill_neg(s: f64, k: f64, n: f64) -> f64 {
    if s <= 0.0 || k <= 0.0 {
        return 1.0;
    }
    let sn = s.powf(n);
    let kn = k.powf(n);
    kn / (kn + sn)
}

/// Hill positive with raw (non-ratio) form: s^n / (k^n + s^n)
/// Used for follicle growth where MATLAB uses explicit powf
#[inline]
fn hill_pos_raw(s: f64, k: f64, n: f64) -> f64 {
    let sn = s.powf(n);
    let kn = k.powf(n);
    sn / (kn + sn)
}

/// Hill negative with raw form: k^n / (k^n + s^n)
#[inline]
fn hill_neg_raw(s: f64, k: f64, n: f64) -> f64 {
    let sn = s.powf(n);
    let kn = k.powf(n);
    kn / (kn + sn)
}
