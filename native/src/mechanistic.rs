use crate::mechanistic_params::{get_default_params, ParaOde, INITIAL_HORMONES};
use crate::zavala::{ZavalaParams, ZavalaState};
use rand::Rng;
use rand_distr::{Distribution, Normal};
use std::f64::consts::PI;

const N_HORMONES: usize = 54;

#[derive(Clone, Debug)]
pub struct Follicle {
    pub size: f64,
    pub fsh_sensitivity: f64,
    pub destiny: i32, 
    pub time_decrease: f64, 
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

    // Zavala 2020 neuroendocrine network
    pub zavala: ZavalaState,

    // Reusable buffers
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
    pub fn new(_seed: u64, rng: &mut impl Rng) -> Self {
        Self {
            params: get_default_params(),
            para_ode: ParaOde::default(),
            hormones: INITIAL_HORMONES,
            follicles: Vec::new(),
            t: 0.0,
            last_ovulation_t: 14.0,
            stress_input: 0.0,
            cumulative_ovulations: Vec::new(),

            zavala: ZavalaState::new(ZavalaParams::default(), rng),

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

    /// Randomize the starting phase of the individual.
    pub fn randomize_initial_state(&mut self, rng: &mut impl Rng) {
        // 1. Randomize last ovulation time (up to 28 days ago)
        let offset: f64 = rng.random_range(0.0..28.0);
        self.last_ovulation_t = -offset;

        // 2. Zavala phases already randomized in ZavalaState::new()

        // 3. Recruit some initial follicles of random sizes
        let n_initial = rng.random_range(2..8);
        for _ in 0..n_initial {
            let size = rng.random_range(4.0..15.0);
            let val: f64 = Normal::new(0.6, 0.33).unwrap().sample(rng);
            let sens = val.max(0.1);
            self.follicles.push(Follicle {
                size,
                fsh_sensitivity: sens,
                destiny: -1,
                time_decrease: 0.0,
                birth_time: -rng.random_range(0.0..10.0),
            });
        }
    }

    pub fn step(&mut self, dt: f64, rng: &mut impl Rng) -> bool {
        // --- 1. Update Zavala neuroendocrine network ---
        // Compute estrous/menstrual phase from last ovulation
        let days_since_ov = self.t - self.last_ovulation_t;
        let cycle_len = if self.cumulative_ovulations.len() >= 2 {
            let n = self.cumulative_ovulations.len();
            let interval = self.cumulative_ovulations[n - 1] - self.cumulative_ovulations[n - 2];
            if interval >= 18.0 && interval <= 45.0 { interval } else { 28.0 }
        } else {
            28.0
        };
        let phi_e = (2.0 * PI * days_since_ov / cycle_len).clamp(0.0, 2.0 * PI - 1e-10);

        self.zavala.step(dt, self.stress_input, phi_e, self.t);
        let f_k = self.zavala.compute_f_k(self.stress_input);

        // --- 2. Recruit new follicles ---
        let y_fsh = self.hormones[7].max(0.0);
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

        // --- 3. Collect active follicle indices and build ODE state ---
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

        self.y_buf.resize(n, 0.0);
        for (j, &idx) in self.active_idx.iter().enumerate() {
            self.y_buf[j] = self.follicles[idx].size;
        }
        self.y_buf[nf..].copy_from_slice(&self.hormones);

        project_algebraic(
            &mut self.y_buf, self.t, nf, self.last_ovulation_t,
            &self.active_idx, &self.follicles,
        );

        // --- 4. RK4 integration ---
        let n_substeps = 64;
        let h = dt / n_substeps as f64;

        self.dy_buf.resize(n, 0.0);
        self.k1.resize(n, 0.0);
        self.k2.resize(n, 0.0);
        self.k3.resize(n, 0.0);
        self.k4.resize(n, 0.0);
        self.y_tmp.resize(n, 0.0);

        for _sub in 0..n_substeps {
            let t_sub = self.t + (_sub as f64) * h;

            rhs(t_sub, &self.y_buf, &mut self.k1, nf, self.last_ovulation_t,
                f_k, &self.para_ode, &self.params, &self.fsh_sens_buf, &self.active_idx, &self.follicles);

            for i in 0..n { self.y_tmp[i] = self.y_buf[i] + 0.5 * h * self.k1[i]; }
            rhs(t_sub + 0.5 * h, &self.y_tmp, &mut self.k2, nf, self.last_ovulation_t,
                f_k, &self.para_ode, &self.params, &self.fsh_sens_buf, &self.active_idx, &self.follicles);

            for i in 0..n { self.y_tmp[i] = self.y_buf[i] + 0.5 * h * self.k2[i]; }
            rhs(t_sub + 0.5 * h, &self.y_tmp, &mut self.k3, nf, self.last_ovulation_t,
                f_k, &self.para_ode, &self.params, &self.fsh_sens_buf, &self.active_idx, &self.follicles);

            for i in 0..n { self.y_tmp[i] = self.y_buf[i] + h * self.k3[i]; }
            rhs(t_sub + h, &self.y_tmp, &mut self.k4, nf, self.last_ovulation_t,
                f_k, &self.para_ode, &self.params, &self.fsh_sens_buf, &self.active_idx, &self.follicles);

            for i in 0..n {
                self.y_buf[i] += h / 6.0 * (self.k1[i] + 2.0 * self.k2[i] + 2.0 * self.k3[i] + self.k4[i]);
                if self.y_buf[i] < 0.0 { self.y_buf[i] = 0.0; }
            }

            project_algebraic(
                &mut self.y_buf, t_sub + h, nf, self.last_ovulation_t,
                &self.active_idx, &self.follicles,
            );
        }

        // --- 5. Write back state ---
        for (j, &idx) in self.active_idx.iter().enumerate() {
            self.follicles[idx].size = self.y_buf[j];
        }
        self.hormones.copy_from_slice(&self.y_buf[nf..]);
        self.t += dt;

        // --- 6. Follicle fate decisions ---
        self.dy_buf.fill(0.0);
        rhs(self.t, &self.y_buf, &mut self.dy_buf, nf, self.last_ovulation_t,
            f_k, &self.para_ode, &self.params, &self.fsh_sens_buf, &self.active_idx, &self.follicles);

        let y_lh = self.hormones[5];
        let mut ovulated = false;

        for j in 0..self.active_idx.len() {
            let idx = self.active_idx[j];
            let f = &mut self.follicles[idx];
            let growth_rate = self.dy_buf[j];
            let age = self.t - f.birth_time;

            match f.destiny {
                -1 => {
                    if growth_rate <= 0.01 || (growth_rate <= 0.1 && age >= 2.0) {
                        f.destiny = -2;
                        f.time_decrease = self.t;
                    } else if f.size >= self.para_ode.folmax {
                        if y_lh < 25.0 { f.destiny = 3; f.time_decrease = self.t; }
                        else { 
                            if self.t - self.last_ovulation_t > 10.0 && !ovulated {
                                f.destiny = 4; 
                                f.time_decrease = self.t; 
                            } else {
                                f.destiny = -2;
                                f.time_decrease = self.t;
                            }
                        }
                    }
                }
                3 => {
                    if y_lh >= 25.0 { 
                        if self.t - self.last_ovulation_t > 10.0 && !ovulated {
                            f.destiny = 4; 
                            f.time_decrease = self.t; 
                        } else {
                            f.destiny = -2;
                            f.time_decrease = self.t;
                        }
                    } else if self.t - f.time_decrease >= 5.0 { 
                        f.destiny = -2; 
                        f.time_decrease = self.t; 
                    }
                }
                4 => {
                    if self.t - f.time_decrease >= 0.5 && !ovulated {
                        f.destiny = 1;
                        self.last_ovulation_t = self.t;
                        self.cumulative_ovulations.push(self.t);
                        ovulated = true;
                    } else if ovulated {
                        // Another follicle already ovulated this step or just before
                        f.destiny = -2;
                        f.time_decrease = self.t;
                    }
                }
                _ => {}
            }
        }

        if ovulated {
            for f in &mut self.follicles {
                if f.destiny == -1 || f.destiny == 3 || f.destiny == 4 {
                    f.destiny = -2;
                    f.time_decrease = self.t;
                }
            }
        }

        for f in &mut self.follicles {
            if f.destiny == -2 && f.size > 0.1 {
                let dt_since = self.t - f.time_decrease;
                f.size *= (-0.05 * dt_since * dt).exp();
                if f.size < 0.5 { f.size = 0.0; }
            }
        }

        if self.follicles.len() > 50 {
            self.follicles.retain(|f| f.destiny == -1 || f.destiny == 3 || f.destiny == 4);
        }

        ovulated
    }
}

fn rhs(
    _t: f64,
    y: &[f64],
    dy: &mut [f64],
    nf: usize,
    _last_ovulation_t: f64,
    f_k: f64,
    po: &ParaOde,
    p: &[f64],
    fsh_sens: &[f64],
    active_idx: &[usize],
    follicles: &[Follicle],
) {
    dy.fill(0.0);

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

    // --- GnRH modulation via Zavala f_K (0 to 1) ---
    // f_K from the Zavala neuroendocrine network gates GnRH pulse frequency.
    // f_K = 1 → normal pulsing, f_K → 0 → suppressed by stress/cortisol.
    let yg_freq = 16.0 / (1.0 + (y_p4 / p[203]).powf(p[204]))
        * (1.0 + hill_pos(y_e2, p[205], p[206]))
        * f_k;

    let yg_mass = hill_pos(y_e2, p[208], p[209]) + hill_neg(y_e2, p[210], p[211]);

    let hp_e2 = hill_pos(y_e2, p[3], p[6]);
    let f_lh_prod = (p[1] + p[2] * hp_e2) / (1.0 + (y_p4 / p[4]).powf(p[7]));
    let f_lh_rel = p[16] + p[5] * hill_pos(y[igrec_a], p[8], p[9]);

    dy[irp_lh] = f_lh_prod - f_lh_rel * y[irp_lh];
    dy[ilh] = (1.0 / p[12]) * f_lh_rel * y[irp_lh] - p[230] * y_lh * y[irlh] - p[231] * y_lh;

    let hm_freq = 1.0 / (1.0 + (yg_freq / p[11]).powf(p[13]));
    let f_fsh_prod = (p[21] / (1.0 + (y_p4 / 5.0).powi(2))) * hm_freq;
    let f_fsh_rel = p[17] + p[28] * hill_pos(y[igrec_a], p[18], p[20]);

    dy[irp_fsh] = f_fsh_prod - f_fsh_rel * y[irp_fsh];
    dy[ifsh] = (1.0 / p[12]) * f_fsh_rel * y[irp_fsh] - (p[243] + p[241]) * y_fsh;

    dy[ifshfoll] = p[243] * y_fsh * p[246] - p[240] * y[ifshfoll] * y[irfsh];
    dy[irfsh] = p[242] * y[irfsh_des] - p[240] * y[ifshfoll] * y[irfsh];
    dy[ifshr] = p[240] * y[ifshfoll] * y[irfsh] - p[244] * y[ifshr];
    dy[irfsh_des] = p[244] * y[ifshr] - p[242] * y[irfsh_des];

    dy[ignrh] = p[301] * yg_mass * yg_freq - p[302] * y[ignrh] * y[irec_ga] + p[303] * y[igrec_a] - p[300] * y[ignrh];
    dy[irec_ga] = p[303] * y[igrec_a] - p[302] * y[ignrh] * y[irec_ga] - p[306] * y[irec_ga] + p[307] * y[irec_gi];
    dy[irec_gi] = p[311] + p[306] * y[irec_ga] - p[307] * y[irec_gi] + p[305] * y[igrec_i] - p[308] * y[irec_gi];
    dy[igrec_a] = p[302] * y[ignrh] * y[irec_ga] - p[303] * y[igrec_a] - p[309] * y[igrec_a] + p[310] * y[igrec_i];
    dy[igrec_i] = p[309] * y[igrec_a] - p[310] * y[igrec_i] - p[305] * y[igrec_i] - p[304] * y[igrec_i];

    let fsh_rez_comp = y[ifshr].max(0.0);
    let sum_v: f64 = (0..nf).map(|j| y[j].max(0.0).powf(po.v)).sum();

    for j in 0..nf {
        let d = follicles[active_idx[j]].destiny;
        if d == -2 || d == 4 { continue; }
        let x = y[j].max(0.0);
        let sens = fsh_sens[j];
        let f_fsh = fsh_rez_comp / (fsh_rez_comp + sens);
        // Stronger P4 inhibition of growth: 1 / (1 + (P4/2)^4)
        let gamma_p4 = po.gamma * (1.0 / (1.0 + (y_p4 / 2.0).powi(4)));
        let kappa_fsh = po.k * hill_neg_raw(fsh_rez_comp, 0.85, 25.0);
        dy[j] = f_fsh * (po.xi - x) * x * (gamma_p4 - kappa_fsh * (sum_v - po.mu * x.powf(po.v)));
    }
}

fn compute_sf(y: &[f64], nf: usize, active_idx: &[usize], follicles: &[Follicle]) -> f64 {
    if nf == 0 { return 0.0; }
    let mut sum = 0.0;
    for j in 0..nf {
        let d = follicles[active_idx[j]].destiny;
        if d == -2 || d == 4 { continue; }
        let x = y[j].max(0.0);
        let x5 = x.powi(5);
        sum += (x5 / (x5 + 15.0f64.powi(5))) * x.powi(2);
    }
    PI * sum
}

fn project_algebraic(y: &mut [f64], t: f64, nf: usize, last_ovulation_t: f64, active_idx: &[usize], follicles: &[Follicle]) {
    let ie2 = nf + 1;
    let ip4 = nf + 3;
    let sf = compute_sf(y, nf, active_idx, follicles);
    
    let dt_ov = t - last_ovulation_t;
    let p4_target = if dt_ov >= 0.0 && dt_ov < 20.0 {
        let sigma = if dt_ov < 7.0 { 3.0 } else { 6.0 };
        18.0 * (-(dt_ov - 7.0).powi(2) / (2.0 * sigma * sigma)).exp()
    } else {
        0.1
    };

    let e2_luteal = if dt_ov >= 0.0 && dt_ov < 20.0 {
        let sigma = if dt_ov < 7.0 { 4.0 } else { 7.0 };
        120.0 * (-(dt_ov - 8.0).powi(2) / (2.0 * sigma * sigma)).exp()
    } else {
        0.0
    };

    y[ie2] = 13.0 * (1.0 + 0.02 * sf) + e2_luteal;
    y[ip4] = p4_target.max(0.1);
}

#[inline]
fn hill_pos(s: f64, k: f64, n: f64) -> f64 {
    if s <= 0.0 || k <= 0.0 { return 0.0; }
    let sn = s.powf(n);
    let kn = k.powf(n);
    sn / (kn + sn)
}

#[inline]
fn hill_neg(s: f64, k: f64, n: f64) -> f64 {
    if s <= 0.0 || k <= 0.0 { return 1.0; }
    let sn = s.powf(n);
    let kn = k.powf(n);
    kn / (kn + sn)
}

#[inline]
fn hill_pos_raw(s: f64, k: f64, n: f64) -> f64 {
    let sn = s.powf(n);
    let kn = k.powf(n);
    sn / (kn + sn)
}

#[inline]
fn hill_neg_raw(s: f64, k: f64, n: f64) -> f64 {
    let sn = s.powf(n);
    let kn = k.powf(n);
    kn / (kn + sn)
}
