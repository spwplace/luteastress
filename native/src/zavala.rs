use rand::Rng;
use std::f64::consts::PI;

/// Parameters from Zavala et al. (2020) Table S2.
/// All rates converted from min⁻¹ to day⁻¹ (×1440) where noted.
pub struct ZavalaParams {
    /// Circadian base frequency: π/(24×60) /min → π/24 /day
    pub omega_h0: f64,
    /// CORT ultradian frequency: π/75 /min → 1440π/75 /day
    pub omega_c0: f64,
    /// Slow genomic CORT rate: 0.005 /min → 7.2 /day
    pub delta: f64,
    /// Circadian coupling strength
    pub alpha: f64,
    /// CORT-circadian coupling
    pub beta: f64,
    /// Estrous phase skew
    pub sigma: f64,
    /// Stress half-saturation (Hill)
    pub k_s: f64,
    /// Cortisol half-saturation (Hill)
    pub k_c: f64,
    /// Sigmoid steepness
    pub h_sig: f64,
    /// Stress weight in g()
    pub m: f64,
    /// Cortisol weight in g()
    pub l: f64,
    /// Hill exponent
    pub n_hill: f64,
    /// Stress amplitude
    pub p_a: f64,
    /// Stress→amplitude coupling
    pub epsilon: f64,
    /// Stress window phase start (fraction of 2π; dark onset = π)
    pub phi_s: f64,
    /// Stress window width (radians)
    pub theta_s: f64,
}

impl Default for ZavalaParams {
    fn default() -> Self {
        Self {
            omega_h0: PI / 24.0,
            omega_c0: 1440.0 * PI / 75.0,
            delta: 7.2,
            alpha: 0.05,
            beta: 0.9,
            sigma: 0.5,
            k_s: 0.6,
            k_c: 1.7,
            h_sig: 50.0,
            m: 13.0,
            l: 6.0,
            n_hill: 4.0,
            p_a: 1.0,
            epsilon: 0.1,
            // Stress window centered around activity period (φ_H ∈ [0, π])
            phi_s: 0.0,
            theta_s: PI,
        }
    }
}

/// State variables for the Zavala neuroendocrine network.
/// 6 dynamical variables tracking the circadian→HPA→KNDy→GnRH chain.
pub struct ZavalaState {
    /// Circadian phase (0–2π), period ≈ 24h
    pub phi_h: f64,
    /// CORT ultradian phase (0–2π), period ≈ 2.5h
    pub phi_c: f64,
    /// CORT amplitude — quasi-steady-state (algebraic)
    pub a_c: f64,
    /// Slow genomic cortisol: dC̃/dt = δ(C − C̃)
    pub c_tilde: f64,
    /// Neurokinin B (KNDy): dN/dt = f_E − N
    pub n_kndy: f64,
    /// Dynorphin (KNDy): dD/dt = f_E − D
    pub d_kndy: f64,
    /// Parameters
    pub params: ZavalaParams,
}

impl ZavalaState {
    /// Create a new ZavalaState with randomized initial phases.
    pub fn new(params: ZavalaParams, rng: &mut impl Rng) -> Self {
        let phi_h = rng.random_range(0.0..(2.0 * PI));
        let phi_c = rng.random_range(0.0..(2.0 * PI));
        let f_h = phi_h.sin().powi(2);
        let a_c = f_h; // initial quasi-steady with no stress
        let c = a_c * phi_c.sin().powi(2);
        Self {
            phi_h,
            phi_c,
            a_c,
            c_tilde: c,
            n_kndy: 0.5,
            d_kndy: 0.5,
            params,
        }
    }

    /// S1: Stress gating by circadian phase.
    /// s(φ_H) = p_A · H(φ_H − φ_s) · H(φ_s + θ_s − φ_H)
    /// where H is the Heaviside step function.
    /// The external `stress_input` modulates the amplitude.
    fn s_gated(&self, stress_input: f64) -> f64 {
        // Wrap phase into [0, 2π)
        let phi = self.phi_h % (2.0 * PI);
        let phi = if phi < 0.0 { phi + 2.0 * PI } else { phi };
        let p = &self.params;
        // Heaviside: active when φ_s ≤ φ_H ≤ φ_s + θ_s
        let in_window = if p.phi_s + p.theta_s <= 2.0 * PI {
            phi >= p.phi_s && phi <= p.phi_s + p.theta_s
        } else {
            // Window wraps around 2π
            phi >= p.phi_s || phi <= (p.phi_s + p.theta_s) % (2.0 * PI)
        };
        if in_window {
            p.p_a * stress_input
        } else {
            0.0
        }
    }

    /// S2: Circadian modulation f_H(φ_H) = sin²(φ_H)
    /// (a=0, b=1 physiological case from Zavala)
    #[inline]
    fn f_h(&self) -> f64 {
        self.phi_h.sin().powi(2)
    }

    /// S3: Estrous function with skew
    /// f_E(φ_E) = sin²(φ_E − σ·sin²(φ_E))
    #[inline]
    fn f_e(phi_e: f64, sigma: f64) -> f64 {
        let arg = phi_e - sigma * phi_e.sin().powi(2);
        arg.sin().powi(2)
    }

    /// S7: Instantaneous cortisol C = A_C · sin²(φ_C)
    #[inline]
    pub fn cortisol(&self) -> f64 {
        self.a_c * self.phi_c.sin().powi(2)
    }

    /// S10: Threshold function g(s, C̃)
    /// g = (1/h) · (1 + m·sⁿ/(K_sⁿ+sⁿ) + l·C̃ⁿ/(K_Cⁿ+C̃ⁿ))
    fn g_threshold(&self, stress_input: f64) -> f64 {
        let p = &self.params;
        let n = p.n_hill;
        let s = self.s_gated(stress_input);
        let s_hill = if s > 0.0 {
            s.powf(n) / (p.k_s.powf(n) + s.powf(n))
        } else {
            0.0
        };
        let c_hill = if self.c_tilde > 0.0 {
            self.c_tilde.powf(n) / (p.k_c.powf(n) + self.c_tilde.powf(n))
        } else {
            0.0
        };
        (1.0 / p.h_sig) * (1.0 + p.m * s_hill + p.l * c_hill)
    }

    /// Sigmoid σ_sig(x) = 1/(1+e^{-x})
    #[inline]
    fn sigma_sig(x: f64) -> f64 {
        1.0 / (1.0 + (-x).exp())
    }

    /// S9: KNDy gating function f_K
    /// f_K = σ_sig(h·(x−g)) · σ_sig(h·(y−g))
    /// where x = D − N/2, y = N − D/2
    /// Returns a value in (0, 1) that modulates GnRH pulse frequency.
    pub fn compute_f_k(&self, stress_input: f64) -> f64 {
        let p = &self.params;
        let g = self.g_threshold(stress_input);
        let x = self.d_kndy - self.n_kndy / 2.0;
        let y = self.n_kndy - self.d_kndy / 2.0;
        Self::sigma_sig(p.h_sig * (x - g)) * Self::sigma_sig(p.h_sig * (y - g))
    }

    /// Advance the Zavala network by one substep.
    ///
    /// # Arguments
    /// * `dt` - timestep in days
    /// * `stress_input` - current stress level (0 = none, >0 = active)
    /// * `phi_e` - estrous/menstrual phase (0–2π), computed from cycle model
    pub fn step(&mut self, dt: f64, stress_input: f64, phi_e: f64) {
        let p = &self.params;
        let s = self.s_gated(stress_input);

        // Phase equations (always stable — just accumulation)
        // dφ_H/dt = ω_H0 + β·(A_C·sin²(φ_C) − ⟨C⟩)
        // Simplified: circadian phase advances at base rate with CORT coupling
        let c = self.cortisol();
        let d_phi_h = p.omega_h0 + p.beta * (c - self.c_tilde);
        self.phi_h = (self.phi_h + d_phi_h * dt) % (2.0 * PI);
        if self.phi_h < 0.0 {
            self.phi_h += 2.0 * PI;
        }

        // dφ_C/dt = ω_C0 − α·s(φ_H)
        let d_phi_c = p.omega_c0 - p.alpha * s;
        self.phi_c = (self.phi_c + d_phi_c * dt) % (2.0 * PI);
        if self.phi_c < 0.0 {
            self.phi_c += 2.0 * PI;
        }

        // Quasi-steady-state for a_c (τ_a = 1 min, way too fast for our dt)
        // Steady state of: dA_C/dt = (1/τ_a)·(f_H·(1 + ε·s) − A_C)
        let f_h = self.f_h();
        self.a_c = f_h * (1.0 + p.epsilon * s);

        // S8: dC̃/dt = δ·(C − C̃) — slow genomic cortisol (Euler, stable)
        let c_new = self.cortisol(); // recompute with updated a_c
        let d_c_tilde = p.delta * (c_new - self.c_tilde);
        self.c_tilde = (self.c_tilde + d_c_tilde * dt).max(0.0);

        // S5–S6: KNDy neuron tracking (Euler, very stable)
        let f_e = Self::f_e(phi_e, p.sigma);
        let d_n = f_e - self.n_kndy;
        let d_d = f_e - self.d_kndy;
        self.n_kndy = (self.n_kndy + d_n * dt).max(0.0);
        self.d_kndy = (self.d_kndy + d_d * dt).max(0.0);
    }
}
