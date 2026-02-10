use rand::Rng;
use rand_distr::{Distribution, StandardNormal};

#[derive(Clone, Debug)]
pub struct StressConfig {
    /// OU mean-reversion rate (1/day). Calibrated from ESM perceived stress
    /// AR(1)=0.277 at 4h intervals (PMC9986209) → θ ≈ 7.7/day (rounded to 7.0), t½ ≈ 2.2h.
    pub theta: f64,
    /// Mean stress level in [0, 1]. The OU processes fluctuate around this.
    pub mu: f64,
    /// Total OU noise amplitude. Stationary SD = sigma_total / sqrt(2*theta).
    /// With theta=7, sigma_total=0.45 gives SD ≈ 0.12.
    pub sigma_total: f64,
    /// Fraction of stress variance shared between treatment individuals (0–1).
    /// Fraction of total stress variance from shared OU component.
    /// Theoretical Pearson r = shared_fraction before clamping;
    /// realized r is lower due to clamp(0,1) truncation.
    /// This is the key experimental variable to sweep.
    /// Empirical anchors: r=0.2-0.3 romantic partners (Saxbe & Repetti 2010),
    /// r=0.18-0.61 family dyads (Elicit review). Unknown for roommates.
    pub shared_fraction: f64,
    /// Temporal resolution: 1 = daily, 4 = 6-hourly, etc.
    pub steps_per_day: usize,
}

impl Default for StressConfig {
    fn default() -> Self {
        Self {
            theta: 7.0,
            mu: 0.25,
            sigma_total: 0.45,
            shared_fraction: 0.3,
            steps_per_day: 1,
        }
    }
}

/// Generate a zero-mean Ornstein-Uhlenbeck process using the exact
/// (analytical) conditional distribution. Unconditionally stable for any dt.
///
/// The process satisfies dS = -theta * S * dt + sigma * dW,
/// with stationary distribution N(0, sigma^2 / (2*theta)).
fn generate_ou_process(
    n_steps: usize,
    theta: f64,
    sigma: f64,
    dt: f64,
    rng: &mut impl Rng,
) -> Vec<f64> {
    let mut timeline = Vec::with_capacity(n_steps);
    if n_steps == 0 || sigma == 0.0 {
        timeline.resize(n_steps, 0.0);
        return timeline;
    }

    let decay = (-theta * dt).exp();
    let var = sigma * sigma * (1.0 - (-2.0 * theta * dt).exp()) / (2.0 * theta);
    let sd = var.sqrt();

    let mut s = 0.0;
    for _ in 0..n_steps {
        timeline.push(s);
        let z: f64 = StandardNormal.sample(rng);
        s = s * decay + sd * z;
    }
    timeline
}

/// Generate stress timelines for the treatment arm (shared stress).
///
/// Each individual's stress = clamp(mu + s_shared + s_individual, 0, 1).
/// The shared OU component is the same for all individuals in the trial;
/// individual OU components are independent per person.
pub fn generate_stress_timeline(
    n_individuals: usize,
    n_days: usize,
    config: &StressConfig,
    rng: &mut impl Rng,
) -> Vec<Vec<f64>> {
    let spd = config.steps_per_day.max(1);
    let n_steps = n_days * spd;
    let dt = 1.0 / spd as f64;

    // Decompose total noise into shared + individual components.
    // Var(shared) + Var(individual) = Var(total) at stationarity.
    // Since OU stationary var = sigma^2/(2*theta), this means
    // sigma_shared^2 + sigma_indiv^2 = sigma_total^2.
    let sf = config.shared_fraction.clamp(0.0, 1.0);
    let sigma_shared = config.sigma_total * sf.sqrt();
    let sigma_indiv = config.sigma_total * (1.0 - sf).sqrt();

    let shared = generate_ou_process(n_steps, config.theta, sigma_shared, dt, rng);

    let mut stress = Vec::with_capacity(n_individuals);
    for _ in 0..n_individuals {
        let indiv = generate_ou_process(n_steps, config.theta, sigma_indiv, dt, rng);
        let timeline: Vec<f64> = shared
            .iter()
            .zip(indiv.iter())
            .map(|(&s, &i)| (config.mu + s + i).clamp(0.0, 1.0))
            .collect();
        stress.push(timeline);
    }
    stress
}

/// Generate stress timelines for the control arm (independent stress).
///
/// Each individual's stress = clamp(mu + s_individual, 0, 1).
/// All OU components are independent. The total sigma is set so that the
/// marginal distribution matches the treatment arm exactly.
pub fn generate_independent_stress(
    n_individuals: usize,
    n_days: usize,
    config: &StressConfig,
    rng: &mut impl Rng,
) -> Vec<Vec<f64>> {
    let spd = config.steps_per_day.max(1);
    let n_steps = n_days * spd;
    let dt = 1.0 / spd as f64;

    // Control uses sigma_total directly (= sqrt(sigma_shared^2 + sigma_indiv^2))
    // so marginal variance matches treatment.
    let mut stress = Vec::with_capacity(n_individuals);
    for _ in 0..n_individuals {
        let indiv = generate_ou_process(n_steps, config.theta, config.sigma_total, dt, rng);
        let timeline: Vec<f64> = indiv
            .iter()
            .map(|&i| (config.mu + i).clamp(0.0, 1.0))
            .collect();
        stress.push(timeline);
    }
    stress
}
