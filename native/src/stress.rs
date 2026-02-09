use rand::Rng;
use rand_distr::{Distribution, LogNormal, Normal, Poisson};

#[derive(Clone, Debug)]
pub struct StressConfig {
    pub shared_event_rate: f64,
    pub individual_event_rate: f64,
    pub magnitude_mean_log: f64,
    pub magnitude_sd_log: f64,
    pub duration_mean: f64,
    pub duration_sd: f64,
    pub shared_exposure_prob: f64,
}

impl Default for StressConfig {
    fn default() -> Self {
        Self {
            shared_event_rate: 0.05,
            individual_event_rate: 0.07,
            magnitude_mean_log: 0.0,
            magnitude_sd_log: 0.7,
            duration_mean: 2.0,
            duration_sd: 1.0,
            shared_exposure_prob: 0.85,
        }
    }
}

/// Generate stress matrix: [n_individuals][n_days]
pub fn generate_stress_timeline(
    n_individuals: usize,
    n_days: usize,
    config: &StressConfig,
    rng: &mut impl Rng,
) -> Vec<Vec<f64>> {
    let mut stress = vec![vec![0.0f64; n_days]; n_individuals];
    let mag_dist = LogNormal::new(config.magnitude_mean_log, config.magnitude_sd_log).unwrap();
    let dur_dist = Normal::new(config.duration_mean, config.duration_sd).unwrap();
    let personal_dist = LogNormal::new(0.0, 0.3).unwrap();

    // Shared events
    let lambda_shared = config.shared_event_rate * n_days as f64;
    if lambda_shared > 0.0 {
        let n_shared: u64 = Poisson::new(lambda_shared).unwrap().sample(rng) as u64;
        for _ in 0..n_shared {
            let onset: usize = rng.random_range(0..n_days);
            let magnitude: f64 = mag_dist.sample(rng);
            let duration = (dur_dist.sample(rng) as usize).max(1);

            for day_offset in 0..duration {
                let day = onset + day_offset;
                if day >= n_days {
                    break;
                }
                let decay = (-0.5 * day_offset as f64).exp();
                for ind in stress.iter_mut() {
                    if rng.random_range(0.0..1.0) < config.shared_exposure_prob {
                        let personal: f64 = personal_dist.sample(rng);
                        ind[day] += magnitude * decay * personal;
                    }
                }
            }
        }
    }

    // Individual events
    let lambda_indiv = config.individual_event_rate * n_days as f64;
    if lambda_indiv > 0.0 {
        for ind in stress.iter_mut() {
            let n_events = Poisson::new(lambda_indiv).unwrap().sample(rng) as u64;
            for _ in 0..n_events {
                let onset: usize = rng.random_range(0..n_days);
                let magnitude: f64 = mag_dist.sample(rng);
                let duration = (dur_dist.sample(rng) as usize).max(1);

                for day_offset in 0..duration {
                    let day = onset + day_offset;
                    if day >= n_days {
                        break;
                    }
                    let decay = (-0.5 * day_offset as f64).exp();
                    let personal: f64 = personal_dist.sample(rng);
                    ind[day] += magnitude * decay * personal;
                }
            }
        }
    }

    stress
}

/// Control: all-individual stress (same total rate, no shared component)
pub fn generate_independent_stress(
    n_individuals: usize,
    n_days: usize,
    config: &StressConfig,
    rng: &mut impl Rng,
) -> Vec<Vec<f64>> {
    let control_config = StressConfig {
        shared_event_rate: 0.0,
        individual_event_rate: config.shared_event_rate + config.individual_event_rate,
        shared_exposure_prob: 0.0,
        ..config.clone()
    };
    generate_stress_timeline(n_individuals, n_days, &control_config, rng)
}
