use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;

use rand_distr::{Distribution, LogNormal, Normal};

use crate::cycle::{CycleParams, CycleState, advance_day, create_population, init_cycle};
use crate::mechanistic::MenstrualMechanisticModel;
use crate::stress::{StressConfig, generate_independent_stress, generate_stress_timeline};
use crate::zavala::ZavalaParams;
use crate::synchrony::{
    compute_phases, mean_resultant_length, onset_synchrony_index, rayleigh_test,
};

#[derive(Clone, Debug)]
pub struct ExperimentConfig {
    pub n_individuals: usize,
    pub n_days: usize,
    pub n_trials: usize,
    pub seed: u64,
    pub cycle_params: CycleParams,
    pub stress_config: StressConfig,
    pub burnin_days: usize,
    pub heterogeneity: f64,
    /// Zeitgeber coupling strength for mechanistic model (default 0.5).
    /// Heuristic; sweep for sensitivity analysis.
    pub k_zeitgeber: f64,
}

impl Default for ExperimentConfig {
    fn default() -> Self {
        Self {
            n_individuals: 6,
            n_days: 365 * 3,
            n_trials: 200,
            seed: 42,
            cycle_params: CycleParams::default(),
            stress_config: StressConfig::default(),
            burnin_days: 90,
            heterogeneity: 0.15,
            k_zeitgeber: 0.5,
        }
    }
}

#[derive(Clone, Debug)]
pub struct TrialResult {
    pub treatment_r: f64,
    pub control_r: f64,
    pub treatment_osi: f64,
    pub control_osi: f64,
    pub treatment_rayleigh_p: f64,
    pub control_rayleigh_p: f64,
}

pub struct ExperimentResult {
    pub trials: Vec<TrialResult>,
}

fn run_trial(config: &ExperimentConfig, trial_seed: u64) -> TrialResult {
    let mut rng = ChaCha8Rng::seed_from_u64(trial_seed);

    let pop_treatment = create_population(
        config.n_individuals,
        &config.cycle_params,
        config.heterogeneity,
        &mut rng,
    );

    let mut rng_ctrl = ChaCha8Rng::seed_from_u64(trial_seed.wrapping_add(1_000_000));
    let pop_control: Vec<(CycleParams, CycleState)> = pop_treatment
        .iter()
        .map(|(p, _)| {
            let s = init_cycle(p, &mut rng_ctrl);
            (p.clone(), s)
        })
        .collect();

    let stress_treatment = generate_stress_timeline(
        config.n_individuals,
        config.n_days,
        &config.stress_config,
        &mut ChaCha8Rng::seed_from_u64(trial_seed.wrapping_add(2_000_000)),
    );
    let stress_control = generate_independent_stress(
        config.n_individuals,
        config.n_days,
        &config.stress_config,
        &mut ChaCha8Rng::seed_from_u64(trial_seed.wrapping_add(3_000_000)),
    );

    let mut pop_t = pop_treatment;
    let mut pop_c = pop_control;
    let mut rng_t = ChaCha8Rng::seed_from_u64(trial_seed.wrapping_add(4_000_000));
    let mut rng_c = ChaCha8Rng::seed_from_u64(trial_seed.wrapping_add(5_000_000));

    for day in 0..config.n_days {
        for (i, (params, state)) in pop_t.iter_mut().enumerate() {
            let s = stress_treatment[i][day];
            if advance_day(state, params, s, &mut rng_t) {
                state.onset_days.push(day as u32);
            }
        }
        for (i, (params, state)) in pop_c.iter_mut().enumerate() {
            let s = stress_control[i][day];
            if advance_day(state, params, s, &mut rng_c) {
                state.onset_days.push(day as u32);
            }
        }
    }

    let phases_t = compute_phases(&pop_t);
    let phases_c = compute_phases(&pop_c);
    let r_t = mean_resultant_length(&phases_t);
    let r_c = mean_resultant_length(&phases_c);
    let (_, p_t) = rayleigh_test(&phases_t);
    let (_, p_c) = rayleigh_test(&phases_c);

    let burnin = config.burnin_days as u32;
    let onsets_t: Vec<Vec<u32>> = pop_t
        .iter()
        .map(|(_, s)| s.onset_days.iter().copied().filter(|&d| d >= burnin).collect())
        .collect();
    let onsets_c: Vec<Vec<u32>> = pop_c
        .iter()
        .map(|(_, s)| s.onset_days.iter().copied().filter(|&d| d >= burnin).collect())
        .collect();

    let osi_t = onset_synchrony_index(&onsets_t, 5.0);
    let osi_c = onset_synchrony_index(&onsets_c, 5.0);

    TrialResult {
        treatment_r: r_t,
        control_r: r_c,
        treatment_osi: osi_t,
        control_osi: osi_c,
        treatment_rayleigh_p: p_t,
        control_rayleigh_p: p_c,
    }
}

#[derive(Clone, Debug)]
pub struct MechanisticTrialResult {
    pub treatment_osi: f64,
    pub control_osi: f64,
    /// Mean pairwise Pearson r of treatment stress timelines.
    /// Validates shared_fraction → realized correlation mapping.
    pub treatment_stress_corr: f64,
}

/// Mean pairwise Pearson correlation across all individual pairs.
fn mean_pairwise_correlation(timelines: &[Vec<f64>]) -> f64 {
    let n = timelines.len();
    if n < 2 {
        return 0.0;
    }
    let mut sum_r = 0.0;
    let mut count = 0;
    for i in 0..n {
        let len_i = timelines[i].len();
        let mean_i: f64 = timelines[i].iter().sum::<f64>() / len_i as f64;
        for j in (i + 1)..n {
            let len = len_i.min(timelines[j].len());
            let mean_j: f64 = timelines[j][..len].iter().sum::<f64>() / len as f64;
            let mut cov = 0.0;
            let mut var_i = 0.0;
            let mut var_j = 0.0;
            for k in 0..len {
                let di = timelines[i][k] - mean_i;
                let dj = timelines[j][k] - mean_j;
                cov += di * dj;
                var_i += di * di;
                var_j += dj * dj;
            }
            let denom = (var_i * var_j).sqrt();
            if denom > 0.0 {
                sum_r += cov / denom;
            }
            count += 1;
        }
    }
    if count > 0 { sum_r / count as f64 } else { 0.0 }
}

pub fn run_mechanistic_trial(config: &ExperimentConfig, trial_seed: u64) -> MechanisticTrialResult {
    let mut rng = ChaCha8Rng::seed_from_u64(trial_seed);

    let steps_per_day: usize = 4;
    let dt = 1.0 / steps_per_day as f64;

    // Sub-day stress resolution for mechanistic model
    let mut mech_stress_config = config.stress_config.clone();
    mech_stress_config.steps_per_day = steps_per_day;

    let stress_treatment = generate_stress_timeline(
        config.n_individuals,
        config.n_days,
        &mech_stress_config,
        &mut ChaCha8Rng::seed_from_u64(trial_seed.wrapping_add(2_000_000)),
    );
    let stress_control = generate_independent_stress(
        config.n_individuals,
        config.n_days,
        &mech_stress_config,
        &mut ChaCha8Rng::seed_from_u64(trial_seed.wrapping_add(3_000_000)),
    );

    // Realized stress correlation for validation
    let treatment_stress_corr = mean_pairwise_correlation(&stress_treatment);

    // Between-person heterogeneity distributions
    let z_normal = Normal::new(0.0, 1.0).unwrap();
    let mut z_defaults = ZavalaParams::default();
    z_defaults.k_zeitgeber = config.k_zeitgeber;
    let epsilon_dist = LogNormal::new(z_defaults.epsilon.ln(), 0.5).unwrap();

    let mut pop_t: Vec<MenstrualMechanisticModel> = (0..config.n_individuals)
        .map(|i| {
            let mut m = MenstrualMechanisticModel::new(trial_seed.wrapping_add(i as u64), &mut rng);
            // Heterogeneity: stress→CORT coupling (LogNormal, right-skewed)
            m.zavala.params.epsilon = epsilon_dist.sample(&mut rng).max(0.01);
            // Heterogeneity: intrinsic period (±0.38h at 95% CI; Duffy 2011, forced desynchrony)
            m.zavala.params.omega_h0 = z_defaults.omega_h0 * (1.0 + z_normal.sample(&mut rng) * 0.008);
            // Heterogeneity: stress→circadian coupling
            m.zavala.params.alpha = (z_defaults.alpha + z_normal.sample(&mut rng) * 0.01).max(0.0);
            // Zeitgeber coupling (from config for sensitivity sweep)
            m.zavala.params.k_zeitgeber = z_defaults.k_zeitgeber;
            m.randomize_initial_state(&mut rng);
            m
        })
        .collect();

    // Treatment: shared zeitgeber phase (cohabitants on same light schedule).
    // Continuous Kuramoto forcing entrains all circadian oscillators to the
    // same external phase, maintaining coherence despite HPA perturbations.
    let shared_phi_z0 = pop_t[0].zavala.phi_z0;
    for m in &mut pop_t {
        m.zavala.phi_z0 = shared_phi_z0;
    }

    let mut pop_c: Vec<MenstrualMechanisticModel> = (0..config.n_individuals)
        .map(|i| {
            let mut m = MenstrualMechanisticModel::new(trial_seed.wrapping_add(i as u64 + 1000), &mut rng);
            m.zavala.params.epsilon = epsilon_dist.sample(&mut rng).max(0.01);
            m.zavala.params.omega_h0 = z_defaults.omega_h0 * (1.0 + z_normal.sample(&mut rng) * 0.008);
            m.zavala.params.alpha = (z_defaults.alpha + z_normal.sample(&mut rng) * 0.01).max(0.0);
            m.zavala.params.k_zeitgeber = z_defaults.k_zeitgeber;
            m.randomize_initial_state(&mut rng);
            m
        })
        .collect();
    // Control: independent circadian phases (already the default)

    let mut rng_t = ChaCha8Rng::seed_from_u64(trial_seed.wrapping_add(4_000_000));
    let mut rng_c = ChaCha8Rng::seed_from_u64(trial_seed.wrapping_add(5_000_000));

    for day in 0..config.n_days {
        for step in 0..steps_per_day {
            let stress_idx = day * steps_per_day + step;
            for (i, model) in pop_t.iter_mut().enumerate() {
                model.stress_input = stress_treatment[i][stress_idx];
                model.step(dt, &mut rng_t);
            }

            for (i, model) in pop_c.iter_mut().enumerate() {
                model.stress_input = stress_control[i][stress_idx];
                model.step(dt, &mut rng_c);
            }
        }
    }

    let burnin = config.burnin_days as f64;
    let onsets_t: Vec<Vec<u32>> = pop_t
        .iter()
        .map(|m| {
            m.cumulative_ovulations
                .iter()
                .filter(|&&t| t >= burnin)
                .map(|&t| t as u32)
                .collect()
        })
        .collect();
    let onsets_c: Vec<Vec<u32>> = pop_c
        .iter()
        .map(|m| {
            m.cumulative_ovulations
                .iter()
                .filter(|&&t| t >= burnin)
                .map(|&t| t as u32)
                .collect()
        })
        .collect();

    MechanisticTrialResult {
        treatment_osi: onset_synchrony_index(&onsets_t, 5.0),
        control_osi: onset_synchrony_index(&onsets_c, 5.0),
        treatment_stress_corr,
    }
}

pub fn run_experiment(config: &ExperimentConfig) -> ExperimentResult {
    let mut seed_rng = ChaCha8Rng::seed_from_u64(config.seed);
    let trial_seeds: Vec<u64> = (0..config.n_trials)
        .map(|_| rand::Rng::random(&mut seed_rng))
        .collect();

    let trials: Vec<TrialResult> = trial_seeds
        .par_iter()
        .map(|&seed| run_trial(config, seed))
        .collect();

    ExperimentResult { trials }
}

pub fn run_mechanistic_experiment(config: &ExperimentConfig) -> Vec<MechanisticTrialResult> {
    let mut seed_rng = ChaCha8Rng::seed_from_u64(config.seed);
    let trial_seeds: Vec<u64> = (0..config.n_trials)
        .map(|_| rand::Rng::random(&mut seed_rng))
        .collect();

    trial_seeds
        .par_iter()
        .map(|&seed| run_mechanistic_trial(config, seed))
        .collect()
}

pub fn run_mechanistic_sweep(
    base: &ExperimentConfig,
    param_name: &str,
    param_values: &[f64],
) -> Vec<Vec<MechanisticTrialResult>> {
    param_values
        .iter()
        .map(|&val| {
            let mut cfg = base.clone();
            match param_name {
                "shared_fraction" => cfg.stress_config.shared_fraction = val,
                "theta" => cfg.stress_config.theta = val,
                "mu" => cfg.stress_config.mu = val,
                "sigma_total" => cfg.stress_config.sigma_total = val,
                "k_zeitgeber" => cfg.k_zeitgeber = val,
                "stress_sensitivity" => {
                    cfg.cycle_params.stress_sensitivity = val;
                }
                _ => {}
            }
            run_mechanistic_experiment(&cfg)
        })
        .collect()
}