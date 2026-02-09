use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use rayon::prelude::*;

use crate::cycle::{CycleParams, CycleState, advance_day, create_population, init_cycle};
use crate::mechanistic::MenstrualMechanisticModel;
use crate::stress::{StressConfig, generate_independent_stress, generate_stress_timeline};
use crate::synchrony::{
    compute_phases, mean_resultant_length, onset_synchrony_index, rayleigh_test,
};

#[derive(Clone, Debug)]
pub struct MechanisticTrialResult {
    pub treatment_osi: f64,
    pub control_osi: f64,
}

pub fn run_mechanistic_trial(config: &ExperimentConfig, trial_seed: u64) -> MechanisticTrialResult {
    let mut rng = ChaCha8Rng::seed_from_u64(trial_seed);

    // Generate stress
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

    // Initialize populations
    let mut pop_t: Vec<MenstrualMechanisticModel> = (0..config.n_individuals)
        .map(|i| MenstrualMechanisticModel::new(trial_seed.wrapping_add(i as u64)))
        .collect();
    let mut pop_c: Vec<MenstrualMechanisticModel> = (0..config.n_individuals)
        .map(|i| MenstrualMechanisticModel::new(trial_seed.wrapping_add(i as u64 + 1000)))
        .collect();

    let mut rng_t = ChaCha8Rng::seed_from_u64(trial_seed.wrapping_add(4_000_000));
    let mut rng_c = ChaCha8Rng::seed_from_u64(trial_seed.wrapping_add(5_000_000));

    // Run simulation
    // We'll use 0.25 day steps
    let steps_per_day = 4;
    let dt = 1.0 / steps_per_day as f64;

    for day in 0..config.n_days {
        for _ in 0..steps_per_day {
            for (i, model) in pop_t.iter_mut().enumerate() {
                model.stress_input = stress_treatment[i][day];
                model.step(dt, &mut rng_t);
            }
            for (i, model) in pop_c.iter_mut().enumerate() {
                model.stress_input = stress_control[i][day];
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
    }
}

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

    // Create population
    let pop_treatment = create_population(
        config.n_individuals,
        &config.cycle_params,
        config.heterogeneity,
        &mut rng,
    );

    // Clone params, re-init states for control
    let mut rng_ctrl = ChaCha8Rng::seed_from_u64(trial_seed.wrapping_add(1_000_000));
    let pop_control: Vec<(CycleParams, CycleState)> = pop_treatment
        .iter()
        .map(|(p, _)| {
            let s = init_cycle(p, &mut rng_ctrl);
            (p.clone(), s)
        })
        .collect();

    // Generate stress
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

    // Run simulation
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

    // Compute metrics
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

pub fn run_experiment(config: &ExperimentConfig) -> ExperimentResult {
    // Generate per-trial seeds deterministically
    let mut seed_rng = ChaCha8Rng::seed_from_u64(config.seed);
    let trial_seeds: Vec<u64> = (0..config.n_trials)
        .map(|_| rand::Rng::random(&mut seed_rng))
        .collect();

    // Parallel execution with rayon
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
