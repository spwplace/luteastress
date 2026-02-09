mod cycle;
mod mechanistic;
mod mechanistic_params;
mod simulation;
mod stress;
mod synchrony;

use pyo3::prelude::*;
use pyo3::types::PyDict;

use cycle::CycleParams;
use simulation::ExperimentConfig;
use stress::StressConfig;

/// Run a single experiment, returning trial-level arrays as a dict.
#[pyfunction]
#[pyo3(signature = (
    n_individuals = 6,
    n_days = 1095,
    n_trials = 200,
    seed = 42,
    burnin_days = 90,
    heterogeneity = 0.15,
    follicular_mean = 14.0,
    luteal_mean = 14.0,
    follicular_sd = 3.0,
    luteal_sd = 1.5,
    stress_sensitivity = 1.0,
    autocorrelation = 0.3,
    shared_event_rate = 0.05,
    individual_event_rate = 0.07,
    magnitude_mean_log = 0.0,
    magnitude_sd_log = 0.7,
    duration_mean = 2.0,
    duration_sd = 1.0,
    shared_exposure_prob = 0.85,
))]
#[allow(clippy::too_many_arguments)]
fn run_experiment(
    py: Python<'_>,
    n_individuals: usize,
    n_days: usize,
    n_trials: usize,
    seed: u64,
    burnin_days: usize,
    heterogeneity: f64,
    follicular_mean: f64,
    luteal_mean: f64,
    follicular_sd: f64,
    luteal_sd: f64,
    stress_sensitivity: f64,
    autocorrelation: f64,
    shared_event_rate: f64,
    individual_event_rate: f64,
    magnitude_mean_log: f64,
    magnitude_sd_log: f64,
    duration_mean: f64,
    duration_sd: f64,
    shared_exposure_prob: f64,
) -> PyResult<Py<PyDict>> {
    let config = ExperimentConfig {
        n_individuals,
        n_days,
        n_trials,
        seed,
        burnin_days,
        heterogeneity,
        cycle_params: CycleParams {
            follicular_mean,
            luteal_mean,
            follicular_sd,
            luteal_sd,
            stress_sensitivity,
            autocorrelation,
        },
        stress_config: StressConfig {
            shared_event_rate,
            individual_event_rate,
            magnitude_mean_log,
            magnitude_sd_log,
            duration_mean,
            duration_sd,
            shared_exposure_prob,
        },
    };

    let result = py.detach(|| simulation::run_experiment(&config));

    let dict = PyDict::new(py);
    let treatment_r: Vec<f64> = result.trials.iter().map(|t| t.treatment_r).collect();
    let control_r: Vec<f64> = result.trials.iter().map(|t| t.control_r).collect();
    let treatment_osi: Vec<f64> = result.trials.iter().map(|t| t.treatment_osi).collect();
    let control_osi: Vec<f64> = result.trials.iter().map(|t| t.control_osi).collect();
    let treatment_rayleigh_p: Vec<f64> =
        result.trials.iter().map(|t| t.treatment_rayleigh_p).collect();
    let control_rayleigh_p: Vec<f64> =
        result.trials.iter().map(|t| t.control_rayleigh_p).collect();

    dict.set_item("treatment_R", treatment_r)?;
    dict.set_item("control_R", control_r)?;
    dict.set_item("treatment_osi", treatment_osi)?;
    dict.set_item("control_osi", control_osi)?;
    dict.set_item("treatment_rayleigh_p", treatment_rayleigh_p)?;
    dict.set_item("control_rayleigh_p", control_rayleigh_p)?;

    Ok(dict.into())
}

/// Run a parameter sweep in parallel, returning a list of experiment results.
#[pyfunction]
#[pyo3(signature = (
    param_name,
    param_values,
    n_individuals = 6,
    n_days = 730,
    n_trials = 500,
    seed = 42,
    burnin_days = 90,
    heterogeneity = 0.15,
    follicular_mean = 14.0,
    luteal_mean = 14.0,
    follicular_sd = 3.0,
    luteal_sd = 1.5,
    stress_sensitivity = 1.0,
    autocorrelation = 0.3,
    shared_event_rate = 0.05,
    individual_event_rate = 0.07,
    magnitude_mean_log = 0.0,
    magnitude_sd_log = 0.7,
    duration_mean = 2.0,
    duration_sd = 1.0,
    shared_exposure_prob = 0.85,
))]
#[allow(clippy::too_many_arguments)]
fn run_sweep(
    py: Python<'_>,
    param_name: &str,
    param_values: Vec<f64>,
    n_individuals: usize,
    n_days: usize,
    n_trials: usize,
    seed: u64,
    burnin_days: usize,
    heterogeneity: f64,
    follicular_mean: f64,
    luteal_mean: f64,
    follicular_sd: f64,
    luteal_sd: f64,
    stress_sensitivity: f64,
    autocorrelation: f64,
    shared_event_rate: f64,
    individual_event_rate: f64,
    magnitude_mean_log: f64,
    magnitude_sd_log: f64,
    duration_mean: f64,
    duration_sd: f64,
    shared_exposure_prob: f64,
) -> PyResult<Vec<Py<PyDict>>> {
    let base = ExperimentConfig {
        n_individuals,
        n_days,
        n_trials,
        seed,
        burnin_days,
        heterogeneity,
        cycle_params: CycleParams {
            follicular_mean,
            luteal_mean,
            follicular_sd,
            luteal_sd,
            stress_sensitivity,
            autocorrelation,
        },
        stress_config: StressConfig {
            shared_event_rate,
            individual_event_rate,
            magnitude_mean_log,
            magnitude_sd_log,
            duration_mean,
            duration_sd,
            shared_exposure_prob,
        },
    };

    // Build configs for each parameter value
    let configs: Vec<ExperimentConfig> = param_values
        .iter()
        .map(|&val| {
            let mut cfg = base.clone();
            match param_name {
                "shared_event_rate" => cfg.stress_config.shared_event_rate = val,
                "individual_event_rate" => cfg.stress_config.individual_event_rate = val,
                "shared_exposure_prob" => cfg.stress_config.shared_exposure_prob = val,
                "stress_sensitivity" => cfg.cycle_params.stress_sensitivity = val,
                "follicular_mean" => cfg.cycle_params.follicular_mean = val,
                "luteal_mean" => cfg.cycle_params.luteal_mean = val,
                "follicular_sd" => cfg.cycle_params.follicular_sd = val,
                "luteal_sd" => cfg.cycle_params.luteal_sd = val,
                "autocorrelation" => cfg.cycle_params.autocorrelation = val,
                "heterogeneity" => cfg.heterogeneity = val,
                "n_individuals" => cfg.n_individuals = val as usize,
                "magnitude_sd_log" => cfg.stress_config.magnitude_sd_log = val,
                "duration_mean" => cfg.stress_config.duration_mean = val,
                _ => {}
            }
            cfg
        })
        .collect();

    // Run all experiments (each experiment is already parallel internally via rayon)
    let results: Vec<_> = py.detach(|| {
        configs
            .iter()
            .map(|cfg| simulation::run_experiment(cfg))
            .collect::<Vec<_>>()
    });

    // Convert to Python dicts
    results
        .into_iter()
        .map(|result| {
            let dict = PyDict::new(py);
            let tr: Vec<f64> = result.trials.iter().map(|t| t.treatment_r).collect();
            let cr: Vec<f64> = result.trials.iter().map(|t| t.control_r).collect();
            let to: Vec<f64> = result.trials.iter().map(|t| t.treatment_osi).collect();
            let co: Vec<f64> = result.trials.iter().map(|t| t.control_osi).collect();
            let tp: Vec<f64> = result.trials.iter().map(|t| t.treatment_rayleigh_p).collect();
            let cp: Vec<f64> = result.trials.iter().map(|t| t.control_rayleigh_p).collect();
            dict.set_item("treatment_R", tr)?;
            dict.set_item("control_R", cr)?;
            dict.set_item("treatment_osi", to)?;
            dict.set_item("control_osi", co)?;
            dict.set_item("treatment_rayleigh_p", tp)?;
            dict.set_item("control_rayleigh_p", cp)?;
            Ok(dict.into())
        })
        .collect()
}

/// Run a 2D parameter sweep — grid of (param_a × param_b), returning flattened results.
#[pyfunction]
#[pyo3(signature = (
    param_a_name, param_a_values,
    param_b_name, param_b_values,
    n_individuals = 6,
    n_days = 730,
    n_trials = 500,
    seed = 42,
    burnin_days = 90,
    heterogeneity = 0.15,
    follicular_mean = 14.0,
    luteal_mean = 14.0,
    follicular_sd = 3.0,
    luteal_sd = 1.5,
    stress_sensitivity = 1.0,
    autocorrelation = 0.3,
    shared_event_rate = 0.05,
    individual_event_rate = 0.07,
    magnitude_mean_log = 0.0,
    magnitude_sd_log = 0.7,
    duration_mean = 2.0,
    duration_sd = 1.0,
    shared_exposure_prob = 0.85,
))]
#[allow(clippy::too_many_arguments)]
fn run_sweep_2d(
    py: Python<'_>,
    param_a_name: &str,
    param_a_values: Vec<f64>,
    param_b_name: &str,
    param_b_values: Vec<f64>,
    n_individuals: usize,
    n_days: usize,
    n_trials: usize,
    seed: u64,
    burnin_days: usize,
    heterogeneity: f64,
    follicular_mean: f64,
    luteal_mean: f64,
    follicular_sd: f64,
    luteal_sd: f64,
    stress_sensitivity: f64,
    autocorrelation: f64,
    shared_event_rate: f64,
    individual_event_rate: f64,
    magnitude_mean_log: f64,
    magnitude_sd_log: f64,
    duration_mean: f64,
    duration_sd: f64,
    shared_exposure_prob: f64,
) -> PyResult<Py<PyDict>> {
    let base = ExperimentConfig {
        n_individuals,
        n_days,
        n_trials,
        seed,
        burnin_days,
        heterogeneity,
        cycle_params: CycleParams {
            follicular_mean,
            luteal_mean,
            follicular_sd,
            luteal_sd,
            stress_sensitivity,
            autocorrelation,
        },
        stress_config: StressConfig {
            shared_event_rate,
            individual_event_rate,
            magnitude_mean_log,
            magnitude_sd_log,
            duration_mean,
            duration_sd,
            shared_exposure_prob,
        },
    };

    fn set_param(cfg: &mut ExperimentConfig, name: &str, val: f64) {
        match name {
            "shared_event_rate" => cfg.stress_config.shared_event_rate = val,
            "individual_event_rate" => cfg.stress_config.individual_event_rate = val,
            "shared_exposure_prob" => cfg.stress_config.shared_exposure_prob = val,
            "stress_sensitivity" => cfg.cycle_params.stress_sensitivity = val,
            "follicular_mean" => cfg.cycle_params.follicular_mean = val,
            "luteal_mean" => cfg.cycle_params.luteal_mean = val,
            "follicular_sd" => cfg.cycle_params.follicular_sd = val,
            "heterogeneity" => cfg.heterogeneity = val,
            "n_individuals" => cfg.n_individuals = val as usize,
            "autocorrelation" => cfg.cycle_params.autocorrelation = val,
            _ => {}
        }
    }

    let na = param_a_values.len();
    let nb = param_b_values.len();

    // Build grid of configs
    let configs: Vec<ExperimentConfig> = param_a_values
        .iter()
        .flat_map(|&va| {
            let base = base.clone();
            param_b_values.iter().map(move |&vb| {
                let mut cfg = base.clone();
                set_param(&mut cfg, param_a_name, va);
                set_param(&mut cfg, param_b_name, vb);
                cfg
            })
        })
        .collect();

    // Run all in parallel (rayon parallelizes within each experiment)
    let results: Vec<_> = py.detach(|| {
        configs
            .iter()
            .map(|cfg| simulation::run_experiment(cfg))
            .collect::<Vec<_>>()
    });

    // Flatten to arrays: shape [na * nb] in row-major order
    let dict = PyDict::new(py);
    let mean_tr: Vec<f64> = results
        .iter()
        .map(|r| r.trials.iter().map(|t| t.treatment_r).sum::<f64>() / r.trials.len() as f64)
        .collect();
    let mean_cr: Vec<f64> = results
        .iter()
        .map(|r| r.trials.iter().map(|t| t.control_r).sum::<f64>() / r.trials.len() as f64)
        .collect();
    let mean_to: Vec<f64> = results
        .iter()
        .map(|r| r.trials.iter().map(|t| t.treatment_osi).sum::<f64>() / r.trials.len() as f64)
        .collect();
    let mean_co: Vec<f64> = results
        .iter()
        .map(|r| r.trials.iter().map(|t| t.control_osi).sum::<f64>() / r.trials.len() as f64)
        .collect();

    dict.set_item("mean_treatment_R", mean_tr)?;
    dict.set_item("mean_control_R", mean_cr)?;
    dict.set_item("mean_treatment_osi", mean_to)?;
    dict.set_item("mean_control_osi", mean_co)?;
    dict.set_item("n_a", na)?;
    dict.set_item("n_b", nb)?;
    dict.set_item("param_a_values", param_a_values)?;
    dict.set_item("param_b_values", param_b_values)?;

    Ok(dict.into())
}

#[pymodule]
fn luteastress_native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(run_experiment, m)?)?;
    m.add_function(wrap_pyfunction!(run_sweep, m)?)?;
    m.add_function(wrap_pyfunction!(run_sweep_2d, m)?)?;
    m.add_function(wrap_pyfunction!(run_mechanistic_experiment, m)?)?;
    Ok(())
}

#[pyfunction]
#[pyo3(signature = (
    n_individuals = 6,
    n_days = 365,
    n_trials = 10,
    seed = 42,
    burnin_days = 60,
    shared_event_rate = 0.05,
    individual_event_rate = 0.07,
    magnitude_mean_log = 0.0,
    magnitude_sd_log = 0.7,
    duration_mean = 2.0,
    duration_sd = 1.0,
    shared_exposure_prob = 0.85,
))]
#[allow(clippy::too_many_arguments)]
fn run_mechanistic_experiment(
    py: Python<'_>,
    n_individuals: usize,
    n_days: usize,
    n_trials: usize,
    seed: u64,
    burnin_days: usize,
    shared_event_rate: f64,
    individual_event_rate: f64,
    magnitude_mean_log: f64,
    magnitude_sd_log: f64,
    duration_mean: f64,
    duration_sd: f64,
    shared_exposure_prob: f64,
) -> PyResult<Py<PyDict>> {
    let config = ExperimentConfig {
        n_individuals,
        n_days,
        n_trials,
        seed,
        burnin_days,
        heterogeneity: 0.0, // Not used in mechanistic yet
        cycle_params: CycleParams::default(),
        stress_config: StressConfig {
            shared_event_rate,
            individual_event_rate,
            magnitude_mean_log,
            magnitude_sd_log,
            duration_mean,
            duration_sd,
            shared_exposure_prob,
        },
    };

    let results = py.detach(|| simulation::run_mechanistic_experiment(&config));

    let dict = PyDict::new(py);
    let treatment_osi: Vec<f64> = results.iter().map(|t| t.treatment_osi).collect();
    let control_osi: Vec<f64> = results.iter().map(|t| t.control_osi).collect();

    dict.set_item("treatment_osi", treatment_osi)?;
    dict.set_item("control_osi", control_osi)?;

    Ok(dict.into())
}
