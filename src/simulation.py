"""
Monte Carlo simulation engine.

Runs paired experiments:
  - Treatment: cohabitants with shared + individual stressors
  - Control: same individuals with fully independent stressors (same total rate)

Supports parameter sweeps over stress rate, shared proportion, group size, etc.
"""

import numpy as np
from dataclasses import dataclass
from tqdm import tqdm

from .cycle_model import CycleParams, CycleState, create_population, advance_day, init_cycle
from .stress_model import StressConfig, generate_stress_timeline, generate_independent_stress
from .synchrony import (compute_phases, mean_resultant_length, rayleigh_test,
                         mean_pairwise_distance, onset_synchrony_index,
                         track_synchrony_over_time)


@dataclass
class ExperimentConfig:
    """Configuration for a single experiment."""
    n_individuals: int = 6          # Group size
    n_days: int = 365 * 3           # Simulation duration (3 years)
    n_trials: int = 200             # Monte Carlo repetitions
    seed: int = 42

    # Cycle parameters
    cycle_params: CycleParams | None = None
    heterogeneity: float = 0.15     # Between-person variability

    # Stress parameters
    stress_config: StressConfig | None = None

    # Measurement
    burnin_days: int = 90           # Discard initial transient
    measure_interval: int = 7       # Days between synchrony snapshots


@dataclass
class TrialResult:
    """Results from a single trial."""
    # Final synchrony metrics (post-burnin average)
    treatment_R: float
    control_R: float
    treatment_mean_dist: float
    control_mean_dist: float
    treatment_onset_sync: float
    control_onset_sync: float
    # Rayleigh test at end
    treatment_rayleigh_p: float
    control_rayleigh_p: float
    # Time series (optional, only for detailed runs)
    treatment_timeseries: dict | None = None
    control_timeseries: dict | None = None


@dataclass
class ExperimentResult:
    """Aggregated results across all trials."""
    config: ExperimentConfig
    trials: list[TrialResult]

    @property
    def treatment_R_mean(self) -> float:
        return float(np.mean([t.treatment_R for t in self.trials]))

    @property
    def control_R_mean(self) -> float:
        return float(np.mean([t.control_R for t in self.trials]))

    @property
    def R_difference(self) -> float:
        return self.treatment_R_mean - self.control_R_mean

    @property
    def treatment_onset_sync_mean(self) -> float:
        return float(np.mean([t.treatment_onset_sync for t in self.trials]))

    @property
    def control_onset_sync_mean(self) -> float:
        return float(np.mean([t.control_onset_sync for t in self.trials]))

    @property
    def fraction_significant(self) -> float:
        """Fraction of treatment trials with Rayleigh p < 0.05."""
        return float(np.mean([t.treatment_rayleigh_p < 0.05 for t in self.trials]))

    @property
    def control_fraction_significant(self) -> float:
        """Fraction of control trials with Rayleigh p < 0.05 (false positive rate)."""
        return float(np.mean([t.control_rayleigh_p < 0.05 for t in self.trials]))


def run_trial(config: ExperimentConfig, trial_seed: int,
              save_timeseries: bool = False) -> TrialResult:
    """Run a single paired trial (treatment + control)."""
    rng = np.random.default_rng(trial_seed)
    stress_config = config.stress_config or StressConfig()
    cycle_params = config.cycle_params or CycleParams()

    # Create two identical populations (same parameters, different starting phases)
    pop_treatment = create_population(
        config.n_individuals, rng, cycle_params, config.heterogeneity)
    # Clone parameters for control, re-randomize starting phases
    rng_ctrl = np.random.default_rng(trial_seed + 1_000_000)
    pop_control = [(p, init_cycle(p, rng_ctrl)) for p, _ in pop_treatment]

    # Generate stress timelines
    stress_treatment = generate_stress_timeline(
        config.n_individuals, config.n_days, stress_config, rng)
    stress_control = generate_independent_stress(
        config.n_individuals, config.n_days, stress_config,
        np.random.default_rng(trial_seed + 2_000_000))

    # Run simulations
    rng_t = np.random.default_rng(trial_seed + 3_000_000)
    rng_c = np.random.default_rng(trial_seed + 4_000_000)

    if save_timeseries:
        ts_treatment = track_synchrony_over_time(
            pop_treatment, config.n_days, stress_treatment, rng_t,
            config.measure_interval)
        ts_control = track_synchrony_over_time(
            pop_control, config.n_days, stress_control, rng_c,
            config.measure_interval)
    else:
        # Just run forward without tracking every interval
        for day in range(config.n_days):
            for i, (params, state) in enumerate(pop_treatment):
                onset = advance_day(state, params, stress_treatment[i, day], rng_t)
                if onset:
                    state.onset_days.append(day)
            for i, (params, state) in enumerate(pop_control):
                onset = advance_day(state, params, stress_control[i, day], rng_c)
                if onset:
                    state.onset_days.append(day)
        ts_treatment = None
        ts_control = None

    # Compute final synchrony metrics (post-burnin)
    burnin = config.burnin_days

    # Phase-based metrics at simulation end
    phases_t = compute_phases(pop_treatment, config.n_days)
    phases_c = compute_phases(pop_control, config.n_days)
    R_t = mean_resultant_length(phases_t)
    R_c = mean_resultant_length(phases_c)
    md_t = mean_pairwise_distance(phases_t)
    md_c = mean_pairwise_distance(phases_c)
    Z_t, p_t = rayleigh_test(phases_t)
    Z_c, p_c = rayleigh_test(phases_c)

    # Onset-based metrics (post-burnin)
    onsets_t = [
        [d for d in state.onset_days if d >= burnin]
        for _, state in pop_treatment
    ]
    onsets_c = [
        [d for d in state.onset_days if d >= burnin]
        for _, state in pop_control
    ]
    osi_t = onset_synchrony_index(onsets_t, window=5.0)
    osi_c = onset_synchrony_index(onsets_c, window=5.0)

    return TrialResult(
        treatment_R=R_t, control_R=R_c,
        treatment_mean_dist=md_t, control_mean_dist=md_c,
        treatment_onset_sync=osi_t, control_onset_sync=osi_c,
        treatment_rayleigh_p=p_t, control_rayleigh_p=p_c,
        treatment_timeseries=ts_treatment,
        control_timeseries=ts_control,
    )


def run_experiment(config: ExperimentConfig,
                   save_timeseries_for: int = 0,
                   show_progress: bool = True) -> ExperimentResult:
    """
    Run a full Monte Carlo experiment.

    Args:
        config: Experiment configuration
        save_timeseries_for: Number of trials to save full timeseries for
        show_progress: Show tqdm progress bar
    """
    base_rng = np.random.default_rng(config.seed)
    trial_seeds = base_rng.integers(0, 2**31, size=config.n_trials)

    trials = []
    iterator = tqdm(range(config.n_trials), desc="Trials",
                    disable=not show_progress)

    for i in iterator:
        save_ts = i < save_timeseries_for
        result = run_trial(config, int(trial_seeds[i]), save_timeseries=save_ts)
        trials.append(result)

    return ExperimentResult(config=config, trials=trials)


def parameter_sweep(base_config: ExperimentConfig,
                    param_name: str,
                    param_values: list,
                    show_progress: bool = True) -> list[ExperimentResult]:
    """
    Sweep over a parameter, running full experiments for each value.

    Args:
        base_config: Base configuration to modify
        param_name: Dot-separated parameter path (e.g., "stress_config.shared_event_rate")
        param_values: List of values to sweep

    Returns:
        List of ExperimentResults, one per parameter value
    """
    results = []
    for val in (tqdm(param_values, desc=f"Sweep {param_name}")
                if show_progress else param_values):
        # Deep-ish copy of config
        import copy
        cfg = copy.deepcopy(base_config)

        # Set nested parameter
        parts = param_name.split(".")
        obj = cfg
        for part in parts[:-1]:
            child = getattr(obj, part)
            if child is None:
                # Initialize default
                if part == "stress_config":
                    setattr(obj, part, StressConfig())
                elif part == "cycle_params":
                    setattr(obj, part, CycleParams())
                child = getattr(obj, part)
            obj = child
        setattr(obj, parts[-1], val)

        result = run_experiment(cfg, save_timeseries_for=1,
                               show_progress=False)
        results.append(result)

    return results
