"""
Stressor generation model for cohabiting individuals.

Stressors are categorized as:
  - Shared: affect all cohabitants (illness spreading, shared deadlines,
    household conflicts, seasonal events, weather)
  - Individual: affect only one person (personal issues, work stress)

The key hypothesis: cohabitants share a significant proportion of stressors,
and these shared perturbations could drive apparent cycle entrainment.
"""

import numpy as np
from dataclasses import dataclass


@dataclass
class StressConfig:
    """Configuration for stress generation."""
    # Rate of stressor events (per person per day)
    # Based on exam/acute stress studies: significant stressors ~1-3x per month
    shared_event_rate: float = 0.05    # ~1.5 shared events/month
    individual_event_rate: float = 0.07  # ~2 individual events/month

    # Magnitude distribution parameters (lognormal)
    # magnitude=1.0 corresponds to a moderate acute stressor
    magnitude_mean_log: float = 0.0   # median magnitude = 1.0
    magnitude_sd_log: float = 0.7     # right-skewed, occasional large stressors

    # Duration of stress effect (days)
    # Acute stressors: 1-3 days; chronic: modeled as repeated events
    duration_mean: float = 2.0
    duration_sd: float = 1.0

    # Proportion of shared stress that each individual actually experiences
    # Models differential exposure (not everyone is equally affected)
    shared_exposure_prob: float = 0.85


def generate_stress_timeline(n_individuals: int, n_days: int,
                             config: StressConfig,
                             rng: np.random.Generator) -> np.ndarray:
    """
    Generate daily stress magnitudes for a group of cohabitants.

    Returns:
        Array of shape (n_individuals, n_days) with stress magnitudes ≥ 0
    """
    stress = np.zeros((n_individuals, n_days))

    # Generate shared stressor events
    n_shared_events = rng.poisson(config.shared_event_rate * n_days)
    for _ in range(n_shared_events):
        # Random onset day
        onset = rng.integers(0, n_days)
        # Random magnitude
        magnitude = rng.lognormal(config.magnitude_mean_log, config.magnitude_sd_log)
        # Random duration
        duration = max(1, int(rng.normal(config.duration_mean, config.duration_sd)))

        for day_offset in range(duration):
            day = onset + day_offset
            if day >= n_days:
                break
            # Decay the magnitude over duration
            decay = np.exp(-0.5 * day_offset)
            for i in range(n_individuals):
                # Each individual may or may not be exposed
                if rng.random() < config.shared_exposure_prob:
                    # Individual variation in experienced magnitude
                    personal_factor = rng.lognormal(0, 0.3)
                    stress[i, day] += magnitude * decay * personal_factor

    # Generate individual stressor events
    for i in range(n_individuals):
        n_indiv_events = rng.poisson(config.individual_event_rate * n_days)
        for _ in range(n_indiv_events):
            onset = rng.integers(0, n_days)
            magnitude = rng.lognormal(config.magnitude_mean_log, config.magnitude_sd_log)
            duration = max(1, int(rng.normal(config.duration_mean, config.duration_sd)))

            for day_offset in range(duration):
                day = onset + day_offset
                if day >= n_days:
                    break
                decay = np.exp(-0.5 * day_offset)
                personal_factor = rng.lognormal(0, 0.3)
                stress[i, day] += magnitude * decay * personal_factor

    return stress


def generate_independent_stress(n_individuals: int, n_days: int,
                                config: StressConfig,
                                rng: np.random.Generator) -> np.ndarray:
    """
    Generate fully independent stress timelines (control condition).

    Same total stress rate but all events are individual — no shared component.
    """
    control_config = StressConfig(
        shared_event_rate=0.0,
        individual_event_rate=(config.shared_event_rate + config.individual_event_rate),
        magnitude_mean_log=config.magnitude_mean_log,
        magnitude_sd_log=config.magnitude_sd_log,
        duration_mean=config.duration_mean,
        duration_sd=config.duration_sd,
        shared_exposure_prob=0.0,
    )
    return generate_stress_timeline(n_individuals, n_days, control_config, rng)
