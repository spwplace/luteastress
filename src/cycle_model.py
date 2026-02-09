"""
Individual menstrual cycle model as a stochastic oscillator with phase response curve.

Based on:
- Apple Women's Health Study (2023): mean 28.7d, SD 6.1d, right-skewed
- Tremin database (Treloar 1967): 16-54d range (5th-95th percentile)
- Derry & Derry (2010): chaotic dynamics, Dc ≈ 5.2
- Xiao et al. (1998/1999): follicular phase stress → 4-20+ day delay;
  luteal phase stress → no timing change

The cycle is modeled as two phases:
  - Follicular phase: variable, stress-sensitive (mean ~14d, SD ~3d)
  - Luteal phase: relatively fixed (mean ~14d, SD ~1.5d)

Stress perturbation is applied via a phase response curve (PRC):
  - Maximal sensitivity during mid-follicular phase
  - Near-zero sensitivity during luteal phase
  - Smooth transition between phases
"""

import numpy as np
from dataclasses import dataclass, field


@dataclass
class CycleParams:
    """Parameters for an individual's cycle model."""
    # Mean phase durations (days)
    follicular_mean: float = 14.0
    luteal_mean: float = 14.0

    # Within-person variability (SD in days)
    follicular_sd: float = 3.0
    luteal_sd: float = 1.5

    # Stress sensitivity (dimensionless scale factor for PRC)
    # Higher = more sensitive to stress-induced delays
    stress_sensitivity: float = 1.0

    # Autocorrelation of cycle-to-cycle variability
    # Models the chaotic attractor structure from Derry & Derry 2010
    # 0 = iid noise, 1 = fully correlated (random walk)
    autocorrelation: float = 0.3


@dataclass
class CycleState:
    """Mutable state for tracking an individual's cycle progression."""
    day_in_cycle: float = 0.0
    phase: str = "follicular"  # "follicular" or "luteal"
    current_follicular_length: float = 14.0
    current_luteal_length: float = 14.0
    cycle_count: int = 0
    # Track onset days for synchrony measurement
    onset_days: list = field(default_factory=list)
    # Running deviation for autocorrelated noise
    _prev_follicular_deviation: float = 0.0
    _prev_luteal_deviation: float = 0.0
    # Accumulated stress delay in current follicular phase
    stress_delay: float = 0.0


def phase_response_curve(day_in_cycle: float, follicular_length: float,
                         total_cycle_length: float) -> float:
    """
    Phase response curve: how much a unit stress delays the cycle.

    Based on Xiao et al. (1998): stress during follicular phase causes
    dramatic delays (follicular phase 10→31 days in Group 1), while
    stress during luteal phase has no timing effect.

    Returns a value in [0, 1]:
      - ~1.0 during early-mid follicular phase (peak sensitivity)
      - Drops to ~0 by ovulation and through luteal phase

    The PRC shape is a smoothed step function (sigmoid-based).
    """
    if total_cycle_length <= 0:
        return 0.0

    # Normalize position: 0 = cycle start, 1 = cycle end
    phase_fraction = day_in_cycle / total_cycle_length
    ovulation_fraction = follicular_length / total_cycle_length

    # Sigmoid centered at ovulation, falling off sharply
    # Width parameter controls transition sharpness
    width = 0.08
    sensitivity = 1.0 / (1.0 + np.exp((phase_fraction - ovulation_fraction) / width))

    # Reduce sensitivity right at cycle start (menses itself is less affected)
    # Ramp up over first ~3 days
    menses_ramp = np.clip(day_in_cycle / 3.0, 0.0, 1.0)

    return float(sensitivity * menses_ramp)


def generate_phase_length(mean: float, sd: float, prev_deviation: float,
                          autocorrelation: float, rng: np.random.Generator,
                          min_len: float = 7.0, max_len: float = 40.0) -> tuple[float, float]:
    """
    Generate a phase length with autocorrelated noise.

    Uses AR(1) process for cycle-to-cycle correlation, then clips to
    biologically plausible range. Returns (length, deviation).
    """
    # AR(1) innovation
    innovation = rng.normal(0, sd * np.sqrt(1 - autocorrelation**2))
    deviation = autocorrelation * prev_deviation + innovation
    length = np.clip(mean + deviation, min_len, max_len)

    return float(length), float(deviation)


def init_cycle(params: CycleParams, rng: np.random.Generator) -> CycleState:
    """Initialize a new individual's cycle state at a random phase."""
    state = CycleState()

    # Generate first cycle lengths
    foll, dev_f = generate_phase_length(
        params.follicular_mean, params.follicular_sd, 0.0,
        params.autocorrelation, rng, min_len=7.0, max_len=40.0)
    lut, dev_l = generate_phase_length(
        params.luteal_mean, params.luteal_sd, 0.0,
        params.autocorrelation, rng, min_len=9.0, max_len=18.0)

    state.current_follicular_length = foll
    state.current_luteal_length = lut
    state._prev_follicular_deviation = dev_f
    state._prev_luteal_deviation = dev_l

    # Start at random point in cycle (uniform phase)
    total = foll + lut
    state.day_in_cycle = rng.uniform(0, total)
    if state.day_in_cycle < foll:
        state.phase = "follicular"
    else:
        state.phase = "luteal"

    return state


def advance_day(state: CycleState, params: CycleParams,
                stress_magnitude: float, rng: np.random.Generator) -> bool:
    """
    Advance the cycle by one day, applying any stress perturbation.

    Args:
        state: Current cycle state (mutated in place)
        params: Individual's cycle parameters
        stress_magnitude: Stress intensity for this day [0, ∞)
        rng: Random number generator

    Returns:
        True if a new cycle onset occurred this day
    """
    new_onset = False
    total_length = state.current_follicular_length + state.current_luteal_length

    # Apply stress via PRC
    if stress_magnitude > 0 and state.phase == "follicular":
        prc = phase_response_curve(
            state.day_in_cycle,
            state.current_follicular_length,
            total_length
        )
        # Delay is PRC * sensitivity * stress magnitude
        # Scaled so stress_magnitude=1 with sensitivity=1 gives ~2-3 day delay per event
        delay = prc * params.stress_sensitivity * stress_magnitude * 2.5
        state.stress_delay += delay
        # Extend the current follicular phase
        state.current_follicular_length += delay
        # Cap at biological maximum (extreme stress → amenorrhea-like)
        state.current_follicular_length = min(state.current_follicular_length, 120)
        total_length = state.current_follicular_length + state.current_luteal_length

    # Advance one day
    state.day_in_cycle += 1.0

    # Check for phase transition: follicular → luteal
    if state.phase == "follicular" and state.day_in_cycle >= state.current_follicular_length:
        state.phase = "luteal"

    # Check for cycle completion: luteal → new follicular
    if state.day_in_cycle >= total_length:
        state.day_in_cycle = 0.0
        state.phase = "follicular"
        state.cycle_count += 1
        state.stress_delay = 0.0
        new_onset = True

        # Generate next cycle's phase lengths
        foll, dev_f = generate_phase_length(
            params.follicular_mean, params.follicular_sd,
            state._prev_follicular_deviation,
            params.autocorrelation, rng,
            min_len=7.0, max_len=40.0)
        lut, dev_l = generate_phase_length(
            params.luteal_mean, params.luteal_sd,
            state._prev_luteal_deviation,
            params.autocorrelation, rng,
            min_len=9.0, max_len=18.0)

        state.current_follicular_length = foll
        state.current_luteal_length = lut
        state._prev_follicular_deviation = dev_f
        state._prev_luteal_deviation = dev_l

    return new_onset


def create_population(n: int, rng: np.random.Generator,
                      params: CycleParams | None = None,
                      heterogeneity: float = 0.15) -> list[tuple[CycleParams, CycleState]]:
    """
    Create a population of individuals with heterogeneous cycle parameters.

    Args:
        n: Number of individuals
        rng: Random number generator
        params: Base parameters (defaults used if None)
        heterogeneity: Coefficient of variation for between-person differences
                       in mean cycle parameters

    Returns:
        List of (params, state) tuples
    """
    base = params or CycleParams()
    population = []

    for _ in range(n):
        # Individual variation in mean phase lengths
        p = CycleParams(
            follicular_mean=max(8, rng.normal(base.follicular_mean,
                                              base.follicular_mean * heterogeneity)),
            luteal_mean=max(9, rng.normal(base.luteal_mean,
                                          base.luteal_mean * heterogeneity * 0.5)),
            follicular_sd=max(0.5, rng.normal(base.follicular_sd,
                                               base.follicular_sd * 0.3)),
            luteal_sd=max(0.3, rng.normal(base.luteal_sd,
                                           base.luteal_sd * 0.3)),
            stress_sensitivity=max(0.1, rng.lognormal(
                np.log(max(1e-10, base.stress_sensitivity)), 0.5)),
            autocorrelation=np.clip(rng.normal(base.autocorrelation, 0.1), 0, 0.8),
        )
        s = init_cycle(p, rng)
        population.append((p, s))

    return population
