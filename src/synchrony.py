"""
Synchrony measurement using circular statistics.

The menstrual cycle is naturally a circular variable — the phase of each
individual can be represented as an angle on [0, 2π). Circular statistics
then provide the natural framework for measuring synchrony.

Key metrics:
  - Mean resultant length (R): 0 = uniform/no synchrony, 1 = perfect synchrony
  - Rayleigh test: statistical test for non-uniformity (i.e., synchrony)
  - Pairwise phase differences: track convergence over time
"""

import numpy as np
from scipy import stats as scipy_stats


def compute_phases(population, day: int) -> np.ndarray:
    """
    Compute the cycle phase (as angle in [0, 2π)) for each individual.

    Phase = 2π * (day_in_cycle / expected_cycle_length)
    """
    phases = []
    for params, state in population:
        total = state.current_follicular_length + state.current_luteal_length
        if total > 0:
            phase = 2 * np.pi * (state.day_in_cycle / total)
        else:
            phase = 0.0
        phases.append(phase)
    return np.array(phases)


def mean_resultant_length(phases: np.ndarray) -> float:
    """
    Compute the mean resultant length R of circular data.

    R = |mean(exp(i*theta))|

    R ∈ [0, 1]:
      0 → phases uniformly distributed (no synchrony)
      1 → all phases identical (perfect synchrony)
    """
    z = np.exp(1j * phases)
    return float(np.abs(np.mean(z)))


def mean_direction(phases: np.ndarray) -> float:
    """Compute the mean direction of circular data."""
    z = np.exp(1j * phases)
    return float(np.angle(np.mean(z)) % (2 * np.pi))


def rayleigh_test(phases: np.ndarray) -> tuple[float, float]:
    """
    Rayleigh test for non-uniformity of circular data.

    Returns (test_statistic, p_value).
    Under H0 (uniform distribution), the test statistic 2nR²
    follows approximately a chi-squared distribution with 2 df.
    """
    n = len(phases)
    R = mean_resultant_length(phases)
    # Test statistic
    Z = n * R**2
    # P-value approximation (Mardia & Jupp, 2000)
    p = np.exp(-Z) * (1 + (2*Z - Z**2) / (4*n) - (24*Z - 132*Z**2 + 76*Z**3 - 9*Z**4) / (288*n**2))
    p = max(0.0, min(1.0, float(p)))
    return float(Z), p


def pairwise_phase_differences(phases: np.ndarray) -> np.ndarray:
    """
    Compute all pairwise circular phase differences.

    Returns array of |phase_i - phase_j| mapped to [0, π].
    """
    n = len(phases)
    diffs = []
    for i in range(n):
        for j in range(i+1, n):
            diff = abs(phases[i] - phases[j])
            diff = min(diff, 2*np.pi - diff)  # circular distance
            diffs.append(diff)
    return np.array(diffs)


def mean_pairwise_distance(phases: np.ndarray) -> float:
    """
    Mean pairwise circular distance.

    Expected value under uniformity: π/2 ≈ 1.571
    Lower values indicate synchrony.
    """
    diffs = pairwise_phase_differences(phases)
    if len(diffs) == 0:
        return 0.0
    return float(np.mean(diffs))


def onset_synchrony_index(onset_days_list: list[list[float]],
                          window: float = 5.0) -> float:
    """
    Compute synchrony based on proximity of cycle onsets.

    For each pair of individuals, compute the fraction of onsets
    that fall within `window` days of an onset from the other person.

    This is a more intuitive metric: "how often do periods start
    within a week of each other?"
    """
    n = len(onset_days_list)
    if n < 2:
        return 0.0

    total_pairs = 0
    close_pairs = 0

    for i in range(n):
        for j in range(i+1, n):
            onsets_i = np.array(onset_days_list[i])
            onsets_j = np.array(onset_days_list[j])
            if len(onsets_i) == 0 or len(onsets_j) == 0:
                continue

            for oi in onsets_i:
                min_dist = np.min(np.abs(onsets_j - oi))
                total_pairs += 1
                if min_dist <= window:
                    close_pairs += 1

    if total_pairs == 0:
        return 0.0
    return close_pairs / total_pairs


def track_synchrony_over_time(population, n_days: int,
                               stress_matrix: np.ndarray,
                               rng: np.random.Generator,
                               measure_interval: int = 7) -> dict:
    """
    Run simulation and track synchrony metrics over time.

    Args:
        population: List of (CycleParams, CycleState) tuples
        n_days: Simulation duration
        stress_matrix: Shape (n_individuals, n_days)
        rng: RNG for cycle model
        measure_interval: Days between synchrony measurements

    Returns:
        Dictionary with time series of synchrony metrics
    """
    from .cycle_model import advance_day

    n = len(population)
    days = []
    R_values = []
    mean_distances = []
    rayleigh_Zs = []
    rayleigh_ps = []

    for day in range(n_days):
        # Advance all individuals
        for i, (params, state) in enumerate(population):
            stress = stress_matrix[i, day] if day < stress_matrix.shape[1] else 0.0
            onset = advance_day(state, params, stress, rng)
            if onset:
                state.onset_days.append(day)

        # Measure synchrony periodically
        if day % measure_interval == 0:
            phases = compute_phases(population, day)
            R = mean_resultant_length(phases)
            md = mean_pairwise_distance(phases)
            Z, p = rayleigh_test(phases)

            days.append(day)
            R_values.append(R)
            mean_distances.append(md)
            rayleigh_Zs.append(Z)
            rayleigh_ps.append(p)

    return {
        "days": np.array(days),
        "R": np.array(R_values),
        "mean_distance": np.array(mean_distances),
        "rayleigh_Z": np.array(rayleigh_Zs),
        "rayleigh_p": np.array(rayleigh_ps),
    }
