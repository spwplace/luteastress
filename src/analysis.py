"""
Statistical analysis and visualization.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
from scipy import stats as scipy_stats

from .simulation import ExperimentResult, TrialResult


def plot_synchrony_timeseries(result: ExperimentResult, trial_idx: int = 0,
                               save_path: Path | None = None):
    """Plot synchrony metrics over time for treatment vs control."""
    trial = result.trials[trial_idx]
    if trial.treatment_timeseries is None:
        raise ValueError("No timeseries data saved for this trial")

    ts_t = trial.treatment_timeseries
    ts_c = trial.control_timeseries

    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

    burnin = result.config.burnin_days

    # Mean resultant length R
    ax = axes[0]
    ax.plot(ts_t["days"], ts_t["R"], alpha=0.7, label="Shared stress (treatment)")
    ax.plot(ts_c["days"], ts_c["R"], alpha=0.7, label="Independent stress (control)")
    ax.axvline(burnin, color="gray", linestyle="--", alpha=0.5, label="Burn-in")
    ax.set_ylabel("Mean resultant length R")
    ax.set_ylim(0, 1)
    ax.legend(loc="upper right")
    ax.set_title("Cycle synchrony over time")

    # Mean pairwise distance
    ax = axes[1]
    ax.plot(ts_t["days"], ts_t["mean_distance"], alpha=0.7, label="Treatment")
    ax.plot(ts_c["days"], ts_c["mean_distance"], alpha=0.7, label="Control")
    ax.axhline(np.pi/2, color="gray", linestyle=":", alpha=0.5, label="Expected (uniform)")
    ax.axvline(burnin, color="gray", linestyle="--", alpha=0.5)
    ax.set_ylabel("Mean pairwise phase distance")
    ax.legend(loc="upper right")

    # Rayleigh p-value
    ax = axes[2]
    ax.semilogy(ts_t["days"], ts_t["rayleigh_p"], alpha=0.7, label="Treatment")
    ax.semilogy(ts_c["days"], ts_c["rayleigh_p"], alpha=0.7, label="Control")
    ax.axhline(0.05, color="red", linestyle="--", alpha=0.5, label="p = 0.05")
    ax.axvline(burnin, color="gray", linestyle="--", alpha=0.5)
    ax.set_ylabel("Rayleigh test p-value")
    ax.set_xlabel("Day")
    ax.legend(loc="upper right")

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()
    return fig


def plot_trial_distributions(result: ExperimentResult,
                              save_path: Path | None = None):
    """Plot distributions of synchrony metrics across all trials."""
    R_t = [t.treatment_R for t in result.trials]
    R_c = [t.control_R for t in result.trials]
    osi_t = [t.treatment_onset_sync for t in result.trials]
    osi_c = [t.control_onset_sync for t in result.trials]

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # R distribution
    ax = axes[0]
    bins = np.linspace(0, 1, 30)
    ax.hist(R_t, bins=bins, alpha=0.6, label=f"Shared stress (mean={np.mean(R_t):.3f})")
    ax.hist(R_c, bins=bins, alpha=0.6, label=f"Independent (mean={np.mean(R_c):.3f})")
    ax.set_xlabel("Mean resultant length R")
    ax.set_ylabel("Count")
    ax.set_title("Phase synchrony distribution")
    ax.legend()

    # Onset synchrony distribution
    ax = axes[1]
    bins = np.linspace(0, 1, 30)
    ax.hist(osi_t, bins=bins, alpha=0.6,
            label=f"Shared stress (mean={np.mean(osi_t):.3f})")
    ax.hist(osi_c, bins=bins, alpha=0.6,
            label=f"Independent (mean={np.mean(osi_c):.3f})")
    ax.set_xlabel("Onset synchrony index (±5 day window)")
    ax.set_ylabel("Count")
    ax.set_title("Onset proximity distribution")
    ax.legend()

    plt.suptitle(f"N={result.config.n_individuals} individuals, "
                 f"{result.config.n_days} days, {result.config.n_trials} trials",
                 fontsize=12)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()
    return fig


def plot_parameter_sweep(results: list[ExperimentResult],
                          param_name: str, param_values: list,
                          save_path: Path | None = None):
    """Plot synchrony metrics as a function of a swept parameter."""
    R_means_t = [r.treatment_R_mean for r in results]
    R_means_c = [r.control_R_mean for r in results]
    R_sds_t = [np.std([t.treatment_R for t in r.trials]) for r in results]
    R_sds_c = [np.std([t.control_R for t in r.trials]) for r in results]

    osi_means_t = [r.treatment_onset_sync_mean for r in results]
    osi_means_c = [r.control_onset_sync_mean for r in results]

    frac_sig_t = [r.fraction_significant for r in results]
    frac_sig_c = [r.control_fraction_significant for r in results]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # R vs parameter
    ax = axes[0]
    ax.errorbar(param_values, R_means_t, yerr=R_sds_t, marker="o",
                capsize=3, label="Shared stress")
    ax.errorbar(param_values, R_means_c, yerr=R_sds_c, marker="s",
                capsize=3, label="Independent")
    ax.set_xlabel(param_name)
    ax.set_ylabel("Mean R")
    ax.set_title("Phase synchrony")
    ax.legend()

    # Onset synchrony vs parameter
    ax = axes[1]
    ax.plot(param_values, osi_means_t, "o-", label="Shared stress")
    ax.plot(param_values, osi_means_c, "s-", label="Independent")
    ax.set_xlabel(param_name)
    ax.set_ylabel("Onset synchrony index")
    ax.set_title("Onset proximity")
    ax.legend()

    # Fraction significant vs parameter
    ax = axes[2]
    ax.plot(param_values, frac_sig_t, "o-", label="Shared stress")
    ax.plot(param_values, frac_sig_c, "s-", label="Independent")
    ax.axhline(0.05, color="red", linestyle="--", alpha=0.5, label="Type I rate")
    ax.set_xlabel(param_name)
    ax.set_ylabel("Fraction Rayleigh p < 0.05")
    ax.set_title("Detection rate")
    ax.legend()

    plt.suptitle(f"Parameter sweep: {param_name}", fontsize=12)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()
    return fig


def summary_statistics(result: ExperimentResult) -> dict:
    """Compute summary statistics for an experiment."""
    R_t = np.array([t.treatment_R for t in result.trials])
    R_c = np.array([t.control_R for t in result.trials])
    osi_t = np.array([t.treatment_onset_sync for t in result.trials])
    osi_c = np.array([t.control_onset_sync for t in result.trials])

    # Paired test: treatment vs control R
    t_stat_R, p_R = scipy_stats.wilcoxon(R_t, R_c, alternative="greater")
    t_stat_osi, p_osi = scipy_stats.wilcoxon(osi_t, osi_c, alternative="greater")

    # Effect size (Cohen's d for paired samples)
    diff_R = R_t - R_c
    d_R = np.mean(diff_R) / np.std(diff_R) if np.std(diff_R) > 0 else 0
    diff_osi = osi_t - osi_c
    d_osi = np.mean(diff_osi) / np.std(diff_osi) if np.std(diff_osi) > 0 else 0

    return {
        "n_trials": len(result.trials),
        "n_individuals": result.config.n_individuals,
        "n_days": result.config.n_days,
        "treatment_R_mean": float(np.mean(R_t)),
        "treatment_R_sd": float(np.std(R_t)),
        "control_R_mean": float(np.mean(R_c)),
        "control_R_sd": float(np.std(R_c)),
        "R_difference": float(np.mean(diff_R)),
        "R_wilcoxon_p": float(p_R),
        "R_cohens_d": float(d_R),
        "treatment_onset_sync_mean": float(np.mean(osi_t)),
        "control_onset_sync_mean": float(np.mean(osi_c)),
        "onset_sync_difference": float(np.mean(diff_osi)),
        "onset_sync_wilcoxon_p": float(p_osi),
        "onset_sync_cohens_d": float(d_osi),
        "treatment_fraction_significant": result.fraction_significant,
        "control_fraction_significant": result.control_fraction_significant,
    }


def print_summary(stats: dict):
    """Print a formatted summary of experiment results."""
    print("=" * 60)
    print("EXPERIMENT SUMMARY")
    print(f"  {stats['n_trials']} trials, {stats['n_individuals']} individuals, "
          f"{stats['n_days']} days")
    print("=" * 60)
    print()
    print("Phase Synchrony (Mean Resultant Length R):")
    print(f"  Treatment (shared stress):   {stats['treatment_R_mean']:.4f} "
          f"(SD {stats['treatment_R_sd']:.4f})")
    print(f"  Control (independent stress): {stats['control_R_mean']:.4f} "
          f"(SD {stats['control_R_sd']:.4f})")
    print(f"  Difference:                   {stats['R_difference']:.4f}")
    print(f"  Wilcoxon signed-rank p:       {stats['R_wilcoxon_p']:.2e}")
    print(f"  Cohen's d:                    {stats['R_cohens_d']:.3f}")
    print()
    print("Onset Synchrony (±5 day window):")
    print(f"  Treatment: {stats['treatment_onset_sync_mean']:.4f}")
    print(f"  Control:   {stats['control_onset_sync_mean']:.4f}")
    print(f"  Difference: {stats['onset_sync_difference']:.4f}")
    print(f"  Wilcoxon p: {stats['onset_sync_wilcoxon_p']:.2e}")
    print(f"  Cohen's d:  {stats['onset_sync_cohens_d']:.3f}")
    print()
    print("Rayleigh Test Detection Rate (p < 0.05):")
    print(f"  Treatment: {stats['treatment_fraction_significant']:.1%}")
    print(f"  Control:   {stats['control_fraction_significant']:.1%}")
    print("=" * 60)
