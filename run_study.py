#!/usr/bin/env python3
"""
Lutea Stress: Computational study of menstrual cycle alignment via shared stressors.

Tests whether shared environmental stressors among cohabitants can produce
apparent menstrual cycle synchrony (entrainment) via phase-dependent delays.

Usage:
    python run_study.py                 # Run default experiment
    python run_study.py --sweep         # Run parameter sweep
    python run_study.py --quick         # Quick validation run
"""

import argparse
import json
from pathlib import Path
import time

import numpy as np
import luteastress_native

from src.cycle_model import CycleParams
from src.stress_model import StressConfig
from src.simulation import ExperimentConfig, run_experiment, parameter_sweep
from src.analysis import (
    plot_synchrony_timeseries, plot_trial_distributions,
    plot_parameter_sweep, summary_statistics, print_summary,
)


OUTPUT_DIR = Path("output")


def run_default_experiment():
    """Run the main experiment with literature-based parameters."""
    OUTPUT_DIR.mkdir(exist_ok=True)

    print("Running default experiment...")
    print("  Cycle model: follicular (14±3d) + luteal (14±1.5d)")
    print("  Stress: shared (0.05/d) + individual (0.07/d)")
    print("  Phase response: follicular-sensitive, luteal-insensitive")
    print()

    config = ExperimentConfig(
        n_individuals=6,
        n_days=365 * 3,
        n_trials=200,
        seed=42,
        cycle_params=CycleParams(
            follicular_mean=14.0,
            luteal_mean=14.0,
            follicular_sd=3.0,
            luteal_sd=1.5,
            stress_sensitivity=1.0,
            autocorrelation=0.3,
        ),
        stress_config=StressConfig(
            shared_event_rate=0.05,
            individual_event_rate=0.07,
            magnitude_mean_log=0.0,
            magnitude_sd_log=0.7,
            duration_mean=2.0,
            duration_sd=1.0,
            shared_exposure_prob=0.85,
        ),
        burnin_days=90,
        measure_interval=7,
    )

    result = run_experiment(config, save_timeseries_for=3)

    stats = summary_statistics(result)
    print_summary(stats)

    # Save results
    with open(OUTPUT_DIR / "default_summary.json", "w") as f:
        json.dump(stats, f, indent=2)

    # Plots
    plot_synchrony_timeseries(result, trial_idx=0,
                               save_path=OUTPUT_DIR / "synchrony_timeseries.png")
    plot_trial_distributions(result,
                              save_path=OUTPUT_DIR / "trial_distributions.png")
    print(f"\nPlots saved to {OUTPUT_DIR}/")


def run_parameter_sweeps():
    """Sweep key parameters to map the entrainment landscape."""
    OUTPUT_DIR.mkdir(exist_ok=True)

    base = ExperimentConfig(
        n_individuals=6,
        n_days=365 * 2,
        n_trials=100,
        seed=42,
        cycle_params=CycleParams(),
        stress_config=StressConfig(),
        burnin_days=90,
        measure_interval=14,
    )

    # Sweep 1: shared stress rate
    print("Sweep 1: Shared stress event rate...")
    rates = [0.0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3]
    results_rate = parameter_sweep(base, "stress_config.shared_event_rate", rates)
    plot_parameter_sweep(results_rate, "shared_event_rate", rates,
                          save_path=OUTPUT_DIR / "sweep_shared_rate.png")

    # Sweep 2: group size
    print("\nSweep 2: Group size...")
    sizes = [2, 3, 4, 6, 8, 10, 15]
    results_size = parameter_sweep(base, "n_individuals", sizes)
    plot_parameter_sweep(results_size, "n_individuals", sizes,
                          save_path=OUTPUT_DIR / "sweep_group_size.png")

    # Sweep 3: stress sensitivity
    print("\nSweep 3: Stress sensitivity...")
    sensitivities = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]
    results_sens = parameter_sweep(base, "cycle_params.stress_sensitivity", sensitivities)
    plot_parameter_sweep(results_sens, "stress_sensitivity", sensitivities,
                          save_path=OUTPUT_DIR / "sweep_sensitivity.png")

    # Sweep 4: shared exposure probability
    print("\nSweep 4: Shared exposure probability...")
    probs = [0.0, 0.2, 0.4, 0.6, 0.8, 0.95, 1.0]
    results_prob = parameter_sweep(base, "stress_config.shared_exposure_prob", probs)
    plot_parameter_sweep(results_prob, "shared_exposure_prob", probs,
                          save_path=OUTPUT_DIR / "sweep_exposure_prob.png")

    print(f"\nAll sweep plots saved to {OUTPUT_DIR}/")

    # Summary table
    print("\n" + "=" * 70)
    print("PARAMETER SWEEP SUMMARY")
    print("=" * 70)
    for name, values, results in [
        ("shared_event_rate", rates, results_rate),
        ("n_individuals", sizes, results_size),
        ("stress_sensitivity", sensitivities, results_sens),
        ("shared_exposure_prob", probs, results_prob),
    ]:
        print(f"\n{name}:")
        print(f"  {'Value':>8}  {'R_treat':>8}  {'R_ctrl':>8}  {'ΔR':>8}  {'%sig_t':>8}  {'%sig_c':>8}")
        for val, res in zip(values, results):
            print(f"  {val:>8.2f}  {res.treatment_R_mean:>8.4f}  "
                  f"{res.control_R_mean:>8.4f}  "
                  f"{res.R_difference:>8.4f}  "
                  f"{res.fraction_significant:>7.1%}  "
                  f"{res.control_fraction_significant:>7.1%}")


def run_quick():
    """Quick validation run with fewer trials."""
    print("Quick validation run...")
    config = ExperimentConfig(
        n_individuals=4,
        n_days=365,
        n_trials=20,
        seed=42,
        cycle_params=CycleParams(),
        stress_config=StressConfig(),
        burnin_days=60,
        measure_interval=7,
    )
    result = run_experiment(config, save_timeseries_for=1)
    stats = summary_statistics(result)
    print_summary(stats)


def run_large_sweep():
    """Run large-scale 2D parameter sweeps using the Rust backend."""
    import time
    import matplotlib.pyplot as plt
    from luteastress_native import run_sweep, run_sweep_2d

    OUTPUT_DIR.mkdir(exist_ok=True)
    n_trials = 1000
    n_days = 365 * 5

    print(f"Large-scale parameter sweeps (Rust backend)")
    print(f"  {n_trials} trials per point, {n_days // 365} years each")
    print()

    # --- 1D sweeps with high power ---
    print("=" * 60)
    print("1D SWEEPS")
    print("=" * 60)

    sweeps_1d = {
        "shared_event_rate": [0.0, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3],
        "stress_sensitivity": [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0],
        "shared_exposure_prob": [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    }

    for param_name, values in sweeps_1d.items():
        t0 = time.perf_counter()
        results = run_sweep(param_name, values, n_trials=n_trials, n_days=n_days)
        elapsed = time.perf_counter() - t0

        print(f"\n{param_name} ({elapsed:.1f}s):")
        print(f"  {'Value':>8}  {'t_R':>7}  {'c_R':>7}  {'ΔR':>7}  {'t_OSI':>7}  {'c_OSI':>7}  {'ΔOSI':>7}")
        for val, res in zip(values, results):
            tr = np.mean(res["treatment_R"])
            cr = np.mean(res["control_R"])
            to = np.mean(res["treatment_osi"])
            co = np.mean(res["control_osi"])
            print(f"  {val:>8.3f}  {tr:>7.4f}  {cr:>7.4f}  {tr-cr:>+7.4f}  "
                  f"{to:>7.4f}  {co:>7.4f}  {to-co:>+7.4f}")

        # Plot 1D sweep
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        t_r = [np.mean(r["treatment_R"]) for r in results]
        c_r = [np.mean(r["control_R"]) for r in results]
        t_o = [np.mean(r["treatment_osi"]) for r in results]
        c_o = [np.mean(r["control_osi"]) for r in results]

        axes[0].plot(values, t_r, "o-", label="Shared stress")
        axes[0].plot(values, c_r, "s-", label="Independent")
        axes[0].set_xlabel(param_name)
        axes[0].set_ylabel("Mean R")
        axes[0].set_title("Phase synchrony")
        axes[0].legend()

        axes[1].plot(values, t_o, "o-", label="Shared stress")
        axes[1].plot(values, c_o, "s-", label="Independent")
        axes[1].set_xlabel(param_name)
        axes[1].set_ylabel("Onset synchrony index")
        axes[1].set_title("Onset proximity")
        axes[1].legend()

        plt.suptitle(f"{param_name} sweep ({n_trials} trials, {n_days//365}yr)")
        plt.tight_layout()
        fig.savefig(OUTPUT_DIR / f"large_sweep_{param_name}.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    # --- 2D sweeps ---
    print("\n" + "=" * 60)
    print("2D SWEEPS (heatmaps)")
    print("=" * 60)

    grids_2d = [
        {
            "a": ("shared_event_rate", [0.0, 0.02, 0.05, 0.08, 0.1, 0.15, 0.2]),
            "b": ("stress_sensitivity", [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0]),
        },
        {
            "a": ("shared_exposure_prob", [0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.95]),
            "b": ("shared_event_rate", [0.0, 0.02, 0.05, 0.08, 0.1, 0.15, 0.2]),
        },
    ]

    for grid_spec in grids_2d:
        a_name, a_vals = grid_spec["a"]
        b_name, b_vals = grid_spec["b"]

        t0 = time.perf_counter()
        result = run_sweep_2d(
            a_name, a_vals,
            b_name, b_vals,
            n_trials=n_trials,
            n_days=n_days,
        )
        elapsed = time.perf_counter() - t0

        na, nb = result["n_a"], result["n_b"]
        print(f"\n{a_name} x {b_name} ({na}x{nb} = {na*nb} points, {elapsed:.1f}s)")

        # Compute delta OSI heatmap
        t_osi = np.array(result["mean_treatment_osi"]).reshape(na, nb)
        c_osi = np.array(result["mean_control_osi"]).reshape(na, nb)
        delta_osi = t_osi - c_osi

        t_r = np.array(result["mean_treatment_R"]).reshape(na, nb)
        c_r = np.array(result["mean_control_R"]).reshape(na, nb)
        delta_r = t_r - c_r

        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        im0 = axes[0].imshow(delta_r, origin="lower", aspect="auto", cmap="RdYlGn")
        axes[0].set_xticks(range(nb))
        axes[0].set_xticklabels([f"{v:.2f}" for v in b_vals], rotation=45)
        axes[0].set_yticks(range(na))
        axes[0].set_yticklabels([f"{v:.2f}" for v in a_vals])
        axes[0].set_xlabel(b_name)
        axes[0].set_ylabel(a_name)
        axes[0].set_title("ΔR (treatment - control)")
        plt.colorbar(im0, ax=axes[0])

        im1 = axes[1].imshow(delta_osi, origin="lower", aspect="auto", cmap="RdYlGn")
        axes[1].set_xticks(range(nb))
        axes[1].set_xticklabels([f"{v:.2f}" for v in b_vals], rotation=45)
        axes[1].set_yticks(range(na))
        axes[1].set_yticklabels([f"{v:.2f}" for v in a_vals])
        axes[1].set_xlabel(b_name)
        axes[1].set_ylabel(a_name)
        axes[1].set_title("ΔOSI (treatment - control)")
        plt.colorbar(im1, ax=axes[1])

        plt.suptitle(f"2D sweep: {a_name} × {b_name} ({n_trials} trials, {n_days//365}yr)")
        plt.tight_layout()
        fname = f"large_sweep_2d_{a_name}_x_{b_name}.png"
        fig.savefig(OUTPUT_DIR / fname, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved {fname}")
        print(f"  ΔR range: [{delta_r.min():.4f}, {delta_r.max():.4f}]")
        print(f"  ΔOSI range: [{delta_osi.min():.4f}, {delta_osi.max():.4f}]")

    print(f"\nAll plots saved to {OUTPUT_DIR}/")


def run_mechanistic_study():
    """Run the study using the high-fidelity mechanistic model in Rust."""
    import matplotlib.pyplot as plt
    from scipy import stats as sp_stats

    OUTPUT_DIR.mkdir(exist_ok=True)

    # Validation settings
    n_trials = 50
    n_days = 365

    print("Running study with High-Fidelity Mechanistic ODE model (Rust)")
    print("  Chain: Shared Schedule -> Circadian -> HPA (Cortisol) -> GnRH")
    print(f"  {n_trials} trials, {n_days} days, 6 individuals")
    print()

    # --- Main experiment ---
    print("Main experiment...")
    start = time.time()
    results = luteastress_native.run_mechanistic_experiment(
        n_individuals=6,
        n_days=n_days,
        n_trials=n_trials,
        seed=42,
    )
    duration = time.time() - start

    t_osi = np.array(results["treatment_osi"])
    c_osi = np.array(results["control_osi"])

    print(f"  Completed in {duration:.1f}s")
    print("=" * 60)
    print("MECHANISTIC MODEL RESULTS")
    print("=" * 60)
    print(f"Onset Synchrony Index (±5 days):")
    print(f"  Treatment: {np.mean(t_osi):.4f} (SD {np.std(t_osi):.4f})")
    print(f"  Control:   {np.mean(c_osi):.4f} (SD {np.std(c_osi):.4f})")
    print(f"  Difference: {np.mean(t_osi - c_osi):+.4f}")

    if n_trials > 1:
        _, p = sp_stats.wilcoxon(t_osi, c_osi, alternative="greater")
        print(f"  Wilcoxon p: {p:.2e}")
    print("=" * 60)

    # --- 1D sweeps ---
    print("\n" + "=" * 60)
    print("MECHANISTIC 1D SWEEPS")
    print("=" * 60)

    sweeps = {
        "shared_event_rate": [0.0, 0.05, 0.1, 0.2], 
        "shared_exposure_prob": [0.0, 0.5, 0.95],
    }

    if hasattr(luteastress_native, "run_mechanistic_sweep"):
        for param_name, values in sweeps.items():
            t0 = time.perf_counter()
            sweep_results = luteastress_native.run_mechanistic_sweep(
                param_name, values,
                n_trials=n_trials,
                n_days=n_days,
            )
            elapsed = time.perf_counter() - t0

            print(f"\n{param_name} ({elapsed:.1f}s):")
            print(f"  {'Value':>8}  {'t_OSI':>7}  {'c_OSI':>7}  {'ΔOSI':>7}")

            t_means = []
            c_means = []
            for val, res in zip(values, sweep_results):
                to = np.mean(res["treatment_osi"])
                co = np.mean(res["control_osi"])
                t_means.append(to)
                c_means.append(co)
                print(f"  {val:>8.3f}  {to:>7.4f}  {co:>7.4f}  {to - co:>+7.4f}")

            # Plot
            fig, ax = plt.subplots(figsize=(7, 4.5))
            ax.plot(values, t_means, "o-", label="Shared stress")
            ax.plot(values, c_means, "s-", label="Independent")
            ax.set_xlabel(param_name)
            ax.set_ylabel("Onset synchrony index")
            ax.set_title(f"Mechanistic model: {param_name} sweep ({n_trials} trials, {n_days}d)")
            ax.legend()
            plt.tight_layout()
            fname = f"mechanistic_sweep_{param_name}.png"
            fig.savefig(OUTPUT_DIR / fname, dpi=150, bbox_inches="tight")
            plt.close(fig)
            print(f"  Saved {fname}")
    
    print(f"\nAll mechanistic results processed.")


def main():
    parser = argparse.ArgumentParser(
        description="Lutea Stress: Menstrual cycle alignment via shared stressors")
    parser.add_argument("--sweep", action="store_true",
                        help="Run parameter sweeps")
    parser.add_argument("--quick", action="store_true",
                        help="Quick validation run")
    parser.add_argument("--large-sweep", action="store_true",
                        help="Large-scale 2D sweeps using Rust backend")
    parser.add_argument("--mechanistic", action="store_true",
                        help="Run high-fidelity mechanistic model (Rust)")
    args = parser.parse_args()

    if args.mechanistic:
        run_mechanistic_study()
    elif args.quick:
        run_quick()
    elif args.large_sweep:
        run_large_sweep()
    elif args.sweep:
        run_parameter_sweeps()
    else:
        run_default_experiment()


if __name__ == "__main__":
    main()
