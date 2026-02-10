# Lutea Stress

Computational study of menstrual cycle alignment via shared stressors.

## Overview

**Lutea Stress** is a research project exploring the hypothesis that apparent menstrual cycle synchrony (the "McClintock effect") among cohabitants can be driven by shared environmental stressors. Rather than pheromonal signaling, this model investigates how shared perturbations can lead to entrainment through phase-dependent delays in the menstrual cycle.

The project implements two distinct modeling approaches:
1.  **Cycle Model:** A stochastic oscillator model with a Phase Response Curve (PRC), where the follicular phase is sensitive to stress-induced delays while the luteal phase remains relatively fixed.
2.  **Mechanistic Model:** A high-fidelity ODE-based model (implemented in Rust) that simulates the neuroendocrine network, including the GnRH pulse generator, FSH/LH dynamics, follicle recruitment/growth/atresia, and HPA axis (Cortisol) modulation of GnRH frequency.

## Key Features

- **Hybrid Python/Rust Architecture:** High-performance Rust backend for heavy ODE simulations and large-scale parameter sweeps, integrated seamlessly with Python for analysis.
- **Stochastic Phase Response:** Models individual variability and chaotic dynamics (AR(1) noise) in cycle lengths.
- **Shared Stress Generation:** Simulates both shared (cohabitant) and individual stressors with lognormal magnitude distributions.
- **Comprehensive Analysis:** Tools for calculating phase synchrony ($R$) and Onset Synchrony Index (OSI), with automated plotting for parameter sweeps and distributions.
- **Large-Scale Sweeps:** High-speed 1D and 2D parameter sweeps (via Rust/Rayon) to map the entrainment landscape.

## Project Structure

- `src/`: Python source code.
    - `cycle_model.py`: Stochastic oscillator implementation.
    - `stress_model.py`: Shared and individual stressor generation.
    - `simulation.py`: Experiment orchestration.
    - `analysis.py`: Statistical metrics and visualization.
- `native/`: Rust implementation of the high-fidelity mechanistic model.
    - `src/mechanistic.rs`: ODE system and follicle dynamics.
    - `src/zavala.rs`: Implementation of the Zavala (2020) neuroendocrine network.
    - `src/lib.rs`: PyO3 bridge to Python.
- `run_study.py`: Main entry point for running experiments and sweeps.

## Installation

### Prerequisites

- Python 3.11+
- Rust toolchain (for building the native extension)
- `uv` (recommended) or `pip`

### Setup

1. Install Python dependencies:
   ```bash
   uv sync
   ```
2. Build the Rust extension:
   ```bash
   # Using maturin (installed via uv/pip)
   maturin develop
   ```

## Usage

The `run_study.py` script provides several modes of operation:

```bash
# Run the default experiment (literature-based parameters)
python run_study.py

# Run a quick validation experiment
python run_study.py --quick

# Run parameter sweeps (Python-based)
python run_study.py --sweep

# Run large-scale 2D sweeps (Rust backend)
python run_study.py --large-sweep

# Run the high-fidelity mechanistic ODE study (Rust)
python run_study.py --mechanistic
```

Results and plots are saved to the `output/` directory.

## Key Findings

Based on computational experiments (as of February 2026):

- **Stochastic Phase-Delay (SPD) Model:** Shared stressors significantly increased the Onset Synchrony Index (OSI) by $+0.0263$ ($p < 0.001$), suggesting that shared environmental timing can produce apparent onset proximity in simplified settings.
- **Mechanistic ODE Model:** In the higher-fidelity model linking shared stress to circadian phase and HPA/GnRH signaling, the synchrony signal was much weaker ($+0.0165$ OSI difference) and did not reach statistical significance ($p = 0.257$).
- **Conclusion:** While shared stress can increase apparent onset proximity in simplified models, robust biological entrainment is not currently supported in higher-fidelity neuroendocrine simulations.

## Models

### Phase Response Model
Based on literature (e.g., Xiao et al. 1998), stress during the follicular phase causes significant delays by suppressing GnRH/LH pulsatility, whereas stress during the luteal phase has negligible effect on timing. This asymmetry is the primary driver of potential entrainment in the model.

### Mechanistic ODE Model
This model builds on the work of Zavala et al. (2020), implementing a detailed neuroendocrine network. Stress acts as an input to the HPA axis, which in turn modulates the frequency and mass of GnRH pulses, affecting the downstream hormonal cascade and follicle maturation.

## License

MIT
