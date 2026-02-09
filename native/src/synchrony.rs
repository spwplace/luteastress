use crate::cycle::{CycleParams, CycleState};

pub fn compute_phases(population: &[(CycleParams, CycleState)]) -> Vec<f64> {
    population
        .iter()
        .map(|(_, state)| {
            let total = state.current_follicular_length + state.current_luteal_length;
            if total > 0.0 {
                2.0 * std::f64::consts::PI * (state.day_in_cycle / total)
            } else {
                0.0
            }
        })
        .collect()
}

pub fn mean_resultant_length(phases: &[f64]) -> f64 {
    let n = phases.len() as f64;
    if n == 0.0 {
        return 0.0;
    }
    let sum_cos: f64 = phases.iter().map(|&p| p.cos()).sum();
    let sum_sin: f64 = phases.iter().map(|&p| p.sin()).sum();
    ((sum_cos / n).powi(2) + (sum_sin / n).powi(2)).sqrt()
}

pub fn rayleigh_test(phases: &[f64]) -> (f64, f64) {
    let n = phases.len() as f64;
    let r = mean_resultant_length(phases);
    let z = n * r * r;
    // P-value approximation (Mardia & Jupp)
    let p = (-z).exp()
        * (1.0 + (2.0 * z - z * z) / (4.0 * n)
            - (24.0 * z - 132.0 * z * z + 76.0 * z.powi(3) - 9.0 * z.powi(4))
                / (288.0 * n * n));
    (z, p.clamp(0.0, 1.0))
}

/// Fraction of onsets within ±window days of any onset from another individual
pub fn onset_synchrony_index(onset_days_list: &[Vec<u32>], window: f64) -> f64 {
    let n = onset_days_list.len();
    if n < 2 {
        return 0.0;
    }

    let mut total = 0u64;
    let mut close = 0u64;

    for i in 0..n {
        for j in (i + 1)..n {
            let oi = &onset_days_list[i];
            let oj = &onset_days_list[j];
            if oi.is_empty() || oj.is_empty() {
                continue;
            }
            for &day_i in oi {
                let min_dist = oj
                    .iter()
                    .map(|&day_j| (day_i as f64 - day_j as f64).abs())
                    .fold(f64::INFINITY, f64::min);
                total += 1;
                if min_dist <= window {
                    close += 1;
                }
            }
        }
    }

    if total == 0 {
        0.0
    } else {
        close as f64 / total as f64
    }
}
