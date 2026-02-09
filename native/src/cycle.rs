use rand::Rng;
use rand_distr::{Distribution, LogNormal, Normal};

#[derive(Clone, Debug)]
pub struct CycleParams {
    pub follicular_mean: f64,
    pub luteal_mean: f64,
    pub follicular_sd: f64,
    pub luteal_sd: f64,
    pub stress_sensitivity: f64,
    pub autocorrelation: f64,
}

impl Default for CycleParams {
    fn default() -> Self {
        Self {
            follicular_mean: 14.0,
            luteal_mean: 14.0,
            follicular_sd: 3.0,
            luteal_sd: 1.5,
            stress_sensitivity: 1.0,
            autocorrelation: 0.3,
        }
    }
}

#[derive(Clone, Debug)]
pub struct CycleState {
    pub day_in_cycle: f64,
    pub in_follicular: bool,
    pub current_follicular_length: f64,
    pub current_luteal_length: f64,
    pub cycle_count: u32,
    pub onset_days: Vec<u32>,
    prev_follicular_deviation: f64,
    prev_luteal_deviation: f64,
    pub stress_delay: f64,
}

/// Phase response curve: how much a unit stress delays the cycle.
///
/// Follicular phase: high sensitivity (sigmoid dropping at ovulation)
/// Luteal phase: near-zero sensitivity
pub fn phase_response_curve(day_in_cycle: f64, follicular_length: f64, total_length: f64) -> f64 {
    if total_length <= 0.0 {
        return 0.0;
    }
    let phase_frac = day_in_cycle / total_length;
    let ov_frac = follicular_length / total_length;
    let width = 0.08;
    let sensitivity = 1.0 / (1.0 + ((phase_frac - ov_frac) / width).exp());
    let menses_ramp = (day_in_cycle / 3.0).clamp(0.0, 1.0);
    sensitivity * menses_ramp
}

fn generate_phase_length(
    mean: f64,
    sd: f64,
    prev_deviation: f64,
    autocorrelation: f64,
    rng: &mut impl Rng,
    min_len: f64,
    max_len: f64,
) -> (f64, f64) {
    let innovation_sd = sd * (1.0 - autocorrelation * autocorrelation).sqrt();
    let innovation: f64 = Normal::new(0.0, innovation_sd).unwrap().sample(rng);
    let deviation = autocorrelation * prev_deviation + innovation;
    let length = (mean + deviation).clamp(min_len, max_len);
    (length, deviation)
}

pub fn init_cycle(params: &CycleParams, rng: &mut impl Rng) -> CycleState {
    let (foll, dev_f) = generate_phase_length(
        params.follicular_mean,
        params.follicular_sd,
        0.0,
        params.autocorrelation,
        rng,
        7.0,
        40.0,
    );
    let (lut, dev_l) = generate_phase_length(
        params.luteal_mean,
        params.luteal_sd,
        0.0,
        params.autocorrelation,
        rng,
        9.0,
        18.0,
    );
    let total = foll + lut;
    let day_in_cycle: f64 = rng.random_range(0.0..total);
    let in_follicular = day_in_cycle < foll;

    CycleState {
        day_in_cycle,
        in_follicular,
        current_follicular_length: foll,
        current_luteal_length: lut,
        cycle_count: 0,
        onset_days: Vec::new(),
        prev_follicular_deviation: dev_f,
        prev_luteal_deviation: dev_l,
        stress_delay: 0.0,
    }
}

/// Advance cycle by one day. Returns true if a new onset occurred.
pub fn advance_day(
    state: &mut CycleState,
    params: &CycleParams,
    stress_magnitude: f64,
    rng: &mut impl Rng,
) -> bool {
    let total = state.current_follicular_length + state.current_luteal_length;

    // Apply stress via PRC (only during follicular phase)
    if stress_magnitude > 0.0 && state.in_follicular {
        let prc = phase_response_curve(
            state.day_in_cycle,
            state.current_follicular_length,
            total,
        );
        let delay = prc * params.stress_sensitivity * stress_magnitude * 2.5;
        state.stress_delay += delay;
        state.current_follicular_length = (state.current_follicular_length + delay).min(120.0);
    }

    state.day_in_cycle += 1.0;

    let total = state.current_follicular_length + state.current_luteal_length;

    // Phase transition: follicular → luteal
    if state.in_follicular && state.day_in_cycle >= state.current_follicular_length {
        state.in_follicular = false;
    }

    // Cycle completion
    if state.day_in_cycle >= total {
        state.day_in_cycle = 0.0;
        state.in_follicular = true;
        state.cycle_count += 1;
        state.stress_delay = 0.0;

        let (foll, dev_f) = generate_phase_length(
            params.follicular_mean,
            params.follicular_sd,
            state.prev_follicular_deviation,
            params.autocorrelation,
            rng,
            7.0,
            40.0,
        );
        let (lut, dev_l) = generate_phase_length(
            params.luteal_mean,
            params.luteal_sd,
            state.prev_luteal_deviation,
            params.autocorrelation,
            rng,
            9.0,
            18.0,
        );
        state.current_follicular_length = foll;
        state.current_luteal_length = lut;
        state.prev_follicular_deviation = dev_f;
        state.prev_luteal_deviation = dev_l;
        return true;
    }
    false
}

pub fn create_population(
    n: usize,
    base: &CycleParams,
    heterogeneity: f64,
    rng: &mut impl Rng,
) -> Vec<(CycleParams, CycleState)> {
    let normal = Normal::new(0.0, 1.0).unwrap();
    (0..n)
        .map(|_| {
            let p = CycleParams {
                follicular_mean: (base.follicular_mean
                    + normal.sample(rng) * base.follicular_mean * heterogeneity)
                    .max(8.0),
                luteal_mean: (base.luteal_mean
                    + normal.sample(rng) * base.luteal_mean * heterogeneity * 0.5)
                    .max(9.0),
                follicular_sd: (base.follicular_sd
                    + normal.sample(rng) * base.follicular_sd * 0.3)
                    .max(0.5),
                luteal_sd: (base.luteal_sd + normal.sample(rng) * base.luteal_sd * 0.3).max(0.3),
                stress_sensitivity: LogNormal::new(
                    base.stress_sensitivity.max(1e-10).ln(),
                    0.5,
                )
                .unwrap()
                .sample(rng)
                .max(0.1),
                autocorrelation: (base.autocorrelation + normal.sample(rng) * 0.1).clamp(0.0, 0.8),
            };
            let s = init_cycle(&p, rng);
            (p, s)
        })
        .collect()
}
