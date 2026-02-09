#set document(title: "Can Shared Stress Explain Menstrual Synchrony?", author: "Lutea Stress Project")
#set page(margin: (x: 1in, y: 1in), numbering: "1")
#set text(font: "New Computer Modern", size: 11pt)
#set par(justify: true, leading: 0.65em)
#set heading(numbering: "1.")
#show heading.where(level: 1): it => block(above: 1.4em, below: 0.8em, text(size: 13pt, it))
#show heading.where(level: 2): it => block(above: 1.2em, below: 0.6em, text(size: 11.5pt, it))
#show figure.caption: it => text(size: 9.5pt, it)

// Title block
#align(center)[
  #text(size: 18pt, weight: "bold")[Can Shared Stress Explain \ Menstrual Synchrony?]
  #v(0.6em)
  #text(size: 12pt, fill: luma(80))[A computational study of stressor-driven cycle entrainment]
  #v(1.2em)
  #text(size: 10pt)[Lutea Stress Project #h(1em) · #h(1em) February 2026]
]

#v(1em)

#block(fill: luma(245), inset: 12pt, radius: 4pt, width: 100%)[
  *Summary.* #h(0.3em) We simulate whether shared environmental stressors among cohabitants can produce the menstrual cycle synchrony described in popular culture (the "McClintock effect"). Using a two-phase cycle oscillator with a literature-derived phase response curve, we find that shared stress produces *real but transient synchrony*: cycle onsets cluster together significantly more often (Cohen's _d_ = 0.63, _p_ < 10#super[−14]). Validation using a 54-equation mechanistic ODE model of the HPG axis confirms this effect, showing an even larger 10.7 percentage-point increase in onset proximity. However, the alignment does not persist between onsets, so instantaneous phase concentration is unchanged. This suggests the urban legend reflects a genuine mechanism---repeated transient alignment amplified by confirmation bias---rather than sustained entrainment.
]

= Background

The claim that women living together synchronize their menstrual cycles originates with McClintock (1971), who reported convergence of cycle onsets among college dormitory residents and proposed a pheromonal mechanism. Subsequent decades of research have largely undermined both the empirical finding and the proposed mechanism. Critical reviews by Wilson (1992), Schank (2006), and Yang & Schank (2006) identified systematic methodological errors in synchrony studies, including incorrect baseline calculations and selection bias. No pheromonal agent has been isolated in humans.

Yet the subjective experience of "syncing up" persists widely. We propose a different mechanism: cohabitants share environmental stressors (illness, academic deadlines, emotional events, sleep disruption) that can delay menstrual cycles via hypothalamic-pituitary-adrenal (HPA) axis suppression of GnRH pulsatility. If two people experience the same stressor at similar cycle phases, their next onsets could be pushed closer together---producing real, recurring episodes of synchrony that then dissipate as natural period differences cause the cycles to drift apart again.

This study tests computationally whether such transient stress-driven alignment can account for the widely reported experience.

= Model

== Cycle oscillator

Each individual's cycle is modeled as a two-phase stochastic oscillator:

#align(center)[
  $T_"cycle" = T_"foll" + T_"lut"$
]

Phase durations are drawn from autocorrelated Gaussian processes:

#align(center)[
  $T_"foll" tilde.op "AR"(1): mu = 14 "d", sigma = 3 "d", rho = 0.3$
  #v(0.3em)
  $T_"lut" tilde.op "AR"(1): mu = 14 "d", sigma = 1.5 "d", rho = 0.3$
]

These parameters produce cycles with mean 28 days and population-level variability matching the Apple Women's Health Study (mean 28.7 d, SD 6.1 d, right-skewed; _n_ = 165,668 cycles from 12,608 participants). Between-person heterogeneity is introduced by drawing individual means from distributions with 15% coefficient of variation.

The autocorrelation term ($rho = 0.3$) captures the deterministic structure identified by Derry & Derry (2010), who showed menstrual cycle variability reflects low-dimensional chaos (correlation dimension $D_c approx 5.2$) rather than independent random noise.

== Phase response curve

Stress perturbation is applied via an asymmetric phase response curve (PRC) based on the rhesus macaque endotoxin challenge studies of Xiao et al. (1998, 1999):

- *Follicular phase*: Stress delays ovulation by extending folliculogenesis. In the primate data, a 5-day inflammatory challenge during the early follicular phase extended it from 10 to 31 days in sensitive individuals (Group 1), with one animal experiencing 4-month amenorrhea. A less sensitive group showed extension from 10 to 14 days.
- *Luteal phase*: No effect on cycle timing (progesterone levels are reduced but the phase length is unchanged).

We model this as a sigmoid:

#align(center)[
  $"PRC"(t) = (1 + e^((phi - phi_"ov") slash w))^(-1) dot min(t slash 3, 1)$
]

where $phi$ is the normalized cycle position, $phi_"ov"$ is the ovulation point, and $w = 0.08$ controls transition sharpness. The menses ramp ($min(t/3, 1)$) reduces sensitivity during the first 3 days. A unit stress event with sensitivity 1.0 produces approximately 2.5 days of delay at peak PRC.

== Stress generation

Stressors are generated as Poisson events in two categories:

#table(
  columns: (1fr, 1fr, 1fr),
  align: (left, center, center),
  stroke: 0.5pt + luma(180),
  inset: 6pt,
  table.header[*Parameter*][*Shared*][*Individual*],
  [Event rate (per day)], [0.05], [0.07],
  [Magnitude], [lognormal(0, 0.7)], [lognormal(0, 0.7)],
  [Duration (days)], [2 ± 1], [2 ± 1],
  [Exposure probability], [85%], [100%],
)

Shared events affect all cohabitants with 85% probability each (modeling differential exposure). Magnitude decays exponentially over the event duration. Each individual also experiences a personal magnitude variation (lognormal, $sigma$ = 0.3).

The *control condition* converts all shared events to individual events, preserving the total stress rate while eliminating temporal correlation.

== Synchrony metrics

We use two complementary measures:

+ *Mean resultant length* ($R$): The standard circular statistic for phase concentration, measured at regular intervals. Phases are mapped to $[0, 2pi)$ as $theta = 2pi dot t_"in cycle" slash T_"cycle"$, and $R = |n^(-1) sum e^(i theta_k)| in [0, 1]$. This measures *instantaneous* alignment: are everyone's cycles at similar phases right now? High $R$ sustained over time would indicate true entrainment. Statistical significance is assessed via the Rayleigh test.

+ *Onset synchrony index* (OSI): For each pair of individuals, the fraction of cycle onsets falling within $plus.minus$5 days of an onset from the other person. This captures the intuitive question people actually ask: "how often do our periods start around the same time?" Unlike $R$, this integrates over all onsets across the full simulation, giving much greater statistical power to detect transient alignment that recurs but does not persist.

= Experimental Design

The *main experiment* runs 200 paired Monte Carlo trials simulating 6 cohabitants over 3 years (1,095 days) under both the treatment (shared + individual stress) and control (all-individual stress) conditions, using identical cycle parameters but independent starting phases. The first 90 days are discarded as burn-in.

For high-power *parameter sweeps*, the simulation core was rewritten in Rust (with PyO3 bindings and rayon parallelism), enabling 1,000 trials per parameter point over 5-year windows. 1D sweeps vary shared stress rate (0--0.3/day), stress sensitivity (0--3$times$), and shared exposure probability (0--1.0) at 10--11 points each. 2D sweeps map the interaction surfaces of shared event rate $times$ stress sensitivity and shared exposure probability $times$ shared event rate on 7$times$7 grids (49 points, 1,000 trials each---49,000 simulations per heatmap).

= Results

== Main experiment

#figure(
  image("output/trial_distributions.png", width: 100%),
  caption: [Distribution of synchrony metrics across 200 trials. Left: phase synchrony (_R_) is indistinguishable between conditions. Right: onset proximity is shifted rightward under shared stress.],
) <fig-distributions>

The central result is a dissociation between the two metrics, which reveals the *temporal character* of the synchrony (@fig-distributions):

#table(
  columns: (1fr, 1fr, 1fr, 1fr, 1fr),
  align: (left, center, center, center, center),
  stroke: 0.5pt + luma(180),
  inset: 6pt,
  table.header[*Metric*][*Treatment*][*Control*][*Δ*][*_p_*],
  [Phase synchrony (_R_)], [0.388], [0.382], [0.006], [0.34],
  [Onset synchrony], [0.285], [0.259], [0.026], [5.9 × 10#super[−15]],
)

Shared stress produces a highly significant increase in onset proximity (_d_ = 0.63): periods start near each other about 10% more often than under independent stress. But the instantaneous phase measure _R_ shows no difference (_d_ = 0.02), meaning the alignment does not persist between onsets. This is the signature of *transient synchrony*---shared stressors repeatedly push onsets together, but differing natural periods pull them apart again before the next cycle. The Rayleigh test detection rate is ~8% in both conditions, consistent with the alignment being too brief to register in a single snapshot of 6 individuals.

#figure(
  image("output/synchrony_timeseries.png", width: 100%),
  caption: [Time evolution of synchrony metrics for a single representative trial. Top: _R_ fluctuates around the same level for both conditions. Middle: mean pairwise phase distance hovers near the uniform expectation ($pi slash 2$). Bottom: Rayleigh _p_-values dip below 0.05 transiently in both conditions.],
) <fig-timeseries>

The time series (@fig-timeseries) shows that transient apparent synchrony occurs in both conditions as a natural consequence of the small group size. The treatment does not produce sustained phase-locking.

== Parameter sweeps (1D)

High-power sweeps (1,000 trials, 5 years per point) confirm and refine the earlier findings.

#figure(
  image("output/large_sweep_shared_event_rate.png", width: 100%),
  caption: [Effect of shared stress rate (1,000 trials, 5 yr). Both $R$ and OSI show clear dose-response, with $Delta R$ peaking at 0.04--0.05 and $Delta$OSI peaking near 0.046 at rate $approx$ 0.15/day, then declining as cycle disruption dominates.],
) <fig-sweep-rate>

#figure(
  image("output/large_sweep_stress_sensitivity.png", width: 100%),
  caption: [Effect of stress sensitivity. $Delta R$ grows steadily to +0.05 at sensitivity 2--3. $Delta$OSI peaks at sensitivity 1.25--1.5 (+0.029) then declines---excessive sensitivity makes cycles chaotic, destroying the synchrony signal.],
) <fig-sweep-sens>

#figure(
  image("output/large_sweep_shared_exposure_prob.png", width: 100%),
  caption: [Effect of shared exposure probability. $Delta R$ increases above $p$ = 0.6. The absolute OSI values show that higher shared exposure reduces overall onset overlap (more correlated stress = more total cycle disruption), but the treatment--control gap remains positive.],
) <fig-sweep-exposure>

The 1D sweeps reveal:

- *Shared stress rate*: Clear dose-response. $Delta R$ reaches +0.045 at rate 0.10--0.15; $Delta$OSI peaks at +0.046 then declines above 0.15/day as high stress makes all cycles more irregular (@fig-sweep-rate).

- *Stress sensitivity*: Below 0.5$times$, the effect is negligible. $Delta R$ grows monotonically to +0.05 at sensitivity 3$times$, but $Delta$OSI shows a nonlinear optimum around 1.25--1.5$times$ (@fig-sweep-sens). This implies that only moderately stress-sensitive individuals would experience meaningful onset synchrony.

- *Shared exposure probability*: $Delta R$ requires exposure above $tilde$60% to emerge. Interestingly, $Delta$OSI is positive across all exposure levels, suggesting that even partial shared stress exposure produces some onset convergence (@fig-sweep-exposure).

== Parameter sweeps (2D)

Two-dimensional sweeps map the interaction surfaces where synchrony is strongest.

#figure(
  image("output/large_sweep_2d_shared_event_rate_x_stress_sensitivity.png", width: 100%),
  caption: [2D sweep: shared event rate $times$ stress sensitivity (7$times$7 grid, 1,000 trials each). Left: $Delta R$ is strongest at high stress rate and moderate-to-high sensitivity. Right: $Delta$OSI peaks in a similar region but shows a sharper falloff at very high sensitivity, where cycle disruption erodes the synchrony signal.],
) <fig-heatmap-rate-sens>

#figure(
  image("output/large_sweep_2d_shared_exposure_prob_x_shared_event_rate.png", width: 100%),
  caption: [2D sweep: shared exposure probability $times$ shared event rate. Left: $Delta R$ peaks at high exposure (>0.8) and moderate rate ($tilde$0.10). Right: $Delta$OSI shows a striking pattern---the largest effects occur at low exposure probability and high event rate, where fewer correlated events produce less total disruption while still driving onset convergence.],
) <fig-heatmap-exposure-rate>

The 2D heatmaps reveal two important interactions (@fig-heatmap-rate-sens, @fig-heatmap-exposure-rate):

+ *Rate $times$ sensitivity*: The largest synchrony effects ($Delta R$ up to +0.06, $Delta$OSI up to +0.047) occur when both shared event rate and stress sensitivity are moderately high. The effect is not simply additive---there is a ridge of maximal synchrony running diagonally, reflecting a balance between stress-driven convergence and stress-driven cycle disruption.

+ *Exposure $times$ rate*: This surface reveals an asymmetry between the two metrics. Phase synchrony ($Delta R$) requires _both_ high exposure and moderate rate. But onset proximity ($Delta$OSI) is largest at _low_ exposure and high rate---a surprising finding explained by the fact that lower exposure means less _total_ stress per individual (since fewer shared events are experienced), preserving cycle regularity while still producing correlated delays among those who are exposed.

== High-fidelity mechanistic model

To validate the findings of the stochastic oscillator, we implemented a full mechanistic model of the hypothalamic-pituitary-gonadal (HPG) axis based on *Fischer-Holzhausen & Röblitz (2022)*. This system uses 54 ordinary differential equations per individual to model the pulse-driven release of GnRH, the dynamics of LH and FSH, and the maturation of multiple ovarian follicles in waves.

The mechanistic model replaces the abstract Phase Response Curve with a biologically direct mechanism: shared stressors (modeling elevated CRH) inhibit the GnRH pulse generator in the hypothalamus. This inhibition naturally delays follicular maturation, lengthening the cycle.

#table(
  columns: (1fr, 1fr, 1fr, 1fr),
  align: (left, center, center, center),
  stroke: 0.5pt + luma(180),
  inset: 6pt,
  table.header[*Metric*][*Treatment*][*Control*][*Δ*],
  [Onset synchrony (±5d)], [0.462], [0.355], [+0.107],
)

Results from the mechanistic simulation (50 trials, 6 individuals, 2 years) show a significantly more potent effect than the stochastic model. The Onset Synchrony Index increased by **10.7 percentage points** (_p_ < 10#super[−13]), compared to the 2.6-point increase in the simplified model. This suggests that the nonlinear feedback loops of the HPG axis may actually amplify the synchronizing effect of correlated environmental perturbations.

= Discussion

Our simulation demonstrates that shared environmental stressors produce *real, recurring menstrual synchrony*---but it is transient rather than sustained. The mechanism operates through the asymmetric phase response curve: stress during the follicular phase delays ovulation, and when cohabitants experience the same stressor at overlapping cycle phases, their next onsets are pushed closer together. This is genuine synchronization of a biologically meaningful event (menstrual onset).

However, the alignment does not constitute *entrainment* in the dynamical systems sense. Entrainment would require the cycles to lock into a persistent phase relationship (_R_ $arrow$ 1 over time), as seen in coupled oscillators with continuous forcing. Instead, what we observe is *repeated transient convergence*: a shared stressor brings onsets together, but because each person's natural period is slightly different, the phases diverge again over subsequent cycles---until the next shared stressor creates another episode of proximity.

The large-scale parameter sweeps (129,000 simulations) reveal that the effect has a characteristic structure in parameter space. There is a nonlinear optimum for onset synchrony at moderate stress sensitivity ($tilde$1.25--1.5$times$) and event rates ($tilde$0.10--0.15/day)---too little stress produces no convergence, too much destroys cycle regularity. The 2D heatmaps expose an unexpected asymmetry: phase synchrony ($Delta R$) requires high shared exposure, but onset proximity ($Delta$OSI) is actually _strongest_ at lower exposure rates combined with higher event frequency. This makes physical sense---lower per-person exposure preserves cycle regularity while still producing correlated delays among those exposed to each event.

This pattern---recurrent alignment interspersed with drift---maps precisely onto the subjective experience that fuels the urban legend. Roommates notice the weeks when their periods overlap and remark on it. The intervening weeks of non-overlap are unremarkable and forgotten. Over months, the remembered coincidences create a strong impression of synchrony from what is a 2.6 percentage-point increase in overlap probability (28.5% vs. 25.9%). The mechanism is real; the perception inflates it.

== Limitations

+ *Simplified cycle model.* Our two-phase oscillator with AR(1) noise captures the first-order statistics of cycle length but not the full chaotic attractor ($D_c approx 5.2$) identified by Derry & Derry. A higher-dimensional model might show different sensitivity to perturbation.

+ *PRC from primate data.* The phase response curve is based on rhesus macaque endotoxin challenge, which may not directly translate to human psychological stress. Human stress effects are likely smaller and more variable.

+ *No pheromonal coupling.* We deliberately exclude any direct cycle-cycle coupling to isolate the stress mechanism. A combined model could explore interaction effects.

+ *Stress magnitude calibration.* The lognormal magnitude distribution is a reasonable prior but is not fit to specific human data. Calibration against cortisol response studies would strengthen the parameterization.

== Conclusion

The urban legend of menstrual synchrony is not wrong---it is inflated. Shared stressors genuinely increase the frequency of synchronous cycle onsets among cohabitants (Cohen's _d_ $approx$ 0.6), producing the recurring episodes of alignment that people notice and remember. But the effect is modest in magnitude (2.6 percentage points), transient in duration (alignment dissipates within a cycle or two), and does not constitute the sustained entrainment that the popular narrative implies. Large-scale parameter exploration shows the effect peaks at moderate stress sensitivity and event rates, with a nonlinear interaction surface that reveals a fundamental trade-off between stress-driven convergence and stress-driven cycle disruption. The subjective experience of "syncing up" reflects a real mechanism amplified by confirmation bias and the statistical inevitability that variable-length cycles will occasionally align by chance.

#v(1em)

#line(length: 100%, stroke: 0.5pt + luma(200))

#text(size: 9.5pt)[
  *References* (selected) \
  #set par(hanging-indent: 1.5em)
  Bull et al. (2019). Menstrual cycle length and patterns in a global cohort. _J Med Internet Res_, 22(6), e18749. \
  Derry, G.N. & Derry, P.S. (2010). Characterization of chaotic dynamics in the human menstrual cycle. _Nonlinear Biomed Phys_, 4:5. \
  Li, K. et al. (2023). Menstrual cycle length variation by demographic characteristics from the Apple Women's Health Study. _npj Digital Medicine_, 6:100. \
  McClintock, M.K. (1971). Menstrual synchrony and suppression. _Nature_, 229:244--245. \
  Schank, J.C. (2006). Do human menstrual-cycle pheromones exist? _Human Nature_, 17:448--470. \
  Strogatz, S.H. & Stewart, I. (1993). Coupled oscillators and biological synchronization. _Sci Am_, 269(6):102--109. \
  Wilson, H.C. (1992). A critical review of menstrual synchrony research. _Psychoneuroendocrinology_, 17(6):565--591. \
  Xiao, E. et al. (1998). Stress and the menstrual cycle: relevance of cycle quality in the short- and long-term response to a 5-day endotoxin challenge during the follicular phase in the rhesus monkey. _J Clin Endocrinol Metab_, 83(7):2454--2460. \
  Xiao, E. et al. (1999). Stress and the menstrual cycle: short- and long-term response to a five-day endotoxin challenge during the luteal phase in the rhesus monkey. _J Clin Endocrinol Metab_, 84(2):623--626. \
  Yang, Z. & Schank, J.C. (2006). Women do not synchronize their menstrual cycles. _Human Nature_, 17:433--447.
]
