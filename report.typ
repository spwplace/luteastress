#set document(title: "Shared Stress and Menstrual Alignment: Updated Computational Results", author: "Lutea Stress Project")
#set page(margin: (x: 1in, y: 1in), numbering: "1")
#set text(font: "New Computer Modern", size: 11pt)
#set par(justify: true, leading: 0.65em)
#set heading(numbering: "1.")
#show heading.where(level: 1): it => block(above: 1.4em, below: 0.8em, text(size: 13pt, it))
#show heading.where(level: 2): it => block(above: 1.2em, below: 0.6em, text(size: 11.5pt, it))
#show figure.caption: it => text(size: 9.5pt, it)

#align(center)[
  #text(size: 18pt, weight: "bold")[Shared Stress and Menstrual Alignment]
  #v(0.6em)
  #text(size: 12pt, fill: luma(80))[Updated computational results from phenomenological and mechanistic models]
  #v(1.2em)
  #text(size: 10pt)[Lutea Stress Project #h(1em) February 10, 2026]
]

#v(1.2em)

= Abstract
We tested whether shared psychosocial stress can drive measurable menstrual alignment among cohabitants. The study compared paired treatment and control groups under two model families: (1) a stochastic phase-delay model (SPD) with follicular stress sensitivity and (2) a mechanistic chain model linking shared stress, circadian phase, HPA signaling, and ovulatory timing. In the SPD experiment (200 trials, 6 individuals, 3 years per trial), shared stress increased onset-window overlap (OSI difference $+0.0263$, 95% CI [$0.0203$, $0.0323$], one-sided Wilcoxon $p=2.34e-14$) but did not increase phase concentration ($R$ difference $-0.0139$, $p=0.75$). In the mechanistic experiment (50 trials, 6 individuals, 1 year per trial), OSI was slightly higher with shared stress (difference $+0.0165$) but not statistically significant (95% CI [$-0.0252$, $0.0582$], $p=0.257$). Realized stress sharing in the mechanistic runs was present (mean pairwise stress correlation $0.299 +/- 0.014$ SD). The updated evidence supports a narrow claim: shared stress can increase apparent onset proximity in simplified settings, but robust biological entrainment is not supported in the higher-fidelity model.

= Introduction
Menstrual synchrony remains disputed. Early work proposed true social entrainment, while later analyses argued that much of the pattern can emerge from stochastic overlap and selective observation. Most debate has focused on pheromonal or social-cue mechanisms. A separate hypothesis is that cohabitants experience correlated stressors and schedules, which could create weak shared forcing on cycle timing.

This project asks a targeted question: does shared stress produce *detectable alignment* relative to independent-stress controls, and does that signal survive when moving from a simple cycle model to an explicit neuroendocrine chain?

= Methods
== Study Design
For each trial we ran a paired comparison:
1. Treatment: cohabitants exposed to shared plus individual stress.
2. Control: matched individuals exposed only to independent stress.

Primary endpoints were:
1. Mean resultant length ($R$): phase concentration at the end of simulation.
2. Onset synchrony index (OSI): fraction of onset pairs within a plus/minus 5 day window after burn-in.

Directional tests used one-sided paired Wilcoxon signed-rank tests (treatment greater than control), consistent with the original entrainment hypothesis.

== Model 1: Stochastic Phase-Delay (SPD)
The SPD model represents cycle timing as follicular plus luteal stages with stochastic variability. Stress delays are applied primarily in the follicular-sensitive regime. Core settings:
1. 6 individuals per group.
2. 3 years simulated per trial.
3. 200 Monte Carlo trials.
4. Burn-in of 90 days.
5. Shared event rate 0.05/day, individual event rate 0.07/day, shared exposure probability 0.85.

== Model 2: Mechanistic Chain
The mechanistic model (Rust backend) simulates a causal chain:
1. Shared stress schedule.
2. Circadian oscillator and zeitgeber forcing.
3. HPA/cortisol dynamics.
4. Stress-modulated GnRH/KNDy gating.
5. Multi-follicular ovarian progression with ovulation tracking.

Mechanistic settings:
1. 6 individuals per group.
2. 365 days per trial.
3. 50 Monte Carlo trials.
4. Burn-in of 60 days.

We also recorded realized pairwise stress correlation in treatment groups to verify that the shared-stress generator produced correlated exposures.

= Results
== SPD Outcomes
Shared stress did *not* increase terminal phase concentration:
1. Treatment $R = 0.3551$.
2. Control $R = 0.3690$.
3. Difference $-0.0139$ (95% CI [$-0.0518$, $0.0240$], Wilcoxon $p=0.750$).

Shared stress did increase onset-window proximity:
1. Treatment OSI $= 0.2838$.
2. Control OSI $= 0.2575$.
3. Difference $+0.0263$ (95% CI [$0.0203$, $0.0323$], Wilcoxon $p=2.34e-14$).

Rayleigh detection rates were similar in treatment and control (both 6.5%), reinforcing that the main detectable signal was onset proximity rather than global phase locking.

== Mechanistic Outcomes
The mechanistic model showed a small positive OSI difference that was not statistically significant:
1. Treatment OSI $= 0.3119$ (SD 0.1048).
2. Control OSI $= 0.2954$ (SD 0.0984).
3. Difference $+0.0165$ (95% CI [$-0.0252$, $0.0582$], Wilcoxon $p=0.257$).

The stress process in treatment cohorts did produce correlated exposure:
1. Mean pairwise stress correlation $=0.299$.
2. SD across trials $=0.014$.
3. Range $[0.269, 0.328]$.

These diagnostics suggest the null/weak mechanistic effect is not explained by absent shared stress input.

= Discussion
Across model families, the updated pattern is consistent:
1. In the simple SPD system, shared stress increases overlap in onset windows.
2. In the mechanistic system, that effect attenuates and becomes statistically uncertain.

This split supports an interpretation in which stress can generate *apparent* alignment in coarse metrics without producing strong entrainment once endocrine and cycle-internal dynamics are explicitly represented. The mechanistic chain likely acts as a low-pass filter: weak shared forcing competes with substantial internal variability in follicle recruitment and hormone feedback.

Limitations remain important:
1. Non-romantic cohabitant stress-correlation targets are uncertain in current literature.
2. Mechanistic runs used 1-year horizons and moderate trial count; small effects may require larger powered runs.
3. OSI (plus/minus 5 day window) captures apparent proximity, not strict phase locking.
4. Parameter sensitivity for mechanistic-specific couplings (for example, zeitgeber strength and stress-to-cortisol gain) should be expanded in future sweeps.

= Conclusion
The revised experiments do not support strong stress-driven menstrual entrainment as a general mechanism. Shared stress can produce modest increases in onset proximity in a simplified phase-delay model, but the higher-fidelity mechanistic model yields only a small, non-significant difference under current calibration. The most defensible conclusion is that shared stress may contribute to apparent synchrony signals, while stochastic overlap and internal cycle variability remain dominant drivers.

#v(1.0em)
#text(size: 9.5pt)[
  *Selected references* \
  #set par(hanging-indent: 1.5em)
  McClintock, M.K. (1971). Menstrual synchrony and suppression. _Nature_, 229:244-245. \
  Yang, Z., and Schank, J.C. (2006). Women do not synchronize their menstrual cycles. _Human Reproduction_, 21(2):466-471. \
  Harris, A.L., Vitzthum, V.J., and Thornburg, J. (2015). Symptoms are associated with symptom severity but not cycle synchrony in large groups. _Human Nature_. \
  Xiao, E. et al. (1998). Stress and the menstrual cycle: follicular vs luteal response. _J Clin Endocrinol Metab_, 83. \
  Zavala, E. et al. (2020). HPA-GnRH crosstalk modeling framework. _Endocrinology_. \
  Phillips, A. et al. (2017). Circadian clock entrainment modeling. _PLoS Comput Biol_. \
  Bull, J.R. et al. (2019). Menstrual cycle length and patterns in a global cohort. _J Med Internet Res_, 22(6), e18749.
]
