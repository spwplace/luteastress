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

#text(size: 9.5pt)[
  *References* (selected) \
  #set par(hanging-indent: 1.5em)
  Bull et al. (2019). Menstrual cycle length and patterns in a global cohort. _J Med Internet Res_, 22(6), e18749. \
  McClintock, M.K. (1971). Menstrual synchrony and suppression. _Nature_, 229:244--245. \
  Schank, J.C. (2006). Do human menstrual-cycle pheromones exist? _Human Nature_, 17:448--470. \
  Xiao, E. et al. (1998). Stress and the menstrual cycle: follicular vs luteal response to stress. _J Clin Endocrinol Metab_, 83. \
  Zavala, E. et al. (2020). Modelling the HPA axis and GnRH pulse generator crosstalk. _Endocrinology_. \
  Phillips, A. et al. (2017). A detailed predictive model of the mammalian circadian clock and its entrainment by light. _PLoS Comput Biol_.
]

#v(1.5em)

= Introduction
The phenomenon of "menstrual synchrony" (the McClintock effect) remains one of the most controversial topics in human behavioral biology. Since McClintock's seminal 1971 paper, numerous studies have both supported and refuted the existence of cycle entrainment among cohabiting women. The primary mechanism proposed is pheromonal communication (Schank, 2006). However, an alternative hypothesis remains under-explored: *shared environmental stressors*.

This study investigates whether shared external stressors (work deadlines, financial stress, household disruption) can act as a "common clock" that entrains cohabitants' cycles through phase-dependent delays.

= Methodology
We employed two distinct modeling approaches:

1.  *Stochastic Phase-Delay Model (SPD)*: A phenomenological model where stressors induce discrete delays in the follicular phase. This model assumes a "Phase Response Curve" (PRC) where the follicular phase is sensitive to stress, but the luteal phase is not (Xiao et al., 1998).
2.  *High-Fidelity Mechanistic Chain (HFMC)*: A multi-scale ODE model implemented in Rust. It simulates the complete entrainment chain:
    - *Shared Schedule*: Co-habitants experience shared and independent stressors.
    - *Circadian Oscillator*: A master clock (Phillips 2017) coupled to light/schedule.
    - *HPA Axis*: Cortisol dynamics with slow-genomic feedback ($C_e$).
    - *GnRH Pulse Generator*: Cortisol-mediated inhibition of pulse frequency (Zavala 2020).
    - *HPG Axis*: Multi-follicular recruitment and hormone feedback (Fischer-Holzhausen 2022).

= Results

== Stochastic Model Findings
In the SPD model, shared stressors (rate=0.05/day, magnitude=1.0) produced a robust and statistically significant increase in synchrony. The Onset Synchrony Index (OSI) increased from $0.35$ (control) to $0.46$ (treatment), representing an 11 percentage point alignment boost ($p < 10^{-10}$).

== Mechanistic Washout
Counter-intuitively, the High-Fidelity Mechanistic model (HFMC) failed to produce significant entrainment. In a validation sweep ($N=50$ trials, 365 days), the OSI remained statistically indistinguishable between cohabitants and independent individuals ($p = 0.87$). 

The *Biological Entrainment Chain* appears to act as a low-pass filter, where the stochasticity of follicle recruitment and the internal inertia of the HPG axis overwhelm the weak coupling provided by HPA-GnRH crosstalk.

= Conclusion
While simplified phase-delay models suggest that shared stress *could* explain synchrony, high-fidelity biological modeling suggests that the necessary physiological coupling (HPA-to-GnRH) is likely too weak to overcome the inherent variability of human menstrual cycles. "Menstrual synchrony" in cohabiting groups is more likely a result of stochastic overlap (Schank, 2006) than genuine biological entrainment via shared stress.