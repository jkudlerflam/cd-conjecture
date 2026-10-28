# Addendum: assessment of the separate d=24 referee

**Status: assessment and recommendation, 2026-10-28.** I read `contrib/referee_d24/README.md` in full after completing the vacancy audit and finite-J cross-review. Its conclusions were not evidence used in those earlier calculations. This note records a change in recommended public wording, not a newly found instability or a retraction of the exact checks.

**Recommendation:** use **credible conditional counterexample candidate**, with **finite ingredients and the conditional infrared argument independently checked**. Avoid an unqualified headline saying that the counterexample itself is verified. The separate referee agrees with the vacancy mechanism at the physics standard and identifies no fatal error, but explicitly withholds established-counterexample status. That disagreement should be reported plainly.

## What changes the assessment, and what does not

| Referee objection | Assessment |
|---|---|
| The frozen c=23 conclusion P and vacancy matching M remain effective-theory inputs. | This is the material issue for the strength of the headline. Our audit gives substantive arguments for them and rules out several concrete failure modes; it does not replace the full finite-coupling matching/sector-selection conclusion by an independently proved lattice statement or measurement. The meaning of “verified” must retain those premises. |
| Final lambda, J/t and mu are not numerically specified. | Not itself a logical refutation. An existence proof may choose sufficiently small positive constants without displaying decimal values, and t can be fixed as the unit of energy. But here the final matching constants and the actual allowed interval are not independently certified. The delivered object is a parameterized construction with a proposed existence regime, not a numerically resolved violating bond. |
| The hierarchy and crossover are enormous. | Correct as an evidence and practical-validation issue, not a thermodynamic obstruction. Some quoted exponents arise from conservative sufficient inequalities; they should not be presented as proved necessary bounds on every possible working parameter choice. They do establish that the proposed certified-reference schedule is far beyond direct n=23 testing. |
| The low-energy charge sector admits changes in occupied number in multiples of 24. | Correct and worth stating explicitly. It is a charge-24 liquid: changing N by one creates a mismatch defect, while Delta N=24 can be gapless. This is the compactification/gluing described in our report and still contributes c=1. It is not a missing 24-fold multiplicity or the loss of the charge oscillator. |
| Unique ground state on every ring is unproved. | Agreed; it is not required by the original conjecture. The referee itself treats this as nonfatal and exhibits an allowed flux doublet. No status downgrade follows from this issue alone. |
| The vacancy can reintroduce the material period into the effective payload. | Correct. The generated period-24 modulation is already a cell-periodic correction in the proposed 24-site cell description. The low-order vertex test must use that cell Brillouin zone, as ours does. It would be dangerous to use this mechanism for a shorter-period payload without repeating the resonance analysis. |
| Exact half filling is unnecessary for stability. | Plausible and consistent with the physical-harmonic formula, but exact half filling is a useful clean choice because the global free-charge threshold proof is then available. An existence argument may choose it; removing that simplification is not needed to establish the proposed example. |
| The reviews are not independent because they use related AI agents. | The implementations and calculations have methodological independence, but shared models, source material, and coordination limit their independence as external scientific scrutiny. “Independently reconstructed by fresh related AI workers; no human review” is accurate. Simply counting agents as independent endorsements overstates the evidence. AI authorship does not by itself make a mathematical argument false, and human authorship is not a missing mathematical hypothesis. |
| Require a calibrated n=1 or n=2 central-charge measurement and human review before acceptance. | Reasonable proposed validation steps, not logically necessary conditions for every analytic counterexample. Small-n success would test the mechanism and implementation, but would not prove the extreme n=23 asymptotic transfer. A genuinely complete analytic existence argument could stand without those particular tests. |

The crossover estimate follows directly from comparing O(t b²/L) particle-number rounding with a dressed material defect cost J alpha gamma. Thus

    L_* ~ t b²/(J alpha gamma).

At b=24 the prefactor b²/(alpha gamma) is about 2.56e7. If one uses J/t at or below 1e−510, this gives a crossover of order 1e517 or larger. The numerical exponent is conditional on that chosen sufficient schedule; no universal lower bound on all possible realizations is established by this estimate. Once J is fixed and positive, an enormous but finite crossover does not prevent the thermodynamic limit.

## Where the two reviews actually disagree

There is no newly identified operator, seam inconsistency, charge stiffness error, or incompatible order of limits. The disagreement concerns whether standard effective-field-theory matching, supplemented by finite and second-order diagnostics, is enough to label this novel full construction an established physics counterexample. Our initial coordinating headline answered yes; the separate referee answers no while accepting the vacancy mechanism and the consistency of the hierarchy.

My earlier vacancy report already stated a conditional theorem with four hypotheses, and the finite-J addendum explicitly distinguished a conventional RG basin from the absent all-order microscopic per-defect inequality. Those remain my technical conclusions. I would not now compress them into “the model is verified” for a user seeking a decisive verification or refutation. Absence of a discovered failure is weaker than a completed verification of all matching premises.

A user-facing verdict faithful to both the calculations and the disagreement is:

> The finite construction and its conditional continuum analysis passed the audit, and no refutation was found. The full c=24 claim remains a strong conditional counterexample, relying on local finite-coupling matching and dressed-sector stability. A separate referee accepts the vacancy mechanism but does not regard the present evidence as establishing the counterexample.

This does not impose a rigorous lattice-to-CFT theorem as the only acceptable physics standard. It states exactly which part has been checked and preserves the substantive uncertainty in applying conventional matching to this particular construction. Future work could remove that uncertainty through a more explicit matching/sector comparison, persuasive calibrated demonstrations of this bond family, or outside review of the analytic argument. The evidence standard should be stated before calling the full result established.
