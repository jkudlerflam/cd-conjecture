# C85 progress

[Proved, algebraic] C84 has three distinct outgoing grammar states. Three vacancy colors give d=14. The first requested count is n=13, p=14 data sites, b=15. The carrier has 29 vertices and the product resolvent has 841 vertices.

[Completed] Regenerated C84 integer and rational pipeline at n=13, built the literal three-color lift, tested nonzero bits and illegal weight, transferred F3/F4, and derived the complete leading grammar self-energy.

[Scope] Work stays in this folder. No network commands, external modifications, or subagents. Project memory will be recorded locally in STATE.md.

[Proved, computer-assisted] The n=13 certificate passes at loop exponent 66. It has 26 simple roots, delta2>0.0156742281, delta4>3.2513114e-6, 351 unordered pair sums, and normalized allowance 2^-35. Material schedule g=2^-272, lambda=2^-608, K=2^2560 passes the physical-unit error budget.

[Numerical] Nonzero-bit two-cell closure dimension 3348 has illegal weight 8.30676e-5, 8.39127e-7, 8.39976e-9 at K=100,1000,10000. Full nonzero-port closures with one and two vacancies have dimensions 368 and 1945 and also show inverse-square suppression.

[Proved, algebraic; numerical audit through five sites] The leading grammar correction is diagonal. Its t^2/(JK) term counts |I(a)|-1 on occupied-vacancy contacts. Counts are 1 on L_0,L_1,D,X,E, 2 on R_0, and 0 on R_1,M. The second-order mixed t/K term vanishes. The J/K term contains the vacuum-tail bit field and a 2 lambda^2 wrong-end port self-energy. Every illegal intermediate has one mismatch. All leading formulas agree on every legal ring of lengths 3,4,5.

[Implementation checks] Two initial utility checks were corrected before acceptance. The alpha lower bound at denominator 16 is equal to 2^-18, so its rational assertion uses equality. Selector output now records outward-rounded 128-bit dyadic bounds, avoiding Python's decimal integer-string limit. Neither issue rejected a microscopic parameter.

[Completed, 2026-10-29] `python3 -B run_all.py` passes all seven scripts. The integer/rational certificate, selector, literal lift, exhaustive virtual-channel checks, resolved bit-dependent energy tests, full frozen sectors and local SW diagnostic all pass. `verification.json` records the run. `RESULT.md` gives the conclusion, `PROOF.md` gives the derivations and scope, `MEMBER.json` gives the exact schedule, and `STATE.md` gives the handoff.

[Open, inherited infrared inputs] No microscopic c measurement or unconditional counterexample was obtained. Finite-J matching, uniform mobile dressed-defect selection and interacting multichannel/conformal control remain assumptions. They are separated from the proved frozen F3 and finite-order local F4 results. No work outside this folder was performed.
