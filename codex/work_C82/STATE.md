# C82 handoff

Updated 2026-10-28. All edits are local to work_C82.

[Proved] Use two vacancies V_0,V_1 and the automaton in PROOF Section 1.
The literal bond is bond17.py. Every bad material pair forces a local
extended mismatch, even through vacancies. Projected hopping is ordinary
vacancy hopping. Finite-K leakage remains, with an actual penalty gap JK.

[Numerical] Small-ring ED reproduces the uncolored 0.0256 plateau and
shows K^-2 suppression for the repair at J=t. Dimensions are at most
4608. Full nonzero-port check dimension 1344. See vacancy_checks.json
and full_checks.json. No central charge was measured.

[Proved] PROOF Section 4 gives the extra finite-order local SW layer.
Section 5 closes frozen nonzero-port marker-sector selection using full
finite-K cells. Neither proof asserts mobile finite-J defect matching.

[Proved, computer-assisted] New n=16 finite certificate: 32 simple roots,
delta_2>0.0145298937, delta_4>5.92937e-7, budget 2^-35. Parameters are
D=1000, beta=1/d_17, t_d=beta 2^-92, g=2^-288, K>=2^1280,
lambda<=2^-704, normalized field 2^-149, b=17 midpoint selector.
Copied scalar routines retain C75 attribution. Reproduce with
certify_n16.py and budget_n16.py. The raw threshold g<=2^-264 is not the
final g choice.

[Assumed] Choose positive J/t in the inherited matching regime, then
increase K to suppress nu^2/(JK) relative to J times the smallest
retained logical tolerance. A numerical J window is not certified.
Finite-J matching, dressed defects, and Gaussian infrared stability
remain assumptions. Conditional candidate d=17,c=17; n=15 alone gives
c=16=d-1. No unconditional counterexample and no d=16 no-go theorem.
