# C82 progress

2026-10-28. Completed the requested local repair, proof-layer write-up,
and frozen-sector adaptation. Work and outputs are confined to work_C82.
Read RESULT.md for the conclusion and PROOF.md for all derivations.

[Proved] F5 is repaired with two vacancy colors, d=17. Vacancy color
transmits the preceding material grammar state through any vacancy run.
Every illegal material pair forces a local bad pair. Projected hopping
is the original vacancy hopping, with a unique allowed color assignment.
The finite-K bond retains all off-grammar transitions.

[Numerical] At J=t, the referee's uncolored leakage plateau is reproduced.
The repaired two-cell n=3 closure has dimension 4608 and illegal weights
1.07551e-4, 1.08213e-6, 1.08279e-8 at K=100,1000,10000. A literal
nonzero-port, nonzero-conversion four-site closure of dimension 1344
also shows K^-2 suppression. All leakage-test residuals are below
2.2e-11. Exhaustive four-site grammar and charge checks passed.

[Proved] F4 is supplied as a local finite-order normal form around
kappa N_bad, with kappa=K frozen and JK mobile. The interaction remainder
bound is independent of volume. It is not a global norm estimate, a
global low-energy-band gap, or an all-order infrared theorem.

[Proved, inherited zero-port estimate] F3 is closed for compatible frozen
rings by applying the boundary-energy lemma to complete finite-K cells.
No bad-pair superselection is used. The original b=16 and new b=17
parameter schedules satisfy the strict wrong-cell and marker-free bounds.

[Proved, computer-assisted] A new n=16 certificate avoids losing the
candidate violation when one vacancy color is added. It has 32 simple
roots, delta_2>0.0145298937, delta_4>5.92937e-7, all 528 pair sums checked,
and finite-g,K errors below one eighth of 2^-35. The final material
schedule uses g=2^-288, K>=2^1280, lambda<=2^-704, loop eta=2^-92 and
normalized field 2^-149. See budget_n16.json, not the raw root routine's
intermediate conversion threshold, for the chosen parameters.

[Proved scale constraint] K must be chosen after J for the repaired
mobile model. Require JK>>nu and nu^2/(JK) below the retained J times
logical tolerance. No finite-J matching window is certified numerically.

[Heuristic, conditional] The new count is d=17, c=16+1=17. It requires
finite-J matching, dressed defect stability, and the inherited Gaussian
infrared assumptions. C75's original n=15 alone would give c=16=d-1.
C81's d=16 lift was not duplicated or modified.

All six reproduction commands in RESULT.md passed. The occupied block
also matches C75's original builder entry by entry at two resolved points.
During development,
one syntax error in check_full.py and an integer-to-string digit limit
in budget_n16.py were corrected before the successful runs. No result
uses a failed run. No network commands or external writes were made.
