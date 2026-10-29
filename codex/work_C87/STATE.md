# C87 handoff

[Completed] Read RESULT.md, then PROOF.md. No <=10-letter microscopic
bond or certified member n>=a. No new central-charge inference.

[Proved] On open fixed-multiset word spaces, uniform pair-dependent
phases of adjacent swaps are removed by the inversion-count gauge.
This includes arbitrary token numbers, interacting real diagonals and
configuration cycles. With the common-scalar C84 coefficient, point
sources and real diagonal conversions in that gauge, the resulting
Bloch spectrum is reciprocal and its pair-sum margin vanishes for n>=2.
Conversions that acquire spectator phases are outside the theorem.

[Proved] Position-only codes with fixed token number have polynomial
dimension. Fixed letter multiplicities make additive charge scalar on
the code. Connected negative-hopping components have unique ground
states. These are three distinct obstructions, not a universal exclusion
of sparse, approximate or fragmented codes.

[Proved] A full-Fock diagonal phase dressing preserves a connected
number-conserving quadratic iff it is exp(i F(N)) times on-mode phases.
Factorized nonconstant amplitude dressings of a scalar pair coefficient
produce mixed density terms. Exact witnesses in checks.py give distances
at least 1 and 9/4 from every Gaussian for the chosen examples.

[Proved] Frozen markers with k charged bits and two injective half-charge
source shifts need at least 2^(k+1)+binomial(k,floor(k/2)) states, allowing
source overlap. k=1 needs five. With four backgrounds and two carriers
this gives a>=11. Dropping payload-preserving port or frozen-projector
assumptions evades this particular count.

[Proved positive example] Four complementary projector source channels
give exactly Gaussian K2 on four logical modes, each channel separately
interacting. Full finite gadget dimension 80, fourth order interacting.
No spatial implementation or alphabet count. General interacting ports
remain open.

[Checked] checks.py passes. Exact gauge graphs through 210 vertices,
36-dimensional Fraction product inverse, all 256 three-mode diagonal
sign dressings, integer Gaussian identities and capacity counts.
Numerical carrier ED through dimension 72; gadget ED dimensions 32,64,80
at three lambda values, all 16 low levels with fourth-order errors.
checks.json stores outputs. input_audit.py reads all sibling prose,
Python and JSON, verifies the 118 inherited obstruction signs and C84
finite budget; it does not rerun the inherited root pipeline. Binary
build arrays were not loaded for orientation.

[Rejected checker, repaired] Fraction(numpy.int64) caused overflow and
an inverse-identity failure. Explicit Python int conversion repairs it.
rejected_checks.json records this. No rejected result is used.

[Scope] All work stays in work_C87. No network commands, subagents,
cluster jobs or writes outside the folder. The project STATE was read
but not modified because it is outside the user's allowed scope.
