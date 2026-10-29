# C86 handoff

[Completed] No at-most-ten-letter certified member was obtained.
Read RESULT.md for the outcome and PROOF.md for all derivations.
The fallback strengthens the restricted obstruction, without claiming
a universal eleven-letter minimum or a new central charge.

[Proved] For the C84 common-scalar two-species coefficient, arbitrary
complex H>0 and pure endpoint sources u,v give
Re C_ij=(1/2) integral f_(u,i)(t) f_(v,j)(t) dt.
A unitary carrier reflection exchanging sources and conversion
projectors gives F=GJ, JAJ=A, and h(-k)=J h(k) J.
Herm(F) has blocks G_++ and -G_--; it cannot be definite on both
parities. At four or more distinct roots, reciprocal pair sums
coincide. Complex reflection-invariant flux and common port phases
do not remove these conclusions. Equal sources also force reciprocity.

[Proved, restricted] The eleven-letter lower bound now covers
single-token product cells with four binary L/R backgrounds, frozen
on-site markers, matched common species, nearest-neighbor token
motion, and payload-independent pair ports. At a<=10 only one
transport type fits. Its graph is a path and its sources point
vectors, so the C84 leading pair certificate is impossible.
This is broader than the faithful-fold count, but contextual storage,
unequal species, data-dependent ports, and different leading orders
remain outside it.

[Proved] A quadratic replacement with intercell hopping rank r has
at most 2r isolated Fermi points. A fixed-charge number-preserving
single species cannot flip one charged data bit. These are separate
restrictions, not a no-go for every single-species interacting model.

[Proved, computer-assisted] certify_obstructions.py gives 20
single-rail and 98 two-vertex-source rejection certificates using
256-bit floor-rounded solves. All decisions are rational. The
largest product dimension is 121. A two-vertex formula with unequal
masses is proved for all positive gapped parameters and verified at
12 rational points. independent_checks.py checks four enclosures by
direct Fraction inversion. consistency_checks.py checks the analytic
filter and reflection identities using complex finite graphs.

[Rejected checker, repaired] The first extended-source checker
incorrectly required positive entries. It was replaced by a signed
interval determinant bound. rejected_checks.json preserves the
failure. Tiny negative numerical Gram eigenvalues are roundoff-scale
diagnostics, never positivity decisions.

[Scope] Candidates failed before the positive LDL/root stages of
C84. No new root certificate, MEMBER, finite-K schedule, or infrared
claim. All files stay in work_C86. The project STATE was read but
not modified because it is outside the permitted scope. No network,
cluster, subagent, or independent referee was used.
