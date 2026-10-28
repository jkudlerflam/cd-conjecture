# C32 progress

## Scope and inputs

[Proved record] Read project STATE, STATUS and the requested spatial, compile,
C13, C20, C25 and C8 padding, helix, vacancy and selector sources. No prior
PROGRESS existed. All writes remain in this folder. Project STATE is outside
the allowed write scope, so a local STATE handoff will be written instead.

## Current findings

[Proved from the supplied constructions] C8 padding transport is a code
identity, not a local parent theorem. Its local vacancy realization requires
an initialized frustration-free core. The supplied 23-letter sine-vacuum
core does not directly meet that hypothesis. C25 also disproves exact
second-order equalization by endpoint compensation alone at finite scales.

[In progress] Construct a positive, flat-vacuum version of the folded bus,
with the same 23 core letters and exact U(1). Apply the initialized C8
wrapper at d0=96*23+48=2256, subject to checking all local parent terms.
Generalize the finite Schur-complement bound to the noncommuting folded-bus
source sectors. All-n crossing and distinct-center assertions will be
separated. No proof of the latter has yet been obtained.

## Next

Write explicit parent and wrapper matrices. Check finite cores and mobile
projected transport, retain scalar length interactions, then carry out
bounded two- and three-cell diagonalizations. Record any missing general-n
lemma rather than promoting a numerical pattern to a proof.

## Positive core and wrapper constructed

[Proved construction, numerical core checks] The final core adds five left
marker letters to the 23-letter alphabet. Both markers have vacuum and
three excited states plus an empty state. The vacuum moves on a positive
path Laplacian. Excited hopping uses positive edge matrices and positive
onsite conversion blocks. There are n data sites and two markers, k=n+2,
with an exact 2^n code. The C8 wrapper therefore has d0=96*28+48=2736.
The preliminary 23-state wrapper idea was replaced by this endpoint design.

[Numerical, core dimensions 30 and 88] Code residual is zero; n=1,2
core gaps are positive. Exact source calculations retain up to 7560 states
per end, including reachable helix excitations. At a very small hopping
signal the n=2 traceless subtraction loses useful relative accuracy. That
run is preserved in small_signal_diagnostic.json. The check is being rerun
at a larger diagnostic conversion, and the scalar errors remain reported.
No small-signal result is being treated as a certified band calculation.

[In progress] The full mobile microscopic Hilbert space is too large for
three-cell n=2 diagonalization within minutes. The bounded ED will retain
all legal core excitations and active mobile lengths in the controlled
fast-wrapper limit. This limitation will be stated in the result.
