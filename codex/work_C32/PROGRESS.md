# C32 final handoff

[Proved record] RESULT.md and PROOFS.md contain the final result. This is a
partial completion of the requested research task, not a closed route-A
construction with proved nonresonance for every n.

[Proved] A positive 28-letter folded-bus core has an exact 2^n code and
exact U(1)/parity. C8's complete wrapper gives d0=2736 and lengths n+3,n+4
with equal payload. The leading length motion is free. The full second-
order bond is equalizable in operator norm, with explicit missing-bond
and fast-wrapper errors, including scalar length interactions.

[Proved] PROOFS Sections 5 and 6 give a finite all-n crossing schedule.
The positive Gram matrix is a rational Lyapunov construction; a Cauchy
resolvent bound replaces the unquantified large-express-energy limit.
The schedule is extremely conservative and all constants are finite at
fixed n. It does not prove distinct centers or Sidon margins.

[Assumed, unresolved] A proof that the loop family avoids every dangerous
pair-sum identity for each n is still missing. The full stable-phase
hierarchy is conditional on that witness and on the stated local massive-
sector elimination/bosonization inputs. No unconditional c=n+1 claim.

[Numerical] Source scans at n=1,2 retain up to 7560 states per endpoint,
including all reachable helix excitations. Direct transport, charge and
raw-inverse checks pass at the errors recorded in RESULT. Small-signal
loss of accuracy is preserved in small_signal_diagnostic.json.

[Numerical] Two- and three-cell ED at n=1,2 uses active lengths and every
legal core excitation in the fast-wrapper limit, in Q=0,1 sectors. The
largest sector has 91872 states; the error is proportional to lambda^4.
This is not full microscopic ED. Initial clustered-eigenvalue ARPACK
failures were repaired by increasing the Krylov space. Final saved runs
have residuals below 2.4e-12. Each current solver call has a 180-second cap.

[Proved record] Python syntax and every JSON output were checked. Only
this folder was written. No network, cluster, subagents, commit or push.
Project STATE is outside the authorized write scope, so STATE.md here
contains the handoff instead.

## Remaining work

[Assumed, unresolved] (i) Prove the all-n nonidentity/Sidon lemma on the
actual small parameter family, or change the local construction to make
it provable. (ii) Extend ED to the full finite-wrapper microscopic mobile
model, rather than interpreting the controlled-limit ED as that test.
(iii) Apply the thermodynamic infrared inputs only after (i) succeeds.

## Earlier chronology

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
