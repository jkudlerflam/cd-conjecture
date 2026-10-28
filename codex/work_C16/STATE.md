# C16 local state

Completed 2026-10-27. Scope is the compactness/singular-neighborhood part of
the d=2 program. Read RESULT.md for the complete tagged derivations and
PROGRESS.md for the work history. No counterexample or global proof.

## Main verdict

Compactness plus local boundedness only gives a finite unspecified global
bound if every point has such a neighborhood. It cannot exclude an isolated
finite-c>1 CFT. Neither exhaustive singular-neighborhood control nor a
classification of reachable c>1 CFTs is available. RESULT.md provides a finite
research triage list, explicitly not a finite certification theorem.

## Derived facts

- Product-eigenstate equations: Jxx=Jyy, Jxy=-Jyx,
  Bx=-(Jxz+Jzx), By=-(Jyz+Jzy). Mixed axial/transverse bonds still leak out
  of the one-magnon sector. The invariant U(1) dispersion has one quadratic
  minimum and no nonflat quartic minimum.
- Free XY+DM velocity collapse has a cubic slow branch and a linear branch.
  A representative is Jxx=1/2,Jyy=Jzz=0,Jxy=-Jyx=1/4,Bz=1/2.
- Positive uniform FF qubit bonds have exponential ground degeneracy only
  for a rank-one forbidden vector a tensor a. Entangled rank-one parents
  have open ground dimension L+1 by an invertible transfer map. Other
  long-chain cases have at most polynomial or bounded degeneracy.
- For H0=U sum n_i n_(i+1), all first-order two-site perturbations reduce to
  constrained hopping/current, chemical potential, and complex source.
  The explicit coefficient map is RESULT.md (17). Higher-order terms
  increase range and may be relevant.
- The number-conserving constrained gas has c=1, K=(1-rho)^2 and
  v=2t sin(pi rho/(1-rho))/(1-rho). This is a controlled slice, not an
  exclusion of every perturbation of the degenerate parent.
- Candidate Ising-product CFT tuning counts depend on translation and
  locality, not just the uniform on-site group. Mixed bosonic/derivative
  couplings can increase them. Bare rank J does not count dressed modes.

## Files and checks

checks.py: magnon/projection identities, cubic edge, constrained hopping
on L=12,15,18,21,24 at rho=1/3 with physical-spin entropy and Casimir energy.
checks_results.json: generated results. Half-size entropy c=1.01048,
with L^-2 correction 0.99750; Casimir c decreases to 1.004516. Within-ring
short-distance entropy fits overestimate c, still 1.07398 at L=24.

ff_checks.py and ff_checks_results.json: representative open and periodic
kernel dimensions on L=3,...,7. All match derived counts.
Both scripts run with python3, NumPy and SciPy only, in seconds.

## Concurrent work and remaining targets

C15 and C17 progress read at final review. C15 independently finds the cubic
edge and strong-exclusion K<1/2. Its new pairing-only factorization line is
an entangled rank-one FF example with L+1 ground states on even rings and a
quadratic defect branch; cited, not rerun. C17 finds no global exhaustion
or emergent-current exclusion either.

Next targets: the source/current constrained family with relevant virtual
corrections; an actual Ising/floating-phase intersection with operator
quantum numbers resolved; a microscopic bound on the number of dressed
fluids; uniform RG-onset control. No computation remains running.
All writes stayed in C16. No network, delegation, DMRG, parent-file edits,
commit, or push. Parent project memory is outside the authorized write scope.
