# C15 progress

2026-10-27. Read project STATE.md, task background, repository STATUS.md, and relevant portions of C6 and contrib/d2_classification/README.md. No preexisting local progress or sibling progress files were present. Work stays in this folder, Python with NumPy/SciPy, no network. No outside project-memory edit is permitted by the task.

## Scope and current conclusions

- [proved, given free theorem] Quadratic axial Z2 chain has c <= 1. Extension by the c-theorem requires a specified relativistic scaling regime; it does not cover a whole connected critical phase diagram, or vanishing-velocity endpoints automatically.
- [standard input / derivation in progress] K >= 1/2 is false globally in the NN XXZ family. At strong repulsion the no-adjacent-particle gas has K=(1-rho)^2, approaching 1/4 at rho -> 1/2 from below. This invalidates a proposed general relevance bound on atom-pair conversion.
- [derivation in progress] For independent fluids, cos(2 theta_1 - theta_2) has dimension 1/K_1 + 1/(4 K_2), with cross-density interactions requiring a full Gaussian covariance matrix. Symmetry alone does not exclude two fluids.
- [derivation in progress] The U(1) NN attractive t-V chain binds low-momentum pairs only at the same threshold as phase separation. Strong attraction favors longer clusters over dimers. The pair band has a single minimum shifted by twice the hopping phase, not two independently tunable minima.

## Next work

1. Derive the free parameter regions, classify controlled strong-coupling corners, and state exclusions precisely.
2. Analyze allowed operators for Ising times a compact boson, retaining symmetry/locality dependence instead of inventing a universal codimension.
3. Write reproducible small ED checks with entropy and spectral/energy calibration, plus a two-/few-particle binding test.
4. Finish RESULT.md and local STATE.md with the residual set and ranked next steps.

## First numerical and analytic checkpoint

- [numerical] `checks.py` completed, L=10,12,14,16 full parity-sector ED. Independent dense Kronecker check at L=4 differs by 1.7e-16. Ising and XX calibrated with entropy size growth and Casimir energy (known velocities). Constrained repulsive gas at (L,N)=(12,4),(15,5),(18,6) gives structure-factor K=.45325,.44929,.44738 and Casimir c=1.01829,1.01164,1.00805, approaching K=4/9,c=1.
- [numerical] Attractive t-V few-particle ED at L=24, N=1..4: below V=-2 the energy per particle decreases with N, supporting clustering rather than an isolated dimer fluid. An axial interacting attractive probe with (Jx,Jy,Jz,D,g)=(-.55,-.45,-.55,.1,0) has entropy estimates tending toward one, not two; this is not a thermodynamic exclusion.
- [proved / numerical] The free Ising/chiral collision (Jx,Jy,Jz,D,g)=(-.75,-.25,0,.25,1) has one cubic soft branch. Its entropy looks like c=1/2 but its dynamics is nonconformal. A nearby untuned interacting probe does not show c=3/2 on these sizes.
- [new unresolved corner] At t=D=0 and large positive V near the empty/CDW boundary, ordinary hopping is absent. Virtual pairing generates hopping by TWO sites at order Delta^2/V, with minima 0 and pi. The effective constrained Hamiltonian has diagonal -sum_{00 bonds}1/(1+n_left+n_right) and allowed two-site hops -1/2, in units Delta^2/V. This corner is not covered by the regular no-adjacent gas. `singular_corner.py` tests its few-particle spectrum by staggered charge Q. Third-order virtual processes create separated opposite-sublattice pairs; the second-order U(1) total particle conservation is accidental. Must keep this as residual unless the dilute instability is resolved.

## Singular-corner resolution at leading order and on an exact line

- [proved] In a fixed-N independent-set configuration let N_3 count cyclic particle separations equal to 3. The second-order pairing-only Hamiltonian is
  H_eff = -L + 2N + graph_Laplacian + (2/3) N_3,
  in units Delta^2/V. Each allowed hop by two sites is a graph edge of weight 1/2. Hence E >= -L+2N. States supported on one sublattice attain equality for every N<=L/2. These are ferromagnetic degeneracies, not two stable fluids. `singular_corner.json` confirms this and shows the mixed-staggered-charge sectors above them at L=24,N<=6.
- [proved] There is an exact pairing-only line t=D=0, 0<2|Delta|<V,
  mu_f=V-sqrt(V^2-4 Delta^2).
  The local bond shifted by mu_f/2 is a positive rank-one projector times V. An invertible staggered on-site filter maps its zero-energy constraints to ferromagnetic symmetric-subspace constraints. On even rings the ground-space dimension is L+1. A closed one-defect sector above a Neel ground state has energies V +/- sqrt((V-mu_f)^2+4 Delta^2 cos^2(k/2)), so its lower band is quadratic. This line is nonconformal, not c=2 or c=3/2.
- [numerical] `factorization.py` verifies ground degeneracy 5,7,9 at L=4,6,8 and the exact defect branch to 3.2e-15. The branch need not be the absolute first excited state.
- [open] Perturbations near this singular line do not inherit a relativistic c-theorem bound. Their full classification is residual. The ordinary strong-repulsion LL argument is restricted to hopping dominating virtual corrections and to a compact interior density window.

## Completed deliverables and final audit

`RESULT.md` is complete, with tagged step-by-step derivations, numerical tables, reproduction commands, a precise operational controlled set and residual set, and ranked next steps. `STATE.md` is the local handoff memory. The parent project's STATE.md was not edited because this task forbids writes outside the folder.

The extra singular-neighborhood probe at (Jx,Jy,Jz,D,h)=(.49,-.51,2,.01,sqrt(15)) has c_chord=.86784,.52258,.38250,.26271 at L=10,12,14,16. Half entropy is nonmonotonic and its last size-step estimate is negative. It does not support c>1 and is not a tuned multicritical exclusion.

`operator_checks.py` verifies the full second-order Schrieffer-Wolff operator at L=6 with zero error and the separated-pair third-order coefficient 2/3 at L=8. `factorization.py`, `corner_probe.py`, and `summarize.py` provide reproducible records for the additional checks. A temporary heredoc terminator typo in corner_probe.py was corrected before its successful run; no results were produced by the failed invocation.

Final caveats retained: small ED is not a global exclusion, the Ising-times-LL codimension is embedding-dependent, K>=1/2 is not a universal bound, and the c-theorem does not control arbitrary unequal-velocity or nonconformal lattice neighborhoods. No theorem excludes every emergent two-fluid state at intermediate coupling. No c>1 counterexample was constructed.

C16 and C17 progress files were reread at completion. They independently cover some shared elementary diagnostics; the new detailed pairing-only effective Hamiltonian and exact degeneracy line are documented here. No delegation or edits to sibling folders were performed. All computations have completed.
