# C127 independent result

2026-11-06. Derivations are in [PROOF.md](PROOF.md), numerical scope and failures in [CHECKS.md](CHECKS.md).

[Assumed or unresolved] The universal interacting scalar majorant is not settled. No lattice counterexample to that bound with arbitrary dimension-dependent K and S has been constructed. Tasks 2 and 4 therefore remain open in their unrestricted form. The results below isolate a false inference in the proposed route, prove the scalar bound for XX, and test the specified interacting chains.

[Proved] The single-insertion map at a signed Matsubara frequency has rank at most d^2-1, after removing the identity. Its range need not be invariant under crossed evolution. The corresponding spectral cyclic space can be infinite-dimensional. A single Sz insertion in the actual XX chain detects infinitely many distinct energies at fixed crossed momentum. Thus d^2-1 counts source directions, not oscillator species or detected eigenstates.

[Proved conditional on a compatible reflected insertion realization] The C125 estimate gives the lower spectral support bound |omega|/(4eM) for 0 < |omega| < 4eM, with correlator prefactor

    2||A||||B|| (4eM)/[|omega|(4eM-|omega|)].

The needed reflected positive spectral measure is additional to a positive crossed eigenvalue representation. At high frequency the available rate is logarithmic, not linear.

[Proved] For exact XX, including both physical spin parity sectors and all beta > 0,

    log Ztilde(L) <= log 3 + (2 beta/pi) L/(L^2-1)
                  <= log 3 + 8 beta/(3 pi L),  even L >= 2.

It follows that the proposed scalar bosonic majorant holds with

    S = 6, V' = 4eM = 2e,
    K = log 3 + 8/(3e) = 2.07962413179.

The proof uses the full exact fermionic product and convex midpoint quadrature. It does not assume a global linear dispersion.

[Proved] A momentum-preserving majorization by linear-energy bosons is false already on XX for every finite V'. Locally detected crossed energies grow logarithmically in momentum, while every state of that reference costs at least |p|/V'. This failure is compatible with the proved scalar XX majorant. A diagonal ferromagnetic d-state lattice chain separately forces K(d) >= log d in any oscillator-only scalar bound.

[Proved] Allowing arbitrary dimension-dependent S and K makes existence of the scalar bosonic product bound equivalent, up to explicit constants, to existence of K0+kappa M beta/L. Consequently the source-rank observation does not reduce the universal target to a counting fact.

[Proved conditional closure] A many-insertion ordered-energy comparison with an additive Fock reference would suffice. Even the softened dispersion

    alpha(omega) = omega/(4eM) below 4eM,
                   log(omega/(4M)) above 4eM

has a reference heat trace bounded by S c_alpha M beta/L for L >= 2, where c_alpha = 3.38828937933. If charge energies admit the stated Gaussian comparison and sum rho_a <= 4M, their entire factor obeys log Z0 <= log g+16 M beta/L, independently of the number of charge directions. What is missing is the interacting spectral comparison, including neutral sectors. Neither positivity, insertion rank, nor the stiffness budget proves it.

[Numerical] All 48 small checkerboard transfer matrices were resolved under discrete time rotation. At m = 4, beta = 2, one Sz frequency source detects 3, 8, 8, and 10 distinct energies for XX, Delta = 0.5, Heisenberg, and Delta = 2. All 192 scalar comparisons with Kstar pass; K = log 3 fails in two coarse m = 1 cases. Of 480 folded-momentum comparisons with Kstar, 56 fail. The original and crossed Trotter contractions agree to 4.03e-12. The exact XX product agrees with the physical trace to 5.59e-9 on the L >= 4 sample; the slower L = 2 tail is bounded explicitly. These computations do not establish a continuum interacting inequality.

[Assumed work record] Concurrent parent-folder edits replaced checks.py and its output schema during verification. They are preserved. This independent_audit subfolder contains this run's stable proof, scripts, outputs and handoff. All work stayed inside work_C127. Larger Trotter calculations are SLURM-ready but were not submitted. No general proof, universal counterexample, or independently verified minimal extra lattice condition is claimed.
