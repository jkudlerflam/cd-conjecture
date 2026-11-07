# C131 numerical record

2026-11-07. All entries in this file are [Numerical] unless explicitly marked otherwise. Multistart local optimization gives lower bounds on the requested suprema, not certified global values.

## Results and what was actually run

Twelve distinct starting directions were optimized separately at L=2,3,4. The pilot allowed 1000 objective evaluations per start. The refined search repeated the same starts with a cap of 4000. The larger-ring multistart computations were not run locally. They are prepared in submit_c131.slurm for Typhon, which this sandbox cannot access.

| L | Parity | Best Q found by multistart | Location | Heisenberg Q at this L | Status |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 2 | even | 4.0000000000 | Heisenberg or its even-ring staggered-rotation image | 4.0000000000 | 12 starts, rings at most 8 sites |
| 4 | even | 2.7912527285 | Heisenberg or its even-ring staggered-rotation image | 2.7912527285 | 12 starts, rings at most 8 sites |
| 6 | even | not computed | not determined | 2.6179242962 | baseline only, multistart deferred to SLURM |
| 8 | even | not computed | not determined | 2.5582473161 | baseline only, multistart deferred to SLURM |
| 3 | odd | 2.0369383345 | DM and field bond listed below | -7.8166538264 | 12 starts, rings at most 8 sites |
| 5 | odd | not computed | not determined | -7.7937837699 | baseline only, multistart deferred to SLURM |
| 7 | odd | not computed | not determined | -7.7446742780 | baseline only, multistart deferred to SLURM |

The exact mathematical suprema remain unknown in this run, including at L=2. [Proved, inherited] Heisenberg attains Q_2=4, so C_even(2)>=4. The search does not prove the reverse inequality.

The normalized even-ring maximizing representative is

    (Jx,Jy,Jz,Dx,Dy,Dz,bx,by,bz)=(2,2,2,0,0,0,0,0,0), M=1.

The representative (2,-2,-2,0,0,0,0,0,0) has the same spectra on even rings by a staggered pi rotation about x. It is not a simultaneous on-site rotation of Heisenberg; the equality uses the even circumference. Both occur among the retained starts.

The best odd-L=3 candidate, in the coefficient convention of PROOF.md, is:

    1.0686113388, 1.0078083333, 1.0055092227, -1.6742951001, 0.3056922939, 0.2926993673, 0.9048978594, -0.1661910178, -0.1674606636

It comes from seed 1. Its gauge bracket before M=1 rescaling is [0.377276521975, 0.377276729581], width 2.08e-07. The dual solver's termination flag is False. Its feasible dual point still supplies a numerical lower bound; the finite bracket width, not a success flag, is the relevant normalization check. The true value at this numerical bond, ignoring floating-point roundoff, lies between 2.036938335 and 2.036939455. No phase or asymptotic stratum is assigned from a six-site spectrum.

## Search and accuracy diagnostics

The nine coefficients specify every two-qubit bond up to scalar/coboundary gauge and simultaneous rotation. Eight angles parametrize their Euclidean sphere; a different orthogonal chart at each start avoids polar singularities at symmetric seeds. Since Q is invariant under positive rescaling, this spans the M-normalized search space. The field sign is +b.S, so b_z=-B in C129 notation.

For each trial bond, Nelder-Mead minimizes the convex half-spread over three coboundary coefficients. Each finalist is recomputed with six gauge starts and checked against the dual problem

    maximize tr(Xh),  tr X=0, tr(XG_a)=0, ||X||_1<=1,
    G_a=S_a tensor I-I tensor S_a.

Any such X obeys tr(Xh)=tr[X(h+g)]<=||h+g|| for every gauge g, so it gives a lower bound on M. The finite-dimensional solver uses the twelve gauge-orthogonal Hermitian basis matrices and retains the best feasible rescaled iterate, rather than trusting the final SLSQP iterate, to satisfy the trace-norm constraint. These are floating-point brackets, not interval-arithmetic certificates. A primal estimate of M from minimization is an upper bound on M; for a positive defect it underestimates Q, so inaccurate normalization can conceal a larger value. The brackets quantify that issue at the saved finalists.

The spin Hamiltonian includes all oriented ring bonds, including both at N=2. Its lowest energy is computed in the full spin space, so Jordan-Wigner parity and charge sectors are not discarded. Dense diagonalization is used through N=6, sparse eigsh above it. Every finalist is recomputed at eigsh tolerance 1e-12 and its eigenvector residual is saved. A residual controls proximity to an eigenvalue, not a certified global lower bound on the spectrum. The tests below also compare against independently constructed full matrices and exact free energies.

| L | Starts | Optimizer termination criterion reached | Budget exhausted | Largest finalist eigenvector residual | Largest M bracket width across finalists |
| :--- | :--- | :--- | :--- | :--- | :--- |
| 2 | 12 | 6 | 6 | 9.33e-16 | 5.63e-05 |
| 3 | 12 | 12 | 0 | 2.68e-15 | 5.02e-07 |
| 4 | 12 | 10 | 2 | 1.32e-13 | 1.02e-05 |

checks.py passed the following checks:

- Independent sparse-matrix construction at N=2,3,4,5,6 with generic DM and field bonds. Maximum matrix-entry error 8.88e-16.
- Gauge minimization compared with the survey's six-start implementation at three generic points. Maximum difference 2.22e-16. Scalar/coboundary invariance also passed.
- Ferromagnetic XX exact energies for every N=2,...,10, including both parities. Maximum energy error 5.33e-15.
- Heisenberg baselines through N=16. Q_2 agrees with 4 within floating-point precision.

The positive finite-size defect of the gapped free example in PROOF.md, Section 7, was also checked directly: E_2=-1.118033988749895 and E_4=-2.135779205069858 agree with the displayed radicals to 1.4e-15.

## Interpretation and unfinished computation

The best even-length values found at L=2 and 4 are at Heisenberg or an even-ring equivalent. The separate Heisenberg benchmarks decrease through L=8. These facts support the sharp conjecture C_even(2)=4, but the data do not establish that the global maximizer remains there at L=6 or 8. There is no measured large-L drift or growth test yet. The odd maximum at L=3 occurs elsewhere and stays below 4; it should not be conflated with the even-ring sequence.

The full array covers even L=2,4,6,8 and odd L=3,5,7 separately, with 32 starts and 4000 objective evaluations per start, sparse rings through 16 sites. It has not been submitted. RUNNING.md gives the commands. No upper bound on the full parameter-space supremum or large-L conclusion is claimed from these local searches.

A numerical fraction of critical surfaces covered by C130 is not reported. Its radii and an ambient perturbation theorem are unavailable, so such a fraction has no defined numerical input. This is separate from the completed finite-ring search.
