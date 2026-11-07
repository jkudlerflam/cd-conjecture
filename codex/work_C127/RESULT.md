# C127 result

2026-11-06. Derivations are in [PROOF.md](PROOF.md), and numerical scope and failures are in [CHECKS.md](CHECKS.md).

[Assumed or unresolved] The general scalar bosonic majorant has neither been proved nor refuted. In particular this run does not establish a finite C(d) for interacting positive-transfer chains. Allowing arbitrary dimension-dependent species count and additive constant makes existence of that majorant equivalent, up to explicit constants, to the original all-temperature scalar trace target. The failed stronger statements below do not settle it.

## The insertion rank is not a spectral multiplicity

[Proved] At a nonzero Matsubara frequency, the map from one-site operators to insertion vectors has rank at most d^2-1. Its range need not be invariant under crossed evolution. The cyclic span obtained by applying spectral functions of the crossed Hamiltonian can be infinite-dimensional.

[Proved lattice example] In the actual XX chain, a single Sz insertion at fixed external omega_k=2pi k/beta detects positive spectral masses at

    E_(n,k)=asinh|(2n+1)pi/beta|
             +asinh|(2n+1-2k)pi/beta|.

For n>=k these are infinitely many distinct energies. At k=1 the number below a fixed E grows as beta sinh(E/2)/(2pi)+O(1), despite d=2. The reflected density correlator at odd spatial separation supplies the positive measure explicitly. This is a lattice failure of the proposed rank-to-species inference, not a violation of bosonic many-body counting.

[Proved conditional local cone] If the insertion correlator has its reflected positive spectral representation, the C125 estimate gives

    Hc >= |omega|/(4eM)

on its detected cyclic subspace for 0<|omega|<4eM. The prefactor is 2||A||||B||V/[|omega|(V-|omega|)], V=4eM. This does not cover undetected states. The high-frequency replacement is log(|omega|/(4M)), not a linear cone.

[Proved lattice failure] Any global momentum-preserving comparison to finitely many bosonic species of linear dispersion fails for XX, even with exp(K) and a finite zero-mode factor. The exact detected branch has energy growing as log k at fixed beta and momentum 2pi k/beta. Its trace in a sufficiently large fixed momentum sector decays more slowly with L than the reference sector trace. This does not refute an unlabelled scalar comparison.

## An all-temperature scalar majorant for XX

[Proved] For h=Sx Sx+Sy Sy, M=1/2, every beta>0 and even L>=2,

    log Ztilde(L)<=log 3+beta[csc(pi/L)-L/pi]
                 <=log 3+(2-4/pi) beta/L.

The proof applies convex midpoint quadrature to the exact crossed fermion product, retaining the parity projection and both offset vacua. The beta coefficient in the first bound equals the exact fixed-L zero-temperature slope.

[Proved] A literal bosonic majorant follows with

    S=2, V'=2e=4eM,
    K_XX=log 3-2log(1-exp(-1))=2.0159625794...,
    Ztilde(L)<=exp(K_XX) product_(k>=1)
                     [1-exp(-2pi k L/(V' beta))]^(-2).

This proves that the failure of a linear global energy-momentum bound is compatible with the desired scalar product bound.

## The constant, zero modes, and the missing input

[Proved lattice necessity] Every proposed universal oscillator-only bound must have K(d)>=log d. The diagonal ferromagnetic d-state chain has positive spatial transfer and

    Ztilde(L)=1+(d-1)[(exp(beta J)-1)/(exp(beta J)+d-1)]^L.

Along even L of order exp(beta J/2), beta/L tends to zero but log Ztilde tends to log d. These discrete zero-frequency sectors are not a U(1) stiffness contribution.

[Proved conditional charge estimate] If the relevant zero-frequency charge sectors are actually dominated by product_a sum_Q exp[-Q^2/(2 rho_a beta/L)], with sum_a rho_a<=4M, their log trace is at most 16M beta/L, even without bounding the number of charge directions. The missing assumption is the uniform quadratic spectral comparison and its normalization. The Kohn budget alone does not supply it or control neutral primaries.

[Assumed sufficient addition] An energy-filtered many-insertion comparison with a complete finite-species occupation space would close the route. A global linear dispersion is unnecessary. Even the lattice rate alpha(omega)=omega/(4eM) below 4eM and log(omega/(4M)) above it gives a convergent reference trace bounded by a constant times S M beta/L for L>=2. The needed all-state energy and multiplicity comparison is not established by one-insertion decay, Gaussian trial states, or time-bandwidth counting. PROOF sections 7 and 8 give the exact missing inequalities and explain where interactions enter.

## Numerical outcome

[Numerical] The exact XX scalar bounds pass 30 ring/aspect-ratio samples. Forty-eight momentum-resolved checkerboard matrices cover XX, XXZ at Delta=0.5 and 2, and Heisenberg, with m=2,3,4 and beta=0.25,1,2,8. All 240 scalar comparisons pass the candidate S=6,V=4eM,K=log 3. Momentum-sector comparisons fail. A coarse m=2,beta=8 Heisenberg matrix also violates the continuum low-frequency cone, which is not a finite-step theorem. Trotter errors and product-tail bounds are reported rather than treated as zero.

[Work record] Wrote PROGRESS.md, RESULT.md, PROOF.md, CHECKS.md, a local STATE.md, and reproducible Python checks. Larger m=5,6 transfer calculations are prepared in run_typhon.slurm but were not run or submitted. No outside file was changed. Tasks 2 and 4 remain incomplete for the general scalar majorant itself, and no lattice counterexample to that statement is claimed.
