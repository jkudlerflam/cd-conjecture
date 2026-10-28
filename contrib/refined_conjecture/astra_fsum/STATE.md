# f-sum bound, closed-book calculation

Completed 2026-10-27. This is the local memory for this calculation. The parent project's STATE.md was read but is outside the writable workspace, so it was not edited. No external sources were consulted. No other research files were changed.

Conventions: lattice spacing and hbar are one; positive-frequency Lehmann structure factor has 1/L and a delta function, without an extra 2 pi. Unless q and -q are related by a symmetry, use the average of their structure factors. The double commutator is the first moment of that average, not automatically of the unsymmetrized structure factor.

For D_a = -[t_1^a,[t_1^a,h]], the exact nearest-neighbor first moment is (1-cos q)<D_a>. The nonchiral affine currents have central term (k/2)m delta in the generator convention tr_fund(t^a t^b)=delta/2. Their static weight per generator is k|q|/(4 pi), and their first moment is k v q^2/(4 pi). Therefore k v dim(g) <= 2 pi <sum_a D_a>.

Writing C=C_2(R), M=half the spectral spread of h, and Phi(X)=sum_a t_1^a X t_1^a, one has D=2 Phi(h)-2 C h <= 4 C M I. Consequently k v <= 8 pi T_R M/d. Sugawara gives c_g v <= 8 pi C M/(k+hdual), and also the stronger expectation bound 2 pi <D>/(k+hdual). These bounds constrain the charged sector, not an arbitrary neutral coset contribution to the total Casimir energy.

Sharpness distinctions: 2 pi in the expectation bound is asymptotically sharp for SU(N) fundamental permutation chains as N tends to infinity. The algebraic coefficient 4 in D<=4CM is attained by an isolated SU(2) singlet with h=P. This does not establish optimality of the combined 8 pi coefficient for critical translation-invariant chains. The permutation family requires a coefficient at least 4 pi, and has ratio exactly 1/2 to the 8 pi bound.

SU(N), h=P: M=1, C=(N^2-1)/(2N), D=I-NP. Recalled thermodynamic Bethe-ansatz input gives e0=1+2[psi(1/N)+EulerGamma]/N and v=2 pi/N. The expectation-bound ratio is (N^2-1)/(N(1-N e0)). Values for N=2,3,4 are 0.846219983903, 0.857549394600, 0.871996548610. It tends to 1. e0=-1+pi^2/(3N^2)-2 zeta(3)/N^3+O(N^-4).

U(1): use central term kappa m, so S(q)=kappa |q|/(2 pi) and kappa v<=pi<D_Q>. For hopping h=-t(c1^dagger c2+h.c.), Q=n-1/2, D_Q=-h, <D_Q>=2t sin(kF)/pi, v=2t sin(kF), and kappa=1. The small-q f-sum bound saturates for every 0<kF<pi. Its norm bound v<=pi t is not saturated, even at half filling. The model is considered in a fixed-filling sector; no chemical-potential term is included in the quoted bond norm.

Checks in check_fsum.py (runtime about 1.3 seconds): SU(N) generator traces and double-commutator identity for N=2,3,4; local SU(2) algebraic saturation; exact six-site SU(2) Lehmann first moments at two momenta (maximum discrepancy 4.5e-16); exact free-fermion particle-hole moment integrals at three fillings and two momenta (maximum discrepancy 5.3e-18). The Bethe-ansatz energy and velocity were recalled and algebraically evaluated, not independently established by thermodynamic numerical extrapolation.

Necessary IR qualification: convergence to a positive current spectral contribution with the charge-fixed coefficient, equal left/right levels k and the stated velocity. Additional orthogonal low-energy weight is allowed and does not weaken the lower bound. Positivity alone does not permit subtracting an arbitrary formal current correlator if its realization as an IR spectral component has not been justified. Marginal logarithms are handled by the q to zero limit. Multiple independent factors obey factorwise bounds; abelian charge normalization and level matrices must be specified separately.
