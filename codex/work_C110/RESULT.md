# C110. Semicontinuity with a vacuum-sector prescription

2026-10-29. Derivations are in [PROOF.md](PROOF.md). Small checks are in [CHECKS.md](CHECKS.md).

[Open] The requested bound remains unresolved. I did not construct a finite upper semicontinuous extension agreeing with physical A, prove finite F(d), or find fixed-d bonds with physical A tending to infinity. Several direct prescriptions fail, and a gauge-optimized cutoff construction gives a precise remaining comparison problem.

## Results

[Proved] Infimizing Q_L over boundary twists has the wrong sign for energy minimization. For the even-length ferromagnetic XX chain with hopping t and flux phi in [-pi,pi],

    E_L(phi) = -2t cos(phi/L)csc(pi/L),
    lim Q_L(phi) = 2t-6t phi^2/pi^2.

Thus the twist infimum tends to -4t, whereas physical A=2t. Its positive part gives zero. The twist supremum recovers 2t in this example. These are conserved U(1) twists on an already vacuum-compatible residue class.

[Proved, using C107] An infimum over all sufficiently long vacuum-compatible lengths is also not generally physical A. C107's critical XY chain at gamma=1/2 approaches A from below. Every finite tail infimum is strictly smaller than A. Recovering the limit requires a supremum over tail cutoffs, which does not preserve upper semicontinuity.

[Proved] Arbitrary seam twists are not gauge invariant. A coboundary leaves the seam operator o-UoU^dagger. For the single gauge orbit

    h_r = -(ZI+IZ)/2+r(ZI-IZ),

all periodic Hamiltonians are the same gapped product chain, with M=1 and A=0, but

    sup_U Q_L^U(h_r) = (6/pi)L[max(1,2|r|)-1].

The apparent divergence is a boundary and gauge artifact, not divergent physical A.

[Proved] Smooth energy coefficients can contain a nonzero gapped contribution even in the canonical Hilbert-Schmidt gauge. The explicit nearest-neighbor bond

    h_s = -(ZI+IZ)/2+s(XZ-ZX),   0<|s|<1/2,

has a unique product ground state, e=-1, gap at least 2-4|s|, and 1<=M<=1+2|s|. Nevertheless, for w_j=f(j/L) with f smooth, positive, periodic, and nonconstant, its weighted ground-energy deficit D_L=e sum_j w_j-E_0(sum_j w_j h_j) obeys

    (s^2/2) I_f <= liminf L D_L
               <= limsup L D_L <= s^2 I_f/(2-4|s|),
    I_f = integral_0^1 f'^2/f > 0.

Therefore its 1/L response cannot be identified with A=0. This is proved by an operator gap bound, a product trial state, and a matching bound on the order of the deficit. No extrapolation is needed. A different gauge removes this example's response, so optimization over gauges remains a possible repair.

## A remaining route

[Proved] For fixed nonnegative open-segment weights with total W, define

    D_w([h]) = inf_o {W e(h)-E_0(sum_j w_j(h+oI-Io)_j)}.

Then D_w is gauge invariant, finite, nonnegative, positively homogeneous, and upper semicontinuous, with

    0<=D_w([h])<=2W M([h]).

An infinite-volume translation-invariant minimizing state proves nonnegativity. Continuity at fixed o and the infimum prove upper semicontinuity. A norm-minimizing representative proves the upper bound.

For any fixed family of weights and positive normalizations, B=inf_j eta_j D_(w_j) has the same properties. A proof of A<=B would give a bound directly, with F(d)<=2 eta_j W_j for any member of the family.

[Open] No suitable family and comparison A<=B were established. Neither the correct conformal normalization nor the zero value for every gapped chain follows from this construction alone. Calling B the physical A without that comparison would change the problem. Finiteness of the canonical local upper envelope of physical A remains the missing local boundedness theorem identified in C107.

## Verification and scope

[Numerical] All checks in checks.py pass. Nine XX flux matrices and six gapped-chain matrices have dimension at most 256. Scalar sums verify the smooth trial-state coefficient. The analytic proofs do not depend on these computations.

[Work record] All files are confined to work_C110. No heavy numerics, cluster jobs, network shell commands, or shared-project edits were used. The local STATE.md records the handoff. These results have not had an independent referee review.
