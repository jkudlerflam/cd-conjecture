# Real history marks in the published moving-cutoff convention

Status: **[proved real C2 covariance extension and primary-source scope
check]**. This note does not identify a quartic-only constant-history theorem
with a fully frozen velocity theorem. The transfer is addressed separately
in [CONSTANT_HISTORY_COMPARISON](../vacancy/CONSTANT_HISTORY_COMPARISON.md).

Author: Codex `infrared_theorem`, 2026-10-28. Sections 1--3 independently
checked by the separate Codex agent `alternative_attack`; no human review.

## 1. The cutoff convention can be matched without complex continuation

[Porta--Scola--Singh, (4.1), (4.59), (4.75), (4.85)](https://arxiv.org/html/2509.08665v1#S4.SS3)
use a smooth radial cutoff with plateau through 1 and support ending at
3/2, with its elliptic norm determined by the current real velocity. This
permits the following precise real-history extension of MARKED_COVARIANCE.
The fixed-cutoff complex-analytic lemma there remains valid in its own
convention; here we make no holomorphic claim for a moving smooth cutoff.

Use the notation D_j,L_j,z_(0,j),z_(1,j),kappa_j of that note. Set

    r_v(p)=sqrt(p_0^2+v^2 p_1^2),       s_j=delta 2^j,
    chi_j(p;v)=chi(r_v(p)/s_j),
    v_j^new=(v_j+z_(1,j))/(1+z_(0,j)),
    f_j(p)=chi_j(p;v_j)-chi_(j-1)(p;v_j^new).             (1)

All histories and increments are real in this note. Velocities range in
a small fixed neighborhood of a nonzero value with fixed sign. Take that
neighborhood and the raw z domain sufficiently small that

    r_(v_j)(p)/r_(v_j^new)(p) < 4/3                     (2)

uniformly for p nonzero. There is a fixed strict margin in (2). On the
support of chi_(j-1)(p;v_j^new), its norm is at most (3/4)s_j. Equation
(2) therefore puts that support in the plateau of chi_j(p;v_j). The
nesting needed for quadratic absorption is exact, with the original
dyadic scale ratio 2; it does not require changing the RG scale convention.

Consequently the exact normalized shell identity is

    G_j(p)=(1+z_(0,j))f_j(p)
             /[D_j(p)+chi_j(p;v_j)L_j(p)].               (3)

The proof is the same subtraction of low covariances as in
MARKED_COVARIANCE, now using the just-proved nesting. Both cutoffs are
one on a uniform inner ball and zero outside a uniform outer ball, so
f_j and all its history derivatives have support in an annulus comparable
to s_j. The next velocity in (1) is a function of the current raw data;
no as-yet-unsolved lower-scale velocity coordinate enters (3).

## 2. Two real cutoff derivatives preserve the shell dimension

For real v bounded away from zero,

    partial_v r_v=v p_1^2/r_v,
    partial_v^2 r_v=p_0^2 p_1^2/r_v^3.                  (4)

On the support of a derivative of chi_j, r_v is comparable to s_j. The
chain rule gives

    partial_v chi_j=chi'(r_v/s_j) (partial_v r_v)/s_j,
    partial_v^2 chi_j=chi''(r_v/s_j)(partial_v r_v/s_j)^2
                         +chi'(r_v/s_j)partial_v^2r_v/s_j.

Both are uniformly bounded: a velocity derivative contributes no inverse
shell scale. More generally, for m=0,1,2 and every fixed momentum multiindex,

    |partial_p^alpha partial_v^m chi_j(p;v)|
                         <= C_alpha s_j^(-|alpha|).      (5)

For m=0 and alpha=0 this includes the plateau; derivatives have annular
support. Smoothness at p=0 follows because the cutoff is constant there.
The rational map v_j^new has bounded first and second derivatives on the
fixed raw domain, so (5) also applies to its composition in f_j.

Thus all new derivatives of the numerator and denominator in (3) are
local to scale j and obey exactly the same symbol dimensions as a local
history mark in MARKED_COVARIANCE. Higher-scale normalization marks still
enter only through kappa_j R_0 and obey the previous curvature estimate.
For one scale-specific raw mark and one common real history mark, the
bounds are therefore the same as (8)--(10) of that note, with enlarged
constants. In particular a higher-scale curvature mark is bounded by

    C s_j^(-1-|alpha|) delta 2^(b j)2^[-a(k-j)]

whenever a+b<1-sigma; direct marks are local. Fourier decay, localization
moments, Gram representations, and the cofactor determinant estimates
then follow verbatim from their symbol proofs. The degree-uniform ratio
bound (16) is unchanged, since cutoffs do not occur in the product Z_p/Z_q.

The statement is real C2, uniformly in the number of scales. It is enough
for the real mean-value estimates, contraction and Brouwer steps of the
current flow construction. It does not analytically continue chi(v).

## 3. Replacing common-history Cauchy estimates

The common-history derivative can be estimated by two explicit marks.
At most two interaction-amplitude marks remove at most two factors of eta.
A velocity, cutoff, covariance, or normalization mark removes none. In
particular an initial irrelevant endpoint with budget eta^4 still has at
least eta^2 after two amplitude marks. Fixed-degree field multiplicities
are the same N2/N4 and field-power factors already identified in the
marked covariance and degree-reserve arguments.

This substitution is important: a real C2 cutoff alone would not justify
the complex Cauchy argument. Instead, use the explicit second-derivative
normalization formula. With z=Phi(x,z), A=D_z Phi, P_l=partial_(x_l)Phi,
and a_l=D_(x_l)z,

    (I-A)a_l=P_l,
    (I-A)D_u a_l=D_u^tot P_l+(D_u^tot A)a_l.             (6)

The real C2 alternative in NORMALIZATION_RESOLVENT proves the needed bounds
by this identity and the weighted inverse. The raw two-mark tree count
provides D_u^tot P_l=O(W_a(h,l)) and D_u^tot A_hk=O(eta W_a(h,k)); the
additional common direction D_u z is O(eta). Equation (6) therefore gives
the required bounded second mark with short memory, without a factor equal
to the number of scales. This is the same endpoint/localization tree task
as before, now using real derivatives rather than complex interpolation.

## 4. What the primary reference statements actually say

[Mastropietro--Porta, (8.60)--(8.61), Theorem 8.5](https://arxiv.org/html/2108.13708v3#S8.SS2)
expresses the beta functions as functions of quartic histories after the
dependent velocities and normalizations have been determined. Its frozen
quartic input includes an explicit decaying velocity-beta estimate. The
proof distinguishes the bounded actual velocity flow in (9.108) from the
frozen-quartic velocity-beta conclusion (9.109). Thus the latter is not
restricted merely to the actual running quartic trajectory. It still does
not mean that the full velocity history has been frozen independently.

For the moving-cutoff convention, [Porta--Scola--Singh, Section 5.4,
the infrared paragraph and Proposition 5.3](https://arxiv.org/html/2509.08665v1#S5.SS4)
explicitly uses the lattice infrared procedure for its reference model,
and states uniformly bounded reference flows. This is a published anchor
for that actual-flow convention. It is stronger evidence than an unsupported
assertion that every cutoff change is an ultraviolet endpoint. Its printed
Proposition 5.3 is an actual-flow statement; obtaining fully frozen beta
bounds still requires the trajectory conversion or the usual
bounded-flow-to-frozen-beta argument in the same coordinates.

The finite-root or limiting-velocity conversion in CONSTANT_HISTORY_COMPARISON
is valid once the frozen-quartic estimate is available in this same
convention. Its real velocity matching uses continuity and Brouwer, and
its derivative step may use real C2 interpolation with half the decay
exponent. Complex velocity Ward identities are unnecessary. Reference UV
endpoint differences retain their initial-scale gain and need not have
the smaller physical residual amplitude.

We do not assert that the numerical anomaly coefficient is independent
of cutoff shape. The beta-bound arguments only require appropriate bounded
Ward correction coefficients, but extending that proof to another shape
must actually be checked if one takes that route. The moving-cutoff real-C2
lemmas above remove the need to impose an unproved complex continuation;
they do not by themselves replace the reference comparison argument.

## 5. A route that does not need the frozen-history transfer

The [actual-reference shadowing construction](../vacancy/REFERENCE_SHADOWING.md)
uses the published bounded reference trajectory directly, in the same
moving-cutoff convention. It seeks a physical trajectory whose difference
from it vanishes toward the infrared. Its forcing is the initial-endpoint
and curvature difference; history derivatives are the real marked estimates
above. Therefore the ordinary flow construction need not wait for the
fully frozen beta comparison in Section 4.

For bare matching it is useful to keep the component sizes distinct.
Velocity and relevant forcing can be first order in the reference quartic
size t, whereas quartic forcing has no one-quartic contribution: one R4
endpoint with no contraction still localizes to zero, a quadratic endpoint
cannot produce four external legs, and a reference endpoint with six or
more legs needs at least two bare quartics. Together with the physical
residual budget this gives quartic forcing O(t^2+epsilon^2). Even after
using the looser eta=sqrt(epsilon) derivative bounds in the shadow equation,
the initial quartic matching error is O(epsilon^(3/2)), while the velocity
error is O(epsilon). An anisotropic Brouwer box with these two radii is
therefore available. A crude O(t) quartic error in an isotropic box would
not by itself establish a near-identity self-map.

The reference ultraviolet map has quartic correction O(t^2), while its
velocity correction is O(t). Finite-scale parameter continuity follows
from the smooth real cutoff formulas and the uniformly convergent GN
expansion. One can first use finite ultraviolet and infrared cutoffs and
then take a diagonal limit under the uniform estimates. No derivative
of the bare-to-infrared velocity map, and no complex velocity Ward identity,
is needed for this matching.
