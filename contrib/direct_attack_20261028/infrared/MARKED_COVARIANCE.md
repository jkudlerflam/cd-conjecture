# Normalized covariances with history marks

Status: **[proved covariance, determinant, and normalization-product lemmas]**.
The formulas below retain the anomalous normalization. They establish the
raw covariance input to [the normalization resolvent](../vacancy/NORMALIZATION_RESOLVENT.md),
with independent small normalization increments. The remaining GN endpoint
and localization accounting is distinguished explicitly in Section 8.

Author: Codex `infrared_theorem`, 2026-10-28. Sections 1--6 were independently
checked by the coordinating Codex agent and the separate `alternative_attack`
agent, with the uniform-degree refinement below incorporated. This is an
original derivation from the normalization identities cited below, not a
theorem quoted from those papers. No human review.

## 1. Choice of cutoff and exact shell identity

Work at fixed physical Fermi centers, with a fixed finite branch set. It is
enough to discuss one branch. Its reference velocity v_* is real and nonzero.
Take smooth even cutoffs chi_j(p)=chi_0(2^(-j)p), j<=0, independent of all
running histories. Arrange that chi_j=1 on supp chi_(j-1), and write
f_j=chi_j-chi_(j-1). Its support is an annulus |p| comparable to s_j=delta 2^j.
The fixed positive delta is small enough to stay within the branch patch.

Fixed cutoffs matter: a smooth compactly supported cutoff evaluated at a
complex running velocity has no specified holomorphic continuation. The
fixed-cutoff scheme here has analytic history dependence directly in its
denominator. A change from the running-cutoff convention is a change of RG
coordinates, not an assertion that a smooth cutoff is holomorphic.

Let R_0(p) be the fixed bare curvature remainder. On the patch assume, for
every finite derivative order used below,

    |partial_p^alpha R_0(p)| <= C_alpha |p|^(2-|alpha|).       (1)

For |alpha|>2 this is just a convenient loose bound. In the linear reference
model R_0=0. Initial interaction-generated irrelevant quadratic kernels stay
among the interaction endpoints, so they need not be absorbed into R_0.

Temporarily let v_j be independent velocities close to v_*, and let
z_j=(z_(0,j),z_(1,j)) be independent time- and space-derivative increments.
Define

    Z_(j-1)=Z_j(1+z_(0,j)),       Z_0=1,
    kappa_j=Z_0/Z_j=product_(r=j+1)^0 (1+z_(0,r))^(-1),
    D_j(p)=-i p_0+v_j p_1+kappa_j R_0(p),
    L_j(p)=-i z_(0,j) p_0+z_(1,j) p_1.                    (2)

Before absorbing the local quadratic term, the low covariance is
chi_j/(Z_j D_j). After absorption it is
chi_j/[Z_j(D_j+chi_j L_j)]. Set

    v_j^new=(v_j+z_(1,j))/(1+z_(0,j)),
    D_(j-1)^new=-ip_0+v_j^new p_1
                        +[kappa_j/(1+z_(0,j))]R_0(p).

Then D_j+L_j=(1+z_(0,j))D_(j-1)^new, so subtracting the next low
covariance gives the exact normalized shell

    G_j(p):=Z_(j-1) g^(j)(p)
          =(1+z_(0,j)) f_j(p)/[D_j(p)+chi_j(p)L_j(p)].      (3)

Indeed chi_j=1 wherever chi_(j-1) is nonzero, so the same denominator can
be used when taking the difference. This proves (3) also in the transition
region of chi_j; it does not replace a momentum-dependent Z by a constant.
On a physical trajectory v_(j-1)=v_j^new. In the raw history construction
each shell uses only its current velocity and derivative increments; both
components of z are eliminated by NORMALIZATION_RESOLVENT, and the velocity
flow enforces this last identity. This avoids dependence on an unspecified
lower-scale velocity and makes the causal history indexing explicit. No
cumulative sum of arbitrary velocity increments is assumed small.

The identities being implemented are the quadratic absorption and field
rescaling in [Mastropietro--Porta, (8.35)--(8.42), (8.49)](https://arxiv.org/html/2108.13708v3#S8.SS2),
and [Porta--Scola--Singh, (4.75), (4.85), (4.89)](https://arxiv.org/html/2509.08665v1#S4.SS3).
Equation (3) also displays the accumulated curvature dependence that would
be missed by treating the normalized shell as a function only of v_j.

## 2. A common complex domain and symbol bounds

Choose a fixed r_0>0 sufficiently small, and restrict

    |z_j|<r_0,       |v_j-v_*|<r_v,
    sigma=-log_2(1-r_0)<1/8,                              (4)

where r_v is a sufficiently small fixed fraction of min(1,|v_*|). All
histories have a slightly larger domain of the same kind available for
Cauchy estimates. The radius r_0 is independent of the interaction size
eta. The solved normalization increments will later be O(eta^2).

For j<=0, |kappa_j|<=(1-r_0)^j=2^(-sigma j). Hence on the j-shell,

    |kappa_j R_0(p)| <= C delta 2^[(1-sigma)j] |p|.        (5)

The linear part of the denominator of (3) is

    -i(1+chi_j z_(0,j))p_0+[v_j+chi_j z_(1,j)]p_1.

It is within a fixed small fraction of -ip_0+v_*p_1. Choose delta small
enough that (5) uses another such fraction. The denominator is then bounded
below in modulus by c|p|, uniformly in every complex history in (4), every
number of scales, and j. Differentiation in real momentum and the inverse
rule give

    sup |partial_p^alpha G_j(p)| <= C_alpha s_j^(-1-|alpha|).
                                                                    (6)

This is a direct ellipticity estimate. It does not rely on positivity of
complex propagators, convergence of Z_j, or a Ward identity. Constants may
depend on the fixed branch velocities, patch, and finite derivative order.

## 3. One coordinate mark and one common mark

Let partial_(v_k) or partial_(z_(mu,k)), mu=0,1, be a raw history
coordinate derivative. A bound written partial_(z_k) applies to either
component, or their fixed finite-dimensional norm.
Let U be a common raw direction,

    U=sum_k [a_k partial_(v_k)
                    +sum_(mu=0,1) b_(mu,k) partial_(z_(mu,k))],
    sup_k (|a_k|+sum_mu |b_(mu,k)|)<=1.                   (7)

At a finite cutoff this is an ordinary finite-dimensional direction.
The following bounds have constants independent of its length. Define

    mu_j=C delta 2^[(1-sigma)j].

For m=0 or 1 the marked symbol bounds are

    |partial_p^alpha U^m partial_(v_k) G_j|
       <=C_alpha s_j^(-1-|alpha|) 1_{k=j},

    |partial_p^alpha partial_(z_k) G_j|
       <=C_alpha s_j^(-1-|alpha|)
            [1_{k=j}+1_{k>j} mu_j],

    |partial_p^alpha U partial_(z_k) G_j|
       <=C_alpha s_j^(-1-|alpha|)
            [1_{k=j}+1_{k>j} mu_j(1+|j|)].               (8)

Coordinates with k>0 are absent; z coordinates with k<j do not affect G_j.
The unmarked and once-common-marked symbols also obey (6), after enlarging
its constant. The same statements hold for any fixed finite number of
branch components.

**Proof.** Velocity dependence is only through v_j. The spatial increment
z_(1,k) occurs only when k=j. For the time increment at k>j,

    partial_(z_(0,k)) kappa_j=-kappa_j/(1+z_(0,k)),
    U kappa_j=-kappa_j sum_(r=j+1)^0 b_(0,r)/(1+z_(0,r)). (9)

Consequently the first derivative costs a constant times |kappa_j| and
the mixed derivative a constant times |kappa_j|(1+|j|). These factors
multiply R_0 in the denominator, so (5) supplies mu_j. The local z_j
derivative costs a bounded factor directly. A common derivative of that
term also contains at worst mu_j |j|, bounded uniformly in j. The same
observation handles a common derivative after a velocity derivative.
Differentiate the reciprocal denominator using its lower bound. Momentum
derivatives of chi_j, f_j and R_0 use s_j^(-1) per derivative. This proves
all of (8). QED.

In particular, for any a,b>0 with a+b<1-sigma and j<k<=0,

    mu_j(1+|j|) <= C_(a,b) delta 2^(b j) 2^[-a(k-j)].     (10)

To verify it, divide by the right-side exponent. The remaining scale factor
is (1+|j|)2^[(1-sigma-a-b)j+a k], uniformly bounded for j<=k<=0.
Thus a higher-history curvature mark has both short memory and a remaining
ultraviolet gain. The initial curvature need not be proportional to eta;
these terms belong to the UV alternative in GN_MARK_CALCULUS.

Two arbitrary coordinate marks can likewise be obtained directly from
(9), with the corresponding intersection indicators. The mixed common
mark in (8) is the stronger useful statement for the desired history
coordinates, since it has already summed over the second mark's scale.

## 4. Decay, moments, and Gram representations

Let M_j be the appropriate dimensionless bracket in (8), or 1 in (6).
For each marked kernel K_j and fixed nonnegative integers q,ell,

    |partial_x^ell K_j(x)|
      <=C_N M_j s_j^(1+ell)/(1+s_j|x|)^N,

    integral dx |x|^q |partial_x^ell K_j(x)|
      <=C_(q,ell) M_j s_j^(ell-1-q).                     (11)

The spatial integral is replaced by a lattice sum when appropriate.
Fourier integration by parts proves the first inequality from (6)/(8);
choose N>q+2 and rescale to prove the second. Difference derivatives in
the lattice coordinate have the same scaling. Thus the marks do not spend
an additional momentum dimension or invalidate localization moments.

There is also a direct Gram representation. For any scalar Fourier
multiplier A_j(p) in these bounds, split it into modulus and phase in
L2(dp): one vector has factor |A_j|^(1/2), and the other the same factor
times its phase. Plane-wave factors supply the two positions. Both squared
norms equal integral |A_j|, and

    integral |A_j(p)| dp <= C M_j s_j.                   (12)

For branches take the direct sum of these spaces. Tensoring with the
unit vectors representing the positive BBF interpolation matrix preserves
the same estimate. Hence each marked covariance itself has a Gram bound
C M_j s_j, not a loss proportional to the number of contracted fields.

At finite volume use the Fourier sum. Bounds (11)--(12) are uniform on
shells with beta s_j and L s_j bounded below by a positive constant; the
periodized real-space estimate uses the torus distance. This includes
the thermodynamic order L to infinity at fixed beta, followed by beta to
infinity, and arbitrary fixed shell cutoffs before removing volume.
No assertion is made for L much smaller than beta on shells below 1/L:
a spatial zero mode can then defeat an unmodified O(s_j) Gram bound.
At the last thermal shell, the Matsubara gap supplies the required annular
lower momentum bound; its low cutoff may replace f_j with no change in
the argument. This cutoff qualification is independent of history marks.

## 5. Differentiating a determinant without losing its cancellations

Let G be an N by N BBF Gram matrix built from an unmarked j-shell kernel.
Every minor has the same Gram bound, so |det G|<=(Cs_j)^N. Suppose a mark D
has entrywise bound C M_D s_j, and two marks have the analogous bounds.
The cofactor identity gives

    |D det G| <= C N^2 M_D (Cs_j)^N,

    |D_1 D_2 det G|
      <=C [N^2 M_(12)+N^4 M_1 M_2](Cs_j)^N.              (13)

For N=0 the derivatives vanish. This uses minors, not an absolute expansion
into N! permutations. Thus one or two determinant marks introduce only
fixed-degree polynomials in the number of fields. Their treatment uses
one fixed field-weight reserve, as in DEGREE_REPARAMETERIZATION; there is
no fresh field-radius loss at each scale. The same calculation applies to
a coordinate mark followed by the common direction U of (7).

## 6. Products of anomalous normalization factors

For p<q set R_(p,q)=Z_p/Z_q and ell=q-p. Directly from (2),

    R_(p,q)=product_(r=p+1)^q(1+z_(0,r)),
    |R_(p,q)|+|R_(p,q)^(-1)| <=2 exp(C r_0 ell),
    partial_(z_(0,k)) log R_(p,q)
                          =1_{p<k<=q}/(1+z_(0,k)).       (14)

The logarithm is the sum of principal logarithms on |z_(0,r)|<r_0. The
z_1 derivatives vanish. For U as
in (7), |U log R|<=C ell, while
|U partial_(z_k) log R|<=C 1_{p<k<=q}. Thus a coordinate mark on R or its
inverse costs C times the path indicator; a mixed coordinate/common mark
costs C(1+ell) times that indicator. Two common marks cost C(1+ell)^2.
For a real or half-integer field power m of R, its unmarked bound is instead
exp(C r_0 |m| ell), and the derivative factors are polynomials of degree at
most two in |m|. Analytic branches are fixed by the above logarithm. The
dependence on m cannot be dropped for arbitrary-degree endpoints.

If an unmarked renormalized path has gain 2^(-d ell), choose r_0 so its
normalization growth leaves d'>0, and reserve any a<d'. For r<=2,

    (1+ell)^r exp(C r_0 ell) 2^(-d ell)
        <=C_(r,d,a) 2^(-a ell).                         (15)

This proves the required polynomial-loss absorption for the ACTUAL
normalization products, rather than treating Z as constant. A coordinate
z_k mark on an edge p to q is at scale k on that edge, so the root-to-k
part of the reserved path supplies its short-memory weight. If it hits
the curvature coefficient of a shell below k, (10) supplies that missing
part of the path.

Here is the uniform-degree version of (15). At an internal source-free
cluster with P>=2 external fermion fields, the canonical positive
renormalized dimension is

    d(P)=1 for P=2,4;       d(P)=P/2-2 for P>=6.

In particular d(P)>=1 and d(P)>=P/6. Products over branches of the field
normalization have total absolute field power m<=P/2. Choose the same
r_0 for every degree so that C r_0<=(log 2)/12. Then

    C r_0 m/(log 2) <= P/24 <= d(P)/4.

For any fixed 0<a<1/4, both the unmarked and at-most-two-mark factors are
therefore bounded by the following common majorant:

    (1+m)^2(1+ell)^2 exp(C r_0 m ell) 2^[-d(P)ell]
      <= C_a (1+m)^2 2^(-a ell) 2^[-d(P)ell/2].          (16)

To prove it, the quotient of the left side by the exponential part of the
right side is at most (1+ell)^2 2^[-(d(P)/4-a)ell], uniformly bounded since
d(P)>=1. Thus a fixed short-memory reserve and half the degree-dependent
dimensional gain remain, simultaneously for all P. The polynomial (1+m)^2
is present only at the at most two marked locations and is paid by the
single field-weight reserve. It is not charged at every scale. Constants
in C account for the fixed finite branch set, so they are chosen before
r_0. This removes any degree-dependent smallness requirement on r_0.

After the normalization equations have been solved, the stronger bounds
from NORMALIZATION_RESOLVENT can be substituted into (14). For example,
if |D_(x_k)z_r|<=C eta 2^[-a(k-r)] for r<=k, then

    |D_(x_k) log R_(p,q)|
      <=C eta 2^[-a(k-q)_+],

with zero contribution if k<=p. A common-history mark of log R costs at
most C eta ell, and its mixed coordinate derivative is bounded by
C 2^[-a(k-q)_+]. The mixed derivative of R itself is therefore bounded
by its unmarked majorant times C(1+eta^2 ell)2^[-a(k-q)_+]. Again (15)
absorbs the polynomial loss; the accumulated Z itself need not converge.

## 7. What this supplies to the raw normalization equation

Use the raw domain (4), not a radius O(eta^2), while constructing the trees.
The denominator bound, (11), and (13) give the same covariance dimensions
as in the usual unmarked GN estimate. A mark on a covariance or a ratio
removes no interaction endpoint. Its coordinate scale is handled by the
path argument following (15), or by the curvature gain (10). Consequently
the covariance/normalization part of the raw tree majorant introduces no
new smallness assumption beyond those already used for unmarked trees.

For the wave-function increment Phi_h, a one-local-quartic endpoint gives
a momentum-independent tadpole; a one-local-relevant endpoint is also
constant in momentum. Their localized momentum derivatives vanish exactly,
including after any raw covariance/ratio mark. All other pure-local
contributions have at least two interaction endpoints. A raw z mark leaves
their O(eta^2) factor intact; an interaction-amplitude mark can remove one
eta. An initial irrelevant endpoint instead retains its eta^4 input budget.
Thus, after the ordinary endpoint/localization tree summations, these
calculations give the raw estimates in NORMALIZATION_RESOLVENT with
partial_z Phi=O(eta^2) W_a and partial_x Phi=O(eta) W_a. No marked bound
for the already-solved Z history has been assumed in reaching this point.

For quartic and relevant betas, differentiating a carried field normalization
can remove an independent z while leaving a local amplitude of size eta.
For the velocity beta, however,

    B_v=(z_1-v z_0)/(1+z_0)

has an O(1) raw z derivative. This is necessary and harmless: the updated
NORMALIZATION_RESOLVENT permits an O(1) beta z-Jacobian and gains the eta
factor from the solved history derivative Dz(x). A direct velocity mark
of B_v, with z fixed, is O(eta^2) on that solved z ball. The lattice
tadpole's linear amplitude term has a curvature gain and belongs to the
explicit UV alternative.

## 8. Scope and the remaining GN work

The lemmas in Sections 1--6 are explicit proofs of the covariance and ratio
claims needed for the marked-tree bridge. Initial retardation and arbitrary
field degree do not change these denominators. The proof is uniform in the
number of shells on the stated complex domain, and fixed-cutoff dependence
is analytic in every history coordinate. The auxiliary normalization
fixed point is handled by NORMALIZATION_RESOLVENT.

There remains a logically separate task when admitting arbitrary initial
endpoints: localization zero routing, moments, and field-label sums must
be bounded in one fixed endpoint norm through the complete GN forest.
Equation (13) identifies precisely the additional degree polynomials that
one and two marks contribute; it does not conceal that endpoint-norm task.
The stronger tree-length endpoint norm supplied by the massive matching
geometry is intended to pay it. These covariance lemmas alone do not prove
the full replacement-endpoint theorem, pressure asymptotics, or CFT
identification.

Primary-source context: the unmarked dimensions and determinant bound are
in [Porta--Scola--Singh, (4.111), (4.116)--(4.117)](https://arxiv.org/html/2509.08665v1#S4.SS5).
The present calculation adds the explicitly differentiated normalized
shells and products on the independent-history domain; it does not cite
those unmarked statements as if they already contained this conclusion.
