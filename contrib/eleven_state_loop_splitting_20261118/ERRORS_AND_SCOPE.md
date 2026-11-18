# Splitting scale, microscopic errors and physical scope

## 1. Three different momentum margins

Fix n and finite q,w,nu,D with simple tree roots. Define delta_2 from
distinct Fermi points and delta_nc from distinct unordered pair sums,
excluding comparisons between two opposite-point pairs. These are
already strictly positive in the previous tree construction. Neither
is the Cooper-transfer margin, which is zero at t_d=0.

For a member with distinct exact coefficients c_j from DERIVATION.md
(5), set kappa=min_(i<j)|c_i-c_j|>0. The analytic result supplies such
members for n=2,3,4. For n=1 the comparison is empty. For general n
the sufficient moment condition is explicit but remains unproved.

Let sigma_j(t)=k_j,+(t)+k_j,-(t). These functions are odd. On a
sufficiently small interval take finite bounds L for all root first
derivatives and R for |sigma_j'''|/6. Taylor's theorem gives

    |sigma_j(t)-c_j t| <= R |t|^3.

Choose t nonzero and small enough that

    4L|t| < min(delta_2,delta_nc)/2,
    2R t^2 < kappa/2.

Then all simple roots persist and the new Cooper margin is at least
kappa|t|/2. The other margins stay at least half their tree values.
The full conservative pair margin, including repeated-point pairs,
therefore has the lower bound

    delta=min(delta_2/2,delta_nc/2,kappa|t|/2)>0.      (1)

This is a finite-parameter existence prescription. The derivative and
remainder bounds can be enclosed from the finite algebraic resolvent
and root equations. The saved high-precision example is a validation,
not an interval certificate for those constants.

With eta=t D^p, the first new margin is of order

    |eta| min_(i<j)|A_i-A_j|/D,

at fixed n,w,q,nu, before the eta^3 remainder. D must be chosen large
enough for distinct finite-D derivatives, then eta small enough for
the linear term to dominate. Taking D to infinity at fixed eta
first is not this schedule. In particular the finite express delay
must stay in the reference symbol; it cannot be discarded as an error
larger than the margin it creates.

## 2. Convert momentum separation into an operator tolerance

Normalize energies by W=lambda^2 g^2 zeta D^(-2p-2), not by the
microscopic unit hop. For a fixed finite reference symbol h(z), let
P(z)=z h(z) and z_j=e^(ik_j). Choose
r<min(1/16,delta/100) and disjoint disks |z-z_j|<r. The minimum m of
|det P(z)| on their boundary circles and on the unit-circle complement
is strictly positive. Let M bound ||P|| on that compact contour.
An analytic Hermitian Laurent perturbation E(z) with annulus norm

    ||E|| < tau,
    tau=min(1/2, m/[8n(M+1)^(n-1)])                  (2)

is sufficient since |z|<2 on the contours.
Multilinearity bounds the determinant change by
n ||zE|| (M+1)^(n-1). Rouche's theorem preserves one simple zero in
each disk, and the complement bound excludes new unit-circle zeros.
Hermitian reciprocal symmetry keeps each isolated zero on the circle.
Angular movement is at most 2r; four-root comparisons change by at
most 8r. Thus every momentum margin remains positive. The construction
also yields a positive speed bound from the determinant derivative.
No uniform-in-n or moderate-sized tau is asserted.

## 3. Conversion, grammar and port errors use this new tolerance

All following choices come **after** the finite loop and its margin.
Tying the loop to g^2 in advance can make its margin shrink at the
same rate as the error and defeats a norm-based argument.

**Finite g.** The source norm and conversion bound inherited from C84
give a two-cell remainder at order lambda^2 bounded by

    lambda^2 zeta g^4/4,   |g|<=1/2.

The modified carrier gap is no smaller than the gap used for this
conservative estimate. A sufficient extraction factor for the scalar,
quadratic and remaining Fock operator norms is C_ext(n)=64n^2.
For example, zero/one-particle matrix elements bound each extracted
quadratic coefficient by at most twice the input norm; summing its
at most 4n^2 entries, the scalar and the remaining operator gives
this conservative factor. The annulus introduces only a factor two
on each of the two nearest-cell hopping directions.
Thus a safe normalized finite-g bound is

    E_g <= 16n^2 D^(2p+2) g^2.                       (3)

Choosing E_g<tau/3 controls the finite-g quadratic projection on the
annulus. The nonquadratic part is retained and has a comparable local
bound. Further smallness relative to the physical velocities and
mismatch logarithms may be required by an interacting matching theorem;
the finite matrix calculation does not provide that theorem.

**Finite K.** The biased legal code is not exactly invariant. Its
lambda=0 physical band splitting is already O(B_p^2/K), with

    B_p=(p+1)(3+D+|mu_total|+|h_L|+24).

For a fixed one- or two-cell cluster the Schur estimates also bound
source-map and coefficient errors. Collect their physical constants
into a finite C_K, independent of sufficiently large K. A conservative
normalized condition is

    E_K <= C_K/(K W) < tau/3.                        (4)

The zeroth-port splitting is not multiplied by lambda^2. Choose K
last, and also impose the previous zero-port selector inequality
q^2/(K-4q)<E_(p+1)/[2(p+2)] and the auxiliary-gap requirements. Neither
the old C84 K schedule nor a finite-cell error bound can be silently
assigned to this modified band theorem.

**Higher port orders.** On fixed open cell clusters, port-parity
symmetry removes odd orders. A contour expansion, including the
order-lambda^2 physical field in the bookkeeping, gives a physical
remainder C_lambda lambda^4 with finite cluster-dependent C_lambda.
Consequently a sufficient condition is

    E_port <= C_lambda lambda^2 D^(2p+2)/(g^2 zeta)
              < tau/3.                              (5)

Choose lambda after g and then K after lambda. These conditions can
all be met at each fixed certified band member, but can be extremely
demanding. They are finite-cluster bounds. Volume-independent bounds
for the thermodynamic chain remain an additional matching assumption.

## 4. Actual term and parameter count

The eleven letters are L0,L1,R0,R1,DA,DB,XA,XB,E,MA,MB. There is no
vacancy and no second express bus. The two restored physical rules are

    D_A E <-> L0 M_A, amplitude -t_d exp(i pi),
    D_B E <-> L0 M_B, amplitude -t_d exp(i pi/2).

They share one real magnitude; no new phase is adjusted. Relative to
the modified tree this adds one independent real coefficient and two
undirected transitions, or four directed complex matrix units. In the
previous contribution's explicit matrix-unit convention, the generic
bond has 24 undirected transitions, 48 directed off-diagonal entries
and 115 diagonal entries. The resulting 163 entries are not 163 free
couplings. With phases fixed and the common positive hop scale removed,
the full loop rule family has nine independent real ratios. Constants
and one-site telescoping terms were removed in the rank audit in
[SIMPLICITY.md](../minimal_d_unbounded_audit_20261117/SIMPLICITY.md).

The numerical derivative check uses n=2, D=10^4 and eta=10^-6, so
t_d=10^-18. Already D/t_d=10^22 in that physical rule basis, before
specifying g,lambda,K or the much smaller field. No moderate coefficient
ratio, wide robust region or practical many-body simulation scale is
claimed. The emergent W and the extra 1/D centre suppression are
reported separately from these microscopic coefficients.

## 5. What the band result does and does not establish

The result is an analytic nonresonant leading-band construction for
n=1,...,4 on eleven physical states, plus an all-n exact derivative
formula and separation of two extreme pair centres. The unproved
moment nonvanishing criterion prevents an arbitrary-n conclusion.
The original unmodified n=11 and regenerated n=13 certificates are
separate finite members and are not invalidated or extended by this work.

Even if the all-n band criterion is proved later, a microscopic result
would still require nonzero-port ground-sector selection against dressed
defects, volume-uniform elimination/matching, interacting multichannel
stability, and identification of the conformal observables. Any claimed
c would mean total nonchiral Virasoro charge in a unitary compact CFT.
No central charge is measured or proved here. No mobile lift is built.
