# Exceptional momentum support is incompatible with connected clustering

Status: **[proved distribution theorem and conditional scaling bridge]**.
This theorem requires a comparison on the full nonexceptional momentum
region. The higher-point comparison currently stated by Porta--Scola--Singh
uses a thinner frequency slice, so the theorem is not an automatic corollary
of that paper. A separate chiral Ward-identity test can exploit the slice;
that is a different argument.

Provenance: the coordinating Codex agent proposed the averaging proof;
the frozen-sector agent checked it and derived the stronger weak-clustering
version below, 2026-10-28. AI assistance only, no human review.

## 1. Reduced momenta and the precise cluster hypothesis

For an m-point translation-invariant distribution remove the overall
momentum delta and write

    E=(R^d)^(m-1),    p_m=-sum_(i=1)^(m-1) p_i.

For every nonempty proper subset I of the m insertions define

    q_I(p)=sum_(i in I) p_i,
    Z_I={p:q_I(p)=0},
    Z=union_I Z_I.

Complementary subsets give the same exceptional subspace, so only one
representative of each pair is needed. Each Z_I is a proper linear subspace.
Let D be a tempered distribution on E. Choose one fixed nonzero spatial
translation direction. Denote the corresponding component of q_I by
`ell_I(p)`. The weak connected-cluster condition used below is

    <D, exp(i t ell_I(p)) psi(p)> -> 0 as |t| -> infinity  (WC)

for every Schwartz test function psi and every nonempty proper I.
This is the distributional statement that a properly smeared connected
correlator vanishes when that subset of insertions is translated far away
in the chosen spatial direction. It does not assert pointwise behavior
of an unsmeared singular correlation function.

For an m>=3 Gaussian theory its connected m-point distribution is zero.
Thus, when comparing a physical connected distribution with the Gaussian
one, D is just that physical connected distribution. For a difference of
two clustering connected distributions, (WC) holds by linearity as well.

## 2. The theorem: no exceptional-supported cluster distribution

**[proved]** If `supp D subset Z` and (WC) holds, then `D=0`.
No positivity, homogeneity or conformal assumption is needed for this
statement. In particular, no additional globally bounded-smearing
assumption is required.

### A distribution on one hyperplane

First suppose T has support in `ell(p)=0` for a nonzero linear form ell,
and obeys weak decay under multiplication by `exp(i t ell)`.
Take a compactly supported smooth test psi. The distribution has finite
order N on a compact neighborhood of that support. In linear coordinates
with y=ell(p), it annihilates every test of the form `y^(N+1) g` there.
For clarity, insert a cutoff chi(y/epsilon) equal to one near zero. The
support property permits this insertion, and the order-N estimate of
`chi(y/epsilon)y^(N+1)g` is O(epsilon). Letting epsilon tend to zero proves
the assertion without a structural decomposition theorem.

Taylor expansion therefore gives the exact polynomial identity

    <T, exp(i t ell) psi>
       =sum_(r=0)^N (it)^r/r! <T, ell^r psi>.            (1)

A polynomial tending to zero at infinity is identically zero. At t=0,
(1) gives <T,psi>=0. Since psi was arbitrary, T=0.

Multiplication by a fixed smooth compactly supported Fourier cutoff
preserves (WC): replace psi by that cutoff times psi. This permits the
same argument locally.

### A finite union of subspaces

Remove any subspaces contained in others from the finite collection.
At a point of one remaining subspace L which belongs to no other member,
choose a neighborhood disjoint from the other members. A smooth cutoff
there localizes D to L. Since L is contained in some original Z_I, it is
also contained in the hyperplane `ell_I=0`; the previous argument kills
the localized distribution.

Thus its surviving support lies in the union of pairwise intersections
of the former collection. Remove contained members again and repeat.
Distinct inclusion-maximal subspaces have intersections of strictly
smaller dimension, so this procedure terminates after finitely many
steps. If only one subspace remains, the one-hyperplane argument kills
it directly, including a possible distribution concentrated at the
origin. This proves the theorem. QED.

This proof handles cancellations between different exceptional components
and derivatives of delta functions. Neither can evade all proper-cluster
limits.

## 3. A direct averaging proof under a stronger, familiar hypothesis

The coordinator's argument gives a complementary proof when every
Schwartz smoothing f of the inverse Fourier transform is bounded,
continuous and satisfies ordinary pointwise proper-cluster limits.
Choose a Schwartz function phi with Fourier transform identically one
near zero, and define

    A_R^I f(x)=integral phi(y) f(x+R translation_I(y)) dy.

One may average only in the fixed spatial direction or in full spacetime.
Its Fourier multiplier is `hat(phi)(R ell_I)` in the first case. Since
`1-hat(phi)(R ell_I)` vanishes on a neighborhood of Z_I,

    product_I (1-A_(R_I)^I) f=0                          (2)

for all independently chosen positive R_I. Smooth multipliers vanishing
near the support also annihilate normal delta derivatives, so (2) is
an exact distribution identity, not a limiting assertion.

For fixed x, boundedness and clustering give `A_R^I f(x)->0` by dominated
convergence after the change of variables displayed above. The convergence
is uniformly bounded by `||phi||_1 ||f||_infinity`. Send the finitely many
R_I to infinity one at a time in (2). Every other averaging operator
has a fixed integrable kernel at that step, so dominated convergence
passes through it. After the last limit, f=0. All Schwartz smoothings
vanish, hence the original distribution vanishes.

This proof explains why a rate of clustering, joint uniformity in the R_I,
or a decomposition into individual exceptional terms is unnecessary.
The theorem in section 2 obtains the same conclusion directly from the
weaker Schwartz-distribution condition (WC).

## 4. Conditional bridge from a nonexceptional scaling comparison

Suppose normalized connected scaling distributions D_rho have a limit D
in the tempered-distribution topology as rho tends to zero. Suppose their
restrictions to every compact subset of `E minus Z` are ordinary functions
and converge uniformly there to zero. Then D vanishes on every compactly
supported test in `E minus Z`, so

    supp D subset Z.

If D also has the connected clustering property (WC), the theorem gives
D=0. For four-point functions this turns such a *full nonexceptional*
comparison into a Gaussian four-point scaling distribution once the
two-point scaling covariance and Wick subtraction are fixed.

This weakens the Gaussian-four-point hypothesis in IR_SOURCE_REVIEW and
MULTIVELOCITY_REVIEW to three explicit inputs:

1. Existence of the connected scaling limit as a tempered distribution.
2. Locally uniform vanishing on the full nonexceptional open complement.
3. Proper-subset clustering of that limit in the Schwartz sense (WC).

These inputs still require proof for the microscopic retarded model.
In particular, pointwise convergence at individual momenta is not a
substitute for item 2, and a pure vector called a vacuum is not by itself
a proof of item 3 if it mixes locally distinguishable ordered sectors.

## 5. Why the presently stated PSS kinematics do not meet item 2

**[read]** [Porta--Scola--Singh, Proposition 6.2](https://arxiv.org/html/2509.08665v1#S6.SS1)
compares higher density-current correlations with the first m-1 external
frequencies all equal to eta_beta; the last frequency is their negative
sum. Spatial momenta vary in the stated bounded region. Its estimate is
uniform there, but that set is not the full nonexceptional open complement
in all frequency and momentum variables. It also concerns its specified
microscopic observables. Extending the theorem to arbitrary branch probes
and broader kinematics remains additional work.

This limitation cannot be bypassed using clustering alone. Here is a
concrete distribution-level counterexample with the correct four-current
scaling degree. Take d=2,m=4 and reduced momenta
`p_i=(omega_i,k_i)`, with p4=-p1-p2-p3, and set

    Q(p)=sum_(i=1)^4 (omega_i^2+k_i^2),
    V(p)=product_(1<=i<j<=4) (omega_i-omega_j)^2,
    F(p)=V(p)/Q(p)^7.                                   (3)

F is nonzero, symmetric under permutations of all four insertions, real,
even, smooth away from the origin, and homogeneous of degree -2, as
required for a connected four-point Fourier scaling distribution of
dimension-one fields in two dimensions. It is locally integrable at the
origin of the six-dimensional reduced momentum space and tempered.
On every permutation of the frequency slice `(eta,eta,eta,-3eta)`, it
vanishes identically for all spatial momenta.

Nevertheless F is nonzero on generic nonexceptional momenta. Its inverse
Fourier transform is a homogeneous distribution of degree -4, smooth
away from the origin. One way to see the latter is to expand the smooth
angular part in spherical harmonics: the degree -2 homogeneous Fourier
transform has finite nonzero gamma-function multipliers at all the
relevant harmonics in dimension six, preserving smoothness on the sphere.
It is locally integrable and decays as |x|^-4 away from the origin.
Schwartz smoothing is therefore bounded and tends to zero under every
nontrivial proper-cluster translation. Thus F obeys clustering but is
invisible on the PSS frequency slice.

No claim is made that (3) is the four-point function of a unitary CFT.
It proves that the kinematics, temperedness, scaling degree, permutation
symmetry and clustering alone cannot establish the missing full-support
statement. Additional CFT current Ward identities may detect a non-Abelian
bracket on that very slice; the coordinator's separate ABAB test does
exactly that under its additional chiral-current hypotheses.
