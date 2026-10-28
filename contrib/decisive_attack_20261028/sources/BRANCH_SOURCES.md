# Finite-range probes of every Fermi branch

Status: **[proved algebraic source construction and conditional analytic
matching lemma; infrared current algebra unproved]**.

Author: coordinating GPT-6-family Codex agent, 2026-10-28. This is an
AI-assisted derivation. A separate agent's
[independent audit](../exact/SOURCE_AUDIT.md) checked the proof and the
Fourier convention. The central charge throughout means the unitary
Virasoro charge.

The [C46 Sugawara argument](../../../codex/work_C46/PROOF.md) already proves
that r positive-level chiral U(1) currents in a unitary CFT imply c >= r.
This note does not repeat that argument. It supplies explicit local probes
of the prospective currents and distinguishes their initial full rank from
the still-required infrared theorem.

## 1. Finite-range branch resolution

Consider a translation-invariant quadratic chain with s complex orbitals
per cell and Hermitian Bloch matrix h(k). Suppose its zero-energy crossings
occur at N distinct momenta k_1,...,k_N modulo 2 pi. At each crossing there
is one simple zero eigenvalue with normalized eigenvector u_a and nonzero
velocity. Distinctness of the momenta is an explicit hypothesis; no exact
separate conservation law for each band is assumed.

Define the real trigonometric polynomial

    p_a(k) = product_(b != a) [1-cos(k-k_b)]
             / product_(b != a) [1-cos(k_a-k_b)].                 (1)

**[proved]** It has degree R=N-1, satisfies p_a(k_b)=delta_ab,
and is nonnegative on real k. Write p_a(k)=sum_(r=-R)^R p_ar exp(i r k).
Then p_a,-r=conjugate(p_ar), and

    sum_r |p_ar| <= D_a^(-1),
    D_a=product_(b != a) sin^2[(k_a-k_b)/2] > 0.                 (2)

Indeed each factor sin^2[(k-k_b)/2] has Fourier coefficient l1 norm
one; convolution is submultiplicative. The denominator is D_a. If the
least circular separation is delta, then D_a >= sin(delta/2)^(2R).
The constants can be large, but are finite at every fixed certified model.

Let c_x be the s-component annihilation column. The Hermitian density

    O_a(x) = (1/2) sum_r p_ar
              [c_x^dagger c_(x+r) + c_(x-r)^dagger c_x]          (3)

is a sum of number-conserving bilinears. Each term spans at most R cell
bonds. The union of the terms defining this one density lies in
[x-R,x+R]; the operator is therefore local, with support of at most
2R+1 cells. Its norm is bounded by s sum_r |p_ar|. The bound follows
by summing the s orbital bilinear norms, each at most one.

With consistent Fourier signs, its small transferred-momentum vertex is

    F_a(k,q) = [p_a(k)+p_a(k+q)] I_s / 2.                       (4)

In particular u_b^dagger F_a(k_b,0) u_b = delta_ab. Thus the N-by-N
matrix of bare source couplings to the N low-energy branch densities
is exactly the identity. Off-diagonal branch bilinears have a nonzero
momentum difference and do not contribute to the smooth q near zero
localization. Orthogonal nonzero-energy bands at the same k are massive.

For real k, |p'_a(k)| <= R/D_a. Consequently shifting all Fermi points
by at most eta changes each source-matrix entry by at most R eta/D_a.
Adding any further localized vertex correction E preserves invertibility
whenever, for example,

    N max_a(R eta/D_a) + ||E||_infinity < 1.                    (5)

Here the matrix infinity norm is the largest absolute row sum. This is
the Neumann-series criterion for I plus a matrix of norm less than one.
No microscopic U(1)^N symmetry is required for this algebraic statement.
It also does not prove emergent U(1)^N currents.

## 2. Lifting the probes into a frozen physical cell

Suppose a cell of fixed physical length b has a parity-preserving code
isometry G:Fock(C^s)->K. Every even logical operator supported on a
finite cell interval I has a physical lift

    O_phys,I = G_I O_I G_I^dagger,                              (6)

extended as zero off the code on that interval. This lift has exactly
the same norm, is even, and is reducing for the all-code projection.
The equality of norms follows from G_I^dagger G_I=1. In an ordinary
tensor-product description a logical bilinear includes the fermion
parity factors on the intermediate cells; these stay inside I. Thus
(3) has a physical local realization, with no infinite Jordan-Wigner
string and no change of the Hamiltonian or of its local dimension.

To retain the smaller support of each summand, lift each bilinear on its
own interval and then sum. The triangle norm bound in Section 1 survives.
Alternatively lifting the whole density on [x-R,x+R] preserves its exact
operator norm but adds code projections throughout that larger interval.
The smaller term-support count below uses the first construction.

For the certified n=23 frozen model, N=46 gives R=45. Each bilinear
term occupies at most 46 cells and the full density at most 91 cells.
At b=24 these are at most 1104 and 2184 physical sites, respectively.
These are fixed supports, independent of the ring length. This is an
existence construction of probes, not a proposal for a practical
measurement or for computing their enormous matrices explicitly.

This paragraph applies in a chosen frozen ordered sector. For mobile
cells, a fixed number of cells need not occupy a bounded number of
physical sites. No claim about physical locality in that setting is
made here.

## 3. Analytic source matching at zero temperature

Use the hypotheses and strengthened coefficient norms of the existing
[finite-cell matching theorem](../../direct_attack_20261028/infrared/FINITE_CELL_MATCHING.md).
Let W be its auxiliary effective-action map, K its block-diagonal
vacuum-annihilating interaction, and R its leakage, with norms k and r.
Suppose a physical source is locally reducing for the code. Its embedded
operator can be written as O_log+D_O, where

    D_O P_vac = P_vac D_O = 0.                                (7)

For example (6) has precisely this property. Bounded time-dependent
test functions multiplying such operators retain it pointwise. Write
D(z)=sum_(j=1)^m z_j D_j and assume coefficient-norm bounds ||D_j||<=s_j.
Fix radii rho_j>0 satisfying

    k+r+sum_j rho_j s_j <= 1/(4 omega).                       (8)

**[proved conditional on that matching theorem]** At zero temperature,
for |z_j|<=rho_j,

    F(z) := W(K+D(z)+R)-W(K+D(z)),
    W(K+D(z))=0,                ||F(z)|| <= 4r.              (9)

The second equality follows from the same rightmost-vertex vacuum
annihilation proof: every term in K+D(z) contains both an auxiliary
creator and an auxiliary annihilator. It holds with arbitrary external
logical fields. The bound is the existing Lipschitz proof with K+D(z)
as its reference. Its complex Cauchy disk remains within the prescribed
analyticity ball by (8). Applying the multivariable Cauchy formula gives
for every multi-index alpha

    ||partial_z^alpha F(0)||
       <= 4 r alpha! product_j rho_j^(-alpha_j).              (10)

This proves that a fixed finite number of matched source derivatives
is proportional to the leakage. The logical source O_log factors out
of auxiliary integration, just as the logical Hamiltonian does. The
statement concerns logical effective-action coefficients with bounded
test functions; it is not yet a bound on their gapless expectation
values, a uniform pointwise current correlator, or a thermal theorem.
Scalar coefficients can be included using the local-marker density
norm from Section 6 of the matching theorem, with its extra fixed
constant and Cauchy-radius reserve.

If the Hamiltonian was put into this frame by a quasi-local unitary U,
the source in the original frame is U^dagger O_phys U. It is quasi-local
with the controlled tails of U. Taking the untransformed finite-range
O_phys instead generates additional source leakage. That leakage must
be estimated; it cannot silently be absorbed into the Hamiltonian's r.

## 4. What remains to reach a current lower bound

For N=2n simple branches, the construction provides independent bare
right- and left-branch probes. A source-dependent RG theorem retaining
a perturbatively invertible N-by-N vertex matrix would preserve their
rank. If it also constructs the Gaussian positive-level chiral algebra
in a unitary conformal limit, the existing Sugawara argument gives
c >= n. This could be enough to refute c<=d-1 at n=23,d=23 without
proving that there are no additional neutral massless sectors.

The following steps remain unproved for the candidate:

1. Uniform marginal-flow control for its high-degree retarded action.
2. Source-dependent infrared comparison for these independent probes,
   including scale-dependent vertex renormalizations and contact terms.
3. The conformal/current-algebra interpretation in the physical bulk.

The published [Porta-Scola-Singh paper](https://arxiv.org/html/2509.08665v1),
Sections 6.1 and 6.2, compares specified density/current observables
for its microscopic class. Our construction supplies more observables;
it does not expand that paper's hypotheses or prove its extension.
It is also invalid to infer the full current rank just by counting
microscopic charge labels, or to infer the thermal scalar from normalized
correlators. This note removes an elementary source-availability issue
and supplies a source derivative bound, while leaving those analytic
obligations explicit.

## 5. Reproduction

Run `python3 -B check_interpolation.py` in this folder. Standard-library
Gaussian-rational arithmetic checks all 2920 cardinal entries for
N=2,4,8,12,24,46, including the Hermitian coefficient identity. The points
are rational points on the unit circle, not the C33 Fermi roots; this is
an independent algebra check, not a replacement for their certificate.
The separate reviewer also checked the real-space/Fourier vertex identity
on an 11-site ring in [source_checks.py](../exact/source_checks.py).
