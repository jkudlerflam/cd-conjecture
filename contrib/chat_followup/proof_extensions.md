# Exact extensions of the controlled classes

Imported from the October 2026 chat "Prove c d conjecture" on 2026-10-27.
The elementary arguments below are written out here; the XXZ, XYZ and
antiferromagnetic permutation-chain CFT identifications remain the standard
inputs already used by C6 and REPORT. This note proves no bound for the generic
interacting chain.

## 1. Symmetry of the Hamiltonian is sufficient for bond averaging [proved]

Let a compact on-site unitary group G act through U(g) on C^d, and suppose the
periodic family H_L = sum_j h_(j,j+1), with fixed h, is invariant under
U(g)^{tensor L} for all sufficiently large compatible rings. The symmetry
reductions in this note concern that family, not accidental symmetries of a
single small ring.
The displayed bond h need not itself be invariant. Define

\[
\bar h=\int_G dg\,(U(g)\otimes U(g))h(U(g)\otimes U(g))^\dagger.
\]

Linearity and invariance give sum_j bar h_(j,j+1) = H_L. Thus one may use an
invariant bond without assuming a particular choice of local energy density.
This also handles telescoping one-site terms. For a nonclosed subgroup of U(d),
use its compact closure.

For qubits, a nontrivial on-site U(1) can be put along Z. Its two-site sectors
have dimensions 1,2,1, so, modulo a constant and a telescoping term,

\[
H=\sum_j\left[J_z Z_jZ_{j+1}
 +t\sigma_j^+\sigma_{j+1}^-+\bar t\sigma_j^-\sigma_{j+1}^+
 +h_zZ_j\right].
\]

A site-dependent rotation removes the phase of t from the bulk and leaves a
twist on a ring. This is the exact reduction to XXZ in a longitudinal field
already analyzed in [C6](../../codex/work_C6_d2/c6.tex). With its stated
exact-solution input the conformal regimes have c = 1. A twist shifts conformal
zero modes, so the untwisted vacuum Casimir formula cannot be applied blindly
to each finite ring. The addition here is the explicit bond-averaging argument,
not a new solution of XXZ.

## 2. The permutation-chain result only needs a two-dimensional commutant [proved, with standard IR input]

Suppose the on-site representation of G satisfies

\[
\operatorname{Comm}\{U(g)\otimes U(g):g\in G\}
 =\operatorname{span}\{I,P\},
\]

where P swaps the two physical sites. Averaging as above forces
bar h = a I + b P, with real a,b. Therefore the entire symmetry class reduces
to the permutation chain, even if G is smaller than SU(d).

For b > 0 its standard infrared identification is SU(d)_1 and c = d-1.
For b < 0, one color flip above a polarized ground state has energy
2|b|(1-cos k), hence a quadratic band edge rather than the conformal premise.
For b = 0 the Hamiltonian is constant. This is a sufficient microscopic
criterion; arbitrary continuous symmetry, or an emergent symmetry, does not
imply the commutant condition. The fundamental SU(d) case and its CFT are
already in [the original paper, Eq. (17)](https://arxiv.org/html/2403.17242v2).

## 3. Spin time reversal plus spatial reflection reduces qubits to XYZ [proved, with standard IR input]

Write the general uniform qubit Hamiltonian, up to constants and telescoping
terms, as

\[
H=\sum_j\left[\sum_{a,b}J_{ab}\sigma_j^a\sigma_{j+1}^b
 +\mathbf B\cdot\boldsymbol\sigma_j\right],\qquad J_{ab}\in\mathbb R.
\]

Physical spin time reversal sends every Pauli matrix to its negative and
forces B = 0. Spatial reflection exchanges the order of the two sites and
forces J = J^T. A common spin rotation diagonalizes this real symmetric
matrix, giving zero-field XYZ. With the exact XYZ gap classification stated
in C6, every conformal regime has c <= 1. Spin time reversal alone permits
antisymmetric (DM) exchange and does not give this reduction. Complex
conjugation in a chosen basis is not the same symmetry as spin time reversal.

This identifies a sufficient symmetry class for C6's existing XYZ slice; it
does not solve C6's antiunitary-only residual.

## 4. Closure under exact finite-depth local circuits [proved]

Assume H has a ground-state conformal limit with total Virasoro central charge
c, and let H' = U H U^dagger. Let U be a unitary circuit with depth D and gate
support at most R adjacent sites, where D and R are independent of L. Assume
H' still meets the intended microscopic requirements (in particular strict
nearest-neighbour range and one-site translation invariance when applying
the original conjecture), on the same on-site dimension d.

The spectra agree exactly. For a spatial interval A, gates contained in A or
its complement leave S(A) unchanged. In each layer only a bounded number of
gates cross its two endpoints. A crossing gate changes S(A) by at most twice
the logarithm of the dimension of its support. Summing over the layers gives

\[
|S_{U\psi}(A)-S_\psi(A)|=O(DR\log d),
\]

independent of the interval length and L. Therefore the coefficient of log|A|,
and hence c for these conformal ground states, is unchanged. This extends any
established c-d subclass to its exact circuit-equivalent Hamiltonians that
remain in scope, including Hamiltonians whose displayed terms are interacting.

The assumption on the transformed interaction range is necessary: a generic
finite-depth circuit makes a nearest-neighbour Hamiltonian longer range.
Approximate simulation or perturbative equivalence does not imply this
statement. The same entropy reasoning is already used for one specific
soluble deformation in [C13](../../codex/work_C13/RESULT.md); here it is
recorded as a general closure rule.

Tensor-product closure is already in the original paper and REPORT: c adds
while d multiplies, and sum_a(d_a-1) <= product_a d_a - 1. It is not a new
result of this import.

## Provenance and checks

The original chat's model identity was not available in the retrieved record.
Codex (GPT-6) imported the statements and rederived the averaging, symmetry,
commutant and entropy arguments on 2026-10-27. These are analytic arguments;
no numerical experiments accompany this note and no new verification of the
standard Bethe-ansatz/CFT inputs is claimed. The complete quadratic proof in
that chat substantially duplicates [C1](../../codex/C1_freefermion_astra.out)
and is therefore linked, not copied as a new result.
