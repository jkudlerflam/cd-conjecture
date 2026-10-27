# A1. Symmetric nearest-neighbour bonds, LSM anomalies and admissible WZW levels

You are working closed-book. Do not consult files, papers or the internet. Show every step of every computation, and check each numerical claim by an independent route (a second formula, a special case, or a short script you write and run). State plainly where an argument is a heuristic rather than a proof.

## Setting

Consider a quantum spin chain on (C^d)^{⊗L} with periodic boundary conditions, Hamiltonian H = Σ_j h_{j,j+1}, where h is one Hermitian d^2 × d^2 matrix (the same on every bond, so H commutes with the one-site translation T). A compact connected Lie group G acts on each site through a unitary irreducible representation R of dimension d, and h commutes with R ⊗ R. Suppose the low-energy limit is a unitary CFT whose chiral algebra contains the affine algebra ĝ at level k (a WZW model G_k, possibly tensored with other sectors), with Virasoro central charge c(G_k) = k dim G / (k + h^∨). Normalize levels so that long roots have length² 2.

The Lieb–Schultz–Mattis (LSM) theorem: if R is a projective representation of the adjoint group G/Z(G) whose projective class is non-trivial, then the pair (translation T, G/Z(G)) carries a mixed 't Hooft anomaly that must be matched in the IR. In the IR, T acts as an internal symmetry (often a center element of G, or a combination of a center element and another internal symmetry). In a WZW model G_k, the center element z acts through a simple current J_z, and the anomaly is diagnosed by the representations of G/Z(G) appearing in the z-twisted sector, which for the diagonal modular invariant is ⊕_λ (λ, J_z λ). For SU(2)_k with spin-½ sites this reproduces the known statement that k must be odd. A useful cross-check for any claim about G_k is to restrict to a subgroup H ⊂ G via a conformal embedding or via the embedding index (for the principal SU(2) in Sp(4), the index is 10, and SU(2)_10 ⊂ Sp(4)_1 is a conformal embedding).

"Parton level bound": writing each site as one fermion in d flavours with a Gutzwiller constraint gives SU(d)_1 ⊃ G_{ind(R)}, where ind(R) is the Dynkin embedding index of G ⊂ SU(d) through R. The conjecture under test says Σ_a c_a ≤ d − 1 for any such chain.

## Tasks

For each pair (G, R) below:

1. Give the decomposition of R ⊗ R into irreducibles, hence the dimension of the commutant of R ⊗ R (the number of linearly independent G-invariant bonds, including the identity), and say which irreducibles lie in the symmetric and in the antisymmetric square.
2. Give the projective class of R as a representation of G/Z(G) (trivial or not, and its order), hence whether an LSM anomaly is present.
3. For levels k = 1, …, 8, compute c(G_k), mark which exceed d − 1, and give ind(R).
4. Decide, with an explicit argument via the twisted sectors (and a subgroup restriction cross-check), for which k the WZW model G_k with the diagonal modular invariant can match the LSM anomaly by itself, with T acting as a center element. If no k can, say so and explain what additional sector would be needed. Treat the possibility that T acts as a center element times an outer automorphism, or combined with another internal symmetry of extra sectors, and state precisely what changes.
5. List the smallest-c candidate IR theories (G_k alone, or G_k times a minimal additional unitary sector) that are compatible with the anomaly and have total c > d − 1. These are the sharpest numerical targets for a counterexample search.

The pairs:

- (a) G = SU(2), R = spin 1, d = 3.
- (b) G = SU(2), R = spin 3/2, d = 4.
- (c) G = Sp(4) ≅ Spin(5), R = the 4 (fundamental of Sp(4), spinor of SO(5)), d = 4.
- (d) G = SU(3), R = 3, d = 3.
- (e) G = SU(2) × SU(2), R = (½, ½), d = 4 (the "spin-orbital" chain; note this is the vector of SO(4)).
- (f) G = SO(N), R = vector, d = N, for N = 5, 6, 7 and general N (careful with Spin(N) versus SO(N) and with the center of SO(N) for even N).
- (g) G = G2, R = 7, d = 7.

Then answer two summary questions.

- (S1) Which of the "level-violation targets" Sp(4)_2 on C^4, SU(2)_k (k ≥ 5) on spin 1, SU(3)_2 on C^3, SO(N)_3 on the SO(N) vector, and G2_k (k ≥ 4) on the 7 are excluded by anomaly matching alone, which survive, and what is the minimal anomaly-compatible theory with c > d − 1 in each case.
- (S2) Is there any general statement relating anomaly matching to the parton level bound k ≤ ind(R)? For example, is G_{ind(R)} always anomaly-compatible when R carries an anomaly?

Deliver a single self-contained report, with a summary table at the top and the derivations below it.
