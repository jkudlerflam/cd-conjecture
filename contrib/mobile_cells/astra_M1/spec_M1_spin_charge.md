# Spin–charge separation for mobile cells on a ring (closed-book physics derivation)

Standard: theoretical physics (bosonization, RG power counting, standard CFT finite-size
formulas). Show the steps. Give honest hedges where an argument is heuristic. No need for
mathematical rigor. Do not run long numerics; small symbolic or numeric checks are fine.

## Setup A: the microscopic toy (d = 6)

A ring of L sites, periodic. Each site carries one of six letters
L, R, #↑, #↓, X↑, X↓. Assign "arrows" (left, right) to letters:
L = (0,0), R = (1,1), # = (1,0), X = (0,1) (the spin label is ignored for arrows).
A nearest-neighbour pair (a, b) is "legal" if right-arrow(a) = left-arrow(b).
The Hamiltonian is H = Σ_i h_{i,i+1} with

h = K Q + (μ/2)(n_# ⊗ 1 + 1 ⊗ n_#)
    − t Σ_σ (|L Xσ⟩⟨Xσ R| + h.c.)
    − u Σ_σ (|#σ L⟩⟨R #σ| + h.c.)
    + (J/2)(S⁺⊗S⁻ + S⁻⊗S⁺) + J_z S^z⊗S^z   acting only on adjacent pairs (#,X) and (X,#),

where Q projects onto illegal pairs, n_# projects onto # letters, S are spin-1/2 operators on
the σ label of # and X letters. Facts you may use without proof: each Q_i commutes with H, and
for K large the ground state is in the legal sector (all pairs legal). In the legal sector the
state is a binary arrow sequence b_i on bonds, letters are determined by (b_{i-1}, b_i),
"walls" are the sites with b_{i-1} ≠ b_i, they alternate between X (0→1) and # (1→0), the
number of walls W is even, and each wall carries a spin-1/2.
The term with t moves an X wall by one site; the term with u moves a # wall by one site.
Walls never pass each other. Call the interval from one # to the next # a "cell"; it contains
exactly one X ("head"). Its length ℓ is the distance between consecutive #'s.

Known facts (take as given): at u = 0, J = 0 each head is a particle in a box with
lowest energy −2t cos(π/ℓ); the cell energy is E_ℓ = μ − 2t cos(π/ℓ); for μ between suitable
crossing values a unique cell length b minimizes E_ℓ/ℓ, and at a crossing value μ_c(b) the
lengths b and b+1 are degenerate in energy density. At small u near μ_c(b), projecting onto
cells of length b or b+1 with heads in their box ground states, the length variables
η_j ∈ {0,1} (ℓ_j = b + η_j) of consecutive cells hop as an XX chain on the cell index with
amplitude J_b = u A_b², A_b the overlap of the box ground states. At u = t the walls are free
hard-core particles (free fermions) when J = J_z = 0.

## Setup B: the general construction this toy stands for

Cells of fluctuating length (a "length fluid", one U(1) charge: the number of cells M) each carry
a payload space C^D. When a marker moves, the payloads ride along (transport acts as a scalar
on the payload, possibly times a length-dependent amplitude). Neighbouring payloads couple
through a bond whose strength may depend on the two cell lengths, schematically
λ r(ℓ_j) r(ℓ_{j+1}) h^p_{j,j+1}, with h^p a fixed payload bond. Alone (fixed lengths), the
payload chain h^p on the cell lattice is assumed critical with central charge c_p and its own
low-energy operator content (primaries O_α of dimension Δ_α carrying cell-lattice momenta q_α).

## Questions

1. Exact kinematics on a ring. In the legal sector at fixed W, give the exact decomposition
   of the Hilbert space in "squeezed" coordinates (wall positions × squeezed spin chain), with the
   identification needed on a ring. Derive how physical translation acts, and derive the
   boundary condition that links the physical momentum P to the squeezed-chain (cell-lattice)
   translation eigenvalue. Compare with Ogata–Shiba's factorization for the U → ∞ Hubbard chain.

2. Effective theory. Bosonize the length fluid (Luttinger parameter K_ℓ, cell-lattice density ν of
   long cells, or in the toy the wall density). Treat the payload as a CFT. Enumerate all
   couplings between the two sectors allowed by cell-lattice translation invariance (the
   effective Hamiltonian is a sum over cells of the same local function of nearby lengths and
   payloads), conservation of M, and the other symmetries you identify. For each, give its
   scaling dimension and the resonance condition on ν (or density) under which it is
   non-oscillating. State under what conditions the infrared theory is (length boson) × (payload
   CFT) with total c = c_p + 1, and list the exceptional densities or parameter regions.
   Treat separately: (a) payload with one site per cell whose low-energy operators sit at
   cell-lattice momentum 0 and π (e.g. an XX or Heisenberg chain on the cells), and (b) the toy,
   where both # and X walls carry a spin, so the squeezed spin chain has two sites per cell.
   For (b), identify the microscopic symmetries of the d = 6 chain (including any reflection,
   possibly combined with a relabeling of letters), and say whether they forbid the dimerization
   of the squeezed spin chain when u ≠ t.

3. Ring lengths. With static markers (u = 0) and a unique preferred length b, rings with
   b ∤ L are known to have ground-state degeneracy growing with L. For u > 0 in the fluid regime,
   predict the finite-size spectrum on a ring of length L (any L): the ground-state energy to
   order 1/L including the Casimir term with two velocities, the zero-mode (twist) contributions
   and their dependence on L, W and the squeezed-chain momentum, and the degeneracies.
   Is the low-energy spectrum on every ring a finite union of conformal towers of a c = c_p + 1
   theory? What changes if the payload chain needs its length to be a multiple of some N
   (e.g. SU(N)_1)?

4. For the toy with XX exchange (J_z = 0, J > 0), predict which estimators of c should give 2
   and with what leading finite-size or finite-entanglement corrections: half-chain entropy on
   rings, finite-entanglement scaling S versus log ξ for a translation-invariant MPS, and the
   Casimir energy with velocities. Also say what happens at J = 0.

Return a concise but complete write-up (equations in plain text or LaTeX).
