# C152: c ≤ 1 on the parity-symmetric qubit class (interacting spinless fermions)

Your folder is `../work_C152/`. Read `specs/_qubit_c_context.md` first.

## Setting

Bonds commuting with the global parity Π = ∏σ^z (equivalently with a π rotation about z) are, after a rotation about z, H = Σ [J_x XX + J_y YY + D(XY − YX) + J_z ZZ + h Z]. Jordan–Wigner maps this exactly to spinless fermions with one orbital per site, nearest-neighbour hopping t e^{iφ} (φ from D), p-wave pairing Δ, chemical potential μ (from h) and nearest-neighbour density interaction V (from J_z): the interacting Kitaev chain with a hopping phase. It is a 4-parameter family after scale. It contains the U(1) class (Δ = 0), zero-field XYZ (μ = 0, φ = 0, Baxter-integrable), the transverse-field Ising and XY chains (V = 0), and the frustration-free Peschel–Emery / Katsura lines.

## Tasks

1. Map the complete phase diagram of the family and identify every gapless phase, critical surface and multicritical locus, with its universality class and c. Use the exact solutions on the integrable submanifolds (XYZ at μ = 0, the free line V = 0, the frustration-free lines), bosonization and RG around the U(1) class (where the pairing operator's scaling dimension depends on the Luttinger parameter K, so pairing is irrelevant for strongly repulsive V and a c = 1 phase can survive Δ ≠ 0), mean-field and strong-coupling expansions, and the literature on the interacting Kitaev chain (Hassler–Schuricht, Katsura–Schuricht–Takahashi, Miao et al., and others on tricritical Ising points at attractive V). Pay particular attention to the hopping phase φ, which breaks inversion and can give asymmetric gapless phases.
2. Derive c ≤ 1 at every gapless point of the family at the theoretical-physics standard, or find a point where c > 1. The mechanism to make precise: one spinless fermion per site gives at most one Dirac fermion (c = 1) at weak coupling; at strong coupling a second gapless fluid needs bound pairs to coexist with single particles, and cluster energies (n − 1)V argue against it; multicritical points where an Ising line meets a Luttinger region (a potential Ising × boson, c = 3/2) need special care, since pairing and the Ising sector involve the same fermion.
3. Calibrate with numerics on the critical path: DMRG (copy `codex/work_C144/dmrg.py`) and ED on candidate multicritical points, entropies with open-boundary corrections and Casimir energies with velocities. Use Typhon for scans.
4. State what is proved, what is derived at physics level, and what is numerical, as a theorem-style summary for the class.
