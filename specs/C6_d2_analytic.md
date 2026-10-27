# C6: Settle the c–d conjecture at d = 2 (nearest-neighbour qubit chains), analytically as far as possible

Closed book: do not read or list any files you did not create. Work only inside the current working directory.
Python (numpy/scipy) is allowed for checks; keep total runtime under 60 minutes. Work hard and be honest: a precisely
located residual set is worth more than an overclaimed proof. Flag every step you assume rather than prove.

## The d = 2 problem
Every translation-invariant NN Hermitian qubit chain has bond
  h = Σ_{a,b ∈ {x,y,z}} J_ab σ^a ⊗ σ^b + Σ_a h_a (σ^a⊗1 + 1⊗σ^a)/2      (up to a constant),
since antisymmetric on-site parts A⊗1 − 1⊗A telescope away on a periodic chain. Write J = J_s + J_a with J_s symmetric
and J_a antisymmetric (the Dzyaloshinskii–Moriya vector D). A global on-site rotation in SO(3) diagonalizes J_s, so up
to rotation and overall scale this is an 8-parameter family (J_x, J_y, J_z, D, h). The conjecture at d = 2 says: if
such a chain is critical with a unitary conformal low-energy theory (possibly several sectors with different velocities),
then the total central charge is at most 1.

## Tasks
1. Canonical forms and equivalences. Besides global rotations, list every transformation that maps this family to
   itself (or to a known solvable chain) while preserving the central charge: staggered π-rotations about an axis
   (flipping signs of two components of J), uniform spiral rotations that remove a DM component along an axis of
   U(1) symmetry (at the cost of a twisted boundary condition), time reversal, Jordan–Wigner, Kramers–Wannier and other
   dualities (say when they preserve nearest-neighbour range). Reduce the parameter space as far as possible.
2. Symmetry classes. Classify the family by its exact symmetries: U(1) about some axis; Z2 × Z2; a single unitary Z2;
   only antiunitary symmetries; none. Give the codimension of each class in the 8-parameter space.
3. Prove c ≤ 1 class by class, with complete arguments or precisely stated standard inputs:
   (a) the free-fermion class reachable by Jordan–Wigner (XY couplings, field along the JW axis, DM along that axis;
       the bound c ≤ rank A(1)/2 ≤ 1 is known for quasi-free NN chains with two Majoranas per site);
   (b) the U(1) class (XXZ plus longitudinal field plus DM along z, which reduces to twisted XXZ in a field): use the
       Bethe ansatz / Luttinger-liquid structure, and say which statements are rigorous and which are standard input;
   (c) Z2 × Z2 (XYZ, Baxter's eight-vertex model);
   (d) any further class you can control (for example via exact mappings, frustration-free points, or a
       "one U(1) at most" argument: rank of commuting on-site charges ≤ d − 1 = 1).
   For each class also say whether multicritical points can have c > 1 (for example Ising × boson with c = 3/2).
   Note that NN SU(2)-invariant qubit chains are only Heisenberg, which excludes SU(2)_k with k > 1.
4. The residual set. State precisely which part of the 8-parameter space is not covered by a proof (for example
   generic single-Z2 or symmetry-free chains, where Ising-type and Luttinger-type criticalities could meet). Argue as far
   as you can whether c > 1 is possible there: what structure would a c > 1 point need (two independent gapless
   sectors built from three Pauli operators per site, emergent U(1) symmetries, anomaly or LSM constraints), and whether
   any lattice obstruction applies. Propose the decisive numerical tests for the residual set and run what you can with
   exact diagonalization (momentum-resolved, L up to about 22), using finite-size conformal spectra (towers) as well as
   ground-state energy and entanglement estimators.

## Deliverable
A classification table (class, codimension, status: proved / standard input / numerical / open), complete proofs for
the proved classes, the precise residual set, and the numerics you ran. Save code as c6.py.

## Background (shared by all campaign specs; everything here is given, you may use it)

The c–d conjecture (J. I. Latorre and G. Sierra, arXiv:2403.17242, J. Stat. Mech. 2024, 113103).
A chain of L sites with Hilbert space (C^d)^{⊗L}, periodic boundary conditions, and
H_L = Σ_{i=1}^{L} h_{i,i+1}, where h is one fixed Hermitian d²×d² matrix acting on sites i, i+1
(so H is translation invariant by one site). If the low-energy spectrum of H_L is described by a
collection of unitary CFTs (possibly with different velocities) with central charges c_a, then
  c_T = Σ_a c_a ≤ d − 1.
Equality holds for the SU(d) Uimin–Lai–Sutherland (ULS) chain h = P (the swap), whose IR is SU(d)_1 with
c = d − 1. So the conjecture is equivalent to c_max(d) = d − 1. Longer-range chains reduce to NN ones by
blocking (a term spanning r+1 sites gives local dimension d^r).

Established facts (you may use them; do not re-derive):
1. c_L = c_R for such chains (Kapustin–Spodyneiko, PRL 123, 060601: no equilibrium energy current).
2. Quasi-free fermion NN chains (any period): c = (#zeros of det h(z) on |z|=1)/4 ≤ rank A(1)/2 ≤ log_2 d.
3. Entropic arguments cannot work: any monotone lattice c-function bounded by K log d at the lattice scale
   gives c ≤ K log d, false for SU(d)_1 at large d (block-entropy Casini–Huerta: K = 3, fails at d ≥ 7).
   Equivalently the conformal regime of the SU(d) chain starts only at l* ≳ (d−1)/(3 log d).
4. Thermodynamic / Lieb–Robinson / f-sum inequalities bound only c times a power of the conformal window:
   e.g. s = πcT/(3v) ≤ log d gives c T*/v ≤ (3/π) log d. ULS has ‖h‖ = 1 and spec(h) = {±1} for every d,
   so no bound in terms of ‖h‖, the Lieb–Robinson velocity, or the eigenvalues of h can exist.
5. LSM-type anomaly arguments give LOWER bounds c ≥ d − 1 in the SU(d)-fundamental and Z_d×Z_d classes.
6. WZW IR with on-site G acting on C^d via R: c(G_k) = k dim G/(k + h^∨). Gutzwiller partons give
   SU(d)_1 ⊃ G_{ind R} (ind = embedding index), so k ≤ ind(R) implies c ≤ d − 1 (GKO coset unitarity).
   For SU(N) fundamental, SO(N) vector, SU(2) spin 1, SU(N) antisymmetric², E6 27 the largest allowed k
   equals ind(R) (1, 2, 4, N−2, 6). Pure SU(2)_k needs no bound for spin ≥ 3/2.
7. Near misses (each drops one hypothesis): d = 2 with three-site terms reaches c = 3/2 (Rahmani–Zhu–Franz–
   Affleck interacting Majorana chain) and c = 2 (Σ S_j·S_{j+2}); SU(2)_k anyon chains violate c ≤ d_q − 1
   if d is read as the quantum dimension; Motzkin chains are not conformal; frustration-free gapless chains
   have gap O(1/L²), so they are never conformal.
8. Known saturating examples: SU(d)_1 (ULS), the qutrit XX chain (d=3, c=2), the t–J model (d=3, c=2).
   Lin–McGreevy exact single-site vector fixed points for d = 2 exist with c_Δ ∈ [1.13, 1.85] but none is
   scale consistent, so they say nothing about c_max.
