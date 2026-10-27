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
