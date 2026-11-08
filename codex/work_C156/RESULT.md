This run was done by a Claude Opus agent (2026-11-07), not by Codex.

# C156 result: literature check for nearest-neighbour one-site qubit chains with c > 1

## Flag

No credible translation-invariant, nearest-neighbour, one-site-cell qubit chain with c > 1 was found in the literature. Every report of c = 3/2 or c = 2 in a spin-1/2 chain that I found uses next-nearest-neighbour or three-site terms, a ladder, or a genuinely two-site cell.

One item stays open. The antiferromagnetic Kitaev spin chain in a [110] field (Sorensen, Gordon, Riddell, Wang and Kee, PRR 5, L012027 (2023), arXiv:2209.06221) lies in B_2 (shown below). The authors report a critical line for 0 < h < 0.511 K but give no central charge. In the language of C144 it is the XZ point with the field along a + b, a direction C144 did not test (it tested B_x and B_z separately, which gap everything). The expected value is c = 1 (heuristic, section 3). Small ED and DMRG here are not converged (the emergent energy scale is about 1e-3 K). Typhon job 17594698 runs DMRG to L = 128 (section 3).

Claims are tagged proved, derived, numerical, heuristic or assumed. Scripts in this folder are ed.py, check_equiv.py, check_kg.py, check_screw.py, kit110.py, kit110_dmrg.py and dmrg.py (copied from C144, with einsum optimization added to expect_bonds).

## 1. A gauge lemma that enlarges the list of B_2 models

[proved] Take a chain with bond matrices J_j and fields B_j. Suppose there is a fixed P in SO(3) with J_{j+1} = P J_j P^T and B_{j+1} = P B_j. This is a "screw" symmetry, a one-site translation combined with a global spin rotation. Then the site rotation R_j = P^j maps the chain to the uniform chain with J' = J_0 P and B' = B_0. On a ring the length must be a multiple of the order of P. The proof is one line, since R_j^T J_j R_{j+1} = P^{-j} P^j J_0 P^{-j} P^{j+1}.

[numerical check] Spectra agree to 1e-14 in four cases. The Kitaev-Gamma chain matches at L = 6 (check_kg.py), and the Kitaev chain in a [110] field at L = 8 (check_equiv.py). The Kitaev-XX-Gamma chain with a [110] field and the angle-theta compass chain match at L = 8 (check_screw.py).

So many "bond-dependent" Kitaev-type chains in the literature are B_2 points in disguise:

- The Kitaev-Gamma chain, with alternating x and y bonds, is first rotated by Yang et al.'s six-sublattice U_6 into a period-3 chain with screw symmetry R_a T_a. Here R_a is a rotation by 2 pi/3 about (1,1,1). The uniform form is J' = [[0, -K, 0], [0, 0, -G], [-G, 0, 0]] (rows: components on site j, columns: site j+1).
- The Kitaev-Heisenberg-Gamma chain in the U_6 frame has the same R_a T_a symmetry (Yang, Nocera and Affleck list it in G_1), so it is also in B_2.
- The Kitaev chain and the 90-degree compass chain become K sum S^x_j S^y_{j+1}, a rank-one J = a b^T with a orthogonal to b. This is the XZ point of C144. A compass chain at a general angle gives J = n_1 n_2^T.
- The Kitaev-XX-Gamma chain becomes J' = J_x P, where P is the pi rotation about (1,1,0). A field is allowed if it lies along (1,1,0).
- A generic bond-alternating chain, with different strengths on the two bonds, is not covered and stays outside B_2.

## 2. Annotated bibliography

Pauli normalization H = sum_j [sigma_j^T J sigma_{j+1} + B . sigma_j], up to overall scale. A "yes" under B_2 means a one-site nearest-neighbour qubit chain, directly or through section 1.

### 2a. Models in B_2 (all reported c <= 1)

1. XXZ in a z field. J = diag(1, 1, Delta), B = (0, 0, b). Bethe ansatz (Yang and Yang 1966). Gapless phases have c = 1. Consistent.
2. Zero-field XYZ (Baxter 1972). J = diag(Jx, Jy, Jz). Critical only where the two largest |J| coincide, with c = 1. Consistent.
3. Free-fermion family. This includes XY in a transverse field, with uniform DM along z and symmetric Gamma_xy (KSEA) terms, for example Abbasi, Mahdavifar and Motamedifar, Sci. Rep. 15, 20469 (2025). It also includes the Kitaev-XX-Gamma chain (arXiv:2509.21901), solved by Jordan-Wigner with c = 1 gapless phases and z = 2 lines. The compass chains, solved by Brzezicki, Dziarmaga and Oles (PRB 75, 134415 (2007)), Sun and Chen (arXiv:0906.1711) and Jafari (PRB 84, 035112 (2011)), have Ising-type transitions and a macroscopic degeneracy at the isotropic point. [derived] The quasiparticle energy D sin k + sqrt((h - J cos k)^2 + gamma^2 sin^2 k) has at most two zeros, because the zero condition is quadratic in cos k and the sign condition keeps one root of each pair. So c <= 1, as in C143. Consistent.
4. Interacting Kitaev chain = XYZ in a z field. J = diag(-(t + Delta)/2, -(t - Delta)/2, U), B = (0, 0, -mu/2). This line of work includes Hassler and Schuricht (NJP 14, 125018 (2012)) and Katsura, Schuricht and Takahashi (PRB 92, 115137 (2015), frustration-free line). The most recent field-theory and DMRG study is Kattel and Giamarchi (arXiv:2610.37761). Its phases are TSC, CDW, an incommensurate LL with c = 1 (fit 1.088) and a trivial phase. The TSC to trivial transition is Ising with c = 1/2, and the TSC to CDW meeting point lies on the eight-vertex line with c = 1. No c > 1 is claimed. The Kitaev-Hubbard chain study (arXiv:1911.03156) finds c = 1/2 on the transition line, rising to c = 1 only at the XY point. Consistent.
5. XYZ in a field along x (Dutta and Sen, cond-mat/0208216). Bosonization and RG give a single gapless line with c = 1 between commensurate and incommensurate gapped phases, despite no continuous symmetry. Consistent.
6. XYZ with uniform DM (Peotta et al., J. Stat. Mech. (2014), arXiv:1403.4568). DMRG to L = 400 finds Luttinger-liquid regions (c = 1) and two Ising transitions. Consistent.
7. XXZ with uniform DM and transverse field. Garate and Affleck (PRB 81, 144419 (2010)) and Chan, Jin, Jiang and Starykh (PRB 96, 214441 (2017), DMRG) find two Ising-ordered phases and one Luttinger liquid. No c > 1. Consistent.
8. XXZ with uniform Gamma_xy and a parameter alpha weighting S^y S^x against S^x S^y (Mahdavifar and Liu, Sci. Rep. 14, 30024 (2024)). iTEBD and Lanczos only, with no central charge reported. [derived] With alpha = 1 a global pi/4 rotation about z turns this into XYZ (Baxter). For general alpha it is XYZ plus D_z (item 6). The extended "gapless composite phases" they report should be read against those exact results.
9. The Kitaev-Gamma family, all in B_2 by section 1:
   - Kitaev-Gamma chain: Yang, Nocera, Tummuru, Kee and Affleck, PRL 124, 147205 (2020), arXiv:1910.14304. An emergent SU(2)_1 gapless phase with c = 1 (ED and DMRG, L up to 144), plus O_h to D_4 and O_h to D_3 ordered phases.
   - The c = 1/2 region near the AFM Kitaev point: Luo et al., arXiv:2012.03382.
   - Kitaev-Heisenberg-Gamma: Yang, Nocera and Affleck, PRR 2, 033268 (2020). Nine to eleven phases, four of them Luttinger liquids (SU(2)_1 or single boson, c = 1). A narrow region near the AFM Kitaev point is left unresolved.
   - Kitaev-Gamma with DM: PRB 111, 174414 (2025), arXiv:2204.13810, emergent SU(2)_1 with c = 1.
   - Gamma chain and the nonsymmorphic cubic groups: Yang, Nocera, Xu and Affleck, SciPost Phys. 17, 097 (2024). c = 1 up to theta = 0.49 pi.
   No c > 1 anywhere. Consistent.
10. XXZ in a staggered z field (Zhang and Song, arXiv:2609.29728). [derived] A pi rotation about x on every other site maps it to the uniform chain J = diag(1, -1, -Delta) (XXZ normalization) with a uniform B_z. This is the staggered-charge class of d2_gap_reduction. The paper has an exact frustration-free boundary Delta = cosh q (with h = sinh q) and a linearly degenerate ground space there. DMRG on 2N = 40 sites finds FM, LL and gapped phases. No central charge is given. This bears on C153 but asserts nothing about c > 1.
11. Kitaev spin chain in a [110] field (Sorensen et al., arXiv:2209.06221). [proved, section 1] The uniform form is J = K e_x e_y^T, i.e. K sum sigma^x_j sigma^y_{j+1}, with B along (1,1,0). The paper's results (ED to N = 36, DMRG to N = 1200):
    - a critical line for 0 < h < h_c1 = 0.511 K at phi = 45 degrees, gapped as soon as phi differs from 45 degrees or theta_z is nonzero;
    - a twofold-degenerate soliton phase with staggered vector chirality for h_c1 < h < h_c2 = 0.726 K;
    - correlations C(r) ~ r^{-0.25} at h = 0, and any field along x alone opens a gap.
    No central charge is reported for the line. Open, see section 3.

### 2b. Near misses: reported c > 1, outside B_2

1. Interacting Majorana chain (Rahmani, Zhu, Franz and Affleck, PRL 115, 166401 and PRB 92, 235123 (2015)), revisited in arXiv:2211.15598. DMRG finds a c = 3/2 Ising plus Luttinger-liquid phase or line, and a c = 7/10 tricritical point. [derived] With two Majoranas per spin, the four-Majorana term maps to Z_j Z_{j+1} + X_j X_{j+2}, so the spin chain is -t sum (Z_j + X_j X_{j+1}) + g sum (Z_j Z_{j+1} + X_j X_{j+2}), with a range-2 XX term. Translation by one Majorana is the Kramers-Wannier duality, not an on-site symmetry. After blocking this is d = 4. Not a counterexample.
2. O'Brien and Fendley (PRL 120, 206403 (2018)). Ising plus three-site terms (XXZ + ZXX type), with a tricritical c = 7/10 point. Three-site, so not in B_2.
3. Extended XY model with three-spin terms (Titvinidze and Japaridze, EPJB 32, 383 (2003)). It is exactly solvable, and its next-nearest hopping produces four Fermi points, giving c = 2 for J* > 2J. This is the cleanest free-fermion illustration of why longer range is needed: a second harmonic in the band.
4. Frustrated ferromagnetic J1-J2 chain in a field (Sudan, Luscher and Lauchli, PRB 80, 140402 (2009)). A c = 2 phase made of a magnon LL coexisting with a quadrupolar LL. Next-nearest-neighbour.
5. XXZ with next-nearest S^z S^z (Fitouchi and Chepiga, arXiv:2510.05988). c = 3/2 transition lines and a c = 2 region with two coexisting LLs, one of them nematic. Next-nearest-neighbour.
6. Zigzag chain with Kitaev interactions (arXiv:2202.09192), with a c = 2 point. Kitaev-Gamma ladders (Sorensen et al., PRX 11, 011013 (2021), and npj Quantum Mater. (2024)). Ladders or next-nearest-neighbour.
7. Bond-alternating Kitaev-Gamma chain (Luo et al. 2021). This is a genuine two-site cell (section 1 does not apply when the two bonds differ in strength). It reports c = 1 and c = 1/2 only.
8. Spin-1 Takhtajan-Babujian point with SU(2)_2 and c = 3/2, and c of about 2 near multicritical spin-1 points (arXiv:1404.4898). This is d = 3.

The pattern is uniform. In the spin-1/2 literature c > 1 always comes from a second species made available by range 2 (two Fermi seas, a bound-pair fluid, two sublattice chains) or by a second leg.

## 3. The open case: Kitaev chain in a [110] field

- Model in B_2: H = K sum S^x_j S^y_{j+1} - (h/sqrt 2) sum (S^x_j + S^y_j), with S = sigma/2 and K = 1.
- Symmetry [derived]. Let P be the pi rotation about a + b, which swaps a and b and sends a x b to its negative. P combined with a site reflection maps the bond a b^T to itself and fixes B along a + b. A field component along a - b breaks this symmetry, and the line h < h_c1 is the locus where it is preserved. That fits Sorensen et al.'s finding that any deviation in phi gaps the chain, and that the soliton phase is the twofold symmetry-broken phase on the same line.
- Expected c [heuristic]. At h = 0 the chain is a critical Majorana chain (c = 1/2) times a flat band of N/2 conserved qubits (the bond operators Y_{2j-1}Y_{2j} and X_{2j}X_{2j+1}, which anticommute in sequence, so they form a second Majorana chain with zero Hamiltonian). A field of order h gives that second chain dynamics at order h^2/K. If the symmetric line keeps both Majorana chains critical, they combine into one Dirac fermion and c = 1. This matches the generic picture of a Gaussian line that splits into two Ising lines around an ordered phase (the soliton phase B), and the magnetization profiles with several crossings (N = 384) that Sorensen et al. show below h_c1, which suggest an incommensurate c = 1 liquid. A c = 3/2 value would need the Majorana sector to stay critical alongside a separate bosonic liquid of the flat-band qubits, which this picture does not produce. The heuristic is not a derivation.
- What the numerics here show [numerical, not converged].
  - ED on rings with L = 8 to 16 at h = 0.2 and 0.35 (kit110.py): gaps of 1e-3 to 3e-3 K, level crossings between sizes, and half-chain entropies of 1.4 to 2.0. Chord fits give 2 to 3.5, but they are dominated by near-degenerate mixing and mean nothing at these sizes.
  - Open-chain DMRG at L = 24 and h = 0.35 (truncation 3e-8): S_mid = 1.08, with window fits ranging from 1.0 to 3.9.
  - The emergent scale v of order h^2/K means L = 16 to 24 lies inside the crossover. These numbers are not evidence for c > 1, and they cannot rule it out.
- Pending: Typhon job 17594698 (c156_k110): h in {0.2, 0.3, 0.4, 0.45, 0.6} and L in {48, 64, 96, 128}, with chi up to 160. Outputs are /scratch/lustre/jkudlerflam/c156/k_h*_L*.jsonl and .out, plus slurm_17594698.out. Each jsonl record holds S(l), parity-resolved fits c_fit for windows 0.05, 0.125 and 0.25, the energy history and the magnetization profiles. The decisive check is the slope of S(L/2) in log L at fixed h (c/6 for open chains), together with the parity-averaged fits as in C144's cfit.avg_fit. The h = 0.6 runs are inside the soliton phase and should give c near 0, as a control. Two local runs (h = 0.3, L = 32 and 48) were started at the end of this session, writing local.jsonl. The open case should be passed to C154 or C155 as an adversarial target, together with the YX and ZY directions of C144.

## 4. Rigorous and exact results that constrain c

- Exact solutions. Every Hermitian integrable nearest-neighbour qubit chain known belongs to the eight-vertex (XYZ), six-vertex (XXZ in a z field) or free-fermion classes, and all have c <= 1. De Leeuw, Pribytok and Ryan (J. Phys. A 52, 505201 (2019), arXiv:1904.12005) classify all regular difference-form R-matrices for qubits and find 14 solutions. Eight are eight-vertex or lower, and the six new ones are non-diagonalizable or nilpotent. De Leeuw, Paletta, Pribytok, Retore and Ryan (PRL 125, 031604 (2020)) add two non-difference eight-vertex-type models, related to the AdS_2 and AdS_3 S-matrices. So no exactly solvable counterexample is available. The classification only covers regular R-matrices, so this is a statement about integrable points, not a theorem about B_2.
- Frustration-free chains cannot be CFTs.
  - Bravyi and Gosset (J. Math. Phys. 56, 061902 (2015), arXiv:1503.04035) classify every frustration-free nearest-neighbour qubit chain as gapped or gapless, with gapless cases having an open-chain gap of O(n^-2).
  - Gosset and Mozgunov (J. Math. Phys. 57, 091901 (2016)) give the local gap threshold 6/(n(n+1)), so any gapless frustration-free 1D chain has an open-chain gap of O(n^-2), which means z >= 2.
  - Masaoka, Soejima and Watanabe prove quadratic or softer dispersion for frustration-free spin-1/2 nearest-neighbour systems (PRB 110, 195140 (2024)), and z >= 2 for frustration-free systems with power-law correlations (PRX 15 (2025), arXiv:2406.06415).
  - So the frustration-free strata of B_2 (C141's family, the Zhang-Song boundary, the Katsura-Schuricht-Takahashi line) carry no conformal point, and c is not defined there.
- Gap theorem. Bravyi, Gosset, Liu and Wong (arXiv:2607.14401) prove gap >= 2 mu for the Suzuki-Fisher Heisenberg-type family in a field mu. This is the input to the rank-one c <= 1/2 result of d2_gap_reduction.
- Rigorous RG near free fermions.
  - Benfatto and Mastropietro (Rev. Math. Phys. 13, 1323 (2001); J. Stat. Phys. 138, 1084 (2010)) and Benfatto, Falco and Mastropietro (CMP 330, 153 and 217 (2014)) construct the ground state of non-solvable XYZ-type and general short-range lattice fermion chains near a two-Fermi-point free point, and prove Haldane's Luttinger-liquid relations there. The infrared is a single Luttinger liquid (one boson), which is c = 1 at the physics level. The central charge itself is not proved in these papers.
  - Giuliani and Mastropietro (CMP 324, 179 (2013)) prove that the finite-size correction defining c is exactly 1/2 for non-solvable Ising models with small finite-range perturbations. The quantum-chain analogue (an XY chain at its Ising point with small nearest-neighbour interactions) is covered at the physics level by the same method, but I found no paper stating it.
  - Together these make c <= 1 rigorous only in small neighbourhoods of the free-fermion locus, with c = 1 itself rigorous only for Ising perturbations. Kozlowski's XXZ asymptotics (already in the campaign) cover the U(1) class.
- No general upper bound on c in terms of the local dimension or the range of interaction exists in the rigorous literature, as far as I could find. The Lieb-Schultz-Mattis and Oshikawa-Yamanaka-Affleck theorems give lower bounds (the system must be gapless or degenerate), not upper ones. Free fermions give the only counting statement: one complex fermion per site with a nearest-neighbour single-harmonic band has at most two Fermi points, so c <= 1. Range r allows up to 2r Fermi points (Titvinidze-Japaridze shows this at r = 2).

## 5. Bottom line

The literature contains no nearest-neighbour one-site qubit chain with c > 1, rigorous or numerical. Every known c > 1 spin-1/2 chain needs range 2 or a second leg. The one place where the literature reports a gapless line in B_2 without a central charge is the Kitaev chain in a [110] field. It is expected at c = 1, and the pending Typhon job should settle it. Section 1's gauge lemma shows that the whole Kitaev-Gamma family, the compass chains and the Kitaev-XX-Gamma chain are B_2 points, and all of them report c <= 1.
