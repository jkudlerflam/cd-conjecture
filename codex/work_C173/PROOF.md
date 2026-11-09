# C173 derivations: two floating bosons at codimension 0, the zigzag-ladder route

This run was done by a Claude Opus 5.5 agent (2026-11-08 to 09), not by Codex.

Standard: theoretical physics. Tags [proved], [derived], [numerical], [heuristic], [assumed], [literature] open the paragraph they apply to. Scripts are in this folder: blockade.py (Schrieffer-Wolff check), eff_ed.py and ring_ed.py (ED of the effective model), twobody.py, block_dmrg.py with dmrg_d.py (C152's DMRG with the local dimension read from the bond).

Conventions. H = Σ_j [σ_j^T J σ_{j+1} + B·σ_j] on a chain, as in C154. For the parity class we use C152's fermion form (C152 PROOF eq. 1.2)

    H = Σ_j [ w c†_{j+1} c_j + h.c. − Δ (c_j c_{j+1} + h.c.) + V n_j n_{j+1} − μ n_j ],      (0.1)
    w = J_x + J_y + 2iD,  Δ = J_x − J_y,  V = 4J_z,  μ = 4J_z + 2h_z.

## 1. The chain as a zigzag ladder

[proved] Call the even sites leg A and the odd sites leg B, A_n = 2n, B_n = 2n + 1. Every bond of H joins A_n to B_n or B_n to A_{n+1}, so H is a zigzag ladder with rung couplings only and no intra-leg term. One-site translation T maps A_n → B_n → A_{n+1}, so it exchanges the legs and T² is the leg translation.

Consequence [proved]. The two "channels" of the ladder are not two species. A pair of leg fields (φ_A, θ_A), (φ_B, θ_B) must be permuted by T. In one-site momentum, the bonding and antibonding combinations of the legs are the momenta k and k + π of a single band. "Both legs are independent Luttinger liquids" therefore means "one band with two occupied valleys", at k and near k + π in the leg-decoupled limit and shifted from there by the rungs. This is the form in which the c = 2 question is posed below: a second floating boson from a second valley of the one-particle band. Composite channels (bound states of the host's quanta, C163 Lemma L, C171 R3) are a different route and are not covered by the ladder argument (section 7).

## 2. First-order rungs: the band has one harmonic

[proved] The quadratic (one-particle) part of any nearest-neighbour bond, around any uniform reference, is a single harmonic in k. For hopping, ε_1(k) = ε_0 + 2|t_1| cos(k − α). With pairing the Bogoliubov form is C154 Lemma 3, ω(k) = a_s sin k + √f(cos k) with f quadratic. Intra-leg motion (j → j + 2) only appears at second order, through a virtual intermediate state: t_2 ~ J²/E_v, with E_v the cost of the intermediate.

Lemma S (strong rungs) [proved]. Let ε(k) = ε_0 + 2|t_1| cos(k − α) + R(k), where R is any smooth periodic function with |R′| ≤ r_1 and |R″| ≤ r_2. If r_1 < |t_1|/2 and r_2 < |t_1|, then ε has exactly one local minimum.

Proof. Put u = k − α. At a critical point 2|t_1| sin u = R′, so |sin u| < 1/4. The critical points therefore lie in two arcs, |u| < arcsin(1/4) and |u − π| < arcsin(1/4). On the first arc ε″ = −2|t_1| cos u + R″ < −1.93|t_1| + |t_1| < 0, so only maxima are there. On the second arc ε″ > 1.93|t_1| − |t_1| > 0, so ε is strictly convex there and ε′ vanishes at most once on it. A minimum exists, so there is exactly one. QED.

For R = 2|t_2| cos(2k + β) + Σ_{r≥3} 2|t_r| cos(rk + β_r) the condition holds when |t_1| > 8|t_2| + 2Σ_{r≥3} r²|t_r|. (The exact J1–J2 threshold for a pure cos 2k term is |t_1| = 4|t_2|; Lemma S is a sufficient condition that tolerates every higher harmonic and any phase.)

Consequence [derived]. Wherever the dressed one-particle band is dominated by its first-order rung harmonic, it has one valley, so a second single-particle species never exists. A second valley needs |t_1| ≲ 4|t_2|, i.e. a first-order rung amplitude no larger than the second-order one it generates. With all couplings of order one there is no separation of scales and no control. A controlled weak-rung regime needs the first-order band amplitude to vanish to leading order while a virtual process survives, i.e. a point where the first-order band is flat.

Where the first-order band is flat [proved, from C143's isolation results plus the large-field limit]:
- (i) large field |B| ≫ |J|: the polarized state with t_1 tuned to O(J²/B). Gapped (gap of order |B|), so no gapless phase arises;
- (ii) the XZ point (flat free band, C143), done in C157: two Majorana chains, induced interaction of the wrong sign, c ≤ 1;
- (iii) the Fibonacci point h_F (J = J_z ZZ, B = 2J_z ẑ, extensive degeneracy, C143), with its first-order hop tuned small. This is new, sections 3 to 5.

(The ferromagnetic Heisenberg point is not of this type: the magnon band there is a nonzero single harmonic. C160 covers it.)

## 3. The Fibonacci cusp: an exact second-order zigzag model

### 3.1 Derivation

[derived; second-order Schrieffer-Wolff, checked numerically in 3.2] Work in the parity class (0.1) with V ≫ |Δ| and |w|, |μ| of order Δ²/V. The low-energy space P has no two adjacent particles (Rydberg blockade), and the cost of a configuration outside P is V times the number of adjacent occupied pairs. Hopping and −μn act inside P at first order. Pairing always leaves P, so it acts at second order. Enumerating the processes:

- NNN hop. From a particle at j, create the pair (j + 1, j + 2) (intermediate cost 2V, three in a row) and annihilate (j, j + 1). With fermion signs this gives +Δ² for the product of matrix elements, hence

      t_2 c†_{j+2} c_j + h.c.,    t_2 = −Δ²/(2V).                                  (3.1)

  It needs j + 3 empty. It never passes a particle (j + 1 is empty by the blockade), so it carries no Jordan-Wigner sign.
- Diagonal. In the vacuum each bond contributes −Δ²/V. A particle blocks two bonds and raises the cost of its two next bonds to 2V, a shift +3Δ²/V, absorbed in μ. Two particles at distance 2 share blocked bonds, which gives V_2 = −Δ²/V = 2t_2. At distance 3 the bond between them costs 3V instead of 2V + 2V − V, giving V_3 = −Δ²/(3V) = (2/3)t_2. Distances ≥ 4 give nothing.
- ΔN = 2 inside P appears only at order wΔ/V and Δ³/V² (create, create, annihilate), which in the cusp is smaller than t_2 by a factor of order Δ/V.

The effective model, in units |t_2| = 1, is therefore

    H_eff = P Σ_j [ (w c†_{j+1} c_j + h.c.) − (c†_{j+2} c_j + h.c.) − 2 n_j n_{j+2} − (2/3) n_j n_{j+3} − μ n_j ] P,   (3.2)

with w now measured in units of Δ²/(2V). It conserves N at this order. On an open chain particles never pass each other, so (3.2) is sign-free (fermions and hard-core bosons coincide).

### 3.2 Check against the full chain

[numerical, blockade.py] Open chains of L = 12, Δ = 1, w = (0.3 + 0.7i)|t_2|, μ near the band. The enumeration-based second-order model (fermion signs computed explicitly, edges exact) agrees with the closed form (3.2) to machine precision (0.0), and the 30 lowest levels of the full 2^12-dimensional chain agree with it to 1.2e-2, 1.5e-3, 1.9e-4 at V = 20, 40, 80. The error falls as V^{-3}, i.e. relative to t_2 as (Δ/V)², the size of the next even order (the odd-order ΔN = 2 terms only enter energies at sixth order). So (3.2) is the correct leading-order model, and it is quantitatively faithful for V/Δ ≳ 40.

### 3.3 Structure of (3.2)

[proved, algebra] Each leg alone (w = 0, no particles on the other leg) is hard-core bosons with leg hopping −1 and leg nearest-neighbour interaction −2, i.e. the spin-1/2 chain 2t_2 S·S: the ferromagnetic Heisenberg point, exactly, at second order. So the w = 0 cusp is two ferromagnetic legs coupled by the inter-leg blockade and V_3. This is the staggered-charge class a(XX − YY) + uZZ + bZ at large b, whose factorization line u = r is C153's z = 2 line with L + 1 ground states, so (3.2) reproduces C153's large-b structure. There the leg numbers N_A, N_B are conserved at second order, N_A − N_B exactly, and the third-order conversion vertex (one A and one B particle created at zero total momentum) is resonant and locks the total density, which is C153's c = 1.

[proved] The one-particle band of (3.2) is

    ε(k) = 2|w| cos(k − φ) − 2 cos 2k,   w = |w| e^{iφ},                          (3.3)

with two local minima, near k = 0 and near k = π, when |w| < 4 (roughly). For φ ∉ {0, π} both minima are incommensurate, at k_a ≈ −s and k_b ≈ π + s with s ≈ |w| sin φ/4 at small |w|. Their energies differ by about 4|w| cos φ and their curvatures by 4|w| cos φ. At φ = π/2 the band has the exact symmetry k → π − k, the two valleys are degenerate with equal masses, and k_a + k_b = π. This symmetry is the image of a lattice symmetry of (0.1): site inversion composed with a π rotation about z on odd sites and a π/2 rotation about z, which is a symmetry of the bond exactly when J_x + J_y = 0.

So the cusp is the first controlled regime in this programme where the zigzag ladder has two legs that are effectively decoupled at first order and coupled only through generated and blockade terms. Its two valleys are incommensurate for generic φ, so every inter-valley conversion vertex carries momentum (single transfer k_b − k_a ≈ π + 2s, pair transfer 4s, inter-valley pairing k_a + k_b ≈ π, intra-valley pairing 2k_a, 2k_b) and is forbidden off resonance hypersurfaces. Kinematically this is C154's pattern P3, and the ΔN = 1 and ΔN = 2 terms present in N (section 5.3) do not change that. Whether c = 2 is realized is a question of dynamics, sections 4 and 5.

## 4. The dilute cusp: the blockade makes the two valleys immiscible

[derived; inputs: the 1D low-energy threshold law, the Girardeau mapping, C164's mass-imbalance result] Take (3.2) at low density with the Fermi level near both valley bottoms. Write m_a, m_b for the valley masses, δ = ε_b − ε_a ≥ 0 for the offset (label the lower valley a) and η for the relative mass difference (δ and η are both of order |w| cos φ, section 3.3).

(a) Two particles from the same valley are identical fermions, so the blockade (r = 1) and V_2, V_3 act in the odd channel and shift the energy only at order n⁴.

(b) Two particles from different valleys are distinguishable at small relative momentum q (their valley label is conserved there, since a transfer carries momentum π + 2s), so they scatter in the even channel. That channel has a hard core at r = 0 and r = ±1 and the attractions V_2, V_3 at r = 2, 3. By the generic 1D threshold law the transmission amplitude vanishes linearly, T(q) = O(q), unless the potential has a zero-energy resonance, which is one fine-tuned condition. So at low density inter-valley pairs are impenetrable, and the order of valley labels along the chain is fixed up to exchange processes of order q², which give the pseudospin couplings J_⊥, J_z of order n³/m (C162's pseudospin chain).

(c) Energy. Impenetrable particles of equal mass have the spectrum of free fermions of the total density, with the valley sequence degenerate up to the pseudospin couplings (Girardeau). So a two-valley state of total density n costs E_F(n) + δ n_b + O(J_⊥, J_z), and the single-valley state costs E_F(n), with E_F(n) = π² n³/(6m). The second valley gains no kinetic energy and pays δ n_b. Mass imbalance adds C164's positive interface cost, of order η² n², which beats the n³ pseudospin couplings as n → 0.

Conclusion. In the dilute cusp the ground state fills one valley (c = 1), except in a thin wedge around the symmetric sub-stratum φ = π/2 of section 3.3, where δ = η = 0 by symmetry. There the valley sequence is a pseudospin XXZ chain, and c = 2 needs it to be easy-plane (|J_z| < J_⊥). The wedge, if that holds, is |δ| ≲ J_⊥(n) of order n³ and |η| ≲ n^{1/2}. It is open, so it has codimension 0 in B_2, but it hangs on a symmetric stratum and needs a dynamical sign. This is the same structure as the two-valley soliton gas of C162 and C164 (assumption A_sol), now for point-like blockaded magnons, where C162 noted that the easy-plane window is not exponentially thin.

Comparison with the J1–J2 chain [literature, recalled, not re-derived here]. In the frustrated J1–J2 chain in a field the two valleys at ±k_0 are exchanged by site inversion, so δ = η = 0 identically, and the near-saturation phase is decided by the same intra- versus inter-valley competition (Kolezhuk and Vekua, Hikihara et al.). The four-Fermi-point (c = 2) regime of the zigzag chain is known at intermediate magnetization below a cusp in the magnetization curve (Okunishi and coworkers). A one-site nearest-neighbour bond has no inversion that exchanges the valleys of (3.3) except on the J_x + J_y = 0 sub-stratum.

## 6. The classical spiral regime: the rungs lock the relative leg phase at first order

[derived; classical, with C160's single-period theorem as input] In the classical (large-S, or semiclassical ħ_eff → 0) regime, translation-invariant ground states are uniform or single-period (C160 M1 in the continuum, C123's 48-bond lattice test). In ladder language a single helix of pitch θ is two leg helices of leg pitch 2θ with a fixed relative phase. Suppose instead the two leg helices had independent phases ϑ_A, ϑ_B, both with leg pitch 2θ (T forces equal pitch). The rung energy per site, averaged along the chain, is

    E_rung(ϑ_A − ϑ_B) = s_⊥ cos(ϑ_A − ϑ_B + θ) + D_m sin(ϑ_A − ϑ_B + θ) + const,

up to the charge-two terms, which depend on ϑ_A + ϑ_B and average to zero for incommensurate pitch. Here s_⊥ and D_m are the symmetric transverse exchange and the DM component about the helix axis, the classical counterparts of the hop t_1. If (s_⊥, D_m) ≠ 0, E_rung locks ϑ_A − ϑ_B with stiffness of the order of the exchange, so a single phason survives, with c = 1 after quantization. If s_⊥ = D_m = 0 (the classical counterpart of the cusp), a nearest-neighbour bond has no intra-leg coupling at all at the classical level, because classical spins generate no virtual terms. Then no leg helix exists to begin with, and the classical ground state is set by the remaining rung terms (Ising and charge-two parts and the field). Either way, in the classical regime two independent leg phasons do not occur [derived]. The second-order intra-leg terms that section 3 relies on are quantum, of relative size ħ_eff, so the cusp is invisible to semiclassics. This is why C160's "one soft pair at most" and the present cusp do not conflict.

## 7. Regimes, and the windows of C171 R3 and R4

Summary of what the ladder picture gives:

| regime | statement | tag |
| :-- | :-- | :-- |
| strong rungs, |t_1| > 8|t_2| + 2Σ r²|t_r| in the dressed one-particle band | one valley, no second single-particle species (Lemma S) | proved for the band, derived for the infrared |
| weak rungs at large field | gapped polarized phase | derived |
| weak rungs at XZ | C157: c ≤ 1 | derived (C157) |
| weak rungs at the Fibonacci cusp, dilute | one valley, except a thin wedge around the J_x + J_y = 0 sub-stratum where the pseudospin sign decides | derived (section 4) |
| Fibonacci cusp, intermediate density | section 5 | numerical |
| classical spiral | rungs lock the relative leg phase, one phason | derived (section 6) |
| generic bonds, all couplings of order one | no control: the dressed band can have two valleys, nothing forbids it | open |

C171 R3 (composites with wavevectors in Λ = Zθ + Zπ: n ≥ 3, neutral n = 0, the 2Q_1 + π branch). The ladder argument is a statement about valleys of the one-particle band. A composite of n host quanta has its own band, built from n-body processes, and Lemma S does not constrain it. So R3 is not covered by anything in this run, in any regime. The only overlap is negative: the Fibonacci-cusp mechanism produces its second valley at the single-particle level, not as a composite.

C171 R4 (a repulsive host, K_1 ≈ 0.9, with a θ-commensurate second boson of small K_2, x_2 = 1/(4K_2) > 2 − K_1). Where the second boson would be a second single-particle valley, the ladder results cover R4 in the same regimes as above: Lemma S (no second valley, whatever K_2 is) in the strong-rung regime, and section 4 (impenetrable, immiscible valleys) in the dilute cusp, because neither argument uses incommensurability. They do not cover R4 when all couplings are of order one, or when the second boson is a composite. Two remarks narrow R4 further.
- In T (B = 0) the Fibonacci point is absent (it needs B = 2J_z ẑ), and the flat-band weak-rung points of T are XZ (C157) and the cyclic stratum (C154). So the ladder route offers no mechanism for R4 in T at all. R4 in T stays under C159's A6_TK as C171 says.
- The commensurate-valley version of the cusp exists: at φ = 0 the valleys of (3.3) sit at exactly 0 and π, so a θ-commensurate second boson with a small K_2 is what a second valley of hard rods would look like (the hard-rod Luttinger parameter is K = (1 − n)², C152 section 7). Its pairing and single-particle vertices, of dimensions 1/K and 1/(4K), become irrelevant for dense rods. This is a concrete microscopic home for the R4 window in the parity class and in N, and it is exactly the regime that section 5 tests numerically.

## 5. Numerics on the effective model (3.2)

All runs are open chains. Units |t_2| = 1, |w| = 1 unless stated. The two-size estimator is c = 6 [S̄(L_2) − S̄(L_1)]/ln(L_2/L_1), with S̄ the average of the entanglement over the three middle block cuts (single-size Calabrese-Cardy fits at L = 48 are useless here: the c = 1 control reads 1.74 and the c = 2 control 1.39).

### 5.1 Checks and controls

- [numerical] block_dmrg.py (3-site blocks, local dimension 5, the cross-block blockade as a penalty U = 40, whose virtual corrections are of relative size |w|²/U) reproduces grand-canonical ED at L = 18 (energy to 6e-5, N = 3 exactly).
- [numerical] Controls, 48 → 96: hard rods with only the hop w (c = 1) read 0.90. Two decoupled XX legs (w = 0, no blockade, NNN hop allowed to pass, c = 2) read 1.88. So the estimator separates 1 from 2 at these sizes.
- [numerical, twobody.py] No two-body bound state lies below 2ε_min for |w| = 0.5, 1, 2 at φ = 1.3 (the lowest two-body level sits 0.002 to 0.008 above threshold), so the dilute gas is a gas, not a pair fluid.

### 5.2 The phase structure seen

| φ | μ | n | L | reading | tag |
| :-- | :-- | :-- | :-- | :-- | :-- |
| 0 (commensurate valleys) | −1.8 | 0.25 | 48 → 96 | c = 1.15 | numerical |
| 1.3 | −2.2 | 0.16 | 48 → 96 | c = 1.11 | numerical |
| 1.5 | −2.2 | 0.10 | 48 | S = 0.82, one valley | numerical |
| 1.3, π/2 | −1.8 | ≈ 0.3 | 48, 72, 96 | S saturates (≈ 0.8 to 0.9, period-9 pattern), gapped | numerical |
| 1.3, π/2 | −1.4, −1.65 | 0.45 | 48, 96 | gapped or c ≈ 1.0 near the CDW | numerical |
| 1.3, |w| = 0.5 | −1.8 | 0.5 | 48, 96 | product CDW at L = 96 (S = 4e-7), first-order jump from the dilute side | numerical |
| π/2 | −2.2 | 0.17 | 48, 72, 96 | S̄ = 1.646, 1.790, 1.869: c = 2.13 (48 → 72), 1.93 (48 → 96); not converged at 72 and 96 | numerical |
| 1.56, 1.565 | −2.2 | 0.17 | 48 | identical to π/2 | numerical |
| 1.54 | −2.2 | 0.14 | 48 | S = 1.56 | numerical |

Everything away from φ ≈ π/2 is c ≤ 1: one valley at low density (as section 4 predicts), gapped commensurate or CDW states at intermediate and high density, and a first-order jump into the CDW at small |w| (the ferromagnetic legs of 3.3). The exception is the dilute symmetric point and its neighbourhood down to φ = 1.54, where the two-size estimator reads c ≈ 2.

### 5.3 The dilute symmetric point

[numerical, valley_ed.py] On rings at φ = π/2 the ground state is valley-balanced: for N = 4 (L = 28, 32) the imbalance is 0 and the fully polarized states lie 0.07 to 0.10 higher, for N = 3 (L = 24) the ground doublet has imbalance ±1 and the fully polarized states lie 0.28 higher. So the valley pseudospin of section 4 is not ferromagnetic. It is easy-plane or antiferromagnetic, and only the former gives c = 2.

[numerical, not converged] The DMRG reads c ≈ 2 between L = 48 and 96, against 1.88 for the c = 2 control, but at L = 72 and 96 the truncation is 5e-7 and 4e-6 and the energy still moves by 1e-3 between χ = 64 and 96 at L = 96. Under-converged DMRG usually underestimates the entanglement of a critical state, so this is not obviously an overestimate, but it does not meet C163 §4.2's first-stage criterion (truncation below 1e-7 at L = 128). An antiferromagnetic pseudospin with a gap below the scale reached here would read the same, so a long crossover cannot be excluded.

[heuristic] Why it could be real. At φ = π/2 the cusp is the frustrated zigzag chain in disguise: the gauge c_j → i^j c_j maps (3.2) to real NN hop |w| and NNN hop +1, so the valleys are the ±k_0 pair of a J1–J2 chain with J2/J1 large, exchanged by a reflection. The frustrated J1–J2 chain is the textbook nearest-plus-next-nearest example with a four-Fermi-point c = 2 phase. The blockade and the attractive V_2, V_3 change the pseudospin couplings, and nothing derived here fixes their sign.

### 5.4 Status as a candidate

If confirmed, this is a c = 2 phase of a nearest-neighbour one-site qubit chain on an open set of B_2:
- the wedge around the J_x + J_y = 0 sub-stratum of the parity class is open inside the parity class (the reading persists to φ = 1.54, a valley offset 0.12|t_2|),
- leaving the parity class adds ΔN = 1 terms of order Δ²/V (and smaller), which couple to single-valley vertices at momenta k_a ≈ −s, k_b ≈ π + s with s ≈ |w| sin φ/4 ≈ 0.25, so they are off resonance and irrelevant away from resonance hypersurfaces, and the c = 2 region would be open in B_2 (codimension 0, pattern P3 of C154).

An explicit bond at the centre of the wedge (V/Δ = 40, |w| = 1, φ = π/2, μ_eff = −2.2), in C152's convention h = J_x XX + J_y YY + D(XY − YX) + J_z ZZ + (h_z/2)(Z⊗1 + 1⊗Z):

    J_x = 0.5, J_y = −0.5, D = 0.00625, J_z = 10, h_z = −19.97625.

The low-energy scale is |t_2| = Δ²/(2V) = 0.0125, and the dilute window in h_z has a width of order 1e-3. The bond is in the parity class and on its J_x + J_y = 0 sub-stratum; generic members of the wedge move J_x + J_y off zero by up to about 0.1|t_2| and add parity-odd couplings of order Δ²/V.

Tests (C163 §4.2 rules, adapted). First stage: two-size and profile c above 1.2 at L = 128 with truncation below 1e-7, and the valley pseudospin gapless (a two-valley momentum distribution and an L^{-1} pseudospin gap on rings). Second stage: FES at χ ≥ 256. Then the full qubit chain at V/Δ = 40 and 80. Typhon job 17603887 (pending) runs the effective model to L = 192 at χ = 192 at φ = π/2, 1.54, 1.52, both controls at L = 144, and the full chain at L = 96, 144. A local full-chain run at L = 48 did not converge within an hour (the scale V/|t_2| = 3200 stalls Lanczos, as the U = 1e4 penalty did in 5.1), so the full-chain test needs either the Typhon job or a code that projects onto the blockaded space.
