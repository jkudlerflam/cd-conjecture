# A refined c–d conjecture: central charge against the conformal window

## Claim

Class. Critical, translation-invariant, nearest-neighbour Hermitian chains H = Σ h_{i,i+1} on (ℂ^d)^⊗L whose low-energy theory is a set of unitary CFTs (c_a, v_a), c = Σ c_a (Virasoro). ‖h‖ is half the spread of spec h (minimized over h → h + o⊗1 − 1⊗o, which leaves H unchanged). With s_CFT(T) = (πT/3) Σ c_a/v_a, the window edge T*_ε is the largest T with |s/s_CFT − 1| ≤ ε below it (ε = 0.1 throughout). Hierarchy ratios are R = ‖h‖/T*, the onset length ℓ* = v_max/T*, and the Casimir coefficient A = Σ c_a v_a, defined by E_0(L) = L e_0 − πA/(6L).

1. Bounds [proved at the physics level, §1].
   - (a) Entropy: c ≤ 3ℓ* log d/(π(1−ε)) ≤ (6e/π) R log d/(1−ε).
   - (b) Energy current: c ≤ 12‖j‖/(πT*²) ≤ 24R²/π, j = i[h_12, h_23].
   - (c) Modular count: C7's c ≤ 1 + (6/π)K log d is (a) on a ring, weaker by a factor 2.
   - (d) Casimir: A ≤ 2ec‖h‖ always. A ≤ π‖h‖ for number-conserving free fermions. For a current algebra g_k of on-site generators in an irrep R, k v dim g ≤ 2π⟨D⟩ ≤ 8πC_2(R)‖h‖, D = −Σ_a [t^a_1,[t^a_1,h]], so a Sugawara-complete sector has A ≤ 8πC_2(R)‖h‖/(k + h^∨).
2. The SU(d) chain h = P has A = 2π(1 − 1/d)‖h‖, while its window closes as T* ≈ 6.1‖h‖/d², ℓ* ≈ d [numerical to d = 8, large-d form heuristic].
3. Conjecture: A ≤ 2π‖h‖ for every chain in the class. With (a) and c² ≤ A Σ c_a/v_a it gives c² ≤ 6R log d/(1−ε), so c ≤ d − 1 whenever T*/‖h‖ ≥ 6 log d/((1−ε)(d−1)²).

## Status

1 proved at the physics level, 2 numerical, 3 heuristic. A proof of the Casimir bound beyond the classes in 1(d), or one chain with A > 2π‖h‖, would settle it.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), sub-session of the coordinating session, 2026-10-27.
- AI assistance: gpt-6-astra (high), closed book via `gpt_solve.sh`, on `spec_fsum_casimir.md` (output `astra_fsum.out`) and `spec_ULS_largeN.md` (output `astra_ULS.out`, `astra_ULS/`: Bethe-ansatz v and e_0, a string TBA solver, a large-N limit).
- Checked by me. I derived the f-sum bound before reading Astra's, and the constants agree. My exact Schur–Weyl diagonalization of the SU(d) ring (all d at once, checked against brute force) reproduces Astra's TBA entropy to 3 or 4 digits wherever the ring has converged (d = 3, L = 12, T ≥ 0.44 and d = 4, L = 12, T ≥ 0.25), its lower window edge (0.775 against 0.778 at d = 3), and the Bethe e_0 to five digits. At d = 2, rings of 12, 14 and 18 sites approach its upper edge (s/s_CFT = 1.054, 1.078, 1.095 at T = 0.475, against 1.100). Other numerics are mine.

## Argument

### 1. The bounds

(a) s ≤ log d at every T, and s(T*) ≥ (1−ε)s_CFT(T*), so Σ c_a/v_a ≤ 3 log d/(π(1−ε)T*). Then c ≤ v_max Σ c_a/v_a and v_max ≤ v_LR = 2e‖h‖. Only the lower half of the window is used. The slack is σ* = s(T*)/log d, the share of the entropy budget spent inside the window.

(b) Right-movers thermal at T* and left-movers empty carry the current (π/12)cT*² at half the thermal energy density, whatever the velocities, while any lattice state has ⟨j⟩ ≤ ‖j‖ ≤ 2‖h‖². The hypothesis is that this partitioning-protocol state is conformal before it relaxes. It never binds (cπT*²/12‖j‖ is 0.005 for SU(5), 0.03 for Ising).

(c) Modular invariance gives the ring at β = K/v entropy at least πc/6 (C7 Lemma 2), and K log d bounds it. On SU(d) rings of 8 to 12 sites, d ≤ 6, the count of states below Δ = c/4 + 2 is within a factor 2.5 of the CFT count and far above C7's lower bound (`modular_count.py`).

(d) Lieb–Robinson gives A ≤ c v_LR, hence A ≤ 2π‖h‖ when c ≤ π/e. For the f-sum (`astra_fsum.out`), bondwise charge conservation makes the first moment of the charge structure factor exactly (1 − cos q)⟨D_a⟩. A g_k sector puts weight k|q|/(4π) per generator at ω = v|q|, and spectral weight is positive, so k v ≤ 2π⟨D_a⟩. No window enters, since both sides are O(q²). Then D = 2Φ(h̃) − 2C_2h̃ ≤ 4C_2‖h‖ with Φ(X) = Σ t^a X t^a, and Sugawara gives c_g v ≤ 8πC_2(R)‖h‖/(k + h^∨), which is below 4π‖h‖ for su(d) in its fundamental. For h = P, ⟨D⟩ = 1 − d e_0 with e_0 = 1 + 2[ψ(1/d) + γ_E]/d, and the CFT saturates the f-sum to 0.846, 0.858, 0.872, 0.885 for d = 2 to 5, tending to 1. The factor 2 between 4π and the SU(d) value 2π is lost entirely in D ≤ 4C_2‖h‖. For free fermions take the total charge, so D = −h_hop and each Fermi point carries unit weight, and the pinching inequality gives ‖h‖ ≥ ‖B‖_1 for the hopping matrix B, so A = ½Σ_FP |v_F| ≤ π‖h‖. A random search over one to four orbitals per site (`ff_casimir_search.py`) finds at most A/‖h‖ = 2.0.

### 2. Known models

T* is the two-sided window edge in the thermodynamic limit. Q = (1−ε)c²/(6R log d) measures saturation of the conjectured inequality.

| chain | d | c | ‖h‖ | v | A/2π‖h‖ | T* | ℓ* | R | σ* | Q |
|---|---|---|---|---|---|---|---|---|---|---|
| Ising | 2 | 1/2 | 1.41 | 2 | 0.11 | 0.725 | 2.8 | 1.95 | 0.30 | 0.03 |
| XX | 2 | 1 | 1 | 2 | 0.32 | 0.36 | 5.6 | 2.8 | 0.29 | 0.08 |
| XXZ, Δ = 1/2 | 2 | 1 | 1/2 | 1.30 | 0.41 | < 0.28 | | | | |
| Potts | 3 | 4/5 | 2.37 | 2.60 | 0.14 | 1.1 | 2.4 | 2.2 | 0.35 | 0.04 |
| Takhtajan–Babujian | 3 | 3/2 | 3 | 6.36 (fit) | 0.51 | | | | | |
| SU(2) | 2 | 1 | 1 | π | 0.50 | 0.475 | 6.6 | 2.1 | 0.25 | 0.10 |
| SU(3) | 3 | 2 | 1 | 2π/3 | 0.67 | 0.213 | 9.8 | 4.7 | 0.21 | 0.12 |
| SU(4) | 4 | 3 | 1 | π/2 | 0.75 | 0.128 | 12.3 | 7.8 | 0.20 | 0.13 |
| SU(5) | 5 | 4 | 1 | 2π/5 | 0.80 | 0.088 | 14.3 | 11.4 | 0.20 | 0.13 |
| SU(8) | 8 | 7 | 1 | π/4 | 0.875 | 0.0995 | 7.9 | 10.1 | 0.40 | 0.35 |
| SU(d), d → ∞ | d | d − 1 | 1 | 2π/d | → 1 | 6.1/d² | 1.03d | 0.16d² | 0.92/log d | 0.92/log d |

Ising and XX are exact. Potts and XXZ use rings of 8 and 9, and 12 and 14, sites, trusted where the two agree to 0.01 (above T = 0.63 and 0.28). Velocities match Casimir fits of E_0(L) to 3 × 10⁻⁴ (6 × 10⁻³ for SU(2) and SU(3), marginal logarithms). SU(d) windows are Astra's TBA, checked as above. For d ≤ 5 the entropy overshoots s_CFT by more than 1 + ε near T ≈ 2/d², which sets T*. For d ≥ 8 the overshoot stays below that and the edge is the undershoot near 6.1/d², which my ring data approach as Td² = 7.0, 6.6, 6.6, 6.5 for d = 3 to 6. The last row is Astra's ideal-Fermi-gas limit (N flavours at density 1/N, dispersion k²). It reproduces e_0, v and the CFT entropy, but its uniformity in N is not derived [heuristic].

So the entropy bound is sharp only up to σ* ≈ 0.2 to 0.4 in every known chain, and for SU(d) the slack grows like log d. The saturating family starts being conformal at ℓ* ≈ d, not at the entropic minimum d/log d (REPORT §3.1). In contrast A/2π‖h‖ climbs to 1 along SU(d).

### 3. The construction

Route A of `contrib/programmable_cell/` §4 is a logical SU(N) chain of cells of length b with exchange J. It has T* ≈ 6.1J/N², v = 2πbJ/N, ℓ* ≈ Nb, R ≈ N²‖h‖/(6J) and A ≈ 2πbJ. Its hierarchy needs J ≪ λ ≪ ρ²Δ and J ≪ t/b⁴ (the selection margin), so A/2π‖h‖ ≲ b⁻³ and R ≳ N²b⁴. The clock gaps J_c/L², the margin t/b⁴ and the visibility ρ² are what make R large. Route B (first-order Majoranas, m channels) has A ~ mbα²λ = bλ/(16mq), using α = 1/(4m√q) from C8 §2. Both satisfy (a) and the conjecture with margins that grow with b.

### 4. The refined conjecture

Conjecture. Σ_a c_a v_a ≤ 2π‖h‖ for every chain in the class. Hence, by Cauchy–Schwarz and (a), c² ≤ (Σ c_a v_a)(Σ c_a/v_a) ≤ 6R log d/(1−ε), c ≤ d − 1 whenever T*/‖h‖ ≥ 6 log d/((1−ε)(d−1)²), and unbounded c at fixed d needs R ≳ c²/log d.

Why this form. It is sharp (approached by SU(d) as d → ∞), contains no d, and is additive under stacking. It is proved for c ≤ π/e, for number-conserving free fermions (with π), and for Sugawara-complete on-site current algebras with a constant at most 4π for classical groups in their defining representations. It holds with room in every chain above. It also says what SU(d) shows, that central charge is bought with velocity. Entropy alone gives R ≳ c/log d, the conjecture R ≳ c²/log d, the scaling SU(d) obeys up to log d.

Why the construction cannot beat it [heuristic]. A chain of emergent cells is connected only if the logical operators induced at a cell's two ports fail to commute. In a cell with gap Δ_c, Lieb–Robinson and exponential clustering then force Δ_c ≲ (2v_LR/b) log(1/η), η the port visibility (`contrib/programmable_cell/` §6). The logical bond is at most ηλ with λ ≪ Δ_c. If the logical chain obeys the bound at its own scale, as free-fermion and SU(N) targets do, A = b × (logical Casimir) ≲ 2πbηλ ≲ 4πv_LR η log(1/η) λ/Δ_c ≤ 8π‖h‖λ/Δ_c, small whenever perturbation theory in the cell is controlled. The machinery that buys c dilutes the velocity by the cell length.

Why not c ≤ f(d) log R. Nothing above produces it. The 3^{−O(N²)} couplings of G1 are a choice (`contrib/unbounded_c/G1/G1_result.md` §8), and the tradeoffs of `contrib/programmable_cell/` §6 allow polynomial ones, which would break a log form while respecting R ≳ c²/log d.

What would falsify it.
1. A chain with Σ c_a v_a > 2π‖h‖, which needs c ≥ 2 at velocities above π‖h‖. A Casimir survey of random critical chains at d = 3 and 4 is the cheap test.
2. A violation of c–d with T*/‖h‖ ≥ 6 log d/((1−ε)(d−1)²). At d = 4 that is a window above 1.0‖h‖, eight times SU(4)'s, so the c = 16/5 target of `contrib/tci_target_d4/` would not falsify it.
3. A gapped cell whose two ports act non-commutingly on one register with Δ_c b ≫ v_LR log(1/η), escaping the clustering suppression.

## Reproducing

Numpy and scipy only.

- `uls_schur_weyl.py`: SU(d) ring spectra for all d (Young's orthogonal form, dense or FTLM per irrep) into `data/uls_L*.npz`, 15 minutes for L ≤ 12. `uls_array.slurm` ran L = 13 to 20 on Typhon. L = 15, 16, 20 (unused here) write to `/scratch/lustre/jkudlerflam/refined_conjecture/data/`. `check_schur_weyl.py` compares with brute force.
- `thermo.py` (entropies, SU(N)_1 torus entropy by modular S) and `analyze_uls.py` (window edges).
- `bench_thermo.py` (7 minutes) and `bench_windows.py`: the other chains.
- `lattice_constants.py` (5 minutes): ‖h‖, ‖j‖, Casimir fits, f-sum ratios, into `data/lattice_constants.json`.
- `modular_count.py` and `ff_casimir_search.py n` (up to 30 minutes per orbital number n), with stdout files.
- `astra_ULS/`: Astra's TBA scripts (`tba_check.py`, `roots_tba.py`, `scan_tba.py`, `ideal_gas.py`) and outputs.

## Relation to the rest of the repository

This makes REPORT §3.1 to §3.2 and C7 Theorem 3 quantitative. Suggested line for "Where things stand": "Refined conjecture: Σ c_a v_a ≤ 2π‖h‖, saturated by SU(d) as d → ∞ and proved for free fermions and on-site current algebras. It implies c² ≤ 6 log d ‖h‖/((1−ε)T*), so c–d can fail only in a conformal window parametrically narrower than SU(d)'s (`contrib/refined_conjecture/`)."
