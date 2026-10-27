# Counterexample search at d ≥ 3: symmetric families, anomalies, and the c = 16/5 target

## Claim

Class: one-site translation-invariant, nearest-neighbour, Hermitian chains H = Σ h_{i,i+1} on (ℂ^d)^⊗L with periodic boundary conditions, d = 3, 4, 5, 7, flowing to unitary CFTs with total Virasoro charge c_T = Σ c_a.

No chain with c_T > d − 1 was found in the places where one seemed most plausible. The concrete statements follow.

1. Pure Sp(4)_k, at any level, cannot be the infrared theory of an Sp(4) chain on ℂ^4, because no center twist of Sp(4)_k carries the spinor LSM anomaly. So the "Sp(4)_2 on ℂ^4" target is excluded, and the smallest anomaly-compatible theory above 3 with Sp(4) currents is Sp(4)_1 × tricritical Ising, glued, with c = 16/5. [Astra, checked]
2. That c = 16/5 theory is a natural continuation of the c = 3 phase of spin-3/2 chains, in which a lone singlet Majorana is protected by translation. One tuning could in principle make it tricritical. The boundary of the c = 3 region in the two-parameter spin-3/2 family shows no such point at L ≤ 12. [numerical]
3. The small-angle Sp(4) region (D3's loose end) is dimerized. Its large entropy slopes come from a finite-size crossover controlled by the ferromagnetic boundary point θ = 0, with scaling variable Lθ. [numerical, partly pending]
4. Dense scans of the spin-orbital SU(2)×SU(2) sphere, the spin-3/2 SU(2) sphere and the SO(5)-vector circle give no point with a calibrated estimate above d − 1. The apparent excesses are incommensurate sizes, rotor towers, dimer doublets and ferromagnetic multiplets, except two spin-3/2 points (c_half 3.2 to 3.4) whose gaps seem to open. [numerical]
5. A direct-sum construction gives c_max(d₁ + d₂) ≥ c_max(d₁) + c_max(d₂) + 1. It saturates d − 1 but cannot exceed it, and it implies that any violation at some d₀ would propagate to every larger d. [heuristic]

## Status

Numerical and heuristic. Item 1 is a proof at the level of the WZW twisted-sector calculation. The rest would be upgraded by larger-size DMRG with symmetric tensors, in particular along the c = 3 boundary of the spin-3/2 family.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), a sub-session of the coordinating session, 2026-10-27.
- AI assistance: one closed-book gpt-6-astra (high) run, A1, on anomalies and levels (`specs/A1_anomaly_levels.md`, output `specs/A1_anomaly_levels_astra_out.md`). All numerics are mine.
- Checks: the ED engine against brute-force diagonalization of the whole Hilbert space (agreement to 1e-15 on random charge-conserving bonds, and on the Sp(4) bond at L = 6, 8). A1's central claim I rederived through the principal SU(2), and I recomputed its tricritical-Ising arithmetic. iDMRG against the Heisenberg energy (2e-6) and against ED energies of dimerized Sp(4) chains (1e-4). An early iDMRG version gave a non-variational energy at long run times. I traced it to the convergence test and the prediction step, fixed both, and added a consistency flag comparing the energy difference with the bond energies. No result from the faulty version is used.

## Targets

Ranked as at the start, refined by A1.

1. Level violations. G-invariant bonds are counted by the commutant of R ⊗ R. SU(N) fundamental has 2 (only a + bP), spin 1 has 3 (BLBQ, fully known), Sp(4) on ℂ^4 has 3 (one angle), SO(N) vector has 3 (one angle), SU(2)×SU(2) on (½,½) has 4 (two angles), spin 3/2 has 4 (two angles, contains Sp(4) and SU(4)), and G2 on the 7 has 4 (two angles). A1 then gives the anomaly-compatible targets: Sp(4)_2 excluded, SU(2)_5 (c = 15/7) for spin 1, SU(3)_2 (16/5) or SU(3)_1 × Ising (5/2) for SU(3), SU(2)_1 × SU(2)_5 (22/7) or SU(2)_3 × SU(2)_3 (18/5) for spin-orbital, SO(N)_3 for vectors, G2_4 (7) for G2.
2. Boundaries of saturating phases, above all the c = 3 phase at d = 4.
3. Small-angle Sp(4).
4. Constructions by design.

## Methods and calibration

`ced.py` is a symmetry-resolved ED engine (any d and bond, Cartan charges, momentum blocks, Schmidt values, Casimir labels), and `idmrg.py` a two-site infinite DMRG. The estimators follow.

- c_half(L₁, L₂) = 3[S(L₂/2) − S(L₁/2)]/ln(L₂/L₁), half-chain entropy of the periodic ground state, at commensurate sizes. Velocity free.
- c_ent, the slope of S(ℓ) against log chord(ℓ) for ℓ ≥ 3. It overshoots at these sizes.
- c_Cas from E₀(L)/L = e − πcv/6L² (three sizes, 1/L⁴ term) with v from the lowest level at k₀ ± 2π/L. With split velocities it measures Σc_a v_a / v.
- z_eff from the lowest gap, and spectra labelled by momentum and Casimir.
- iDMRG dimerization D = e_AB − e_BA and correlation length ξ against χ. Saturating D and ξ means dimerized, decreasing D with ξ ∝ χ^κ means critical.

| Control | c | Sizes | c_half | c_ent at L_max | c_Cas |
|---|---|---|---|---|---|
| SU(3)_1, h = P, d = 3 | 2 | 6 to 18 | 2.06 (15→18), 2.00 (12→18) | 2.11 | 1.985 |
| SU(4)_1, d = 4 | 3 | 8, 12, 16 | 3.03 (12→16), 3.05 (8→12) | 3.34 | 2.98 |
| Takhtajan–Babujian spin 1 | 3/2 | 6 to 18 | 1.506 (14→18) | 1.54 | 1.480 |
| Heisenberg A + 3B, decoupled | 2 | 6 to 14 | 1.98 (10→14) | 2.05 | 3.96 (artifact, Σcv/v_min = 4) |
| Sp(4) θ = π, SO(5)_1 × Ising | 3 | 8, 12 | 2.98 (8→12) | 3.52 | 4.96 (single-velocity artifact) |
| SU(5)_1, d = 5 | 4 | 10 | n/a | 4.99 | n/a |

c_half is within 0.07 of the exact value on every control, including the split-velocity ones, so I quote it with a ±0.1 systematic band. c_ent at L = 12 overshoots SU(4)_1 by 0.47, so it is used only relative to the SU(d)_1 value at the same size. Single-velocity Casimir estimates are unusable whenever velocities split, which reproduces D3's c ≈ 4.72 artifact. iDMRG slopes of S against log ξ give 2.09 for the decoupled control and 1.35 to 1.69 for Takhtajan–Babujian at χ ≤ 32. Its two-site cell cannot represent period-3 or period-4 critical states (SU(3)_1, SU(4)_1, the c = 3 phase), and the consistency flag catches those, so I use iDMRG only for dimerization.

## Results

| Model | d | Parameters | Verdict | c estimates | Sizes | Where |
|---|---|---|---|---|---|---|
| Spin-orbital (AB, A, B) on S² | 4 | 302 points (S↔T halved) | no candidate | c_half(8→12) 2.99 to 3.09 near SU(4); max 3.20 at an AB-dominated point with near-degenerate towers, L·gap/2π ≈ 0.02 | 8, 10, 12 | Typhon 17511274 |
| Spin 3/2, Σa_J P_J on S² | 4 | 400 points | no candidate, two residuals | c = 3 region 3.0 to 3.1; of the 24 highest c_half, 8 are rotor towers (S ≥ 2 at k₀ with L·gap/2π ≤ 0.02), 6 dimer doublets, 5 near-ferromagnetic; points 188, 235 have c_half 3.40, 3.20 with z_eff 0.81, 0.65 (gap seemingly opening) | 8, 10, 12; labels at 8, 12 | Typhon 17511274, laptop |
| Spin-3/2 chart, cos θ P₂ + (sin θ ± η)P₁,₃ | 4 | 129 points around the Sp(4) line | c = 3 phase, no TCI | c_half 3.00 to 3.15 inside; 3.48 to 3.79 only where the dimer doublet is forming; iDMRG D = 0.107, 0.106 at χ = 16, 32 at (0.45π, 0.6) | 8, 12; iDMRG χ ≤ 32 | Typhon 17511373, 17511428, 17511429 |
| Sp(4) line, small θ | 4 | 48 angles, 9 at L = 14 | dimerized | doublet ratio → 10⁻⁴ for Lθ ≳ 2.5; iDMRG D = 0.387, 0.385, 0.385 (π/8) and 0.463, 0.455, 0.450 (π/4) at χ = 16, 32, 64 | 6 to 14; χ ≤ 64 | laptop, Typhon 17511258, 17511424 |
| SO(5) vector (P, K) | 5 | 96 angles | no candidate | extended phase c_ent(10) = 4.99 to 5.05, equal to SU(5)_1 (4.99); Reshetikhin point ≈ 2.5 = SO(5)_1; ferromagnetic side z ≈ 1.9 | 6, 8, 10 (12 for 51 angles) | Typhon 17511274 |
| G2 (P₇, P₁₄, P₂₇) on S² | 7 | 300 points | inconclusive | slopes up to 16 with non-scaling gaps; SU(7)_1 itself needs 7 \| L | 5 to 8 | Typhon 17511274 |

The near-SU(4) spin-orbital region behaves as SU(2)_2 × SU(2)_2 (3 + 3 Majoranas, c = 3). I did not add d = 3 scans beyond calibrations, since the spin-1 BLBQ circle is known, the SU(3) bond is a + bP, and D3 covered the U(1) grid.

## Small-angle Sp(4)

At θ = 0 the bond is P₅. The kernel of Σ P₅ is Sym^L(ℂ⁴) plus one Sp(4) singlet at k = π for even L (L ≤ 7, `sp4_theta0.py`), so θ = 0 is the ferromagnetic boundary point, the analog of the SU(3) ferromagnet at θ = −3π/4 of the spin-1 BLBQ chain, next to which that chain dimerizes. The only length scale near θ = 0 is about 1/θ. D3's 0/π crossing obeys θ_c(L)·L ≈ 1.37 at every size. The doublet ratio (k = 0/π singlet splitting over the next gap) collapses in Lθ, giving 0.0092, 0.0109 and 0.0097 at L = 8, 10 and 12 for Lθ ≈ 2.3. For Lθ ≳ 2.5 a clean doublet appears at every angle I could reach (π/16 at L = 12, π/32 already 0.026 at L = 14, π/12 and 3π/8 at 2×10⁻⁴ and 3×10⁻³), with an L-independent gap above it (0.19 at π/12, made of Sp(4) adjoints). The steep entropy slopes (about 8) and z_eff ≈ 2.4 all sit at Lθ < 1.5, inside the crossover window. iDMRG energies match ED (0.20129 against 0.2012 at π/12), and the dimerization is saturated in χ at π/8 and π/4 and large at π/12 (0.30) and 3π/8 (0.49), where only χ = 16 had finished. At π/24, χ ≤ 32 is not yet converged. The runs at π/64, π/48, π/32 and π/16 (job 17511424) were still converging when I stopped, so below π/32 the conclusion rests on the scaling collapse. The whole circle then reads ferromagnet [5π/4, 2π], dimerized (0, π/2), SU(4)_1 at π/2, and c = 3 on (π/2, 5π/4). Nothing exceeds 3.

## Anomalies and the sharpest remaining target

A1 shows that the simple current of Sp(4)_k maps (a, b) to (a, k − a − b), so it preserves the spinor class, and every center-twisted sector is spinor-neutral. I checked this independently. The principal SU(2) gives Sp(4)_k ⊃ SU(2)_{10k}, an even level, which cannot carry the spin-3/2 LSM anomaly. A neutral product factor cannot repair it. The minimal fix with Sp(4) currents is the extension of Sp(4)_1 × M(4,5) by (v, f), with h = ½ + 3/2 = 2 and c = 16/5. The fixed point X = (s, σ′) of weight 3/4 splits into two simple currents carrying the spinor class, and no unitary coset charge lies in (½, 7/10), so 16/5 is the minimum. Replacing tricritical Ising by Ising gives SU(4)_1 back.

This suggests a mechanism. In the c = 3 phase the singlet Majorana's mass term is the Ising energy operator at momentum k₀ + π, which I see numerically as a spin-0 level at k₀ + π with x ≈ 1. At the glued TCI point, the monodromy charge under X is ½ for ε (x = 1/5) and 0 for ε′ (x = 6/5), so ε stays forbidden and ε′ is the single symmetric relevant operator. As in Rahmani et al.'s interacting Majorana chain, one tuning would suffice. With SU(2) only, Sp(4)_1 becomes SU(2)_10 at the same c, so the target is a TCI line in the two-parameter spin-3/2 family with c = 16/5, signalled by the k₀ + π singlet levelling off near x = 1/5 instead of falling to zero. Along the boundary I scanned, this singlet falls continuously to zero (0.001 of the current level at (0.44π, 0.6)), and iDMRG finds stable dimerization, so the transition there is a dimerization onset. The rest of the c = 3 boundary (towards the ferromagnet near θ = 5π/4, and towards lower c at η < −0.25) and larger sizes remain open. This is where I would look next.

## Analytic constructions

- Direct sums. Split ℂ^{d₁+d₂} = V₁ ⊕ V₂ and use a swap-hopping term that exchanges a V₁ and a V₂ state between neighbours, plus bonds h₁, h₂ acting on V₁V₁ and V₂V₂ pairs. At small h_i the Ogata–Shiba factorization gives a free charge mode times two squeezed chains, so c = 1 + c₁ + c₂ with three velocities. The t–J model is the case (2, 1). With c_i ≤ d_i − 1 this saturates d − 1 and never exceeds it, and δ(d) = c_max(d) − (d − 1) is superadditive. It fails to exceed because the squeezed chains stay nearest-neighbour: same-species particles never pass each other, so effective next-nearest squeezed couplings appear only at O(h²/t) and cannot dominate.
- Tensor products give c₁ + c₂ < d₁d₂ − 1. Blocking the known c = 3/2 and 2 range-3 qubit chains lands at d = 4.
- Gauging keeps d per cell by the Gauss law. Kramers–Wannier turns σᶻ_jσᶻ_{j+2} into a nearest-neighbour term but σˣ_jσˣ_{j+1} into a range-3 one, and the c = 3/2 Rahmani chain is self-dual (Majorana translation), so its range cannot shrink.
- Integrable chains in small representations stay far below d − 1: SO(N)_1 (N/2), Sp(2n)_1, G2_1 (14/5 at d = 7), E6_1 (6 at d = 27). Only SU(N)_1 and finite-density graded chains reach d − 1, and the non-compact regimes are non-Hermitian (D1).

## Reproducing

numpy and scipy only. Set `OPENBLAS_NUM_THREADS=1` for ED.

- `ced.py`, `models.py` (Sp(4), spin s, spin-orbital, SO(N) vector, G2 via the octonionic 3-form), `est.py` (estimators).
- `run.py` runs one model at one point, e.g. `python3 run.py --model uls --d 4 --sizes 8 12 --out x.jsonl` (about a minute). L = 14 to 18 needs about 30 GB and an hour on Typhon.
- `scan.py --family {orbital,spin32,so5,g2,sp4}` for the sphere and circle scans, with blocks cached per size. `scan_analyze.py` ranks points.
- `tci_scan.py` and `tci_analyze.py` for the spin-3/2 chart with spin labels. `sp4_theta0.py` and `sp4_analyze.py` for small-angle Sp(4).
- `idmrg.py --model sp4 --theta 0.2618 --chis 16 32 64 --out x.jsonl`.
- `slurm_run.sh`, `slurm_scan.sh`, `slurm_pack.sh` and `cluster/tasks*.txt` give the exact Typhon tasks (partition typhon, 1 to 28 CPUs, 24 to 96 GB). Job IDs: 17511258, 17511274, 17511373, 17511424, 17511425, 17511428, 17511429. 17511322, 17511348 and 17511349 were cancelled (old iDMRG).
- Results are in `data/` as JSON lines, one record per point and size (`data/scan/`, `data/tci/`, `data/idmrg/`).
- A1 ran through `~/.claude/bin/gpt_solve.sh --model gpt-6-astra --effort high` in an empty directory (transcript 20261027T160801-70730, about 22 minutes).
- Still running at hand-off, writing to `/scratch/lustre/jkudlerflam/cdsearch/` on Typhon: SO(5) at L = 12 (17511274, tasks 40 to 43), small-angle Sp(4) iDMRG (17511424, logs `logs/pack_tasks_idmrg1_[0-4].log`) and the spin-3/2 iDMRG at larger χ (17511429).

## Relation to the rest of the repository

This addresses open problem 3 and the Sp(4) loose end of open problem 4 in the top-level README, and REPORT §3.4 and §7. It suggests replacing "Sp(4)_2" by "c = 16/5 (Sp(4)_1 × tricritical Ising)" as the sharp d = 4 target. It is complementary to `contrib/unbounded_c/`, which builds large-d candidates from spontaneously selected cells. The direct-sum argument here is the symmetric, cell-free version and saturates the bound.
