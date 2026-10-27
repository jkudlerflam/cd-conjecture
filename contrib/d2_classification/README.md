# The D2 qubit candidates: none survives as c > 1

## Claim

Class: periodic, one-site translation-invariant, nearest-neighbour Hermitian qubit chains, H = Σ_j [Σ_ab J_ab σ^a_j σ^b_{j+1} + h·σ_j]. J is a general real 3×3 matrix (XYZ anisotropy plus a uniform DM vector) and h is a uniform field, giving 8 parameters up to rotation and scale. The central charge is the total Virasoro c of a unitary CFT.

Tested: the 16 points that D2 (`codex/work_C4_D2`) followed to L = 22 with entropy slopes 1.0 to 1.7, and the points where Ising lines meet c = 1 regions in the structured sub-families.

Conclusion: none is a CFT with c > 1.

- Ten are c = 1 incommensurate Luttinger liquids whose small-L entropy slopes overshoot.
- Three show no logarithmic entanglement growth from L = 48 to 96 and are most likely gapped.
- Three sit near a canted-ferromagnetic boundary, with modulated, non-logarithmic entanglement.
- The Ising/c = 1 junction of the XXZ chain in a transverse field has c ≈ 1.
- The literature has no strictly NN, one-site-cell qubit chain with c > 1.

## Status

Numerical, plus literature. It would be upgraded by finite-entanglement scaling (iDMRG/VUMPS, χ ~ 10³) on the six points in classes (c) and (d).

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), as a subagent of the coordinating session, 2026-10-27. The code here is new and independent of `d2.py`. No Codex/GPT was used.
- Literature: read by a second Claude Opus 5.5 subagent. Tags: [read] means it opened the paper at the place cited, [recalled] means memory or a search snippet.
- `ed.py selftest` matches full-space diagonalization in every momentum sector to 10⁻¹³, for random two- and three-site terms.
- DMRG reproduces the exact open XX ground energy to 6·10⁻¹⁰.
- Both estimators were calibrated on the controls below. Nothing was rerun independently.

## Method

All 16 candidates have a large DM vector (|D| ≈ 0.53 to 0.64, against symmetric couplings of about 0.4). Classical minimization (`classical.py`) gives incommensurate spirals (q ≈ 0.19 to 0.36) for twelve of them. For 269, 1289, 871 and b107 it gives canted ferromagnets (q = 0). The quantum energy lies only 3% to 29% below the classical one.

An incommensurate spiral has one sliding mode, so the generic expectation is a c = 1 floating phase with a large Luttinger parameter. On finite rings this shows up as a dense "rotor" tower at incommensurate momenta, with ground-state momenta that jump irregularly with L. Single-gap, velocity and Casimir estimators then fail. The Casimir estimator gives −5.2, −3.9 and −5.5 on the twisted-XXZ and cluster-Ising controls, so it is not used for verdicts. `frame.py` makes the spiral explicit by a site-dependent rotation (Δ_eff ≈ −0.3 to −0.8), but the oscillating remainder is 15% to 60% of the bond, so this is only a heuristic.

Estimators (`analyze.py`, `make_table.py`):

1. c_S: periodic-ring fit S(l) = (c/3) ln[(L/π) sin(πl/L)] + b, using blocks l ≥ L/4 only. The data come from momentum-resolved ED at L ≤ 24 (25 for three points), with all L sectors and the four lowest levels in each.
2. c_half: slope of S(L/2) against (1/3) ln L over the last four even L.
3. z10, z20: from ln ε_n ∝ −z ln L over L ≥ 12, where ε_n is the n-th excitation energy among all sectors. This is robust to ground-state level crossings.
4. c_mid = 6[S̄(96) − S̄(48)]/ln 2, from open-chain two-site DMRG (`dmrg.py`, χ ≤ 128 or 192, lowest energy over seeds). S̄ is the mean of S(l) over l ∈ [L/4, 3L/4]. For a CFT its L-dependence is exactly (c/6) ln L, and the window averages out incommensurate oscillations.
5. c_fit(L): open-chain fit over the same window, with the bond-energy correction of Laflorencie et al.

## Calibration

| control (true c) | c_S, L = 20 to 24 | c_half | z20 | c_mid, c_fit(96) |
|---|---:|---:|---:|---:|
| Ising (1/2) | 0.50 | 0.50 | 1.0 | 0.52, 0.51 |
| XX (1) | 1.00 | 0.99 | 0.9 | 0.99, 0.98 |
| Heisenberg (1) | 1.01 | 1.00 | 0.8 | - |
| XXZ Δ = 0.5, −0.5, −0.8 (1) | 1.00 to 1.01 | 1.00 | 0.9 to 1.0 | - |
| XXZ Δ = −0.95 (1) | 1.06 | 1.06 | 0.8 | 1.01, 0.99 |
| XXZ Δ = −0.99, L ≤ 20 (1) | 1.20 | 1.28 | 1.4 | c_fit(48) 1.05 |
| XXZ Δ = −0.999, L ≤ 20 (1) | 1.26 | 1.39 | 1.9 | - |
| XXZ + D_z, twisted (1) | 1.00 | 0.99 | 0.9 | - |
| FM-XY + D_z + h_z, twisted (1) | 1.01 | 1.00 | 1.0 | 0.98, 0.96 |
| cluster-Ising, three-site (3/2) | 1.52 ± 0.03 | 1.48 | 0.9 | not NN |

c_S separates 3/2 from 1 away from the ferromagnet. Near it, c_S overshoots by 0.2 to 0.4 at these sizes while c is exactly 1, and z20 rises toward 2. That is the von Neumann ½ ln m law at Δ = −1. So c_S ≈ 1.2 to 1.4 at L ≤ 24 cannot tell c = 3/2 from a near-ferromagnetic c = 1 chain or a gapped chain with ξ ~ L. The open chain at L = 48 to 96 can.

## Results

The ± on c_S is the spread over the three largest even L. Drift is c_S(Lmax) − c_S(Lmax − 4). z errors are regression errors. "-" means the run was not done or not copied back (see Reproducing).

| point | ED Lmax | c_S | drift | c_half | z10 | z20 | c_mid | c_fit 96 [48] | verdict |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| 1256 | 25 | 1.33 ± 0.12 | −0.11 | 1.08 | 1.2 ± 0.1 | 1.7 ± 0.1 | 0.96 | 0.97 [1.18] | (b) c = 1 |
| opt16_121 | 25 | 1.39 ± 0.41 | −0.41 | 1.21 | 0.8 ± 0.1 | 1.1 ± 0.1 | 1.07 | 1.10 [1.09] | (b) c = 1 |
| 830 | 25 | 1.24 ± 0.14 | −0.14 | 0.99 | 1.1 ± 0.3 | 1.8 ± 0.2 | - | 0.76 [-] | (b) c = 1 |
| 1282 | 24 | 1.22 ± 0.08 | −0.08 | 1.19 | 0.9 ± 0.1 | 1.1 ± 0.1 | - | - [1.07] | (b) c = 1 |
| 1144 | 24 | 1.17 ± 0.19 | −0.19 | 0.91 | 0.8 ± 0.2 | 1.0 ± 0.1 | - | - [1.05] | (b) c = 1 |
| 854 | 24 | 1.16 ± 0.16 | −0.04 | 1.28 | 0.8 ± 0.3 | 0.7 ± 0.3 | 1.12 | 0.90 [0.86] | (b) c = 1 |
| 220 | 24 | 1.12 ± 0.07 | −0.06 | 1.10 | 1.1 ± 0.3 | 1.3 ± 0.1 | 0.89 | 0.99 [1.02] | (b) c = 1 |
| b107 | 24 | 1.10 ± 0.12 | +0.09 | 1.15 | 1.1 ± 0.1 | 1.5 ± 0.1 | 0.99 | 0.99 [0.99] | (b) c = 1 |
| 1057 | 24 | 1.06 ± 0.11 | +0.10 | 1.22 | 0.9 ± 0.1 | 1.0 ± 0.0 | 0.76 | 0.86 [0.93] | (b) c = 1 |
| 1004 | 24 | 1.01 ± 0.01 | 0.00 | 1.01 | 0.9 ± 0.1 | 0.9 ± 0.2 | 0.82 | 0.90 [0.95] | (b) c = 1 |
| 323 | 24 | 1.42 ± 0.05 | −0.05 | 0.22 | 1.5 ± 0.1 | 1.5 ± 0.0 | 0.32 | 0.62 [0.94] | (d) gapped |
| opt68 | 24 | 1.07 ± 0.02 | +0.01 | 0.38 | 0.7 ± 0.2 | 0.9 ± 0.2 | −1.66 | 0.17 [0.93] | (d) gapped |
| 628 | 24 | 1.05 ± 0.03 | −0.02 | 1.04 | 0.7 ± 0.2 | 0.7 ± 0.2 | −0.74 | 0.47 [0.31] | (d) gapped |
| 269 | 24 | 1.34 ± 0.02 | +0.01 | 0.58 | 2.5 ± 0.1 | 2.8 ± 0.1 | 0.86 | modulated | (c) |
| 1289 | 24 | 1.24 ± 0.19 | −0.18 | 1.37 | 1.9 ± 0.2 | 2.6 ± 0.1 | - | modulated | (c) |
| 871 | 24 | 1.12 ± 0.08 | −0.05 | 1.00 | 0.8 ± 0.1 | 1.3 ± 0.1 | - | modulated | (c) or (b) |

**(b) c = 1, misestimated at small L.** The ring slope drifts down with L, and the open chain gives c between 0.8 and 1.1. For 1256 the L = 96 run is converged at χ = 192, with two seeds agreeing to 10⁻¹². Values a little below 1 reflect χ ≤ 128 and residual oscillations, as for the controls. The D2 raw maximum (830, c_eff = 2.70 at L = 10) gives c ≈ 0.8 on the open chain at L = 96.

**(d) Gapped.** At L = 96 the open-chain entropy is flat over about 60 central sites, and its window average does not grow from L = 48 (it drops). Truncation errors at χ = 128 are 10⁻¹⁰ to 10⁻¹⁵. For 323 and opt68 the ring half-chain entropy also saturates (c_half 0.22 and 0.38). The large in-ring c_S of 323 is the shape of S(l) at ξ ~ L, not logarithmic growth. For opt68 the boundary energies at L = 48 and 96 agree (0.41, 0.37), so the L = 96 state does not look like a metastable branch. Confidence is moderate: a c = 1 phase with a pinned incommensurate boundary cannot be fully excluded at L = 96.

**(c) Near a canted-ferromagnetic edge.** These three are classically canted ferromagnets (|m| = 0.5 to 1). Their open-chain S(l) is a non-logarithmic modulated pattern between 0.1 and 1.0. In 269 and 1289, ε_n collapses faster than 1/L (z ≈ 2 to 2.8), which points to a commensurate-incommensurate (z = 2) neighbourhood. Whether they are critical is unresolved. None shows growth compatible with c > 1 (269: c_mid = 0.86).

## Multicritical points

**XXZ in a transverse field.** H = Σ(SˣSˣ + SʸSʸ + ΔSᶻSᶻ) − hΣSˣ, coded as J = diag(1, 1, Δ) with field (−2h, 0, 0) (`junction.py`, `results/junction.json`). h_c comes from crossings of L·[E0(π) − E0(0)] at L = 12 and 16. The c values are at L = 20 (c_S) and from L = 12 to 20 fits (c_E).

| Δ | h_c | c_S | c_E | L·gap |
|---:|---:|---:|---:|---:|
| −0.50 | 1.126 | 0.502 | 0.497 | 1.02 |
| −0.80 | 0.875 | 0.503 | 0.502 | 1.15 |
| −0.90 | 0.765 | 0.515 | 0.481 | 1.17 |
| −0.95 | 0.696 | 0.570 | 0.481 | 1.23 |
| −0.99 | 0.609 | 0.846 | 0.774 | 1.22 |
| −1.00 | 0.502 | 1.028 | 0.962 | 0.79 |

The Ising line ends on the U(1) segment Δ = −1 near the KT point h_0 ≈ 0.52 of Dmitriev et al. The apparent c crosses over from 1/2 to 1, and the junction has c ≈ 1 with z = 1. [numerical]

**XYZ in a longitudinal field.** The Ising lines end on the U(1) line at the free-fermion band edge, which is the saturation transition of an XXZ chain, z = 2 exactly. The c = 1 floating lobes meet only at an eight-vertex c = 1 point. (Chepiga–Mila 2206.11754, Eq. 3, Fig. 1 [read].)

**XXZ + D ∥ z + transverse h.** The triple point of the Luttinger liquid and two Néel phases is bounded by KT separatrices of a single boson, so c = 1 (Chan–Jin–Jiang–Starykh 1708.03670, Eqs. 15 and 16, Fig. 4 [read]). Not recomputed, since D/J = 0.1 needs L ~ 10³.

**Generic large-DM region.** This is where all 16 candidates live, and it sorts into a floating spiral phase, gapped pockets, and canted-ferro edges. I found no point where an Ising line runs inside the floating phase, which is the mechanism behind c = 3/2 in the three-site models below.

## Literature

- No strictly NN, one-site-cell spin-1/2 chain with c > 1 was found. Every c = 3/2 example has a three-site term.
  - Rahmani–Zhu–Franz–Affleck 1505.03966, Eq. 7, Figs. 1 and 5: Ising + Luttinger for −2.86 < g/t < −0.285 [read].
  - Chepiga–Laflorencie 2211.15598, Eq. 1, Fig. 2: c = 3/2 only where a self-dual Ising line runs inside the floating phase [read]. Dropping σˣ_jσˣ_{j+2} makes the model NN but removes the self-duality (our inference).
  - Cluster-Ising has c = 3/2 [recalled, Lahtinen–Ardonne 1504.07044]. 2608.24248 finds that an NN term splits it into Ising and Gaussian lines [read, abstract].
  - O'Brien–Fendley 1712.06662 makes no c = 3/2 claim [read].
- XYZ in a field.
  - Chepiga–Mila 2206.11754 and Verresen–Vishwanath–Pollmann 1903.09179 (c = 1.000 ± 0.002) [read].
  - Dmitriev–Krivnov–Ovchinnikov cond-mat/0112145: Ising line, and c = 1 on Δ = −1 up to h_0 ≈ 0.52 [read].
- DM + anisotropy + field.
  - Garate–Affleck 1001.4290 (abstract, Sec. I) and Chan et al. 1708.03670 find a Luttinger liquid and two Néel phases, with nothing above c = 1 [read].
  - Free XY + D_z + h_z has a chiral c = 1 phase [recalled].
  - Free NN fermion chains have c ≤ 1 (lit.md [own]).
- KHΓ. The uniform spin-1/2 KHΓ chain is rotation-equivalent to XYZ. Yang–Nocera–Affleck (PRR 2, 033268) study the bond-alternating chain, which keeps a three-site pattern even after U6, and find only c = 1 phases [read].
- Ferromagnetic boundary.
  - Castro-Alvaredo–Doyon 1011.4706, Prop. 1, Eq. 5.73: S = ½ ln m + O(1) as Δ → −1⁺ [read].
  - Popkov–Salerno quant-ph/0404026 [read, abstract].
  - This explains the calibration overshoot.
- Latorre–Sierra 2403.17242 checked only Ising and XXZ at d = 2 [read].

## Reproducing

numba, numpy and scipy are needed. On Typhon, use `PYTHONNOUSERSITE=1` with the default anaconda.

- `python3 ed.py` runs the self-test.
- `python3 run_point.py NAME --sizes 8 ... 24 [--procs P]` runs ED. NAME is a key of `candidates_params.json` or `controls.py`, and the output is `results/NAME.json` with the spectra of all sectors and the ground-state S(l). L = 24 takes 5 to 15 min on 8 cores.
- `python3 dmrg.py NAME --L 96 --chi 32 64 128 192` runs DMRG, 5 to 25 min on 4 cores.
- `python3 junction.py -0.5 -0.8 -0.9 -0.95 -0.99 -1.0` runs the junction scan, about 3 min.
- `classical.py` and `frame.py` do the spiral analysis.
- `merge_results.py` then `make_table.py` rebuild `results/summary.json`. The laptop DMRG is in `results_dmrg_local/`.

Typhon jobs (partition typhon, /scratch/lustre/jkudlerflam/d2_classification), of the order of 100 core-hours in total:

- 17511297: ED to L = 26, 24 cpus and 180 GB. Three tasks ran.
- 17511301: ED to L = 24, 8 cpus and 30 GB, 24 tasks.
- 17511318, 17511341, 17511376, 17511403: DMRG on 4 cpus, the last two under qos short.

The laptop used about 45 min of wall time.

When this was written the jump host was refusing connections. The last outputs are still on Typhon: L = 96 DMRG for 1144, 1282, 1289 and 871, L = 48 for 830, and ED at L = 25 and 26 for 1256, opt16_121 and 830. Rsync them into `results_cluster/` and rerun the two scripts to fill the "-" entries.

## Relation to the rest of the repository

This addresses the qubit half of open problem 4 (README; REPORT.md §7, D2). The slopes of 1.35 to 1.7 are finite-size effects of near-classical spirals, gapped pockets with ξ ~ L, and canted-ferromagnetic edges. They are not c > 1. That supports "no counterexample at d = 2" without proving it. The heuristic that emerges is that a spin-1/2 spiral has a single magnon branch, and every known c = 3/2 at d = 2 needs a three-site term to put an Ising line inside a floating phase.

## Coordinator's note (2026-10-27, after this agent finished)

The coordinating session retrieved the Typhon outputs that this agent could not copy back. It used the `typhon` ssh alias via the IAS bastion; the `ias-jump` alias does not exist on this machine. They are now in `results_cluster/`. Windowed open-chain fits on them agree with the table to two digits.

A second estimator flags two points. The change of the central-cut entropy from L = 48 to L = 96, 6ΔS/ln 2, calibrates to 0.51 (Ising), 0.96 (XX) and 1.01 to 1.08 (XXZ), but gives 1.93 for point 220 and 2.16 for point 1282. Their profiles are not clean logarithms:
- 220 oscillates with a period of about L/8 at L = 48 and rises almost linearly at L = 96.
- 1282 has a plateau at L = 48 and a rising profile at L = 96.

Both are consistent with an incommensurate spiral pinned by the open boundaries and a crossover length comparable to L, which would make the doubling estimate unreliable. But the table's "c = 1" verdict for these two points should be read as provisional. A follow-up using finite-entanglement scaling with iDMRG (S against log ξ as χ grows) is being run. Its result will be added here.
