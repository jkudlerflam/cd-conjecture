This run was done by Claude Opus 5.5 (claude-opus-5-5), agent C188, not by Codex. Astra was not used (out of credits).

# C188 result: the cyclic stratum's SU(2)_1/ordered boundary

(Registration first. Results and the verdict are appended below it once the data are in.)

## Registration (written 2026-11-09, before any new numerical output was produced or looked at)

Already seen before this registration, and not counted as a test: C154 RESULT/PROOF §3 and C184's summary of the 53-point grid (job 17594765, `work_C184/c154_cyc_copy/cyc_grid.out`), which says the SU(2)_1 region covers γ/β ≥ −0.3, the triplet collapses along γ/β = −0.6 as β grows, γ/β = −0.9 is ordered, and no grid point shows an SU(2)_3 tower.

Setting. h = α XZ + β YX + γ ZY, mapped exactly (C154 Prop. 2) to the period-3 XYZ chain with bond couplings cycling through (α, β, γ). Fundamental triangle α = 1 ≥ β ≥ |γ|. The sign class of αβγ decides the "side": by site-dependent π rotations (Q_j = Z, X, Y cycling, consistent on rings with 3 | L) (1, β, γ) with γ < 0 is equivalent to −(1, β, |γ|), a ferromagnetic period-3 XYZ chain whose strong axis rotates from bond to bond. The ferromagnetic Heisenberg point is (1, 1, −1). I parametrize the ferromagnetic half of the triangle by β ∈ (0, 1] and t = −γ/β ∈ [0, 1].

Predictions.
- Q1 (structure). The SU(2)_1 region and the ordered (ferromagnet-like, A_4- and Θ-breaking) region are separated by a single curve t_c(β) running from near the XZ corner to the edge β = 1, with t_c(1) between 0.5 and 0.6 and t_c decreasing towards small β (the ordered region widens towards the XZ corner).
- Q2 (order of the transition). Along its whole length away from its endpoints the boundary is a first-order level crossing. On rings L = 12, 18, 24 (and 30 where run) the ground state switches from the SU(2)_1 singlet (momentum 0 or π) to the lowest member of an A_4 multiplet of ordered states, the crossing point drifts by an amount that shrinks with L, and the energy-density slope ∂(E_0/L)/∂t jumps across the crossing by an amount that does not shrink with L. On open chains the half-chain entropy jumps.
- Q3 (no c > 1 tower). No point on or near the boundary shows a tower with r_c = x_current/x_triplet near 10/3 together with a Casimir c above 1.3 at two consecutive sizes (the SU(2)_3 signature). Casimir c and ring-entropy c on the SU(2)_1 side of the boundary stay at or below 1.15 at (18, 24) except where the velocity is very small, where any excess falls with L.
- Q4 (velocity). Approaching the boundary from the SU(2)_1 side the velocity v (from the current level) decreases. Whether it reaches zero at the boundary (a z ≠ 1 point) or stays finite (pure level crossing) I do not predict, but in neither case is there a c > 1 CFT.
- Q5 (dimerized pocket). A dimerized region (r_s < 1, marginally relevant isotropic coupling) exists only near the XZ corner. Its boundary, the g = 0 line of the SU(2)_1 marginal coupling (c = 1), meets the ordered region at or near the XZ corner. Along the rest of the SU(2)_1/ordered boundary r_s stays above 1.
- Q6 (field theory, derived before numerics, see PROOF §2). At an SU(2)_1 point of the cyclic stratum the ferromagnetic order parameter is the conserved current J + J̄, so a direct SU(2)_1-to-ferromagnet transition is either first order or has v → 0 (z ≠ 1). An SU(2)_k point with odd k ≥ 3 has one relevant A_4-invariant, centre-even scalar (the j = 1 singlet, dimension 4/(k+2)), and its two sides are SU(2)_1 and a two-fold dimerized phase, not the ferromagnet. So SU(2)_3 can only sit on an SU(2)_1/dimer boundary, and it can touch the SU(2)_1/ordered boundary only at an isolated junction point.
- Q7 (DMRG). Open chains at the boundary and on the SU(2)_1 side give c_pair(48, 96) at most 1.15, not rising with L. Heisenberg control reads 1.0 within 0.1.

Pass criterion for a c > 1 candidate: the spec's adaptation of C163 PROOF §4.2 (converged open-chain fits above 1.2 not falling with L at truncation below 1e-7, ring-spectroscopy confirmation, controls at c = 1 at the same sizes).

## Verdict

No c > 1 point on the SU(2)_1/ordered boundary, and no candidate bond. The boundary is not a first-order level crossing either, which I had registered (Q2). It is the ferromagnetic instability of the SU(2)_1 liquid. Along every cut tested the SU(2)_1 velocity falls linearly to zero at the boundary, the energy density has a continuous first derivative there (the per-site slope jump at the ring crossing falls like L^-2 or faster), and the gaps at the crossing close like 1/L² (z = 2). That is a non-conformal point of the XXZ-at-Δ = −1 type, with c = 1 on the SU(2)_1 side all the way to it. Field theory (PROOF §2) says why. In this stratum A_4 forbids every relevant and every anisotropic marginal perturbation of SU(2)_1, so the liquid can only end by the isotropic marginal coupling changing sign (dimerization, c = 1 at the line), by v → 0 (the ferromagnetic instability, since the ferromagnetic order parameter is the conserved current), or by a level crossing. A continuous z = 1 transition would itself have to be a c > 1 CFT with one relevant invariant flowing to SU(2)_1 and to the A_4-breaking ordered phase. The only odd-level WZW candidate, SU(2)_3, has a dimerized massive side, so it cannot sit on this boundary. The numerics find the v → 0 route on the whole boundary.

## Key claims

- [proved] Lemma 1 (PROOF §1): a uniform π rotation about σ^x maps (α, β, γ) to (−α, −β, γ), so spec H(1, β, −γ) = −spec H(1, β, γ). Checked to 7e-14 at L = 9, 10.
- [derived] PROOF §2.1 to 2.2: within the cyclic stratum the SU(2)_1 phase can end in only four ways (λ sign change to dimerization; v → 0; level crossing; a z = 1 transition, which must be a CFT with c > 1 by Zamolodchikov).
- [derived; heuristic for the last step] PROOF §2.3: SU(2)_3 is the only odd-k WZW point with a single relevant invariant here, and its massive side is dimerized, so it can lie only on an SU(2)_1/dimer boundary; decoupled SU(2)_1 × X products have a second relevant invariant σ_X·N. I know of no other candidate.
- [numerical, scan_table.out, L = 12, 18 all momenta] Boundary location t_c(β) = −γ/β: 0.565 (β = 1), 0.62 (0.8), 0.66 (0.7), 0.72 (0.5), 0.79 (0.4), 0.82 (0.3), 0.87 (0.2), 0.92 (0.1), from the XZ corner to the AF–FM edge.
- [numerical] v falls linearly to zero at the boundary on every cut (β = 1: 1.66 to 0.18 in equal steps for t = 0.30 to 0.55, zero at t ≈ 0.58).
- [numerical, fine_cross_table.out, local24.jsonl] Ring crossing t_x stable in L (β = 1: 0.5655, 0.5628, 0.5664 at L = 12, 18, 24); slope jump per site 0.114, 0.048, 0.015 to 0.024 (β = 1), 0.082, 0.033, 0.018 (β = 0.7), 0.054, 0.017, 0.0035 (β = 0.4).
- [numerical, gap_at_tx.jsonl, local24.jsonl] L² × current gap at t_x: 14.5, 12.2, 11.5 (β = 1); 12.1, 10.3, 8.3 (β = 0.7); 7.5, 7.3, 5.4 (β = 0.4). z ≥ 2.
- [numerical] No c > 1 tower. c_cas(18, 24) on the SU(2)_1 side near the boundary: 0.92, 0.83, 0.43 (β = 1, t = 0.45, 0.5, 0.55), 0.80 (β = 0.7, t = 0.6), 0.81 (β = 0.4, t = 0.7); Heisenberg control 1.004. L = 30 (converged, job 17605308): β = 1, t = 0.5: E_S = −28.418643, x_t = 1.130, v = 0.445, c_cas(24, 30) = 0.90; t = 0.45: E_S = −29.074229, x_t = 2.051 (2.039, 2.045 at L = 18, 24), c_cas(24, 30) = 0.93 with v(24). r_c crosses 10/3 near β = 1, t ≈ 0.54 while c_cas falls, which is the triplet collapse, not SU(2)_3.
- [numerical, ring_cS.out, xxz_control.out] Ring c_S(12, 18) rises above 1 only within about 0.05 of the boundary (up to 1.37); the exactly c = 1 XXZ chain near Δ = −1 shows the same rise (1.04, 1.11 at Δ = −0.9, −0.95) with c_cas below 1, so this is the z = 2 crossover.
- [numerical, local24_table.out] XZ corner (β ≤ 0.2): c_cas > 1 at small L comes from the gapped-but-slow A Majorana of C154 §3.4 ((c_cas − 1)v constant along cuts) and falls with L: 1.25 → 1.11, 1.52 → 1.17, 3.32 → 0.97 (β = 0.2, t = 0.3, 0.6, 0.8; (12, 18) → (18, 24)); 1.97 → 1.38, 4.13 → 2.13 (β = 0.1). The dimerized pocket (r_s < 1, falling with L) reaches into the ferromagnetic half at β = 0.1.
- [numerical, dmrg_summary.out] Open chains (real τ-frame DMRG, truncation below 1e-8 on the SU(2)_1 side). c_pair(48, 96), (96, 192): Heisenberg control 0.916, 0.934; β = 1, t = 0.5: 0.928, 0.931; t = 0.54: 0.827, 0.862; t = 0.56: 0.52; β = 0.7, t = 0.55: 0.956; t = 0.62: 0.862; β = 0.4, t = 0.65: 1.002; t = 0.75: 0.761, 0.793; β = 0.7, t = 0.62 at (96, 192): 0.884. XZ corner (β = 0.2): t = 0.6: 1.144, 1.057; t = 0.8: 1.322, 0.939, falling with L (the slow A-Majorana crossover; ring c_cas(18, 24) = 0.97 at the same point). Ordered side: bulk |⟨τ⟩| = 0.92 (β = 1, t = 0.7) and 0.84 (β = 0.7, t = 0.8) at L = 96, no logarithmic growth.

## Against the registration

- Q1 (single curve from near XZ to β = 1, t_c(1) in 0.5 to 0.6, t_c decreasing toward small β): the curve is as predicted in shape, but t_c increases toward small β (0.565 at β = 1, 0.92 at β = 0.1), the opposite of what I wrote. The ordered region shrinks toward the XZ corner, hugging the XZ–FM edge. Partly wrong.
- Q2 (first-order level crossing with a slope jump that does not shrink with L): wrong. The slope jump falls like L^-2 or faster, and the gaps close as 1/L². The boundary is a z = 2 point.
- Q3 (no SU(2)_3 tower, c_cas ≤ 1.15 on the SU(2)_1 side except where v is small): holds. c_cas is at or below 1 near the boundary at all sizes up to L = 30. The only excess is at the XZ corner, from the slow A Majorana, and it falls with L.
- Q4 (v decreases toward the boundary, no prediction on whether it reaches zero): v falls linearly to zero at the boundary, the z ≠ 1 option.
- Q5 (dimerized pocket only near the XZ corner, its λ = 0 line meeting the ordered region near the corner, r_s > 1 elsewhere): holds. At β = 0.1 the pocket reaches into the ferromagnetic half (r_s(24) = 0.80 at t = 0.6), and r_s > 1 on every cut with β ≥ 0.2.
- Q6 (field theory): the argument stands, and its case 2 (v → 0) is the one realized. My registered expectation that case 3 (first order) would dominate was wrong.
- Q7 (DMRG c_pair ≤ 1.15, not rising, Heisenberg control near 1.0 within 0.1): holds away from the XZ corner (largest 1.002; control 0.916 and 0.934). At the corner, β = 0.2, c_pair(48, 96) reaches 1.32 at t = 0.8 but falls to 0.94 at (96, 192), so the registered bound is exceeded at the smallest pair and the reading falls with L, as the A-Majorana crossover requires.

## What is left open

- The structure of the ordered side right after the boundary. A narrow window has an L-dependent ground-state momentum before the τ'-ferromagnet takes over. It may be an incommensurate sliver next to the z = 2 point. It does not affect c at the boundary.
- The junction near the XZ corner where the dimer line (λ = 0, c = 1) meets the v = 0 line. Both the B velocity and the A-sector gap are small there, so the sizes used here do not resolve it. It is a codimension-2 point of the stratum with v = 0, so not a CFT, but I did not pin down its structure.
- L = 30 at t = 0.45 and 0.55 (β = 1) and the L = 24 grid were still running in job 17605308 at write-up (collection below). They refine the numbers and are not needed for the verdict.
- The heuristic step of PROOF §2.3: I know of no CFT other than SU(2)_3 that could make the boundary continuous with z = 1. The numerics show that none does along the boundary, which is the question asked.

## Typhon jobs

- Job 17604564 (c188, 24 cores): cancelled after 3 h 20 min. Its DMRG was slowed by an unoptimized einsum in the post-processing, and its L = 30 Lanczos (130 steps) was unconverged. Its only output, /scratch/lustre/jkudlerflam/c18x/c188/out/ring30.jsonl, is saved as typhon_out/ring30_job1_unconverged.jsonl and not used.
- Job 17605308 (c188, 24 cores, 11 h limit, started 2026-11-09 15:14 on typhon-node6). Directory /scratch/lustre/jkudlerflam/c18x/c188/, inputs c188b.slurm, run_task2.sh, tasks_a2.txt, ring30_queue.txt, ringk.c/ringk, ringk_run.py, dmrg3.py, dmrg_tau.py (patched with einsum optimize=True during the run). Outputs in out2/:
  - dmrg_L96.jsonl, dmrg_L48.jsonl, dmrg.log (20 points, three cuts, two corner points, Heisenberg control);
  - ring24.jsonl, ring24.log (8 × 12 grid of L = 24 rings, all 13 momenta);
  - ring30.jsonl, ring30.log (L = 30 at β = 1, t = 0.5 done, then t = 0.45, 0.55, 0.4, sectors 15, 10, 14, 300 steps);
  - slurm_17605308.out in the parent directory.
- Collection, from the cd_conjecture directory:

      mkdir -p "codex/work_C188/typhon_out2"
      scp -r -o LogLevel=ERROR typhon:/scratch/lustre/jkudlerflam/c18x/c188/out2 "codex/work_C188/typhon_out2/"
      scp -o LogLevel=ERROR typhon:/scratch/lustre/jkudlerflam/c18x/c188/slurm_17605308.out "codex/work_C188/typhon_out2/"
      cd codex/work_C188
      python3 summarize_dmrg.py typhon_out2/out2 dmrg_local_L48.jsonl dmrg_local_L96.jsonl dmrg_local_L192.jsonl
      python3 analyze_scan.py scan_L12_18.jsonl scan_L12_18_b.jsonl typhon_out2/out2/ring24.jsonl

  (the scp -r puts out2 inside typhon_out2/, so the paths become typhon_out2/out2/...). The files at the top level of typhon_out2/ are an interim snapshot taken at about 19:30, which is what the numbers above use. Check the ring24 grid's c_cas(18, 24) against Q3 and the corner numbers of PROOF §5, and add the L = 30 points to PROOF §3.7 (c_cas(24, 30) and x_t at t = 0.45 and 0.55).

Files in this folder: RESULT.md, PROOF.md, PROGRESS.md; code ringk.c, ringk_run.py, analyze_scan.py, analyze_cross.py, ring_cS.py, xxz_control.py, xxz_open.py, dmrg3.py (C154's dmrg.py plus a prebuilt-MPO option), dmrg_tau.py, summarize_dmrg.py, scan_local*.sh, fine_cross.sh, local24*.sh, local96.sh, typhon/; data scan_L12_18*.jsonl, scan_table.out, fine_cross.jsonl, fine_cross_table.out, gap_at_tx.jsonl, local24.jsonl, local24_table.out, ring_cS.jsonl/.out, xxz_control.jsonl/.out, xxz_open.jsonl, dmrg_local_L48/L96/L192.jsonl, dmrg_summary.out, typhon_out/, typhon_out2/.
