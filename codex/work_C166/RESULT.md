This run was done by a Claude Opus 5.5 agent (C166, 2026-11-08), not by Codex. Astra was not used (out of credits).

# C166 result: open-chain DMRG against infinite-chain VUMPS in the XZ neighbourhood

## Pass criterion (fixed before any new data was looked at)

Written at the start of the run, before the first local check. It is the criterion of the spec, task 6, made numerical.

A bond is a c > 1 candidate only if all three of the following hold.
- (C1) Open chains. c_prof > 1.2 at L = 256 and at L = 512, with c_prof(512) ≥ c_prof(256) − 0.03, the two-size slope 6 [S_mid(512) − S_mid(256)] / ln 2 above 1.2, and maximal discarded weight below 1e-8 in the last sweep at both lengths.
- (C2) VUMPS. Take the unit cell n (among all cells and seeds tried) whose variational energy per site at the largest χ is lowest. Its FES slope (least squares 6 dS/d ln ξ over the last three χ values, χ_max ≥ 256) must exceed 1.2.
- (C3) Energies. The open-chain bulk bond energy (mean of ⟨h_{j,j+1}⟩ over the central L/4 bonds, which spans several incommensurate periods) at L = 512 agrees with the lowest VUMPS e₀ within the open-chain finite-size allowance. The allowance is set from the data as |e_mid(512) − e_mid(256)|, i.e. the change from the previous size; agreement means |e_mid(512) − e₀| ≤ 2 |e_mid(512) − e_mid(256)| + 1e-6.

If C3 fails with e_mid above e₀, the open chain misses the VUMPS state. If it fails with e_mid below e₀ (beyond the χ-error of e₀, taken as |e₀(χ_max) − e₀(χ_max/2)|), VUMPS is trapped.

Decision rules for the hypotheses of the spec.
- (iv) bug or convention mismatch: excluded if open-chain DMRG matches ED to 1e-8 at L ≤ 16 for every bond tested, and the VUMPS e₀ matches the L → ∞ extrapolation of periodic ED within the extrapolation error.
- (i) weak or strong time-reversal breaking with a symmetric cat on the open chain: supported if the lowest VUMPS state has ⟨σ⟩ ≠ 0 (more than 1e-3 per site), and the open chain shows long-range order in ⟨σ^a_{L/4} σ^b_{3L/4}⟩ at the VUMPS value m_a m_b within 20 percent, and a pinning field of order 0.1 at the edges removes the growth of S_mid with L (S_mid(256) − S_mid(128) below 0.05 with the field, against 0.3 without).
- (ii) VUMPS trapped: supported if some cell n > 2 gives e₀ below the two-site e₀ by more than the χ-error, with an FES slope above 0.5.
- (iii) crossover beyond 256 sites: supported if none of (i), (ii), (iv) holds, VUMPS stays gapped in every cell with ξ converging in χ, and the open-chain slopes fall between 256 and 512.

## Prediction registered before the Typhon data (added after the local checks, before job 17602270 started)

The local checks below point to hypothesis (i) in a specific form, a Néel cat whose ln 2 builds up slowly. If that is right, the unpinned S_mid(L) approaches S_VUMPS + ln 2 from below, and the existing sizes 64 to 256 fit a deficit ∝ exp(−L/ℓ) with ℓ ≈ 187 (hi_T_40) and ℓ ≈ 134 (dist_T_50) [numerical, fit to `inputs/ref_open.json`]. The prediction for the open chains of job 17602270 is then
- hi_T_40: S_mid(384) ≈ 1.31, S_mid(512) ≈ 1.35, two-size slope 256 → 512 ≈ 1.08, bounded by S_VUMPS + ln 2 = 1.394;
- dist_T_50: S_mid(384) ≈ 1.34, S_mid(512) ≈ 1.36, slope 256 → 512 ≈ 0.73, bounded by 1.377.
A critical bulk with c ≈ 1.3 instead predicts S_mid(512) ≈ S_mid(256) + (1.3/6) ln 2 = 1.38 and 1.43, the second above the cat bound. S_mid(512) above S_VUMPS + ln 2 + 0.02 for either bond falsifies the cat explanation.

## What is decided by the local checks (before Typhon)

Files: ed_check.out, local/scan.log, local/vc/*.json, local/oc/*.json, local/summary_local.txt (summarize.py on local/).

Task 1, code consistency [numerical, ed_check.out]. For all seven bonds (hi_T_40, dist_T_50, cand_60, cand_61, KI22, KI26, N34) the two-site operator encoded by fesx.mpo_cell equals est.bond_from_JB(J, f) plus the other half of the field to 1e-16. C163's fes/T40.json and fes/T50.json used J = J_f(v) of the same bonds (to 5e-13, B = 0). Open-chain DMRG equals sparse ED to 7e-13 at L = 8, 12, 16 for every bond. vumps_cells.py reproduces C163's χ = 16 VUMPS energies to 1e-10. Periodic ED at even L = 8 to 20 oscillates by 1.5e-2 (incommensurate correlations on short rings) and cannot pin e₀ better than that, so the energy check on e₀ is the open-chain one (C163's two-size energy slopes approach e₀ from above, and the pinned chain below reproduces e₀ to 4e-6 in its middle). Hypothesis (iv) is excluded: both codes solve the same Hamiltonian.

Task 2 at small χ [numerical, local/scan.log]. VUMPS at χ = 8 to 32, two seeds, cells n = 1 to 8, 10 to 14, 20, 21:
- Every even cell (2, 4, 6, 8, 10, 12, 14, 20) converges to one state, with energy equal to the two-site value to 1e-12 (hi_T_40 −0.959308071976, dist_T_50 −0.837123110661 at χ = 32), the same S on every bond (0.6904, 0.6764), and ξ (22.7, 18.4). Larger or incommensurate-friendly cells find nothing lower.
- That state is Néel ordered and breaks time reversal: ⟨σ⟩ alternates between sites with |m| = 0.345 (hi_T_40, along (0.020, 0.001, 0.345)) and 0.366 (dist_T_50, along (0.000, −0.072, −0.359)). The TR mixed-transfer radius is below 1 (0.971 and 0.964 per two sites). Bond energies are uniform within the cell (spread 4e-12). The same staggered moment is already in C163's χ = 192 files (fes/T40.json, `mag`), where it had not been looked at.
- Odd cells cannot hold a Néel state. They return non-injective cats (ξ → ∞, TR radius 1) at higher energy, for instance n = 11 at χ = 32 gives the n = 2 energy at χ = 16. They are excluded from the energy comparison.
So hypothesis (ii) gets no support at χ ≤ 32. Typhon extends it to χ = 512 (n = 2), 384 (n = 4) and 256 (n = 6, 8), plus n = 10 and 20 up to χ = 192.

Tasks 3 and 4 at L ≤ 128 [numerical, local/oc/, local/summary_local.txt]. hi_T_40:
- No field. S(1) = ln 2 and ⟨σ_j⟩ = 0 at every site, so the open chain is the TR-symmetric state. The order-parameter correlation (n·C·n)(−1)^r between L/4 and 3L/4 is 0.0758 at L = 64 (r = 32) and 0.0731 at L = 128 (r = 64), against m² = 0.119. It is nearly flat in r. A c ≈ 1.3 critical state would make it fall as a power of L. The unpinned chain is therefore an ordered superposition of the two Néel sectors with a reduced order amplitude. Its middle bond energy sits above e₀ by 1.3e-5 (L = 64) and 1.6e-5 (L = 128). That excess does not fall like 1/L² (the CFT form would drop it by about 4 between these sizes), so the middle of the unpinned chain is not yet in the bulk ground state.
- Néel-compatible fields 0.2 at both ends (LRst). S_mid = 0.671 (L = 64) and 0.690 (L = 128), approaching the VUMPS S = 0.7007, with c_prof 0.48 → 0.29. The central staggered moment is (−0.0203, −0.0011, −0.3450), the VUMPS vector to 1e-3, and emid − e₀ = −4e-6. Pinning removes the growth of S (0.019 per doubling, against 0.140 without the field) and turns the open chain into the VUMPS state.
- A field 0.2 at the left end along x, nearly perpendicular to the order parameter (Lx), leaves the profile almost unchanged (S_mid 0.918 and 1.039 against 0.926 and 1.066). It pins only through the small x part of the order vector.
- dist_T_50 at L = 64 behaves the same way: S_mid 0.957 without field and 0.653 with LRst, central moment (−0.0005, 0.0737, 0.3674) against the VUMPS vector (0.0005, −0.0724, −0.3591) up to the global sign, and order correlation 0.080 against m² = 0.134.

Reading [numerical for each statement above, heuristic for the mechanism]. The bulk is a gapped, TR-breaking Néel state, and VUMPS describes it correctly. The open chain without edge fields sits in the TR-symmetric superposition of the two Néel sectors. Its extra entanglement over the Néel value builds up towards ln 2 only over a length ℓ much larger than the VUMPS ξ at χ = 192 (the fit above gives ℓ ≈ 130 to 190). The rising c_prof and two-size slopes are this crossover, not bulk criticality. The mechanism behind the long ℓ is heuristic. Domain walls between the two Néel sectors appear to be cheap (incommensurate transfer eigenvalues, ξ still growing in χ while S has converged), and a free edge can emit a single wall, which a pinned edge cannot. This is hypothesis (i), with order that is strong rather than weak.

## Waiting on Typhon (job 17602270, PENDING at submission)

Directory /scratch/lustre/jkudlerflam/c166/. Inputs typhon/c166.slurm, typhon/tasks_heavy.txt, typhon/tasks_light.txt. The job is one 48-core node with a 40 h limit. Expected outputs:
- oc/{hi_T_40,dist_T_50,cand_60,cand_61}_{128,256,512}.json and oc/{hi_T_40,dist_T_50}_384.json: open chains with observables. L = 512 uses χ up to 320 (XZ bonds) or 384 (C155 candidates).
- oc/*_{LRst,Lst,Lx}_{128,256}.json: pinning for hi_T_40 and dist_T_50. oc/cand_6{0,1}_LRst_{128,256}.json for the candidates.
- vc/*.json: VUMPS for hi_T_40, dist_T_50, cand_60, cand_61 (n = 2, χ ≤ 512), for hi_T_40 and dist_T_50 (n = 4, χ ≤ 384; n = 6, 8, χ ≤ 256; n = 2 seed 1, χ ≤ 256), hi_T_40 n = 10 and dist_T_50 n = 20 (χ ≤ 192), cand_60 and cand_61 n = 4 (χ ≤ 256), and KI22, KI26, N34 (n = 2, χ ≤ 384; n = 4, χ ≤ 256; task 5).
- summary.txt, written by summarize.py at the end of the job.

To collect: `scp -r typhon:/scratch/lustre/jkudlerflam/c166/{oc,vc,logs,summary.txt} typhon_out/` (inside codex/work_C166, leaving out vc/*.npz). Then run `python3 summarize.py typhon_out > typhon_out/summary.txt`. Section (4) of the output evaluates C1 to C3 for each XZ bond and the FES second stage for KI22, KI26 and N34. Compare S_mid(384) and S_mid(512) with the prediction above. Check that the pinned L = 256 chains keep S_mid near S_VUMPS, and that the unpinned order correlation stays flat in L.

Local hints for task 5 [numerical, χ = 32 only, local/vc/*_quick.json]. N34 is a near product state at χ = 32 (S = 0.013, ξ = 1.1, |m| = 0.78). This suggests that its first-stage c_prof(128) = 1.68 is a level-crossing artefact of the open chain, which fits C163's remark that all three passes sit on a jump in S_mid(λ). VUMPS at χ = 32 could still be trapped, so this waits on the χ ≤ 384 runs. KI22 (S = 1.14, ξ = 52) and KI26 (S = 1.17, ξ = 99, ξ tripling from χ = 16 to 32) are still growing at χ = 32 and need the Typhon FES. cand_60 and cand_61 at χ = 32 are dimerized (bond energies −0.742 / −0.808 and −0.722 / −0.818) with a uniform x moment and a staggered z moment (0.076, 0.14), TR radius 0.85 and 0.79, ξ 22 and 13.

## Verdict so far

No c > 1 bond [numerical, pending the L = 512 and χ = 512 tests]. In the XZ neighbourhood the conflict is explained by hypothesis (i). VUMPS is right about the bulk (gapped Néel, TR broken, FES slope about 0.01 at χ ≥ 96, the same state in every even cell). The open-chain c ≈ 1.3 is the slow approach of the TR-symmetric open chain to the Néel cat, and an edge field along the VUMPS order parameter removes it (local L = 64 and 128). C1 to C3 cannot be evaluated until the Typhon outputs arrive. The pass criterion and the prediction above were fixed before them.

