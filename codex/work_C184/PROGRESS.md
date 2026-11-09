# C184 progress log (Claude Opus 5.5)

2026-11-09
- Read the spec in full, STATUS, _qubit_c_context, REPORT §9.126 to §9.139, C169 RESULT, C154 PROOF 1, 4 to 9, C174 REDTEAM, C152 RESULT.
- Found that C152's Typhon job 17597010 had completed on 2026-11-08 and was never collected. Copied its outputs (read only, from /scratch/lustre/jkudlerflam/c152/j2/) into c152_j2_copy/ and summarized them with summ_c152j2.py, giving c152j2_summary.out. The XZ window falls toward 1 (1.46, 1.24, 1.13 at the worst point). The Peotta region and strong attraction show nothing above 1.15. See PROOF §6.
- Literature opened: arXiv:1403.4568 (abstract), arXiv:1705.00635 (abstract). Search snippets only, not opened: 2211.15598, 1012.4518, 2510.05988 (J1-J2, range 2).
- Wrote the inventory, PROOF §1 to §5. Registered predictions in RESULT.md before any new numerics.
- Ring ED map at φ = 1, h_z = 0, L = 12 and 16 (ring_map.py, ring_phi1.jsonl, ring_phi1_table.out). Phases found: z-FM at V ≲ −2 to −2.6 (growing with Δ), a Luttinger region at small Δ for −2 < V < −1, an x-ordered phase at Δ ≥ 0.5 and V ≥ −1.4, and a high-entropy band along V ≈ −2(cos 1 + Δ), which is the |J_x| = |J_z| line where Peotta et al.'s DM-split Ising pair runs. On rings the band reads ΔS(12 → 16) ≈ 0.09, close to the c = 1 value 0.096.
- Estimator controls from exact free fermions (ff_controls.py, ff_controls.out). c_pair(48,96), (96,192): XX 0.989, 0.993; TFI 0.515, 0.507; XX×XX 1.977, 1.987; TFI×XX 1.504, 1.501. The d = 4 DMRG ladder controls were too slow locally, so they were replaced by these exact values. The qubit XX and TFI DMRG controls are in the Typhon job.
- Typhon job 17604462 (c184), submitted 2026-11-09. It has 418 DMRG tasks: φ = 1 at Δ = 0.3 to 0.7 and V = −2.8 to −1.0 at L = 48 and 96, L = 192 on the coarse grid at Δ = 0.4 to 0.7, φ = 1.2 at L = 48 and 96, and XX and TFI controls at 48, 96 and 192. Directory /scratch/lustre/jkudlerflam/c18x/c184/, outputs out/dmrg_{48,96,192}.jsonl, out/ctrl_{48,96,192}.jsonl, out/summary.txt (summ.py), out/err.log, slurm_17604462.out.
- Local preview at L = 32 and 64 (local_preview.sh, local_out/) on six points of the band.
- Local preview stopped after one L = 64 point (c_pair(32, 64) = 1.08 at Δ = 0.5, V = −2.4), because the shared CPU made it exceed the 20-minute local limit.
- Collected C154's uncollected job 17594765 (the 53-point cyclic grid) into c154_cyc_copy/. Read it in PROOF §6c.
- Ring L = 20 at eight points (ring20.jsonl, ring_cS.py, ring_cS.out), about 25 minutes on the loaded machine.
- Added the K > 15/8 kinematic window to PROOF 4.1(a) and §5. Wrote PROOF §7 and §8 and RESULT. Job 17604462 still PENDING with estimated start 15:29.
