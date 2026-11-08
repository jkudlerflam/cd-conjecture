This run was done by Claude Opus agents (2026-11-07), not by Codex (interrupted and resumed by a second agent).

# C153 progress log: staggered-charge box

## 2026-11-07, start

- Read specs/C153_staggered_box.md, specs/_qubit_c_context.md, contrib/d2_gap_reduction/staggered_slice.md and additive_charge.md, C15 RESULT sections 2.1 to 2.5, C143 PROOF sections 2.5 to 2.10, C144 RESULT, and contrib/d2_staggered/README.md (prior charge-resolved DMRG at 22 points, b = 0.5 to 4; it found c = 1 everywhere it was critical and left open the near-boundary continuum theory and b > 4).
- Copied udmrg.py (U(1) DMRG, validated in contrib/d2_staggered) and codex/work_C144/dmrg.py (as dmrg_c144.py) into this folder.
- Conventions. Original frame H = sum [XX - YY + u ZZ + b Z], a = 1, b > 0, r = sqrt(1 + b^2/4). Rotated frame (pi about x on odd sites) H' = sum [XX + YY - u ZZ] + b sum (-1)^j Z_j. Q = sum (-1)^j Z_j becomes sum Z_j, i.e. 2 S^z_tot.

## Plan

1. Analytic: bosonization at small b, exact kernel expansion near u = r (first order in r - u), spin-wave stiffness, strong coupling at large b.
2. ED (rings, Sz sectors, twist) for K, v, c and neutral levels; LSWT for the optical branch.
3. DMRG grid on Typhon (udmrg.py), including b = 6, 8 not covered before.

## 2026-11-07, first results

- ed.py: ring ED in Sz sectors with twist. Calibration on XXZ (b = 0): K and v reproduced to 1e-3 and 1 percent at L = 12 to 20 (ed_calib.out).
- Near-boundary expansion (derived, PROOF section 4): e''(m) = 8 r (r - u) to first order in r - u; classical (diamagnetic) stiffness D_cl = 2/(1+r); K0 = (pi/2)/sqrt(r(1+r)(r-u)), v0 = 4 sqrt(r(r-u)/(1+r)). ED (dsk.out) confirms e''/(r-u) -> 8r at every b, but the quantum stiffness at u = r is below D_cl for b > 0: D_s = 1.00, 0.98, 0.84 to 0.88, 0.59, 0.28, about 0.08 to 0.1 at b = 0, 0.5, 1, 2, 4, 8 (D_cl = 1, 0.985, 0.944, 0.828, 0.618, 0.390). So K0 with D_cl is an upper bound on the leading behaviour; kcheck.out shows K from ED approaches pi sqrt(D_s/(8 r (r-u))) as r - u -> 0.
- lswt.py: classical phases (rotated Neel for u < b/2 - 1, canted for b/2 - 1 < u < r, FM for u > r, the classical boundaries; the canted lower edge is shifted upward by quantum fluctuations to the BKT line). LSWT: one Goldstone branch with v = sqrt(e''_cl D_cl), optical branch gapped (minimum 4.0 in the whole canted region at every b tested).
- Typhon job 17594700 (c153, one 48-core node, 95 tasks, 2-day limit) submitted. Contents: FES (VUMPS) at b = 3, 6, 8, 12 near u = r plus xx and isxx controls; OBC DMRG strips at b = 3, 6, 8 with L = 64, 128; L = 64 grid over the box at b = 0.25, 0.75, 1.5, 3, 6. Outputs /scratch/lustre/jkudlerflam/c153/{fes_data,strip_data,grid_data}/*.json, logs in logs/, slurm_17594700.out.

## 2026-11-07, resumed (second Claude Opus agent)

- Job 17594700 was cancelled before running. Its scratch copy also lacked vumps.py, so every FES task would have failed on import; fixed in the resubmission.
- Trimmed resubmission: Typhon job 17596442 (c153b, one 48-core node, 20 h limit, 45 tasks, typhon/tasks2.txt, typhon/c153b.slurm). Contents: FES (chi 8 to 128) for xx and isxx controls, b = 3 (r - u = 0.003, 0.008, 0.015, 0.03), b = 6 (0.0005 to 0.005), b = 8 (0.0002 to 0.0025), b = 12 (0.0001, 0.0003, 0.0006); open-chain strips L = 64, 128 at b = 3 (r - u = 0.01), b = 6 (0.001, 0.002, 0.003), b = 8 (0.0005, 0.001, 0.0015); L = 64 grid at b = 0.25 and b = 6 only (b = 0.5 to 4 is covered by contrib/d2_staggered). Outputs /scratch/lustre/jkudlerflam/c153/{fes_data,strip_data,grid_data}/*.json, logs/, slurm_17596442.out. Analysis: fesfit.py of C158 on fes_data/*.json.
- Strips at b = 6, 8 switched to two neutral excited states (--nexc 2) on scratch before the job started.
- Local work after resumption: kscan.json analysed for BKT bounds (K(L) = 2 crossings); ring spectra at b = 6, 8 (spec_ed_b6_b8.out, 170 s); fesfit.py written for fes.py output and applied to fes_local/ and test_fes/.
- Added PROOF 5.1 (large-b scaling variable g = b^3 (r - u), with corrections to the t_2 factor and to the reading of the 4.0 optical branch), PROOF 7 (numerics) and PROOF 8 (theorem summary). Wrote RESULT.md.
- End of run: job 17596442 PENDING, estimated start 2026-11-12 13:59. Nothing collected from it.
