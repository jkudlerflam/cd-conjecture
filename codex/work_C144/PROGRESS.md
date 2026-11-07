# C144 progress log

Run by a Claude Opus agent, 2026-11-07 (not by Codex). Spec: specs/C144_xz_cone.md.

## Reading done

- specs/_Cd_target.md, specs/_bound_background.md, C143 RESULT and PROOF (2.2, 4.2, 5, 6), scripts ed.py, xz_at_sym.py and outputs, C141 RESULT and PROOF 5 and 6, C131 RESULT (cover lemma).

## Plan

1. Majorana dictionary for all 12 bond terms at XZ, first-order projection onto the degenerate manifold, tuned locus.
2. Exact free-fermion solution of the mixing directions (XX + ZZ, B_y) with ZX.
3. Second-order effective B Hamiltonian on the Ashkin-Teller cone: kernel K(r) from the critical A chain.
4. ED of the B-only effective model at long lengths.
5. ED of Q_L along all directions (L up to 10).
6. NumPy DMRG of the full qubit chain on the cone (open chains), validated on free fermions.

## Log
- kernel.py (kernel.out): second-order kernel K(r) of the critical A chain. Free-fermion formula equals ED of the full XZ ring (L = 8, 10) to 2e-16; K real symmetric (so nearest-neighbour b_j b_{j+1} terms cancel); tail K(r) (-1)^r 2 pi r -> 1; K~(pi) = (1/pi) log L - 0.0577.
- bmodel.py, validate_bmodel.py (validate_bmodel*.out): effective B chain h_B(lam) on L/2 spins reproduces the full qubit spectrum on the cone at beta = 0.004 (L = 8, 10, 12) to the expected O(beta) accuracy. Physical sector (c, parity) = (+1, 0) for L = 0 mod 4, (-1, 1) for L = 2 mod 4.
- bscan.py, bdiag.py (bscan_*.out, bdiag_lam0.out): B chain to L = 40. lam = 0: twofold ground state with exponentially small splitting, <b_j> = +-(-1)^j 0.638 (= 2/pi), gapless tower with L*gap growing like log L. So B self-generates a staggered-sign hopping and is a critical Majorana chain in a period-4 gauge, not a gapped dimer. lam = 0.5: gapped twofold. Larger lam: Ising-like up to a size that grows with lam (marginal flow).
- dmrg.py validated: open L = 10 vs ED (6e-15), free fermions L = 40, 80 (XZ + 0.5 ZX; XZ + 0.3 ZX + 0.2(XX+ZZ) + 0.1 By) energies to 1e-8 or better and S(l) to 1e-10 (test_dmrg.out; ff_open.py needed half fields at open ends).
- Typhon job 17580166 submitted (Q_L ED along 22 direction families, L <= 10, 2L <= 20). Output: /scratch/lustre/jkudlerflam/c144/q_<label>.jsonl, slurm_17580166.out.
- Typhon job 17580176 submitted: DMRG grid L=128 (beta 0.15, 0.3, 0.5; alpha/beta in 15 values), outputs /scratch/lustre/jkudlerflam/c144/d128_b*_r*.jsonl
- Local DMRG L = 64, beta = 0.3 (cone_b03_L64.out, fine_b03_L64.out): c ~ 1 at alpha = 0, 0.1 (edge), 0.3, 0.45; gapped (c_fit < 0.1) at alpha = 0.15, 0.191, 0.23 with strong A dimerization (a_stag ~ -0.25) and b_stag ~ 0.6; b_mid changes sign between alpha = 0.191 and 0.197 (B critical line near alpha* ~ 0.196); at alpha = 0.197 S_mid = 0.81, parity-split fits 1.1 / 0.9, not conclusive at L = 64.
- Directions at L = 64 (dirs_L64.out): XZ + 0.1 B_x gapped (S_mid 0.26); XZ + 0.1 YX critical-looking (c_fit 1.37 / 0.84, b_mid = 0), unexpected, to be checked at L = 128, 256.
- Typhon job 17580192 submitted: tuned-line finite-size scaling (beta = 0.3, alpha in 0.194..0.1985, L = 64, 128, 192, 256), YX (L = 128, 256), ZY, cone + YX, cone + mixing, cone + B_x at L = 128. Outputs /scratch/lustre/jkudlerflam/c144/t_L*_a*.jsonl, yx_*.jsonl, zy_*.jsonl, cone*_L128.jsonl, slurm_17580192.out.
- Typhon Q scan (typhon_out/q_*.jsonl, qscan_typhon_summary.out; 63 of the planned rows were in when copied): global max of Q_L over all directions and amplitudes is 2.3476 at L = 2 (direction r4, s = 0.05), against Q_2(XZ) = 2.3431; the max over L >= 4 is 1.847 (YY, s = 0.4). Gapped directions decay to 0.
- YX direction: ring ED (yx_ed.out) shows L x gap ~ 0.6 (L = 8, 10, 14, 16) and a triplet at ~1.6 to 1.75 at L = 12, 18 (period-3 commensuration), Q_L drops to 0.31 at L = 10, so A is gapped and a slow c ~ 1 sector remains.
- B-critical line, beta = 0.3: L = 96 DMRG at alpha = 0.1965 gives c fits ~1 (S_mid 0.845) but at alpha = 0.198 gives 0.47 / 0.44 with larger S_mid 0.895; the local uniform part of b drifts from +0.22 at the edges to -0.045 in the middle, so the B mass is inhomogeneous and profile fits are biased. Correlations at L = 64, alpha = 0.197 (corr.out): a-a connected decays exponentially (xi_A ~ 1.7), Y-Y fast, b-b as r^{-1.7 to -1.9} (Ising eps). A is gapped, B is Ising-critical: c_tot = 1/2 expected on the line. Decisive test is the scaling of max_alpha S(L/2) with L.
- Typhon job 3 (c144_d3, alpha 0.199..0.205 at L = 128, 192, 256) was sent while the login node was slow; check whether it was submitted.
- Local machine overloaded (load 46) from several BLAS-threaded runs; killed the local L = 128 run.
- Typhon jobs 17580166 (Q scan) and 17580176 (L = 128 grid) finished; outputs copied to typhon_out/. d128_summary.out: for beta = 0.15, 0.3, 0.5 the AT plane at L = 128 shows c ~ 1 for alpha/beta <= 0.3 and >= 1.0, a gapped wedge (c_fit < 0.1, a_stag up to -0.35, b_stag ~ 0.6) of width 0.05, 0.11, 0.17 in alpha, i.e. 1.17, 1.02, 0.98 times beta^2 log(1/beta), and inside it a single line (b_mid = 0) near alpha/beta = 0.64 to 0.68. On that line at beta = 0.5 (alpha = 0.34) the parity fits give 0.46 / 0.44. No fit anywhere exceeds 1.09.
- Jobs still running: 17580192 (t_L{128,192,256}_a*, yx, zy, cone + X at L = 128) and 17580432 (u_L{128,192,256}_a{0.199..0.205}); outputs in /scratch/lustre/jkudlerflam/c144/.
- peak64.out: S(L/2) ~ 0.87 flat for alpha = 0.199..0.203 then 0.17 at 0.205 (L = 64): cat log 2 from B edge zero modes on the ordered side; entropy estimators near the line unreliable, correlations used instead. RESULT.md and PROOF.md written; Typhon jobs 17580192 and 17580432 still running at close.
