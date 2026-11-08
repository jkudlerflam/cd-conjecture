# C157 progress log

This run was done by a Claude Opus agent (2026-11-07), not by Codex.

## Setup (2026-11-07, evening)

- Read specs/C157_xz_counterexample.md, specs/_qubit_c_context.md, C144 RESULT.md and PROOF.md, C143 PROOF.md section 4.2, C154 PROGRESS.md (YX/ZY directions are a period-3 twisted XYZ chain, C154's finding).
- Copied from codex/work_C144 (unmodified): xzlib.py, kernel.py, dmrg.py, dmrg_run.py, cfit.py, bmodel.py.
- Typhon at start: C144 jobs 17580192 (c144_d2) and 17580432 (c144_d3) still RUNNING (6.0 h and 5.2 h elapsed). In /scratch/lustre/jkudlerflam/c144/ only the L = 64 tuned-line files exist; all t_L128+, u_L*, yx_*, zy_* outputs are empty. Nothing to collect yet.
- C144 timing: L = 128 at chi = 128 takes about 2 h per point on Typhon (2 threads). L = 256 is out of reach for a scan in this run, so Typhon will be used for targeted L = 64 to 128 scans only.

## Literature (RZFA)

- Rahmani, Zhu, Franz, Affleck, PRB 92, 235123 (2015), arXiv:1505.03966, Eq. (4): H = sum_j [ i t g_j g_{j+1} + g g_j g_{j+1} g_{j+2} g_{j+3} ], t > 0, with g < 0 called repulsive. Spin form (their Eq. 7): t sum Z - t sum X X - g sum Z Z - g sum X_j X_{j+2}.
- Phases (DMRG): g/t < -2.86 gapped fourfold (CDW); -2.86 < g/t < -0.285 Ising + Luttinger, c = 3/2; Lifshitz point (z = 3) at -0.285 (mean field -0.29); -0.285 < g/t < about 250 Ising, c = 1/2; tricritical Ising at g/t about 250; beyond it gapped twofold (dimerized), first-order line at g/t = infinity. A later study (arXiv:2211.15598) disputes the C-IC endpoint but not the sign of the c = 3/2 region.
- The c = 3/2 phase is on the repulsive side (g < 0) only.

## Kernel and mean field (local, done)

- kernel_q.py -> kernel_q.out. [numerical] The YY kernel K(r) of the A sector is positive semidefinite (min eigenvalue of the L = 60 kernel matrix about 1e-17), its Fourier transform Kt(q) rises monotonically from Kt(0) = 0 (critical A) to its maximum at q = pi, and K(2) > 0 for critical A (0.0748) and for dimerized A with d = 0.05, 0.2, 0.5 (0.0559, 0.0247, 0.0047; tails then exponential). Matched to RZFA: g = 2 beta^2 K(2) > 0, the attractive side.
- hfk.py (uniform Hartree-Fock in k space). Calibration on RZFA (hf_rzfa.out): mean-field Lifshitz point between g/t = -0.29 and -0.30 (RZFA: -0.29 mean field, -0.285 DMRG), c_MF = 3/2 beyond. Good.
- hf_bmodel.py -> hf_bmodel.out. Our B chain (critical or dimerized A kernel), lam = t/beta^2 from -4 to 8: h3/h1 > 0 at every lam, the slope v = sum r h_r never changes sign (|v| >= 2 for critical A), zeros = 2, c_MF = 1/2 throughout. No Lifshitz in the uniform ansatz. The dimerized ansatz wins near lam = 0 (C144).

## Coordinator request (received mid-run): Kitaev chain in a [110] field = XZ + bisector field

- C156 RESULT: H = K sum S^x_j S^y_{j+1} - (h/sqrt 2) sum (S^x + S^y). In the XZ frame (Pauli): H = sum X_j Z_{j+1} + B sum (X_j + Z_j)/sqrt 2, with B = 2h/K. Sorensen et al. (arXiv:2209.06221, read from the PDF): "critical line" for h < h_c1 = 0.511 K (B < 1.022), "difficult to establish numerically", no central charge; soliton phase with staggered vector chirality for 0.511 < h < 0.726. C156 job 17594698 (h = 0.2, 0.3, 0.4, 0.45, 0.6; L = 48 to 128) is pending (Priority).

## Typhon job 17595396 (c157), submitted 2026-11-07 evening

- One packed 48-core node job, 24 parallel tasks of 2 threads, 27 tasks (typhon/tasks.txt, typhon/run.sh, typhon/c157.slurm; runner dmrg_run.py from C144, unmodified).
- Tasks: (i) bisector field (Kitaev [110]) h = 0.3, 0.45 (B = 2h, Bx = Bz = B/sqrt2) at L = 192 and 256, chi up to 160/192; (ii) Ashkin-Teller wedge boundary at beta = 0.5, alpha = 0.21 to 0.27 at L = 64, 128 and alpha = 0.225, 0.24, 0.25 at L = 192; (iii) mixing + YY: XZ + alpha ZX + 0.3 YY + 0.1 (XX + ZZ), alpha = -0.15, -0.05, 0.05, 0.19, at L = 64, 128.
- Outputs: /scratch/lustre/jkudlerflam/c157/{bis_h*_L*, at_L*_a*, mix_L*_a*}.jsonl and .out, slurm_17595396.out. Analysis: python3 cfit.py <files> (and cfit.avg_fit).
- Local sanity check: dmrg_run.py on the bisector chain at open L = 10 equals dense ED (E0 = -6.8535968), L = 16 runs.

## Bisector: exact rewriting and local numerics

- block_check.py -> block_check.out [proved, operator identity to 1e-14 on rings L = 8, 10]: Hadamard on even sites, blocks (2i-1, 2i), and the variables sigma^z_i = tau^z_i nu^z_i tau^z_{i+1}, sigma^x_i = X_{2i}, g^z_i = G_{i+1} (Gauss operator of the XZ point), g^x_i = Z_{2i+1} turn the bisector chain into
  H = sum_i [ sigma^z_i + g^z_i sigma^x_i sigma^x_{i+1} + (B/sqrt2) ( sigma^x_i (1 + g^z_i) + g^x_i (1 + sigma^z_i) ) ].
  (A = critical TFIM sigma with bond signs g^z; the field couples the A order parameter to (1 + g^z) and flips g only where sigma^z = +1.)
- bis_ed.py -> bis_ed.out: rings L = 8 to 18 at B = 0.4 (h = 0.2): (E - E0) L of the 7 lowest levels is 0.003 to 0.3 and decreasing with L. The slow sector is not resolved by ED (killed before L = 20).
- Local DMRG (local/bis_h*_L*.jsonl): incommensurate oscillations of S(l), <a_j>, <b_j> with period about 8 to 10 sites at L = 32. Profile fits are meaningless at these sizes (0.3 to 2.4).

## Calibration: does the pipeline see c = 3/2 when it is there? (local, rzfa_dmrg.py, dmrg_d.py)

- dmrg_d.py = C144 dmrg.py with d inferred from h and an optional uniform on-site MPO term; validated against dense ED with d = 4, n = 5 (random bond and on-site terms): equal to 1e-13.
- RZFA spin chain blocked in pairs. g = 0 (TFIM), L = 64: fits 0.51 to 0.52. g = +1 (Ising phase), L = 64: 0.52. g = -1 (c = 3/2 phase): L = 64 profile fits 0.75 / 2.1 / 6.0 (windows 0.1 / 0.2 / 0.3, incommensurate), smoothed fit 0.91; L = 128 fits 1.50 / 1.58, smoothed 1.69, c_pair(64, 128) = 1.24. So a true c = 3/2 incommensurate phase reads near 1.5 only from L = 128 on, and the pair estimate still lags. Any candidate must be judged at L >= 128 with the same estimators.

## B chain at small lam (bsmall.py -> bsmall.out)

- lam = 0.1 and 0.25 (where uniform HF comes closest to a Lifshitz point): ring ED, L = 24, 32, 40. Twofold ground state (splitting 0.009/L and below), then a gap: L (E - E0) = 29.6, 36.5, 43.0 (lam = 0.1) and 53, 70, 87 (lam = 0.25), growing linearly in L, so a finite gap of about 1.1 and 2.2 (units beta^2). Gapped twofold (dimerized, cat eigenstates), no tower of a c = 3/2 liquid.

## Bisector: C156 outputs so far and a Casimir test of the fast A Majorana (c156_copy/, cscale.py, casimir.py)

- Copied C156's finished jsonl (read only) into c156_copy/. Central-entropy pair estimates: h = 0.3, c_pair(48, 64) = 1.05; h = 0.45, c_pair(48, 96) = 1.04 (profile fits at L = 96: 0.94 to 1.17); h = 0.4, c_pair(48, 64) = 0.44; h = 0.6 (soliton phase) c_pair 0.43, 0.30, 0.21 falling (gapped).
- casimir.py: E0(L) = e L + e_b + f/L on open chains. Calibration, XZ point (A critical, v = 4): f = -0.245 from L = 32, 48, 64 (exact -pi/12 = -0.262 up to 1/L^2). Bisector: h = 0.3 (L = 32, 48, 64) f = +0.146; h = 0.45 (L = 32, 48, 96) f = -0.012. No -0.26/L term, so no fast critical A Majorana on the line [numerical, three sizes, indicative]: A is gapped by the field, and c_tot is the c of the slow incommensurate sector, which reads about 1.

## Bisector: C156 outputs to L = 128 (copied at about 22:45 local)

- cscale.py (central entropy S-bar over l/L in [3/8, 5/8]; pair estimate c = 6 dS-bar / d log L):
  h = 0.4: c_pair = 0.44 (48, 64), 0.94 (64, 96), 0.98 (96, 128); smoothed profile fit 0.97 at L = 128.
  h = 0.45: 1.26 (48, 64), 0.89 (64, 96), 0.97 (96, 128); smoothed fit 1.07 at L = 128.
  h = 0.3: 1.05 (48, 64), 0.84 (64, 96). h = 0.2: 1.16 (48, 64). h = 0.6 (soliton phase): 0.43, 0.30, 0.21 (gapped).
  Against the RZFA calibration (true c = 3/2: c_pair(64, 128) = 1.24, smoothed fit 1.69 at L = 128), the bisector line reads c = 1 at L = 128, not 3/2.
- Correction to the Casimir entry above: with more sizes the triple fits of f jump between -0.31 and +0.22 (h = 0.45: +0.09, -0.31, +0.22; h = 0.4: +0.22, -0.14). The incommensurate oscillations of E0(L) dominate the 1/L term, so the Casimir test does not decide whether A is critical. Withdrawn as evidence.
- corr_bis.py (local, L = 48, h = 0.45): a-a connected correlator oscillates with a slowly decaying envelope (6.8e-2, 1.0e-2, 1.4e-4, -7.6e-4, 6.3e-4, 4.6e-6, 7.0e-5 at r = 2, 4, 6, 8, 12, 16, 24), incommensurate. Not decisive either.

## Close-out (end of session)

- PROOF.md and RESULT.md written. Verdict: no counterexample. The induced B interaction is positive semidefinite (attractive), on the attractive side of RZFA (g = 0.150 beta^2 > 0). The Lifshitz and tricritical routes are closed at second order. The bisector line reads c = 1 at L = 128.
- Pending at close: C157 job 17595396 (PENDING, Priority); C144 jobs 17580192 and 17580432 (RUNNING, no L >= 128 outputs yet); C156 job 17594698 (RUNNING). Output paths in RESULT.md.
