This run was done by a Claude Opus 5.5 agent (2026-11-08), not by Codex. Astra was not used (out of credits).

# C165 result: a staggered Ising transition inside a floating liquid

Derivations are in PROOF.md, the log in PROGRESS.md. Scripts here: classical.py, symcheck.py, scan_cl1..4.py, hess_cl.py (classical); runq.py, anaq.py, dmrgq.py, dmrg_d.py (quantum drivers). Copied unmodified: dmrg.py, dmrg_c.py, est.py (C163 copies of C144/C155), vumps.py and fesx.py (C153/C158 via C163), fesfit_c158.py, fes_c153.py and fesfit_c153.py (C153), est_c155.py. dmrg_d.py is dmrg.py with d inferred from the bond (needed for the d = 4 Ising x XX control). Local outputs are in loc/ and cl*.out, the Typhon package in typhon/.

## Verdict

No c = 3/2 bond was found. The construction exists classically and fails quantum mechanically, for a reason that can be stated precisely. Classically, the obstruction named in the spec is only a property of the rigid ansatz. A charge-2 harmonic (XY anisotropy about the spiral axis) makes the Landau N^4 coefficient positive, and a symmetry-free bond then has a minimizing incommensurate spiral whose staggered Bloch-pi mode softens continuously into a coexisting phase. At spin 1/2 this coexistence is absent. The spiral/staggered switch of every nearest-neighbour bond's U(1) parent is, exactly, the spin-1/2 XXZ chain at Delta = 1 (or its spin-flop line in a field), so the staggered moment and the spiral are the three components of one SU(2)_1 primary, and the classical curvature becomes an operator of dimension at least 8. This is derived to all orders in weak harmonics and confirmed numerically at harmonic amplitude 0.3. What remains is the narrowest assumption A_S (harmonics comparable to the exchange, near the hidden Heisenberg point), which the pending Typhon screen tests.

## Task 1: mechanism

- [proved] For every bond and incommensurate q the rigid staggered cone has energy sin^2 theta e_s + cos^2 theta e_N (PROOF 2.1), so it never coexists.
- [proved] The parent (XXZ + DM along the spiral axis + field) maps by exact site rotations to the spin-1/2 antiferromagnetic XXZ chain with Delta = Dz/|J_perp| in a longitudinal field. The switch is its isotropic point (field zero) or its spin-flop line, and the cone family is the SU(2) orbit there or the classical biconical degeneracy (PROOF 1.2).
- [derived, checked to 4 digits] Second-order relaxation gives the N^4 coefficient c_2 - c_1, with c_2 = g^2 cot^2 q / 4 from XX - YY and XY + YX (convex, favours coexistence) and c_1 = G^2 (1 - cos q)^2/(2(1 + cos q)) from Gamma terms XZ + ZX (concave, first order). DM and fields give no quartic term. The terms that allow coexistence are the charge-2 harmonics (PROOF 2.2, 2.3).
- [derived] Spin 1/2: near the switch the infrared is SU(2)_1, N^a = tr(g sigma^a) contains both orders, the quadratic part of the Landau potential is the marginal XXZ anisotropy and the quartic part has dimension >= 8. To all orders in weak harmonics the outcome is the XXZ dichotomy, a c = 1 liquid or gapped Neel order (PROOF 4.1, 4.2). The liquid's boson is the staggered order here, in the sense the spec asks to check.
- [derived] The classical window is O(g^2) wide and sits at K = 1/2 to 1. The semiclassical regime (large K, near h_F, C160) is disjoint from it, and there the staggered mode is far from soft (PROOF 4.4). In N, C159 Lemma A gives both fields incommensurate, so Ising x boson is kinematically stable at any K. The obstruction is dynamical, not kinematic.

## Task 2: classical lattice search

- [numerical] Ring L = 100. g = 0.3: coexistence for 1.00 < Dz < 1.047, |N|^2 linear in Dz (continuous at both ends). In N (g = 0.3, hx = 0.2, hz = 0.1, symcheck residuals 0.067 and 0.040): coexistence for about 0.99 < Dz < 1.05. G = 0.3 or hx = 0.3 alone: first-order switch to tilted or canted Neel.
- [numerical] Hessian (hess_cl.py): at g = 0.3 the staggered mode of the minimizing spiral goes 0.0995, 0.0605, 0.0214, 0.0018 at Dz = 0.95 to 1.00 and crosses zero at Dz ~ 1.001 while the spiral is still the global minimum. This is a continuous classical transition at Bloch momentum pi. C160 M2 does not cover the lattice.
- [numerical] The window does not exceed a few percent of the exchange: larger g turns the spiral into soliton lattices or a collinear ferromagnet (PROOF 2.4).

## Task 3: quantum test

Criteria, fixed before any Typhon output (typhon/summarize.py, stageB.py): entanglement fits above 1.3 at L = 128 not falling to 256, truncation below 1e-7, staggered order vanishing continuously at the point, incommensurate peak on both sides, FES slope above 1.3 at chi >= 256. Floating + SSB flag: N(64) >= 0.05, N(64)/N(32) >= 0.97, c_pair(48 to 64) >= 0.8, xi(64)/xi(48) >= 1.2.

- [numerical, local, chi <= 40, two-site VUMPS] Best candidate: line Q2, a = 0.9, g = 0.3, hx = 0.2, hz = 0.05, Dx = 0.03 (in N), through its classical coexistence window. At Dz = 1.0, inside the window, the staggered moment decays with chi (0.082 to 0.050) and c_pair = 0.96 to 1.07: a c = 1 floating liquid. At Dz = 1.25 to 1.4 the staggered moment saturates (0.26 to 0.41) while c_pair falls (0.45 to 0.09): gapped Neel. The quantum Neel onset is near Dz = 1.2, outside the classical window. No point is both ordered and gapless. Largest c_pair seen, 1.07 at chi = 40 (Dz = 1.0), consistent with c = 1.
- [numerical, local] Neel plateau edge at Dz = 2.0 with the same harmonics: gapped Neel up to hz = 2.2, then a disordered state. Finer edge scans are in the Typhon screen.
- Typhon job 17601089 (c165, one 48-core node, 14 h, partition typhon), submitted 2026-11-08, PENDING (Priority) with estimated start 2026-11-08 13:07. Scratch /scratch/lustre/jkudlerflam/c165/. Contents:
  - long tasks (tasks_long.txt): DMRG L = 128 (chi 128) and 256 (chi 160) on Q2 at Dz = 1.15, 1.2, 1.25, 1.3, the floating control and XX; TFI L = 128; Ising x XX ladder L = 64, 128; FES to chi = 256 for Ising x XX, XX, TFI (to 160), the floating control and Q2 at Dz = 1.2, 1.25, 1.3. Outputs dm/*.json, fes_cal/*.json, fes/nfloat.json, fes/Q2big_*.json.
  - stage A (tasks_A.txt, 252 points, chi <= 64): A1 charge-2 family (a = 0.7, 0.9, 1.1; g = 0.2, 0.35; hz = 0.05, 0.3; Dz = 1.0 to 1.8), A2 Neel plateau edges (Dz = 1.6, 2.0, three harmonic sets, eight fields each), A3 sixteen random symmetry-free harmonic vectors of norm 0.35 at Dz = 1.0 to 1.5. Outputs fes/A*.json.
  - stage B: stageB.py selects up to 24 points (floating + SSB flags first, then c_last >= 1.25, then the ordering transition of each line), continues them to chi = 256 and runs DMRG L = 128, 256 on the top four.
  - Digest: summary.txt. Collect with scp -r 'typhon:/scratch/lustre/jkudlerflam/c165/{summary.txt,stageB_selection.json,dm,fes_cal}' typhon_out/.

## Task 4: calibration

- [numerical, local] TFI L = 32: c_prof 0.53. Ising x XX ladder (c = 3/2, outside B_2) L = 24, chi = 32: c_prof 1.29, c_sm 1.43. contrib/d2_floating_ising's FES with the same VUMPS code read 1.536 on Ising x XX and 1.000 on XX. The c = 1 floating point in B_2 used for the quantum line (Q2, Dz = 1.0) reads c_pair 0.96 to 1.07.
- The calibrations at the decisive sizes (Ising x XX FES to chi = 256 and DMRG L = 128; the N-stratum floating control at chi = 256 and L = 256) are in job 17601089.

## Task 5: summary

- [proved] Rigid-cone linearity for every bond, and the exact map of every parent switch to the spin-1/2 XXZ isotropic point or spin-flop line.
- [derived] The classical N^4 coefficient is c_2 - c_1 (charge-2 against Gamma). Classical P1 exists in N.
- [derived] At spin 1/2, near the switch and to all orders in weak harmonics, the staggered and spiral orders are one SU(2)_1 multiplet, so the coexistence transition does not exist and c <= 1. Scope: the perturbative neighbourhood of the U(1) class in the hidden-Heisenberg region. Near h_F it is C160 M2. Elsewhere in N, C154 (A6) stands.
- [assumed] Narrowest remaining assumption A_S: in N, with harmonics comparable to the exchange near the hidden Heisenberg point (parent Delta near 1, fields up to the Neel plateau edge), no lattice-scale excitation independent of the liquid's boson orders at momentum pi while the liquid survives. The Typhon screen tests amplitudes 0.2 to 0.6.
