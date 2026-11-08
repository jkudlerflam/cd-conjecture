# C165 progress log

Model: Claude Opus 5.5 agent (2026-11-08). Astra not used (out of credits).

- Read spec C165, shared context, C154 RESULT/PROOF (sections 2, 4, 5, 6, 7), C159 PROOF 2.2 and 4, C160 PROOF 3, C161 REFEREE 4-5, C163 RESULT, contrib/d2_floating_ising README, C152 RESULT, C155 RESULT.
- Copied unmodified: dmrg.py (C144), vumps.py and fesx.py (C163 copies of C153/C158), fesfit_c158.py (C158), est_c155.py (C155 est.py).
- classical.py (energy, gradient checked to 1e-13 relative, multi-start L-BFGS, Fourier diagnostics, tangent Hessian). symcheck.py (residual O(3) stabilizer of (J, B); identity exclusion modulo 2 pi).
- Classical scans (scan_cl1..4.py, outputs cl2_*.out, cl3_*.out, cl4_*.out). Ring L = 100, pitch a = 2 pi 13/100.
  - g = 0 (pure U(1) parent): spiral for Dz < 1, Neel for Dz > 1, degenerate at Dz = 1 (rigid-cone linearity).
  - g = 0.3 (charge-2 harmonic XX - YY): coexistence window 1.00 < Dz < 1.047, |N|^2 linear in Dz (continuous at both ends); relaxed spiral gain c = 0.024 (prediction g^2 cot^2 q / 4 = 0.020 at q = 0.817, 0.033 at the relaxed pitch 0.691).
  - G = 0.3 (Gamma, XZ + ZX): tilted Neel wins, no coexistence. hx = 0.3: canted Neel wins, no coexistence.
  - g = 0.3, hx = 0.2, hz = 0.1 (in N by symcheck): coexistence window about 0.99 < Dz < 1.05 with uniform canting |M| ~ 0.13.
  - Larger g (0.45 at a = 0.817; 0.5-0.9 at a = 1.19): window does not widen, the spiral turns into soliton lattices or collapses to collinear FM.
- hess_cl.py: at g = 0.3 the lowest staggered (Bloch pi) Hessian mode of the minimizing spiral falls linearly, w = 0.0995, 0.0605, 0.0214, 0.0018 at Dz = 0.95, 0.97, 0.99, 1.0, and crosses zero at Dz ~ 1.001 while the spiral is still the global minimum (e = -1.024 < -Dz). g = 0: w = 2(1 - Dz), zero only at the degeneracy.
- Quantum (two-site VUMPS via runq.py -> fesx.py, anaq.py; outputs loc/Q1_*.json (a = 0.8) and loc/Q2_*.json (a = 0.9), g = 0.3, hx = 0.2, hz = 0.05, Dx = 0.03, chi <= 40):
  - Dz = 1.0 to 1.15: staggered moment N decays with chi (0.082 -> 0.050 at Dz = 1.0), c_pair 0.9 to 1.07: c = 1 liquid with no staggered order, inside the classical coexistence window.
  - Dz = 1.25 to 1.4: N saturates (0.256, 0.311, 0.365, 0.413 at chi = 40), c_pair falls with chi (0.45, 0.27, 0.17, 0.09): gapped Neel_z. No point with saturated N and c_pair near 1.
  - a = 0.8 is close to the period-4 lock (transfer phase 0.494-0.499 per cell), so a = 0.9 is used for the quantum line.
- dmrgq.py (open chains, d-generic copy dmrg_d.py, structure factors). Local tests: TFI L = 32 c_prof 0.53; Ising x XX ladder L = 24 chi 32 c_prof 1.29, c_sm 1.43; Q2 Dz = 1.0 L = 32: c_prof 0.84, S_zz peak at pi (4.16), S_xx and S_yy peaks at k = 0.86-0.88.
- Typhon package in typhon/ (make_tasks.py, stageB.py, summarize.py, c165.slurm; 22 long tasks + 252 screen tasks + stage B). Stage B criteria fixed in stageB.py before submission. First ssh attempt: banner-exchange timeout (rate limit).
- Second-order Landau coefficients derived (PROOF 2.2) and checked: c_2 = g^2 cot^2 q/4 matches the spiral gain to 4 digits at g = 0.05, 0.1; Gamma Neel gain matches with the stiffness factor (1 + cos q)/(Dz + cos q).
- Local plateau-edge scan (loc/P_hz*.json, Dz = 2.0): gapped Neel N = 0.73 up to hz = 2.2, unconverged disordered state at hz = 2.6. A2 grid adjusted to Dz = 1.6, 2.0 accordingly.
- Typhon: files copied to /scratch/lustre/jkudlerflam/c165/, job 17601089 submitted (c165, 48 cores, 14 h). Last check: PENDING (Priority), estimated start 2026-11-08 13:07. Not waited on.
- Wrote PROOF.md and RESULT.md.
