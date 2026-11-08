This run was done by a Claude Opus agent (2026-11-08), not by Codex. Astra was not used (out of credits).

# C161 progress log (referee of C160)

- Read specs/C161_referee_C160.md, specs/_qubit_c_context.md, specs/C160_fm_neighbourhood.md, codex/work_C160/PROOF.md (all), RESULT.md, PROGRESS.md, csl_tilt.py and its outputs, codex/work_C159/RESULT.md (Proposition C, open regions) and the h_F paragraphs of work_C159/PROOF.md (3.4, 4), codex/work_C154/RESULT.md (Lemma 3, P1 to P3).
- First pass, analytic. Points to test:
  (a) scope. The C159 window is K > 1 (two bosons) or K > 15/8 (Ising x boson), which reaches hbar_eff of order one. C160 has eps0 unquantified.
  (b) the singular set through h_F is a curve, J = -R(phi) (DM-twisted Heisenberg, exactly h_F after the spiral rotation), and it leaves every eps-ball.
  (c) M2 reduction: multipliers e^{ikT} != 1 lie in the energy-preserving subspace, so the reduction is valid there. Looks sound.
  (d) SC fails at dilute onsets (K -> 1). For solitons the analogue of Lemma 3' is missing. A tilt-unstable isolated soliton has a two-valley band E(P) at P = pi +- P0 (Haldane: P = s x solid angle), which is the J1-J2 two-minimum mechanism for solitons. Test classically whether the ingredients coexist in the K stratum.
- twisted_hF.py (seconds): J = -R(phi) has the h_F open-chain spectrum to 5e-14 for phi = 0.3, 1, 2, pi. C159 thresholds: K = 15/8 at Delta_AF = -0.669 (pure-DM eps = 1.11, hbar_eff = 2.2), K = 1 at Delta_AF = 0. Output twisted_hF.out.
- soliton_valley.py (seconds per point): K-stratum functional with in-plane anisotropy mu. At kappa = 0.05, mu = 2.3 the planar soliton is tilt-unstable at its onset over a stable Z2-symmetric background (lam0 = -0.068, h + kappa > 0). Tilted solitons: E below planar, max n_z = 0.55, solid angle 2 pi (1 -+ 0.443), so P = pi -+ 0.443 pi. Onset of the tilted soliton h_c between 0.02 and 0.03.
- soliton_pair.py: like tilts attract (U_same < 0), opposite tilts repel, |U| ~ exp(-T sqrt(h + kappa)). Onset first order (E_sol + U_same(10) < 0 while E_sol > 0).
- csl_phase.py: tilted (conical-like) CSL, <n_z> = 0.44, is the classical ground state for h in [0, 0.0265] at kappa = 0.05, mu = 2.3. Floating + uniform Z2 SSB exists in K near h_F (C160's statement is for mu = 0 only).
- bloch_hessian.py: Bloch Hessian of the tilted CSL (D.n up to 0.60, two-component gyroscopic). lambda_0 rises monotonically from 0 at k = 0, lambda_1 > 0. M2 confirmed on a non-planar texture. Planar CSL at kappa = 0.15: tilt band minimum at k = 0 (0.015). Outputs saved as console captures (*.out).
- Writing REFEREE.md.
- REFEREE.md written. Done.
