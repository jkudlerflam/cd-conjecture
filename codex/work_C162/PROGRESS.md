This run was done by a Claude Opus agent (2026-11-08), not by Codex. Astra was not used (out of credits).

# C162 progress log (two-valley soliton onsets)

- Read specs/C162_soliton_onset.md, specs/_qubit_c_context.md, codex/work_C161/REFEREE.md (all), its scripts soliton_valley.py, soliton_pair.py, csl_phase.py, bloch_hessian.py and outputs, codex/work_C160/PROOF.md sections 1 to 4, codex/work_C159/RESULT.md.
- First pass, analytic. Plan:
  (a) dilute two-valley gas of hard-core solitons maps to a pseudospin XXZ chain on the squeezed lattice, J_z from the scattering lengths (n^3 scale), J_perp from exchange (pass-through plus flip-flop, exponentially small in 1/hbar_eff). c = 2 needs J_z > -|J_perp| (Delta > -1).
  (b) K stratum: the |n_z| rearrangement (W even in n_z, DM term (n x n')_z independent of n_z) gives V_{++}(R) < V_{+-}(R) at every separation. With the monotonicity identity da/dlambda = (M/hbar^2) int psi^2 dV/dlambda this gives a_{++} < a_{+-}^even, i.e. Delta < -1 for bosons, and Delta -> -infinity semiclassically for any statistics. Classical shadow of Eisenberg-Lieb.
  (c) general tail calculus: slowest channel gives additive (hard-rod) scattering lengths, no pseudospin coupling. Pseudospin coupling comes from the subleading channel, J_z ~ -c_2 (X_+ - X_-)(Y_+ - Y_-). Fix its sign in each symmetric class.
  (d) classify symmetries that protect a two-valley soliton band (P -> -P within one winding sector, fixed background): K (Theta R_pi(a), a perp n_bg) and U = I R_pi(n_bg). Check U by the tail calculus.
  (e) numerics: find a K bond with slow in-plane tail (m_w > m_u) and tilt-unstable soliton, so that the onset could be continuous; check V_{++} < V_{+-} and no binding there.
- kfun.py: general K functional W = s(a n_x^2/2 + b n_y^2/2 - hx n_x - hy n_y), S_zz = 0 gauge, field allowed off the principal axes. Validated against C161: lam0 = -0.0670 at kappa = 0.05, mu = 2.3, h = 0.0238 (C161: -0.067).
- Analytic: for a planar soliton the DM energy is -2 pi (topological), so E_gp scales as sqrt(s) and the onset scale is s_c = (2 pi/E_gp(1))^2. Exact identity H_perp = H_par + (m_w^2 - m_u^2) + 2 phi' + 2(L - L_bg), so lam0 >= m_w^2 - m_u^2 + min[2 phi' + 2(L - L_bg)].
- search_caseB.py, probe_caseB.py, opt_caseB.py: no tilt-unstable planar soliton with m_w >= m_u found (targeted Nelder-Mead: best lam0 = +0.62 at the constraint boundary).
- bound_scan.py (57 s, bound_scan.out): 4719 admissible random shapes, 21 tilt-unstable, all with m_u^2 > m_w^2 (ratio >= 4.6); lam0 - (m_w^2 - m_u^2) >= +0.088 in every shape. So at a K two-valley onset the tilt tail is always the slow one.
- lemma_check.py (lemma_check.out): off-axis field psi = 0.3, kappa = 0.05, mu = 2.3. Tilted soliton n_z single-signed (min/max = 5e-4), tails B_R = 0.36, B_L = 1.03 (asymmetric, both > 0, decay rate 0.2785 = m_w). U_alt - U_same > 0 at T = 10 to 22, U_same < 0 (like binding, first order).
- scatlen.py (scatlen.out): K_0 formula for a single exponential to 1e-6, so single-channel scattering lengths are additive (hard rods). Two-channel sign rule sign J = -sign((x1 - x2)(y1 - y2)) confirmed in K-like, U-like and no-symmetry cases, with J ~ k2^{-(lam2/lam1 - 1)}. Monotonicity da/dlam > 0.
- No miscible region in K or U, so Task 4 (DMRG on a c = 2 bond) does not apply. No Typhon job submitted.
- Writing PROOF.md and RESULT.md.
- Coordinator priority shift (construction side) addressed in PROOF section 10 and RESULT. No candidate, no Typhon job. RESULT.md and PROOF.md written. Done.
