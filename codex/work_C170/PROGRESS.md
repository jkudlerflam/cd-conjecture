C170 progress log. Model: Claude Opus 5.5 agent (2026-11-08). Astra not used (out of credits). No Typhon use.

- Read specs/C169_C174_swarm.md, specs/_qubit_c_context.md, STATUS.md, REPORT.md §9.110 to §9.126, codex/work_C154/RESULT.md, codex/work_C165/RESULT.md, PROOF.md, PROGRESS.md, classical.py (family definition), runq.py, anaq.py, typhon/make_tasks.py, loc/P_*.json.
- Re-derived §1.2 by hand. frame.py: open-chain spectra parent vs XXZ equal to 5.5e-14. Ring identity [H_μ, T_μ e^{−i(a−π)S^z}] = 4.4e-15 with random harmonics.
- Re-derived §2.1 and §2.2 by hand. Found that DM gives a quartic term (C165 used L^T A L = 0, which ignores the frame rotation). Derived the general charge-1 coefficient c_1 = sin²(q/2)|tan(q/2)Γ_c − iD_c|².
- lr2.py: exact linear response on a ring L = 100. It confirms c_2 and Γ exactly, DM −D² sin²(q/2), Γ–DM_x interference, and that fields are linear. check_c1.py: eight random harmonic vectors match the corrected formula to all digits, and C165's formula has the wrong sign in two.
- clmin.py (background, about 10 min): g = 0.2 window 1.000 to 1.019; with D_x = 0.12, 1.004 to 1.019; with D_x = 0.25, first-order jump at 1.020 to 1.0225. Outputs cl_g02.out, cl_g02_D012.out, cl_g02_D025.out.
- heff.py: a uniform field is generated at second order by hx·G and hx·D_x (m = 1.6e-3, quadratic scaling), and not by hx, G or g·hx alone.
- plateau.py: exact plateau edge hz_c1(Δ) (des Cloizeaux–Gaudin), with ring ED in S^z sectors consistent. Values 0.022 (Δ = 1.25), 0.27 (1.6), 0.78 (2.0). Néel onset at hz = 0.05 is Δ = 1.321.
- plateau2.py and pfull.py: C165's plateau-edge VUMPS states at Dz = 2, hz ≥ 1.4 lie about 0.1 per site above the ground-state energy, so they are metastable. Estimator error about 0.02, read from hz = 0.
- Wrote PROOF.md, REFEREE.md, RESULT.md.
