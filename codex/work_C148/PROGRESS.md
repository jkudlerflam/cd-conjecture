# C148 progress log

Run by a Claude Opus agent (2026-11-07), not by Codex. J=1 throughout.

## 2026-11-07

- Read the spec, _Cd_target.md, C129 (RESULT, PROOF), C139 (RESULT, PROOF), C140 (RESULT, PROOF, REFEREE), C142 (RESULT, PROOF), C143 RESULT.
- Plan. For Delta = cosh(lambda) > 1 the Bethe equations live on the circle x in R/(pi Z) with a1, K = a2 positive Poisson kernels. The structure is that of the repulsive side of C139/C140/C142 in the limit gamma -> 0 (q = ||K||_1 = 1, q_J = ||J||_1 = 1/2), on a compact rapidity domain. The plan is to transfer each lemma to the circle, with constants uniform in lambda in (0, infinity), handling the two places where the circle differs from the line: the constant Fourier mode (q = 1), and the wrap-around of the root set when Q > pi/4.
- Copied bethe_scan.py and check_c129.py from work_C142 (they already contain the periodic Delta > 1 solver of C129).
- Built circle.py (kernels on T, Nystrom restricted continuum, r0-quantile seed for the Bethe solver). quick_test.py: Bethe = ED to 4e-14 at L <= 10, closed-form e_bulk(1/2) reproduced.
- scan_defect.py (L to 64) and scan_large.py (L to 512, 2L to 1024): max defect 2.79 (L=4), all intermediate inequalities hold with large margin, sigma_L > 0 on T everywhere.
- Positivity. Plain C140 transfer fails beyond lambda = 0.18 because of constant Fourier modes and the two-sided out arc. Fixed by (i) sum_{m>=1} J^{*m} = K exactly, (ii) location bound via dn/am, L r0(x_n) >= 4 sin(pi/8) kK/pi^2, (iii) direct bound min a1 > max K/2 for lambda >= 1.36. pos_window.py, pos_check.py (exact, sup ratio 0.74), pos_rigorous.py (bracketing, 0.939 on [1,1.36], 0.804 analytic on (0,1]).
- Identification: C142 branch continued through Delta = 1 in XXX rapidities (analytic in s = gamma^2), then along lambda by open-closed; identity theorem.
- ed_checks.py: 324 sectors Bethe = ED (9.2e-14); 5805 end-to-end field/DM samples, max 3.9993; h_c1 = 2 a r0(pi/2) confirmed; odd-L half-filled defect positive, even-L exponentially small.
- constants.py: analytic majorant K <= 218.4 ab <= 4311 M_0, C = 10781; exact grid 4153 M_0, C = 10386.
- Wrote PROOF.md and RESULT.md.
