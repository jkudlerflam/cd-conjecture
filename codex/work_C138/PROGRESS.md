# C138 progress log

Run by a Claude Opus agent (not Astra). J=1 throughout, eta = pi - gamma, epsilon = 1+Delta = 2 sin^2(eta/2).

## Done

1. Read the spec, C136 (PROOF, RESULT, REFEREE), C129 Sections 8 to 10.
2. New solver `bethe_ferro.py` in the step-free form L p(x_j) = sum_{k!=j} phi(x_j-x_k), phi(y)=2 atan(tan eta/tanh y). Agrees with ED to 5.8e-15 on all sectors L<=12 at eta in {0.3, 0.05, 0.01}.
3. `scan_defect.py K`: doubling defect Q = L[F_2L(2n)-2F_L(n)]/M_0 for epsilon=10^-K, K=2..6, L=4..4096 (paired rings to 8192), densities 1/L to 1/2. Output scan_K.json, energies cached in energies_K.json. [Numerical] 0 < Q <= pi v/2 (the Casimir limit) in every case; no band or dilute violation.
4. Ladder estimate of L(F_L - L e_bulk) = -(1/2) sum_j Lambda_j/4^j along doublings. [Numerical] |L(F_L-L e_bulk)| <= pi v/6, about eta/12.
5. `continuum.py`: restricted dressing equations on a sinh-graded mesh with product integration of the eta-narrow kernel; Richardson in the mesh. Checked against the half-filled closed form and against the ladder (agreement to 4e-7 in L(F-Le)).
6. `diag_half.py`: at half filling the cumulative root error L sup|G| is 0.67 to 0.76 for all eta and L (C136 bound B/2 diverges), and eta L sup|sigma_L - r0| is about 1/pi (C136 bound diverges like 1/eta).

7. `diag_identity.py`, `conv_check.py`, `dens_check.py`: at a finite Fermi boundary the cumulative error and the mismatch dN grow slowly as eta->0 (dN about 4.6 at epsilon=1e-6), but eta L sup|sigma_L - r_Q| stays below 0.28. The pointwise density error, not the cumulative error, is the uniform quantity.
8. `cont_diag.py`, `cont_diag_wide.py`: Var(v_Q) <= 8 eta, e_bulk'' <= (4/pi) eta^2, r_Q/r_0 >= 0.76 at distance >= 1 from the Fermi point; edge ratio scales as eta^(1/2).
9. PROOF.md written: Task 1 (normalization, step-free equations, eta->0 log-gas and sum rule), Theorem A (unconditional |F-Le| <= (5pi/4)eta), Lemma 1 (three-region quadrature), Theorem B (conditional (X) with K_* <= 125, C <= 628).

## Running / next

10. `dens_check2.py`: (D) max 0.317 over 120 cases, epsilon 0.9 to 1e-6; (D1) at most 0.22 (a large-epsilon artifact at x=0 identified and removed, corrected max 0.019).
11. RESULT.md written. Run complete. Open: a proof of (D) and (D1), and of (W1)-(W3).
