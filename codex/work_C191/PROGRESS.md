This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

# C191 progress log

1. Read the spec (C191 section), C175 RESULT/PROOF/smatrix.out/scan_w.out, REPORT 9.139, 9.140, 9.158, C179 REFEREE/PROOF, C181 PROOF, C173 PROOF. Did not open codex/work_C189/RESULT.md.
2. Wrote scat191.py (independent two-body S-matrix). Checked: unitarity 1e-13, C175's a_±(q) reproduced, hard-rod phase −q exact (check1.py).
3. Fermi-surface data at n = 1/16 to 1/4 (fs_data.py). Found |r(k_F)| small at the candidate densities, transmission close to +1 (boson frame), like phase δ_s ≈ 0.5 to 0.6.
4. Gaussian dimension code gauss191.py. First version used W instead of S W S in the field Hamiltonian; the one-valley g₂ test caught it (gave 1/K instead of K); fixed and retested (free, one-valley, drift, Gaudin–Yang tests pass).
5. lowE.py: first-order Landau parameters, g₁, K's, exponents, three f₄ variants. R–R kinematics checked separately (rr_check.py): near-total reflection, threshold a₊(K) = 12.0, 9.7, 8.0.
6. aba_valley.py: one valley in the asymptotic Bethe ansatz. Hard rods exact to 1e-4. K_∥ = 1.36, 1.46, 1.49. Check at φ = 1.3 against C179's rings (aba_phi13.py): comparable, ABA slightly low.
7. hybrid.py, predict.py: hybrid estimate and the twelve-variant band; scan in n (scan.out) gives n_KT ≈ 0.07 and the zero of |r| near n ≈ 0.167.
8. ring2_frame.py: two-boson ring L = 600 confirms the boson frame and the sign of the phase (2.7e-6 against 9e-3).
9. gy_check.py: Gaudin–Yang separatrix and gap exponent reproduced.
10. ktflow.py: full KT flow, drift of K_s between L = 128 and 256 below 3e-5 at the candidate densities; ξ_s below n_KT.
11. Registered the prediction in RESULT.md. Wrote PROOF.md. Then read codex/work_C189/RESULT.md for the comparison (appended to RESULT.md).
12. Read C189 RESULT.md: registration only, no numbers. Computed the single-particle wall reflection (wall.py) and the stiffnesses (stiff.py); appended the frame and boundary notes to RESULT.md and PROOF 7.6, 7.7.
