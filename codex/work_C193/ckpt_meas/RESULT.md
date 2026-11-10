This run was done by a Claude Opus 5.5 agent (C193b, 2026-11-10), not by Codex. Astra was not used.

# C193b result: the missing correlators, measured on the C190 checkpoints

Derivations in PROOF.md, log in PROGRESS.md. Brief: `specs/C193_C194_sign_tests.md`, section C193b.

## Registration (copied from the spec before any candidate output of the c193 job was opened)

- (B0) Calibration. On the two-leg control η_b is within 0.03 of 1/2, the pair exponent within 0.08 of 1, and K_c from the structure factor within 5% of its exact value, at L = 128 and 256. If not, everything below is uncalibrated.
- (B1) Tracks the effective model. At L = 256 on the ground branch η_b is within 0.05 of 0.326 (n = 3/16) or 0.351 (n = 1/8), and the pair exponent lies in 0.50 to 0.75, at every key.
- (B2) Gapless sign. With 1/K_c measured on the qubit chain, K_s ≥ 1.1 at L = 256 for every key in both sectors, the Kosterlitz–Thouless bound from the 128 → 256 change gives a lower limit on y_∥ above the upper limit on |y_1| at every key, and the change of K_s from 128 to 256 is no more negative than the change from 64 to 128 (a flow to the gapless line slows down).
- (B3) Gap sign. K_s(256) < 1 at any key, or the bound allows |y_1| ≥ y_∥ at any key, or the downward change of K_s from 128 to 256 exceeds the one from 64 to 128 by more than the combined error.
- Anything between B2 and B3 is reported as undecided, with the numbers.

Conventions fixed with the registration, before any output. Ground branch at each key and L is the parity sector with the lower stored final energy (as in stage one). Fit ranges are frozen on the two controls only (xx2 two-leg XX, c = 2, and xx1 single XX chain, c = 1), at fixed r/L, and then applied unchanged to the candidates and to the effective model's C176 data. The exact control values are η_b = 1/2 (both controls), pair exponent 1 (two-leg, inter-leg pair b_A b_B at d odd), K_c = 1 (two-leg, π times the small-q slope of S(q)) and 1/2 (single XX chain, same normalization, PROOF §4).

How the registration is read, fixed now (before any candidate output was opened).
- "Pair exponent" is EP at d = 3 (pair members three sites apart). On the two-leg control only odd d gives the inter-leg pair (exponent 1); at even d the pair sits on one leg and has exponent 2. For the candidates both d carry the inter-valley pair at momentum near π (PROOF §2), and d = 2 is reported alongside.
- K_s is computed with both estimates of 1/K_c, the pair exponent (K_s^pair) and 1/K_c from the structure factor (K_s^sf). Each clause of B2 and B3 is evaluated with both. B2 is met only if it holds with both, B3 is met if it holds with either, and the letter is stated per estimate when they differ.
- The KT bound is drift.py's: Δl = ln 2 per doubling, |y_1| upper limit sqrt(max(0, −2(ΔK − σ)/(K̄² ln 2))) from the 128 → 256 change, y_∥ lower limit 2(K_s(256) − σ_K(256) − 1). "Combined error" in B3 is the quadrature sum of the errors of the two changes.

## Calibration and frozen ranges (controls only, written before any candidate output was opened)

Commands, from this folder: `python3 -E -s ana193.py controls typhon_out/out --eff` (trial ranges, also on C176's c = 2 control) and the frozen estimators in `frozen.json`.

Checks of the loaded control states (job c193, outputs in typhon_out/out). Energy against the stored value 1e-14 to 1.3e-12, S̄₆ to below 2e-14, the full entropy profile to 6e-14, the stored 48-site Gb to below 7e-15. So the checkpoints are the final states of the runs.

Trial ranges (fractions of L). EB on (0.04, 0.12), (0.06, 0.20), (0.10, 0.30), (0.15, 0.40), (0.20, 0.60), (0.30, 0.80). EP on (0.04, 0.12), (0.06, 0.20), (0.10, 0.30), (0.15, 0.40), (0.12, 0.34), (0.15, 0.38). ES on R/L = 0.25, 0.35, 0.45 with q windows (3, 8) and (4, 10) in units of π/R.

| control | L | EB trial ranges 1 to 6 | EP d = 3 ranges 1 to 4 | ES (0.25,3-8) (0.25,4-10) (0.35,3-8) (0.35,4-10) (0.45,3-8) (0.45,4-10) |
| :-- | :-- | :-- | :-- | :-- |
| two-leg XX | 64 | 0.299 0.490 0.489 0.531 0.496 0.492 | −0.74 0.44 1.058 – | 0.31 0.04 1.21 0.37 2.00 1.24 |
| two-leg XX | 128 | 0.547 0.497 0.503 0.504 0.501 0.495 | 1.043 1.051 1.017 0.998 | 1.98 2.41 1.014 1.43 1.043 1.022 |
| two-leg XX | 256 | 0.492 0.500 0.502 0.503 0.499 0.501 | 1.030 1.010 1.002 1.000 | 1.024 1.016 1.035 1.024 1.047 1.034 |
| XX chain | 64 | 0.487 0.539 0.512 0.505 0.505 0.517 | (exact 2) 0.03 1.48 2.20 – | 1.00 1.15 0.504 0.754 0.511 0.505 |
| XX chain | 128 | 0.522 0.509 0.503 0.500 0.502 0.503 | 2.26 2.13 2.04 2.02 | 0.508 0.505 0.511 0.508 0.515 0.511 |
| XX chain | 256 | 0.506 0.502 0.500 0.501 0.501 0.501 | 2.07 2.03 2.01 2.00 | 0.492 0.494 0.489 0.491 0.485 0.489 |

Exact values: η_b = 1/2 on both; EP at d = 3 is 1 on the two-leg control (inter-leg pair) and 2 on the XX chain (e^{2iθ}); K_c = 1 (two-leg) and 1/2 (XX chain).

Frozen (frozen.json).
- EB: ranges (0.10, 0.30) and (0.15, 0.40), both valley peaks, η_b = mean, error = half spread ⊕ 0.005 systematic (largest control bias at L ≥ 128 is 0.004). The upper end 0.40 L is the longest range that the effective model's stored window (L/2 sites) also reaches, so the same ranges apply to both.
- EP: ranges (0.10, 0.30) and (0.15, 0.38), error = half spread ⊕ 0.015 (largest two-leg bias at L ≥ 128 is 0.014).
- ES: (R/L, q window) = (0.35, 3-8), (0.45, 3-8), (0.45, 4-10), K_c = mean, error = half spread ⊕ 4% (largest control bias at L ≥ 128 is 3.9%).

Frozen values on the controls.

| control | L | η_b | EP d = 3 | K_c (ES) |
| :-- | :-- | :-- | :-- | :-- |
| two-leg XX | 64 | 0.510(22) | 1.043(21) | 1.48(40), fails |
| two-leg XX | 128 | 0.504(5) | 1.014(15) | 1.026(44) |
| two-leg XX | 256 | 0.502(5) | 1.001(15) | 1.039(42) |
| XX chain | 64 | 0.509(6) | 2.145(59) | 0.507(21) |
| XX chain | 128 | 0.501(5) | 2.032(19) | 0.512(21) |
| XX chain | 256 | 0.500(5) | 2.008(15) | 0.488(20) |
| C176 two-leg (effective model) | 64 / 128 / 256 | 0.501 / 0.502 / 0.501 | 1.052 / 1.010 / 1.002 | 1.96 (fails) / 0.968 / 0.968 |

B0 (registered). Two-leg control at L = 128 and 256: η_b = 0.504 and 0.502 (within 0.03 of 1/2), pair exponent 1.014 and 1.001 (within 0.08 of 1), K_c = 1.026 and 1.039 (2.6% and 3.9%, within 5%). B0 holds.

Limits found on the controls, stated before the candidates. At L = 64 the structure-factor estimator fails on the two-leg control (1.48, and 1.96 on the effective model's two-leg control), because R ≤ 0.45 L ≈ 29 sites holds too few wavelengths of the leg density. So K_s^sf at L = 64, and therefore the 64 → 128 change of K_s^sf, is not calibrated. At L = 64 η_b and the pair exponent are calibrated only to about 0.02 to 0.04 (two-leg 0.510, 1.043), against 0.005 and 0.015 at L ≥ 128.

Frozen. Nothing above is changed after this point.

## Results (candidate outputs opened after the section above was written)

Job 17611534 (c193, apollo02, one A100, 4 cores, 64 GB), 07:34 to 07:46 EST, all 42 states done (the three smoke-test controls of job 17611514 were reused). Outputs in typhon_out/out (one `<tag>_meas.json` per state), job logs in typhon_out/logs. Full analysis printout: analysis.txt (`python3 -E -s ana193.py all typhon_out/out --eff`), numbers in analysis.json.

### State checks [numerical]

Every candidate and control state was measured from its final checkpoint (no state is missing, nothing was rerun). For all 42 states: |E − E_stored| ≤ 2.7e-12 (L = 256, against the last sweep's energy) and ≤ 1.1e-13 elsewhere, |ΔS̄₆| ≤ 1.8e-14, whole entropy profile ≤ 5.7e-14, stored 48-site Gb reproduced to ≤ 1.2e-14, and the checkpoint sweep count equals the run's history length.

### R1. Frozen estimators on the candidates, ground branch (*) [numerical]

| key | n (256) | L | η_b | pair (d = 3) | K_c (ES) | K_s^pair | K_s^sf |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| p153c | 0.215 | 128 | 0.274(20) | 0.631(18) | 1.17(15) | 2.14(37) | 4.1(2.3) |
| | | 256 | 0.282(16) | 0.638(16) | 1.703(73) | 2.05(27) | 1.85(23) |
| sym40_n316 | 0.192 | 128 | 0.278(16) | 0.632(18) | 1.22(14) | 2.09(28) | 3.5(1.4) |
| | | 256 | 0.285(12) | 0.638(16) | 1.716(72) | 2.00(20) | 1.80(18) |
| rat_n316 | 0.200 | 128 | 0.277(20) | 0.626(18) | 1.20(15) | 2.08(36) | 3.7(1.8) |
| | | 256 | 0.284(16) | 0.635(15) | 1.709(73) | 2.00(26) | 1.82(22) |
| sym40_n18 | 0.136 | 128 | 0.290(18) | 0.641(19) | 1.21(16) | 1.93(27) | 3.0(1.2) |
| | | 256 | 0.297(14) | 0.643(16) | 1.732(73) | 1.84(20) | 1.64(17) |
| p153_n316 | 0.196 | 128 | 0.277(20) | 0.626(18) | 1.20(15) | 2.08(36) | 3.7(1.8) |
| | | 256 | 0.284(16) | 0.634(15) | 1.710(72) | 1.99(26) | 1.81(22) |
| sym80_n316 | 0.195 | 128 | 0.278(14) | 0.629(18) | 1.22(15) | 2.08(26) | 3.4(1.4) |
| | | 256 | 0.286(12) | 0.640(16) | 1.734(74) | 1.99(19) | 1.77(16) |

Both sectors and L = 64 are in analysis.txt. The other sector agrees with the ground branch within the errors everywhere except sym40_n18 at L = 256 (P0: η_b = 0.311, pair 0.664, K_s^pair 1.73, K_s^sf 1.50). At L = 64, η_b = 0.252 to 0.271, pair 0.66 to 0.76, and K_c^sf = −0.02 to 0.49 (the estimator that failed on the two-leg control at L = 64), so K_s^sf(64) is negative and meaningless. The pair (d = 3) momentum is Q = π (symmetric line) and π − 0.010 to π − 0.014 (φ = 1.53 bonds), and d = 2 gives the same exponent within 0.003 at L ≥ 128. The valley momenta are 0.236 and 2.901 (n = 3/16) and 0.240, 2.901 (1/8) at L = 256.

### R2. The effective model with the same frozen estimators (C176 files) [numerical]

| run | L | η_b | pair (d = 3) | K_c (ES) | K_s^pair | K_s^sf |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| c = 2 control (c2) | 128 / 256 | 0.502 / 0.501 | 1.010 / 1.002 | 0.968 / 0.968 | 1.00 / 1.00 | 1.03 / 1.03 |
| n = 3/16 (sym316) | 64 | 0.315(10) | 0.577(38) | 0.41(20) | 1.47(12) | −0.8 |
| | 128 | 0.320(5) | 0.614(15) | 1.27(15) | 1.50(6) | 2.0(4) |
| | 256 | 0.326(5) | 0.623(15) | 1.789(74) | 1.47(5) | 1.34(6) |
| n = 1/8 (sym18) | 64 | 0.311(16) | 0.610(46) | 0.07(14) | 1.57(19) | −0.1 |
| | 128 | 0.347(5) | 0.613(23) | 1.28(20) | 1.29(5) | 1.65(34) |
| | 256 | 0.352(5) | 0.617(15) | 1.807(76) | 1.26(4) | 1.17(4) |

The frozen estimators reproduce C192's numbers on the effective model (η_b 0.326 and 0.351, pair 0.62, K_c ≈ 1.8 at L = 256). K_c from the structure factor at fixed R/L rises steeply from L = 128 to 256 in the effective model too (1.27 → 1.79), so that L dependence belongs to the estimator in a K_c ≈ 1.8 liquid and not to the qubit chain.

Estimator for estimator at L = 256: pair exponent qubit 0.634 to 0.643 against 0.623 (3/16) and 0.643 against 0.617 (1/8), K_c^sf 1.70 to 1.73 against 1.79 and 1.81, η_b 0.282 to 0.286 against 0.326 (3/16) and 0.297 against 0.352 (1/8). So the pair and the charge stiffness agree to 2 to 4%, while η_b at fixed r/L reads 0.04 to 0.055 lower on the qubit chain. Stage one (r = 8 to 34) saw only a 0.011 to 0.018 difference. R4 below traces this to the open ends.

### R3. Changes of K_s and the KT bound, ground branch, Δl = ln 2 [numerical]

| key | estimate | K_s 64 / 128 / 256 | change 64 → 128 | change 128 → 256 | |y_1| central / upper | y_∥ central / lower |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| p153c | pair | 2.55 / 2.14 / 2.05 | −0.41(69) | −0.09(46) | 0.25 / 0.60 | 2.09 / 1.55 |
| sym40_n316 | pair | 2.45 / 2.09 / 2.00 | −0.36(42) | −0.09(35) | 0.26 / 0.55 | 1.99 / 1.58 |
| rat_n316 | pair | 2.50 / 2.08 / 2.00 | −0.42(65) | −0.08(44) | 0.23 / 0.60 | 2.00 / 1.48 |
| sym40_n18 | pair | 2.54 / 1.93 / 1.84 | −0.62(44) | −0.09(34) | 0.27 / 0.59 | 1.67 / 1.27 |
| p153_n316 | pair | 2.50 / 2.08 / 1.99 | −0.42(66) | −0.09(44) | 0.25 / 0.61 | 1.97 / 1.46 |
| sym80_n316 | pair | 2.40 / 2.08 / 1.99 | −0.32(40) | −0.09(32) | 0.25 / 0.53 | 1.98 / 1.59 |
| p153c | sf | −1.0 / 4.1 / 1.85 | +5.1(2.4) | −2.3(2.3) | 0.86 / 1.21 | 1.71 / 1.24 |
| sym40_n316 | sf | −0.6 / 3.5 / 1.80 | +4.1(1.5) | −1.7(1.4) | 0.83 / 1.13 | 1.60 / 1.24 |
| rat_n316 | sf | −0.8 / 3.7 / 1.82 | +4.5(1.9) | −1.9(1.8) | 0.84 / 1.18 | 1.64 / 1.19 |
| sym40_n18 | sf | −0.1 / 3.0 / 1.64 | +3.1(1.2) | −1.4(1.2) | 0.86 / 1.18 | 1.28 / 0.94 |
| p153_n316 | sf | −0.8 / 3.7 / 1.81 | +4.5(1.9) | −1.9(1.8) | 0.85 / 1.18 | 1.62 / 1.18 |
| sym80_n316 | sf | −0.6 / 3.4 / 1.77 | +4.1(1.5) | −1.7(1.4) | 0.84 / 1.14 | 1.54 / 1.21 |

Per sector (analysis.txt) the pair rows look the same. The tightest pair-estimate case is sym40_n18 P0 (|y_1| ≤ 0.66 against y_∥ ≥ 1.17). In the sf rows the 64 → 128 change is a rise from a negative, uncalibrated K_s^sf(64), and the 128 → 256 errors are dominated by K_s^sf(128) ± 1.2 to 2.3.

### Registered outcome

- B0 holds (two-leg control η_b = 0.504, 0.502, pair 1.014, 1.001, K_c 1.026, 1.039 at L = 128, 256).
- B1 fails, at one key. The five n = 3/16 keys have η_b within 0.041 to 0.044 of 0.326, but sym40_n18 (n = 1/8) has η_b = 0.297 against 0.351, a difference of 0.054 > 0.05. Every pair exponent lies in 0.50 to 0.75 (0.634 to 0.643).
- B2. With 1/K_c from the pair exponent all three clauses hold at every key and in both sectors (K_s^pair(256) ≥ 1.73, y_∥ lower limit 1.17 to 1.64 above the |y_1| upper limit 0.53 to 0.66, and the change 128 → 256 of −0.06 to −0.20 is less negative than the change 64 → 128 of −0.32 to −2.2). With 1/K_c from the structure factor the first clause holds (K_s^sf(256) ≥ 1.50) but the second fails at sym40_n18 (y_∥ ≥ 0.94 against |y_1| ≤ 1.18) and p153_n316 (1.18 against 1.18) on the ground branch, and also in the other sector of rat_n316 (1.18 against 1.20) and sym40_n18 (0.77 against 1.24), and the third fails at every key (+3 to +5 then −1.4 to −2.3). By the rule fixed before the candidates, B2 needs both, so B2 is not met.
- B3. With the pair estimate no clause holds. With the structure-factor estimate the second clause holds (the bound allows |y_1| ≥ y_∥ at sym40_n18 and p153_n316) and the third holds at every key (the 128 → 256 fall exceeds the 64 → 128 change by 4.5 to 7.4, against combined errors of 1.7 to 3.3). By the rule fixed before the candidates (B3 if either estimate meets it), B3 is met.

The registered letter is B3, carried entirely by the structure-factor estimate of 1/K_c. With the pair estimate alone the letter would be B2.

### What I think it means, stated separately

1. The B3 clauses rest on K_s^sf at L = 64 and 128. The calibration section recorded, before any candidate output was opened, that the structure-factor estimator fails on the two-leg control at L = 64 (1.48 for exact 1). On the candidates it gives K_c^sf(64) ≈ 0.4 and a negative K_s^sf(64), which is not a physical value, and K_s^sf(128) = 3 to 4 with errors of 1.2 to 2.3. The same estimator on the N-conserving effective model gives the same sequence (K_c^sf 0.41, 1.27, 1.79; K_s^sf −0.8, 2.0, 1.34), where no pseudospin gap is at stake. So the sf "acceleration" is the finite-size behaviour of the K_c estimator at fixed R/L in a K_c ≈ 1.8 liquid, shared by both models. At L = 256 the sf route gives K_s^sf = 1.64 to 1.85, above 1 everywhere.
2. B1's failure and the large K_s^pair ≈ 2.0 come from the same source, the qubit chain's η_b at fixed r/L. That η_b is biased by the open ends (R4). With the ends treated correctly the qubit chain's η is 0.325 to 0.330 (3/16) and 0.353 (1/8), equal to the effective model's 0.326 and 0.352, and K_s ≈ 1.46 to 1.52 and 1.30 (R4).
3. With the pair estimate, which is calibrated at all three sizes, the qubit chain's K_s falls by about 0.4 from 64 to 128 and by 0.09 from 128 to 256, a flow that slows by a factor of about 4 per doubling, and the bound puts |y_1| ≤ 0.6 against y_∥ ≥ 1.2. This reads as the gapless side. I do not think the registered B3 is evidence for a gap. It is the uncalibrated L = 64 structure factor feeding a drift test.

### R4. The open ends of the qubit chain pin the charge phase (unregistered, after the candidates were opened) [numerical, heuristic]

Diagnostics: diag193.py, outputs diag_eta.txt, diag_mixed.txt, diag_mixed_table.txt, diag_pp.txt.

- The local η_b of the qubit chain falls with r/L at fixed L: at L = 256 (sym40_n316) 0.322, 0.296, 0.274, 0.219, 0.131 on (0.05, 0.10), (0.10, 0.30), (0.15, 0.40), (0.20, 0.60), (0.30, 0.80) L, and at L = 128 the ranges from 0.10 L on give the same numbers within 0.01 (0.292, 0.263, 0.208, 0.125), so the effect is a function of r/L. The effective model rises instead (0.319, 0.326, 0.327), and both controls are flat at 0.500 ± 0.003 out to 0.8 L. The fit rms is 0.013 to 0.025 on the qubit chain against 0.002 on the effective model.
- Mechanism (hypothesis). The qubit chain conserves only parity. Its pair-creation term Δ(c†c† + h.c.) is removed in the bulk by the Schrieffer–Wolff rotation and leaves no relevant bulk pair field (the inter-valley pair carries momentum π). At an open end momentum is not conserved, so a boundary pair field e^{i√2θ_c} survives, and with θ free at the end its boundary dimension is 1/K_c ≈ 0.6 < 1, so it is relevant and pins θ_c (Dirichlet). The pseudospin field θ_s stays free (no boundary term changes it, by parity). Then the boson correlator for pairs symmetric about the centre is tan(u)^{−η} cos(u)^{−2η_c} with u = πr/(2(L + 1)), η_c = 1/(4K_c), instead of tan(u)^{−η} (PROOF §3 with the image sign of θ_c flipped). The effective model conserves N and keeps θ_c free at its ends.
- Test. A two-exponent fit of that form on (0.10, 0.80) L gives, at L = 256, η = 0.325 to 0.330 at the five n = 3/16 keys and 0.353 at n = 1/8, equal to the effective model's 0.326 and 0.352, with η_c = 0.16 to 0.17 (3/16) and 0.21 (1/8), against 1/(4K_c) = 0.14 to 0.16 from K_c^sf or the pair exponent. The same fit returns η_c = −0.002 to +0.006 on both controls and on the effective model (L = 256). The charge part is larger than predicted at n = 1/8 (0.21 against 0.15), and the residual rms stays at 0.015, so the form captures most but not all of the end effect.
- Supporting. ⟨b_j b_{j+3}⟩ has a q = 0 part of 0.009 (V/Δ = 40) that halves to 0.0045 at V/Δ = 80, the O(Δ/V) Schrieffer–Wolff dressing of PROOF §2, and a momentum-π part of 0.002 to 0.016 in the bulk. Its profile in x is too noisy to test the predicted [L sin(πx/L)]^{−1/(2K_c)}.
- With this η and the frozen pair exponent, K_s at L = 64 / 128 / 256 on the ground branch is 4.6(3.2) / 1.63(10) / 1.52(8) (p153c), 3.0(1.0) / 1.55(9) / 1.47(7) (sym40_n316), 3.6(1.7) / 1.54(9) / 1.47(7) (rat_n316), 2.3(6) / 1.32(10) / 1.30(6) (sym40_n18), 3.6(1.7) / 1.54(9) / 1.46(7) (p153_n316), 2.9(9) / 1.55(8) / 1.47(8) (sym80_n316). The change 128 → 256 is −0.02 to −0.12 (σ ≈ 0.11), and the bound gives |y_1| ≤ 0.47 to 0.53 against y_∥ ≥ 0.78 to 0.88 at n = 3/16. At n = 1/8 it is a tie (|y_1| ≤ 0.477 against y_∥ ≥ 0.478). These K_s agree with the effective model's 1.47 (3/16) and 1.26 (1/8) at L = 256. L = 64 is not usable in this fit (the end region is a large part of the chain).

## Rerunning on new states (for example the χ = 1024 states of job 17609630)

Measurement, on Typhon, from `/scratch/lustre/jkudlerflam/c19x/c193/` (needs meas193.py, c193.slurm and lib/ there, as now). Write a list with one `<run.json> <ckpt.npz>` per line and submit it with its own output folder, for example

    cd /scratch/lustre/jkudlerflam/c19x/c193
    for t in p153c_L256_P1 sym40_n316_L256_P1; do f=$(ls /scratch/lustre/jkudlerflam/c18x/c190g/outX1024/${t}_*.json | head -1); echo "$f /scratch/lustre/jkudlerflam/c18x/c190g/outX1024/ckpt/$(basename $f .json).npz"; done > states_x1024.txt
    sbatch c193.slurm states_x1024.txt out_x1024

The tag names in outX1024 were not checked (the job had not started), so check the list before submitting. Use only runs whose JSON has status done. The script reads the checkpoint and the run JSON only, skips finished outputs, and writes `<outdir>/<tag>_meas.json` with the state checks in its `check` field. At χ = 1024 the state is about 8.6 GB on the device, which fits the 40 GB A100.

Collection and analysis, locally from this folder:

    mkdir -p typhon_out/out_x1024 && scp 'typhon:/scratch/lustre/jkudlerflam/c19x/c193/out_x1024/*_meas.json' typhon_out/out_x1024/
    python3 -E -s ana193.py states typhon_out/out_x1024/*_meas.json

`states` applies the frozen estimators (frozen.json) to any measured state and prints η_b, both pair exponents, K_c and both K_s with errors and the state checks. For the boundary-corrected η of R4: `python3 -E -s diag193.py <tag> ...` with the tag as it appears in typhon_out/out (edit the path in diag193.py for another folder).
