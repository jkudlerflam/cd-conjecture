This run was done by a Claude Opus 5.5 agent (C193, 2026-11-10), not by Codex. Astra was not used.

# C193 result: K_s on the qubit chain from the stored C190 correlators

Derivations are in PROOF.md, the log in PROGRESS.md, the script is c193.py.

## Registration (copied from specs/C193_C194_sign_tests.md before any candidate correlator was opened)

- (Q0) Calibration. With the frozen ranges the two-leg control reads η_b within 0.03 of 1/2 and the pair exponent within 0.08 of 1 at L = 128 and 256. If not, everything below is uncalibrated.
- (Q1) Tracks the effective model. At L = 256 on the ground branch, every n = 3/16 key has η_b within 0.05 of 0.326 and the n = 1/8 key within 0.05 of 0.351, and every pair exponent lies in 0.50 to 0.75.
- (Q2) Gapless sign. K_s ≥ 1.1 at L = 128 and 256 for every key in both sectors, K_s(256) − K_s(128) ≥ −0.05, and the drift bound gives y_∥ > |y_1| at every key.
- (Q3) Gap sign. Any key has K_s(256) < 1, or K_s falls by more than 0.1 from 128 to 256, or the drift bound allows |y_1| ≥ y_∥.
- Anything between Q2 and Q3 is reported as undecided, with the numbers.
- If the 48-site window is too short for a stable η_b on the controls, say so, stop the candidate analysis for that estimator, and give the exact command that would recompute the correlators from the Typhon checkpoints under `/scratch/lustre/jkudlerflam/c18x/c190g/` (read only) in a SLURM job named c193.

## Addendum to the registration, written after reading the measurement code and before opening any candidate correlator

Step 1 of the spec (PROOF.md §1 and §2) found that the stored C_P is not an inter-valley pair function. `z2dmrg.measure` stores the connected correlator of the bond operator P_j = X_jX_{j+1} + Y_jY_{j+1} = 2(σ⁺_jσ⁻_{j+1} + h.c.), the nearest-neighbour hop, from site i₀ = L/4. In the effective model this is the valley-polarization operator (2 cos k in momentum space). In C191's frame its components are (PROOF §2):
- the smooth valley polarization ∂φ_s, exponent 2, amplitude ∝ K_s (no exponent information);
- the inter-valley exchange e^{±i(θ_A − θ_B)} = e^{±i√2 θ_s} at momentum k_B − k_A, exponent 1/K_s, with a lattice amplitude ∝ cos((k_A + k_B)/2), which vanishes on the symmetric line (k_A + k_B = π) and is of order (Q − π)/2 ≈ 0.005 off it;
- density harmonics, exponents above 1.5.
So the inter-valley pair exponent 1/K_c of the registration cannot be measured from the stored data. On the two-leg control, C_P is the inter-leg exchange b†_A b_B, whose exponent is 1 (= 1/K_s = 1/K_c there), so the control calibration of the C_P estimator is still meaningful as a fitting check.

How the registered quantities are evaluated, fixed now.
- η_b: from Gb, estimator E1 (PROOF §3), the frozen ranges of the calibration section below.
- "Pair exponent" for Q0: the dominant power-law component of C_P on the two-leg control (estimator E2).
- 1/K_c for K_s and Q1: not measurable here. The primary K_s uses C192's effective-model value at matched density, 1/K_c = 0.62, with the band 0.55 to 0.64 (structure factor and pair estimates, §9.161) as a systematic error, K_s = 1/(4η_b − 1/K_c). Q1's pair-exponent clause is reported as not evaluable, not as passed.
- E2 on the candidates is run with the same frozen ranges, and its result is read by this rule, fixed before looking. If the dominant slow component of C_P sits within 0.1 of the boson valley-momentum difference |k_B − k_A| (mod 2π) and its exponent is below 1.2, it is read as the exchange vertex, and K_s^ex = 1/η_ex is reported as an independent K_s, with an inferred 1/K_c = 4η_b − η_ex. Otherwise C_P is reported as carrying no exponent information on K_s.
- Drift. The registered drift is K_s(256) − K_s(128) at the frozen fit ranges. Because the 48-site window probes the same separations r ≲ 40 at both sizes, this is not a doubling of the RG scale (PROOF §5). As a secondary drift, registered now, the local exponent at L = 256 is compared between the shortest and longest frozen ranges, a scale change of Δl = ln(r̄_long/r̄_short).
- Ground branch: at each L, the parity sector with the lower final energy (last history record), and where two initializations exist the lower-energy one.

(Calibration section and results follow below; nothing above was changed after candidate correlators were opened.)

## Calibration and frozen ranges (controls only; written before any candidate correlator was opened)

`python3 -E -s c193.py controls` (output controls.json). Exact values: two-leg XX (c = 2) η_b = 1/2, C_P exponent 1 (inter-leg exchange, PROOF §2). Single XX chain (c = 1, K = 1): hard-core boson η_b = 1/2 at k = 0, Jordan–Wigner fermion exponent 1 at ±k_F, C_P (kinetic-energy density) exponent 2 in both its q = 0 and 2k_F parts.

E1 (boson, 48-site window, centres averaged over ±4, four binomial passes), η_b on the four trial ranges R1 = [4, 12], R2 = [6, 20], R3 = [8, 28], R4 = [12, 34]:

| control | L | R1 | R2 | R3 | R4 |
| :-- | :-- | :-- | :-- | :-- | :-- |
| two-leg XX | 128 | 0.529 | 0.523 | 0.501 | 0.502 |
| two-leg XX | 256 | 0.522 | 0.525 | 0.502 | 0.501 |
| XX chain | 128 | 0.537 | 0.514 | 0.507 | 0.504 |
| XX chain | 256 | 0.537 | 0.514 | 0.507 | 0.504 |
| two-leg XX | 64 (W = 32) | 0.465 | 0.534 | 0.511 | 0.492 |

The short ranges carry a lattice bias of +0.015 to +0.037. Frozen: η_b = mean over R3 and R4 (and over the two largest valley peaks), error = half the spread, plus a systematic 0.005 (the largest R3/R4 control bias). At L = 64 (W = 32) the analogous ranges [6, 16] and [8, 18] are used, as a check only.

E2 (C_P from i₀ = L/4), trial ranges [6, L/8], [10, L/5], [16, L/3], [24, 0.45L]:

| control | L | q found | 1 | 2 | 3 | 4 |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| two-leg XX | 128 | 0 | 1.057 | 0.990 | 1.001 | 1.005 |
| two-leg XX | 256 | 0 | 1.037 | 1.005 | 1.003 | 1.004 |
| two-leg XX | 64 | 0 | 1.240 | 0.767 | 0.862 | 1.517 |
| XX chain | 128 | 0 | 1.94 | 1.60 | 1.43 | 1.15 |
| XX chain | 256 | 0 | 2.27 | 2.26 | 2.27 | 2.19 |

Frozen: E2 = mean over ranges 2 to 4, error = half spread. The XX chain shows that E2 fails when two components of equal exponent (q = 0 and 2k_F) compete, since the binomial filter suppresses 2k_F ≈ 1.24 only to 0.44 (rms residual 0.17 to 0.28 in log). E2 is calibrated only for a single dominant component, and L = 64 is not usable for E2.

The Jordan–Wigner fermion (Gf) reads 0.96 to 1.04 on the XX chain at R3/R4 but with rms 0.13 in log, for the same filter reason. It is used only as a qualitative frame check.

Q0 (registered). Two-leg control, frozen ranges: η_b = 0.502 (L = 128) and 0.501 (L = 256), within 0.03 of 1/2. C_P exponent 0.999 (L = 128) and 1.004 (L = 256), within 0.08 of 1. Q0 holds, for the C_P estimator in the restricted sense of the addendum.

Window-length check on the effective model (not candidate data; `effcheck.py`, C176 runs cut to their central 48 sites). E1 at R3, R4 reads 0.312, 0.318 (L = 128) and 0.314, 0.320 (L = 256) at n = 3/16, against C192's full-window 0.325 and 0.326, and 0.322, 0.341 and 0.327, 0.340 at n = 1/8 against 0.347 and 0.351. The c = 2 control there reads 0.503 to 0.505. So the 48-site window reads η_b low by about 0.010 at 3/16 and 0.011 to 0.024 at 1/8 in the effective model, and the local exponent rises with r. This is a property of the candidate liquid, not of the window (the controls are flat), and it is relevant for the drift (PROOF §5). For a matched comparison with the qubit chain the 48-site effective-model values are used alongside C192's numbers.

Frozen. Nothing above is changed after this point.

## Results (candidate correlators opened after the section above was written)

Commands, from this folder: `python3 -E -s c193.py candidates` (candidates.json), `python3 -E -s drift.py`, `python3 -E -s effcheck.py <C176 tags>`. Inputs are the 51 GPU runs in `codex/work_C190/gpu/typhon_out/outG/`, all with status done. Where two initializations exist, they agree in energy to below 1e-9 and the lower is used.

### R1. Exponents per key [numerical]

η_b from E1 (frozen R3, R4, both valleys), error = half spread with 0.005 systematic. K_s = 1/(4η_b − 0.62), with the band for 1/K_c = 0.55 to 0.64 in brackets. * marks the ground branch.

| key | n_mid (256) | L = 128 P0 | L = 128 P1 | L = 256 P0 | L = 256 P1 |
| :-- | :-- | :-- | :-- | :-- | :-- |
| p153c | 0.219 | 0.2934(81)*, K_s 1.81 [1.60, 1.87] | 0.2947(78), 1.79 | 0.3098(101), 1.62 | 0.3088(102)*, 1.63 [1.46, 1.68] |
| sym40_n316 | 0.196 | 0.3027(50), 1.69 | 0.2994(50)*, 1.73 [1.54, 1.79] | 0.3151(60), 1.56 | 0.3152(60)*, 1.56 [1.41, 1.61] |
| rat_n316 | 0.196 | 0.2975(73), 1.75 | 0.2981(73)*, 1.75 [1.56, 1.81] | 0.3134(93), 1.58 | 0.3124(95)*, 1.59 [1.43, 1.64] |
| sym40_n18 | 0.133 | 0.3188(50)*, 1.53 [1.38, 1.57] | 0.3177(52), 1.54 | 0.3386(70), 1.36 | 0.3327(70)*, 1.41 [1.28, 1.45] |
| p153_n316 | 0.196 | 0.2983(73), 1.74 | 0.2980(74)*, 1.75 [1.56, 1.81] | 0.3132(94)*, 1.58 [1.42, 1.63] | 0.3126(95), 1.59 |
| sym80_n316 | 0.195 | 0.2999(51), 1.73 | 0.2974(51)*, 1.76 [1.56, 1.82] | 0.3130(56)*, 1.58 [1.42, 1.63] | 0.3139(56), 1.57 |

Error on K_s from η_b alone: 0.05 to 0.11. L = 64 (W = 32, check only): η_b = 0.262 to 0.270 at every key, K_s ≈ 2.2 to 2.3.

- Valleys. Boson peaks at |k| = 0.236 and 2.896 (n = 3/16) and 0.241, 2.901 (1/8) at L = 256, |k_B − k_A| = 2.66 at every key, equal valley exponents on the symmetric line, and a 0.008 to 0.010 difference between the valleys for the φ = 1.53 bonds.
- Jordan–Wigner fermion, frame check. Gf reads 1.7 to 2.0 at L = 128 and 256 (rms large, PROOF §3), above 1, while η_b ≈ 0.3. The boson frame of C191 holds on the qubit chain as well.
- C_P (E2). The dominant component is at q = 0 at every key and size, with exponent 2.04 to 2.09 at L = 128 and 256 and a negative sign, the smooth valley polarization of PROOF §2. By the registered rule (q not within 0.1 of |k_B − k_A| = 2.66), C_P carries no exponent information on K_s, and no exchange exponent or inferred 1/K_c is reported.

### R2. Against the effective model [numerical]

At L = 256 on the ground branch, η_b = 0.309, 0.315, 0.312, 0.313, 0.313 at the five n = 3/16 keys against C192's 0.326 (differences −0.011 to −0.017), and 0.333 at n = 1/8 against 0.351 (−0.018). With the same 48-site window the effective model itself reads 0.314 to 0.320 (3/16) and 0.327 to 0.340 (1/8), so at L = 256 the qubit chain and the effective model agree to 0.005 at equal window.

At fixed window the effective model does not move with L (W = 32: 0.298, 0.299, 0.300 at L = 64, 128, 256 for n = 3/16; W = 48: 0.312/0.318 at 128 and 0.314/0.320 at 256). The qubit chain does: η_b ≈ 0.265 (64), 0.298 (128), 0.313 (256) at 3/16. So the qubit chain approaches the effective-model value from below as L grows, and at L = 256 it has reached it within the error. I do not know the mechanism. χ is not it at 64 → 128 (both converged, C4 passes), but it may contribute at 256, where C4 fails (§9.178).

Q1 (registered). The η_b clause holds at every key (all within 0.05). The pair-exponent clause cannot be evaluated, since the stored data have no pair correlator (addendum). Q1 is therefore met only in its η_b part.

### R3. Drift and the KT bound [numerical; PROOF §5]

Ground branch, Δl = ln 2. |y₁| is from −ΔK_s (central, and with ΔK − σ), y_∥ = 2(K_s(256) − 1) (central, and lowered by its error).

| key | K_s(128) | K_s(256) | ΔK_s (σ) | |y₁| central / upper | y_∥ central / lower |
| :-- | :-- | :-- | :-- | :-- | :-- |
| p153c | 1.81 | 1.63 | −0.18 (0.15) | 0.42 / 0.57 | 1.25 / 1.03 |
| sym40_n316 | 1.73 | 1.56 | −0.17 (0.08) | 0.43 / 0.52 | 1.12 / 1.00 |
| rat_n316 | 1.75 | 1.59 | −0.16 (0.13) | 0.41 / 0.55 | 1.18 / 0.99 |
| sym40_n18 | 1.53 | 1.41 | −0.12 (0.07) | 0.40 / 0.51 | 0.81 / 0.70 |
| p153_n316 | 1.75 | 1.58 | −0.17 (0.13) | 0.42 / 0.56 | 1.16 / 0.97 |
| sym80_n316 | 1.76 | 1.58 | −0.17 (0.08) | 0.43 / 0.52 | 1.16 / 1.05 |

Per sector the drifts are −0.13 to −0.19, and the tightest case is sym40_n18 P0 (|y₁| ≤ 0.57 against y_∥ ≥ 0.62). The errors are dominated by the fit-range spread, which is largely common to both sizes, so σ is an overestimate of the error on ΔK_s.

Secondary drift (registered in the addendum), L = 256, R3 → R4 (Δl = 0.30): η_b rises by 0.005 to 0.010 at every key, dK_s/dl = −0.17 to −0.26, |y₁| ≈ 0.36 to 0.51. The effective model at the same window shows the same rise (0.0055 at 3/16, 0.0135 at 1/8), and the controls are flat (0.001).

Effective-model cross-check with C192's full-window numbers. n = 3/16: K_s = 1.47 → 1.46 at fixed 1/K_c (|y₁| ≈ 0.11, at most 0.21), or 1.50 → 1.46 with C192's measured pair exponents 0.635 and 0.619 (|y₁| ≈ 0.24). n = 1/8: 1.30 → 1.28 (|y₁| ≈ 0.22). y_∥ = 0.92 and 0.55. Gapless by the same criterion.

Unregistered secondary (PROOF §6). The smooth C_P amplitude at fixed exponent 2 gives K_s(256)/K_s(128) = 1.06, 1.05, 1.03, 1.04, 1.03 at the five n = 3/16 keys and 0.95 at n = 1/8. So this estimator sees no fall at 3/16, opposite in sign to the η_b drift, and a 5% fall at 1/8. It has no control calibration.

### Registered outcome

- Q0 holds (two-leg control η_b = 0.502, 0.501, C_P exponent 0.999, 1.004).
- Q1 holds in its η_b part only. The pair clause is not evaluable.
- Q2 fails, on its drift clause alone: K_s ≥ 1.1 everywhere (lowest 1.36, sym40_n18 P0 at L = 256) and y_∥ > |y₁| at every key, but K_s(256) − K_s(128) = −0.12 to −0.18 on the ground branch, below −0.05.
- Q3 holds, on its second clause only: K_s falls by more than 0.1 from 128 to 256 at every key on the ground branch (−0.12 to −0.18, with σ = 0.07 to 0.15, so the fall exceeds 0.1 by less than one σ at every key). Its other two clauses fail: K_s(256) ≥ 1.36 > 1 everywhere, and the drift bound gives y_∥ > |y₁| everywhere, also at the error edges.

The registered letter is Q3.

### Why Q3 here does not read as a gap, stated separately

1. The registered drift compares the same separations (r = 8 to 34) at two L, so it is not a flow over ln 2 (PROOF §5). Applied to the effective model, the same estimator at fixed window gives no L drift at all (0.312 → 0.314), while its long-window value is flat too (C192). The qubit chain's L drift is an approach to the effective-model value, which it reaches within 0.005 at L = 256.
2. Read as a KT flow anyway, the drift implies |y₁| ≈ 0.4 against y_∥ ≈ 0.8 to 1.25, inside the gapless basin, and with 2K_s − 2 ≈ 0.8 to 1.2 the remaining drift is about −0.1 (−0.16 at 1/8), ending near K_s ≈ 1.45 to 1.5 (3/16) and 1.25 (1/8).
3. The independent C_P amplitude estimator (uncalibrated) sees K_s rise by 3 to 6% at 3/16.

What would change this. If the L drift continued at −0.17 per doubling without slowing, K_s would reach 1 after about 2 to 3.5 more doublings, near L ≈ 1000 to 3000, at both densities. Neither these data nor the 48-site window can exclude that. The test is η_b at L = 256 with a longer window (r up to 100), and whether the qubit value stays at the effective-model value from L = 256 on.

### What is left

- No inter-valley pair correlator exists for the qubit chain, so 1/K_c is borrowed from the effective model. A measurement needs new code (a four-point Dpair, as in C176's measurement) run on the C190 checkpoints. Command sketch, not written as a job and not submitted: a c193 SLURM job that loads the L = 256 checkpoints under `/scratch/lustre/jkudlerflam/c18x/c190g/` (read only; file names not checked), measures the full-chain Gb on a 160-site window and Dpair at separations 2 and 3, writes to `/scratch/lustre/jkudlerflam/c19x/c193/`. The fallback clause of the spec does not apply (the window gives a stable η_b on both controls), so I did not prepare or submit it.
- The cause of the qubit-chain L dependence of η_b at fixed r (64 → 128 → 256) is not identified. The χ = 1024 stage at L = 256 (§9.174) will show whether part of the 128 → 256 step is χ.
- The C_P amplitude estimator needs a single-component control.
