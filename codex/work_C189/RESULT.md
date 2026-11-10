This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

# C189 result: is the pseudospin backscattering relevant at the candidate densities?

Derivations are in PROOF.md, the log in PROGRESS.md. Copied unmodified from C176: fibdmrg.py, analyze.py.

## Registration (written 2026-11-09 before any estimator below was evaluated on the candidate series, and before the infinite-chain code existed)

What I had seen before writing this. REPORT 9.139, 9.140, 9.156 to 9.158, C176 RESULT and the C176 summary file (c_two, n_mid, L·Δ_P, and the envelope exponents −2.0 for |G_f| and −0.78 for |D| quoted in REPORT 9.158). I had not fitted any exponent, structure factor or susceptibility myself on any candidate run.

### Conventions

Two-component Luttinger liquid in Giamarchi's normalization, ρ_c = −(√2/π)∂φ_c, ρ_s = −(√2/π)∂φ_s, free fermions at K_c = K_s = 1. The pseudospin is the valley index (d frame valleys at ±k₀, cos k₀ = −1/4 at |w| = 1, φ = π/2). The pseudospin backscattering g₁ cos(√8 φ_s) has scaling dimension 2K_s. At small |g₁| it is relevant iff K_s < 1. The relations used below are derived in PROOF.md, including the point that a hard wall reflects one valley into the other, so the open chain pins φ_c and θ_s (not φ_s).

### Part 1, estimators on the converged C176 data (and C181, C176 off-line runs for the susceptibility)

All from the last line of each jsonl (the final sweep).
- E_c (K_c). Density structure factor S(q) = C(0) + 2 Σ_{r≥1} C(r) cos qr from the four C_nn rows (averaged), slope at small q. S(q) = K_c q/π (two components), K q/(2π) for the spinless hard-rod control.
- E_D (pair). Envelope exponent η_D of |D(r)| (separation 2, from the centre). η_D = K_s + 1/K_c, so K_s^D = η_D − 1/K_c.
- E_F (Friedel). Exponent η_F of the boundary density oscillation at 2(k₀ ± k_F) (d frame). η_F = (K_c + 1/K_s)/2, so K_s^F = 1/(2η_F − K_c).
- E_G (one-body). Envelope exponent η_G of |G_f(r)| (c-frame fermion with the full Jordan–Wigner string). η_G = (K_c + 1/K_c + K_s + 1/K_s)/4. It fixes K_s + 1/K_s only, so it measures |ln K_s| and cannot tell K_s from 1/K_s.
- E_χ (susceptibility against the parity gap). (2/π) L·Δ_P = v_c/K_c + v_s K_s (open chain, θ_s pinned), v_c/K_c = (2/π) ∂²e/∂n², and the bulk pseudospin susceptibility χ_s = ∂(n_a − n_b)/∂(δ/2) = 2K_s/(π v_s). The product gives K_s² = (v_s K_s)(K_s/v_s), free of v_s. Caveat registered now: n_a − n_b is read from the n(k) valley weights, which equal the low-energy N_s only up to a factor Z_N that is 1 for free fermions.

Controls, same files and same code. c1 (hard rods, spinless, K = (1 − n)² = 0.5625 at n = 1/4), c2 (two XX legs, K_c = K_s = 1, ordinary walls), and a third control I add because it has the candidate's valley structure and boundary: ff2v, free fermions with the same band at φ = π/2 (no blockade, no interactions), exact Slater-determinant open chains at the candidate's L and N, K_c = K_s = 1. An estimator counts only if it reads every applicable control within 0.10 of the exact K (absolute). Fit windows (ranges of q, r and boundary distance) are fixed on the controls only, and written below as an amendment before the candidate series is evaluated.

Uncertainty of each K_s value: half the spread over three fit windows (the registered one and the two neighbouring ones), added in quadrature to the control bias of that estimator.

Decision rule for relevance at n = 3/16 (L = 256 primary, 5/32 and 1/8 as checks):
- Relevant: the passing estimators that measure K_s itself (E_D, E_F, E_χ) give K_s + 2σ < 1 in at least two of them, and none gives K_s − 2σ > 1.
- Irrelevant: at least two give K_s − 2σ > 1, none gives K_s + 2σ < 1, and the drift from L = 128 to 256 does not carry the estimate below 1 when extrapolated linearly in 1/ln L.
- Marginal or undecided otherwise. Then the drift decides the reading: K_s(L) approaching 1 from above like 1/ln L is the marginally irrelevant (gapless) side, from below the marginally relevant (exponentially small gap) side.
- E_G is a consistency check: it must agree with K_c and the K_s found, through K_s + 1/K_s, within its 2σ.

### Part 2, infinite chain (registered before the code was written)

A blockade-exact, U(1) VUMPS or iDMRG for C173's model at fixed density n = p/q (q-site unit cell), validated by: (i) the exact energy density of the two controls and of ff2v at χ = 256 to better than 1e-4 relative; (ii) agreement with C176's open-chain bulk energy density at n = 3/16 (from E(256) − E(128) at fixed n, which has an O(1/L²) bias I will state).

Runs: n = 3/16 and 1/8 at φ = π/2, |w| = 1, plus c1, c2 and ff2v, over a χ sequence (64 up to the largest that finishes). For each χ: the half-chain entropy S at a unit-cell boundary, and the transfer-matrix spectrum resolved by the charge ΔQ carried across the cut (ΔQ = 0, 1 with and without the fermion string, 2), giving per-site correlation lengths ξ₀, ξ₁ᵇ, ξ₁ᶠ, ξ₂, and ξ = max.

- c_fit = 6 × slope of S against ln ξ over the largest four χ.
- Calibration: c1 within 0.15 of 1, c2 and ff2v within 0.25 of 2. If not, everything below is uncalibrated.
- c = 2 liquid reading: c_fit ≥ 1.6, and ξ₁ᶠ keeps growing at the rate of ξ₂ (ξ₁ᶠ/ξ₂ changes by less than a factor 1.5 over the last χ doubling).
- Pseudospin gap reading: ξ₁ᶠ grows by less than 10 percent over the last χ doubling while ξ₂ or ξ₀ grows by more than 30 percent, or c_fit < 1.4 with the local slope falling toward 1/6.
- Anything else is inconclusive.
- Stated limit. A pseudospin correlation length longer than the largest ξ reached is invisible to this test. C175's estimate (gap about 1e-4, ξ_s of order 1e4) is far beyond it, so a c = 2 reading here bounds the gap from above at about v_s/ξ_max, and only Part 1 can speak to smaller gaps.

### Amendment A1 to Part 2 (2026-11-09, written before the Typhon job 17606815 started)

What I had seen. Local runs of the controls c1 (χ = 16 to 48) and ff2v (χ = 16 to 32), and of the candidate at n = 1/8 for χ = 16 and 24 only (local/out/), all far below the registered χ range.
- The code is an infinite DMRG (McCulloch growth with a two-cell window, idmrg.py), not VUMPS. A single-site VUMPS with C176's mixed labels stalls at a finite error, because the left-physical and right-physical label dimensions of a bond cannot be fixed independently (PROGRESS, ivumps.py, kept as a record). The two-site SVD of iDMRG sets both consistently.
- Entropy estimator. The converged finite-χ state of ff2v alternates between two values of S at consecutive iterations (it is 2q-periodic, a finite-χ density wave). S̄ is therefore the mean of the last two iterations' S, which samples both inequivalent cell-boundary bonds. ξ comes from the two-cell transfer matrix T[A_t] M[C_t] T[B_t] M[C_{t−1}⁺], which closes on one bond basis.
- The ΔQ = 1 sector with and without the string has identical |eigenvalues| (the string is a diagonal gauge in the charge label, it only shifts the momentum by πn), so ξ₁ᵇ = ξ₁ᶠ ≡ ξ₁.
- c2 cannot be run in the infinite code. With the inter-leg hop exactly zero the leg numbers are separately conserved, and the window never leaves the leg sector of its seed (energy frozen at the first window's value, checked with and without a random Lanczos admixture). The c = 2 calibration of Part 2 is ff2v alone (free fermions with the candidate's band and two valleys, exactly c = 2, K_c = K_s = 1). c1 stays the c = 1 calibration.
- Validation done (validate_inf.out). The first iDMRG window, swept by idmrg.py's own two-site code, reproduces fibdmrg's finite-chain energy to 1e-11 on 8, 16 and 32 sites (real and complex), and S at the middle to 1e-10 (2e-5 in the complex case). The infinite c1 energy density approaches the exact −0.41349667 from above, error 4.5e-6 at χ = 48, and c1 reads c = 1.03, 0.98, 1.00 between consecutive χ from 16 to 48.
- Iteration rule. At each χ the window is iterated until |S_t − S_{t−2}| < 2e-5 and |e_t − e_{t−2}| < 1e-9, with 30 to 300 iterations.

### Amendment A2 to Part 1 (2026-11-09, written after the control calibration and before any estimator was evaluated on a candidate run)

The controls changed two of the registered estimators. Calibration output: kest_controls.out (made by kest_controls.py).
- Boundary condition, correction to the registration and to my first derivation. With NN and NNN hops a hard wall mixes the two valleys only partly. In the exact ff2v open chain the Friedel oscillation is dominated by 2k_F = πn (amplitude 0.24 against 0.03 at 2k₀, L = 1024), i.e. the wall mostly reflects within a valley, which is the ordinary (φ_c, φ_s pinned) boundary. At K_s < 1 that boundary is the stable one (the valley-flip boundary operator e^{i√2θ_s} has boundary dimension 1/K_s > 1), at K_s > 1 the valley-flip boundary is. PROOF 2 is written accordingly.
- E_F is replaced by E_F': the 2k_F = π n_mid Friedel component, fitted with C176-style BCFT chord distance, d(x) = (2(L+1)/π) sin(πx/(L+1)), together with a 2k₀ term, over 8 ≤ x ≤ L/4 from both walls. Exponent η_F = (K_c + K_s)/2, so K_s^F = 2η_F − K_c. Controls at L = 128 to 256: ff2v 0.94 to 0.95 (exact 1), c1 0.54 (exact 0.5625, spinless, 2k_F = 2πn), c2 0.96 to 0.98 (exact 1, k = π/4, 3π/4). Branch rule fixed now: if the fitted 2(k₀ ± k_F) amplitudes in the candidate exceed the 2k_F amplitude, the boundary is the valley-flip one and the original E_F and E_χ formulas apply instead.
- E_χ. With the valley-preserving boundary, L·Δ_P = (π/2)(v_c/K_c + v_s/K_s), and the susceptibility also measures v_s/K_s. The product then carries no information on K_s. E_χ becomes a consistency test of the two-component description: (2/π)L·Δ_P − (2/π)∂²e/∂n² should equal 2/(π χ_s). It is reported, not used for K_s.
- E_D window: 16 ≤ r ≤ L/4, chord distance, envelope over 8 sites. ff2v reads 2.07 at L = 256 (exact 2) but 2.14 at 192 and 2.28 at 128, so E_D passes its control only at L = 256. c1 and c2 have different pair operators and are not applicable. The L-drift of E_D is reported with the ff2v bias at the same L subtracted, as a secondary reading.
- E_G fails its controls (ff2v reads 1.15 to 1.47 against exact 1, c1 1.31 to 1.47 against 1.17), so it is excluded by the registered rule and reported only as an uncalibrated consistency number.
- E_c window: 0.05 ≤ q ≤ 0.3 (below 2k_F at every candidate density). Controls at L = 128 to 256: c1 0.568 to 0.571 (exact 0.5625), c2 0.961 to 0.982, ff2v 0.948 to 0.978.
- Uncertainties: neighbouring windows are (0.05, 0.2) and (0.1, 0.3) for E_c, (8, L/2) and (12, L/4) for E_F', (24, L/4) and (16, 3L/8) for E_D; the control bias is ff2v's at the same L. K_s^F error: √(4σ_F² + σ_c²). K_s^D error: √(σ_D² + σ_c²/K_c⁴).
- With E_χ no longer a K_s estimator, the registered decision rule "at least two of E_D, E_F, E_χ" becomes "both E_D and E_F'". If they disagree beyond their 2σ, the reading is undecided.

### Amendment A3, the frame (2026-11-09, post hoc: written after the registered estimators had been evaluated on the candidates, prompted by the coordinator's frame check)

Everything registered above assumed the fermion frame, two Fermi seas of one fermion species. The registered numbers below show that frame is wrong for this model, and the coordinator independently pointed to the same thing. I re-derived the relations in the boson frame (PROOF 3), where the two valleys are two mutually bosonic quasi-condensates of the hard-core bosons b_j (the effective model's particles never pass, so b and the Jordan–Wigner fermion c describe the same Hilbert space and Hamiltonian, and the frame is a question of which fields are free at the fixed point). The test that separates the frames: in the fermion frame every charge-1 correlation decays at least as fast as 1/r (exponent (K_c + 1/K_c + K_s + 1/K_s)/4 ≥ 1), while the lattice boson in the candidate decays with exponent 0.32 to 0.34 (calibrated on c2, which reads 0.50). The boson-frame estimators below were chosen after seeing the candidate data and are not a registered test. The decision rule (K_s against 1 at 2σ, drift from 128 to 256) is kept as registered.

Frame used by each estimator:
- E_c (S(q)) and K_s^P (pseudospin density, below) involve densities only, so they are frame independent.
- E_D, E_F', E_G as registered are fermion-frame relations. They are reported but superseded.
- K_c^P (pair exponent 1/K_c) and K_s^B (boson exponent 1/(4K_c) + 1/(4K_s)) are boson-frame relations.
- The infinite-chain observables (S, ξ by charge sector) are frame independent; their interpretation through scaling dimensions uses the boson frame.

## Results, Part 1 (C176's converged open chains)

Files: kest_controls.out, kest_run.out (registered estimators), kest_boson.out, kest_cp.out (frame-corrected and pseudospin-density estimators). n = 3/16 unless stated, φ = π/2, |w| = 1.

Registered, fermion frame [numerical].
- E_c, K_c from S(q) on 0.05 ≤ q ≤ 0.3: 1.005, 1.026, 1.020 ± 0.17 at L = 128, 192, 256. Controls pass (c1 0.571 against 0.5625, c2 0.982, ff2v 0.978 at L = 256). But πS(q)/q is not flat in the candidate: it rises from 1.0 at q = 0.4 to about 1.55 at q = 0.06 to 0.1 (flat at 1.0 for ff2v and c2 down to q = 0.15). The long-wavelength K_c is about 1.5 to 1.65. The registered window averages over this curvature and underestimates it.
- E_D: η_D = 0.771 ± 0.074 at L = 256, so K_s^D = η_D − 1/K_c = −0.21 ± 0.18 (5/32: −0.24 ± 0.20; 1/8: −0.32 ± 0.25). Negative, so impossible in the fermion frame.
- E_F': the 2k_F wall oscillation is short ranged in the candidate (amplitude falls by 30 between 6 and 24 sites from the wall), and the wall is dominated instead by a slowly decaying 2k₀ density modulation (0.018 at 16 sites, 0.009 at 96 sites, L = 256). The 2k_F fits are unstable (K_s^F = −0.52 ± 0.21 at 3/16, 3 to 7 ± 2 to 3.6 at the other densities). Uninformative.
- E_G (excluded, fails its controls): fermion envelope exponent 1.93 at L = 256.
- Under the registered rule, E_D and E_F' both give K_s + 2σ < 1 at n = 3/16, which reads "relevant". Both values are negative, so what the registered estimators show is that the fermion-frame two-component Luttinger liquid does not describe the data. I do not take the registered "relevant" reading as a statement about K_s.

Frame-corrected and frame-independent [numerical; relations derived in PROOF 3 and 4].

| n | L | K_c^P = 1/η_P | η_B | K_s^B | K_s^P |
| :-- | :-- | :-- | :-- | :-- | :-- |
| 3/16 | 128 | 1.60 ± 0.09 | 0.324 ± 0.063 | 1.49 | 1.26 ± 0.10 |
| 3/16 | 192 | 1.61 ± 0.07 | 0.331 ± 0.040 | 1.43 | 1.27 ± 0.10 |
| 3/16 | 256 | 1.61 ± 0.05 | 0.318 ± 0.029 | 1.54 (2σ: 1.13 to 2.42) | 1.27 ± 0.11 |
| 5/32 | 256 | 1.61 ± 0.05 | 0.327 ± 0.029 | 1.45 (2σ: 1.08 to 2.22) | 1.28 ± 0.09 |
| 1/8 | 256 | 1.65 ± 0.06 | 0.340 ± 0.030 | 1.33 (2σ: 1.01 to 1.96) | 1.25 ± 0.09 |

- Controls for the boson-frame fits: c2 reads 0.504 (boson, exact 0.5), 1.020 (inter-leg pair, exact 1) and 2.041 (intra-leg pair, exact 2), c1 0.918 (boson, exact 0.889), all at L = 256 with the image-corrected chord form. The uncertainties are the window spread in quadrature with the larger control bias.
- K_s^P = π²A/Z_P², where −A/r² is the smooth tail of the valley-polarization correlation C_P (window L/8 to L/4) and Z_P = ⟨P⟩/(p n) = 1.869 ± 0.003 is measured from the off-line runs (φ = 1.50 to 1.54, L = 128 to 256). The same procedure on ff2v reads 0.94 to 1.04 (exact 1). This estimator uses densities only. It assumes that the n(k) valley weight measures the low-energy N_s (Z_N = 1). If Z_N < 1, K_s^P is larger.
- Consistency checks. K_c^P = 1.6 agrees with the long-wavelength S(q) (1.5 to 1.65). The boson-frame prediction for the Jordan–Wigner fermion exponent, K_c + 1/(4K_c) + 1/(4K_s) = 1.93 with (1.6, 1.5), matches the measured 1.9 to 2.0 (REPORT 9.158), which the fermion frame cannot produce together with η_P. The slowly decaying 2k₀ wall modulation is the inter-valley coherence b_b†b_a ∝ e^{i√2θ_s}, bulk exponent 1/K_s ≈ 0.7, which is slow in the boson frame and impossible to make slow in the fermion frame.
- Drift. K_s^P is flat from L = 128 to 256 at every density (1.256, 1.271, 1.266 at 3/16). K_s^B scatters within its error (1.49, 1.43, 1.54). There is no sign of the logarithmic approach to 1 that a marginal flow would give, and nothing moves toward or below 1.
- Not consistent [numerical, open]. The parity gap: v_c/K_c = (2/π)∂²e/∂n² = 0.50 (L → ∞ about 0.49), and v_s/K_s = 2/(πχ_s) = 1.59 from the susceptibility χ_s = 2pn/δ = 0.40 ± 0.04 (φ = 1.50 to 1.555, L ≥ 128). These predict L·Δ_P = (π/2)(v_c/K_c + v_s/K_s) = 3.3 for a valley-preserving wall, against the measured 2.38. The valley-flip wall predicts more, not less. So the open-chain parity gap is about 30 percent below the two-component prediction, and I have no explanation for that. The onset of p with δ (p ≈ 0 at δ = 0.022, L = 64 and 128) agrees with the finite-size quantization δ_c = π(v_s/K_s)/L of a valley-preserving wall, which needs no pseudospin gap.

- Not done: bipartite fluctuations of the valley imbalance. The block valley number is not a local observable (the valley is a momentum label), and C176 stored only the one-body matrix of an interacting state, which does not fix two-body block fluctuations. The C_P tail measures the same quantity, ⟨ρ_s ρ_s⟩, and is used instead.

## Results, Part 2 (infinite chain)

Local runs (this machine, local/out/, same code and rule as the Typhon jobs) [numerical]:
- c1 control, χ = 16 to 96: c_fit = 1.015 (registered tolerance 0.15). The energy density approaches the exact value from above, error 7.6e-7 at χ = 96.
- ff2v control, χ = 16 to 96: c_fit = 1.989 (tolerance 0.25), local c between consecutive χ 1.92 to 2.03. The energy error is 3.5e-5 at χ = 96. ξ₁/ξ₂ = 2.13 and ξ₁/ξ₀ = 2.17, against the free-fermion scaling-dimension ratios 2 and 2.
- n = 1/8 candidate, χ = 16 to 96 (ξ up to 102 sites): c_fit = 2.008, local c 1.90 to 2.03 at every χ step, ξ₁, ξ₂ and ξ₀ all grow by about 1.9 over the last χ doubling, and ξ₁/ξ₂ changes by a factor 1.026. Registered reading: c = 2 liquid, no pseudospin gap up to ξ ≈ 100. [heuristic, unregistered] The transfer-matrix ratios agree with the boson frame (PROOF 6.5). With K_c = 1.6, ξ₂/ξ₁ = 0.486 = (1 + K_c/K_s)/4 gives K_s = 1.7, and ξ₁/ξ₀ = 2.17, read as Δ₀/Δ₁ with Δ₀ = 1/(2K_s) from the lowest neutral operator e^{i√2θ_s}, gives K_s = 1.35. Both are above 1, at ξ ≈ 100.
- n = 1/8 extended to χ = 256 (ξ up to 236 sites): c_fit = 2.011 over χ = 96 to 256, and the local c is 1.99 to 2.03 on every step from χ = 32 up (the last one, 192 → 256, is 2.00). Over the doubling 128 → 256, ξ₁, ξ₂ and ξ₀ grow by 1.813, 1.782 and 1.849, and ξ₁/ξ₂ changes by a factor 1.017. Registered reading: c = 2 liquid, no pseudospin gap up to ξ ≈ 236. The energy density is −0.2809474 at χ = 256. The ratios are ξ₂/ξ₁ = 0.479 and ξ₁/ξ₀ = 2.13, so the unregistered K_s readings barely move (1.75 and 1.41 with K_c = 1.6).
- n = 3/16, the main density, χ = 16 to 96 (ξ up to 89 sites): c_fit = 1.980 over χ = 32 to 96, local c 1.94 to 2.10 on every step. Over the doubling 48 → 96, ξ₁, ξ₂ and ξ₀ grow by 1.891, 1.855 and 1.954, and ξ₁/ξ₂ changes by a factor 1.019. Registered reading: c = 2 liquid, no pseudospin gap up to ξ ≈ 89. The energy density is −0.4182007 at χ = 96. [heuristic, unregistered] ξ₂/ξ₁ = 0.441 and ξ₁/ξ₀ = 1.97 read as K_s = 2.1 and 1.66 with K_c = 1.61. These two readings of the same number disagree by about 25 percent at both densities, so they only say that K_s is above 1, not what it is.
- The n = 3/16 local run was still going at write-up (χ = 128, then 192 and 256, a detached process writing to local/out/loc_s316.*). To read it, run `python3 summarize_inf.py local/out loc_ > local/summary_local.txt` from codex/work_C189/. If the process has died, `python3 run_inf.py local/cfg/loc_s316.json local/out` resumes it from its checkpoint. The summary as of write-up is in local/summary_local.txt. The loc_c2 block there is the frozen run that amendment A1 excludes, not a result.

The Typhon jobs have not started (see below). The registered c_fit for the χ range up to 1024 is therefore still open. The stated limit stands: this test resolves pseudospin correlation lengths up to the largest ξ reached, now 236 sites at n = 1/8 and 89 sites at n = 3/16.

## Verdict on relevance

On the evidence available, the pseudospin backscattering is irrelevant at n = 1/8, 5/32 and 3/16. The frame-correct estimators put K_s at 1.25 to 1.55 at all three densities. The density-only K_s^P exceeds 1 at 2σ at every density and size except 5/32 at L = 192, where its ff2v calibration is noisy (σ = 0.29). K_s^B exceeds 1 at 2σ at L = 256. Neither drifts from L = 128 to 256. The backscattering cos(√8φ_s) has scaling dimension 2K_s ≈ 2.5 to 3, so it flows to zero. The c = 2 liquid is then stable. The exponentially small pseudospin gap of the attractive fermion analogy does not apply, because the valleys are not two Fermi seas starting at K_s = 1. They are two boson quasi-condensates, which start above 1, and the inter-valley attraction lowers K_s but at these densities leaves it above 1. The infinite-chain runs agree as far as they reach (c_fit = 2.01 at n = 1/8 up to ξ ≈ 236 and 1.98 at n = 3/16 up to ξ ≈ 89, all sector lengths growing together).

Limits of this verdict:
- The boson-frame estimators were chosen after the data were seen (amendment A3). The registered fermion-frame estimators gave impossible values rather than a reading.
- The open-chain parity gap is 30 percent below the two-component prediction, unexplained.
- K_s^P assumes Z_N = 1.
- All of this is C173's effective model at second order, not the qubit chain.

## Typhon jobs (pending at write-up)

The first submission, a single 48-core job (17606815), had a scheduler start estimate of 2026-11-14. I cancelled it and resubmitted as five small jobs, which can start on partly used nodes. All five were PENDING (Priority) at write-up, and at the last check the scheduler's start estimate for them was also 2026-11-14 (14:52), so the split did not buy an earlier start. Each runs one config from typhon/cfg/ with run_inf.py and writes a summary when it ends (3-day limit).
- 17606915 c189_s316: n = 3/16, χ = 32 to 1024, 10 threads.
- 17606916 c189_s18: n = 1/8, χ = 32 to 1024, 8 threads.
- 17606917 c189_s16: n = 1/6, χ = 32 to 1024, 6 threads.
- 17606918 c189_ff2v: control, free fermions with the candidate's band, c = 2, χ = 32 to 1024, 10 threads.
- 17606919 c189_c1: control, hard rods, c = 1, χ = 32 to 768, 4 threads.

Scratch: /scratch/lustre/jkudlerflam/c18x/c189/ (code at the top level, configs in typhon/cfg/, outputs in out/, logs in logs/, checkpoints in out/ckpt/, which need not be copied). Each finished χ appends one line to out/<tag>.jsonl, and out/<tag>.iter.jsonl logs every iteration.

Collection, run from codex/work_C189/:

    mkdir -p typhon_out/out typhon_out/logs
    scp -o LogLevel=ERROR 'typhon:/scratch/lustre/jkudlerflam/c18x/c189/out/*.jsonl' typhon_out/out/
    scp -o LogLevel=ERROR 'typhon:/scratch/lustre/jkudlerflam/c18x/c189/logs/*' typhon_out/logs/
    scp -o LogLevel=ERROR 'typhon:/scratch/lustre/jkudlerflam/c18x/c189/summary_inf_*.txt' typhon_out/
    python3 summarize_inf.py typhon_out/out > typhon_out/summary_local.txt

What to read on the outputs:
- Calibration under the registered rule: inf_c1 must give c_fit within 0.15 of 1, inf_ff2v within 0.25 of 2.
- The registered c = 2 and gap readings for inf_s316, inf_s18 and inf_s16 (printed by summarize_inf.py).
- Secondary, unregistered: the boson-frame ratios ξ₂/ξ₁ = (1 + K_c/K_s)/4 and ξ₁/ξ₀ = Δ₀/Δ₁ at the largest χ (PROOF 6.5), against K_c ≈ 1.6 and K_s ≈ 1.3 to 1.5 from Part 1.
- At n = 1/6, the energy density against C176's ring ED (−0.3722 to −0.3749 at L = 24 to 36) and open-chain −0.37327.

## Desktop package

codex/work_C189/desktop_run/ (README.txt, PROMPT.txt, run_desktop.sh, tasks.txt, cfg/, code). It runs n = 3/16 to χ = 768 and n = 5/32 (the density whose bulk value is closest to C175's 0.17) to χ = 512, two processes of 4 threads. The same registered rule applies, and a smoke test passed. It covers the main density in case the Typhon jobs keep waiting, plus a density Typhon does not run.

## Files

- Registration and results: RESULT.md. Derivations: PROOF.md. Log: PROGRESS.md.
- Part 1: kest.py (estimator library and the exact ff2v generator), kest_controls.py/.out, kest_run.py/.out/.json (registered, fermion frame), kest_boson.py/.out (boson frame), kest_cp.py/.out (pseudospin density).
- Part 2: idmrg.py (infinite DMRG), run_inf.py (driver), summarize_inf.py (registered readings, optional second argument for the tag prefix), validate_inf.py/.out, local/ (local runs: cfg/, out/ with jsonl and logs, summary_local.txt), typhon/ (configs, task list, slurm scripts).
- ivumps.py: the abandoned single-site VUMPS, kept as the record of why it fails with mixed labels (PROOF 6.2).

## Coordinator note (2026-11-09, after the agent finished)

- The local n = 3/16 run has since reached χ = 256. The refreshed `local/summary_local.txt` reads c_fit = 2.011 over χ = 96 to 256, ξ up to 205 sites, growth of ξ₁, ξ₂, ξ₀ by 1.811, 1.793, 1.828 over χ = 128 → 256, and ξ₁/ξ₂ changing by a factor 1.010. This replaces the χ ≤ 96 numbers quoted above for that density.
- `kest_cp.out` was regenerated with the final `kest_cp.py`. The saved copy came from an earlier version of the script and lacked the K_s^P table. The rerun reproduces the K_s^P values quoted above.
- `run_inf.py` now writes its checkpoints atomically (temporary file, then rename). The copy in `desktop_run/` is unchanged.
- The Typhon runs were started as extra steps inside other running allocations of ours. See `../work_C192/TYPHON_OVERLAP.md` and `typhon/overlap_inf.sh`.
- Write-up: REPORT §9.163.

## Typhon runs to χ = 768 (coordinator, Claude Fable 5, collected 2026-11-10 06:40 EST)

The five infinite-chain runs went as extra steps inside other allocations (`../work_C192/TYPHON_OVERLAP.md`), not as the queued jobs listed above. inf_c1 has finished its whole sequence (χ = 32 to 768). The other four have completed χ = 768 and are running χ = 1024, so their numbers here are interim and c_fit will shift once the last χ lands. Collected files are typhon_out/out/<tag>.jsonl (finished χ only), typhon_out/logs/ and typhon_out/summary_local.txt, made with `python3 summarize_inf.py typhon_out/out > typhon_out/summary_local.txt`. Everything below is from that summary [numerical]. Only what goes beyond the local runs (χ ≤ 256) is recorded.

| tag | χ | S̄ | ξ₁, ξ₂, ξ₀ (sites) | local c from the previous χ | iterations |
|---|---|---|---|---|---|
| inf_s316 (n = 3/16) | 384 | 2.3527 | 287.8, 125.3, 149.7 | 1.991 | 40 |
| | 512 | 2.4315 | 364.3, 157.0, 193.2 | 2.005 | 45 |
| | 768 | 2.5435 | 510.6, 220.8, 269.1 | 1.991 | 62 |
| inf_s18 (n = 1/8) | 384 | 2.2514 | 333.6, 161.1, 154.5 | 2.006 | 85 |
| | 512 | 2.3294 | 422.3, 202.3, 196.7 | 1.983 | 30 |
| | 768 | 2.4413 | 590.9, 282.9, 276.4 | 2.000 | 138 |
| inf_s16 (n = 1/6) | 384 | 2.3230 | 299.9, 134.4, 150.4 | 1.990 | 91 |
| | 512 | 2.4013 | 379.4, 168.8, 192.8 | 1.996 | 105 |
| | 768 | 2.5136 | 532.1, 237.1, 269.8 | 1.992 | 126 |
| inf_ff2v (c = 2 control) | 384 | 2.2590 | 301.8, 141.0, 139.2 | 1.975 | 58 |
| | 512 | 2.3393 | 384.5, 177.7, 176.0 | 1.989 | 69 |
| | 768 | 2.4500 | 538.9, 250.5, 248.6 | 1.967 | 92 |
| inf_c1 (c = 1 control) | 384 | 1.5220 | 1434.5, 225.6, 1242.1 | 0.808 | 300 (cap) |
| | 512 | 1.5821 | 2051.3, 320.9, 1776.6 | 1.008 | 300 (cap) |
| | 768 | 1.6200 | 2881.3, 469.1, 2504.3 | 0.669 | 300 (cap) |

Calibration under the registered rule (c_fit over the largest four χ).
- inf_ff2v gives c_fit = 1.977, inside the 0.25 tolerance. Its energy error against the exact value is 1.2e-6 at χ = 768.
- inf_c1 gives c_fit = 0.847, which is 0.153 below 1 and so misses the 0.15 tolerance by 0.003. By the registration the readings below are therefore uncalibrated on the c = 1 side, and I record them with that label.
- What went wrong with the control is visible. It hit the iteration cap (300) at every χ ≥ 64, while the candidates stop on tolerance after 30 to 138 iterations. Its correlation length at χ = 768 is 2881 sites, five times the candidates'. Its local c is 0.93 to 1.07 on every step up to χ = 256 (ξ ≤ 940) and scatters only beyond (0.81, 1.01, 0.67 at ξ = 1400 to 2900). Its energy error is 5.9e-9 at χ = 768. [post hoc, unregistered] Fitted over the ξ range that the candidates cover, the control reads 0.957 (χ = 64 to 192, ξ = 167 to 663) and 0.978 (χ = 32 to 256). The earlier local control (χ ≤ 96) read 1.015.

Readings for the candidates (registered rule, interim at χ = 768):

| tag | c_fit (χ = 256 to 768) | growth of ξ₁, ξ₂, ξ₀ over 384 → 768 | change of ξ₁/ξ₂ | c = 2 liquid | pseudospin gap | e at χ = 768 | discarded weight |
|---|---|---|---|---|---|---|---|
| inf_s316 | 1.995 | 1.774, 1.762, 1.797 | ×1.007 | yes | no | −0.4182637706 | 2.9e-7 |
| inf_s18 | 1.997 | 1.771, 1.756, 1.789 | ×1.009 | yes | no | −0.2809524996 | 1.4e-7 |
| inf_s16 | 1.993 | 1.775, 1.764, 1.794 | ×1.006 | yes | no | −0.3728294392 | 2.5e-7 |

- Beyond χ = 256. The local runs reached ξ = 205 sites at n = 3/16 and 236 at n = 1/8. The Typhon runs reach 511, 591 and 532 sites at n = 3/16, 1/8 and 1/6, and the local slope stays between 1.98 and 2.01 on every step from χ = 256 to 768 at all three densities. All three correlation lengths keep growing at the same rate, a factor 1.76 to 1.80 over the last doubling of χ, and no sector saturates.
- n = 1/6 is new. Its energy density −0.37283 sits where C176's open-chain value (−0.37327) and ring ED (−0.3722 to −0.3749) put it.
- [heuristic, unregistered] The ratios at χ = 768 are ξ₂/ξ₁ = 0.432, 0.479, 0.446 and ξ₁/ξ₀ = 1.90, 2.14, 1.97 at n = 3/16, 1/8, 1/6, against 0.465 and 2.17 for the free control. They have not moved from the local values (0.441 and 1.97 at n = 3/16, 0.479 and 2.13 at n = 1/8).
- Stated limit, unchanged in kind. The test now resolves a pseudospin correlation length up to about 500 to 600 sites. C175's estimate (ξ_s of order 1e4) is still far beyond it.
- Left to do. Collect χ = 1024 for the four running tags and rerun the summary. A rerun of the c = 1 control with a higher iteration cap is needed before the registered calibration can be called met or failed on converged states.
