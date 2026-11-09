This run was done by a Claude Opus 5.5 agent (C183, 2026-11-09), not by Codex. Astra was not used (out of credits).

# C183 result: the full qubit chain near the Fibonacci cusp, a targeted search

Derivations are in PROOF.md and the log in PROGRESS.md.

## Verdict (interim, written 2026-11-09 with both Typhon jobs still running or pending)

No c > 1 candidate is established by this run, and none is excluded. The local work gives four results.
1. [numerical] The far-field two-body map of the full chain (stage T) has an unbound corner r < 10/|w| with an easy-plane dilute pseudospin, but the corner is an artefact of measuring the interactions far below the band. Refitting with the field raised toward the band moves the pair resonance to smaller V/Δ (r ≈ 10 → 8 → 6 → below 4 as μ_eff goes −12 → −8 → −5 → −4, |w| = 1), so close to the gas onset (μ_eff = −4 to −5, the highest field where the two-body extraction is clean) the parity-class pair is bound for every r ≥ 4 tested, and the trend keeps it bound at the onset itself. This confirms the coordinator's correction from C186 with an independent full-chain fit.
2. [derived, checked against the data] The bound pair and the easy-plane two-component gas are the two sides of one zero-energy resonance (Δ_σ → +1 from the unbound side, → −1 from the bound side, PROOF §2.2). The cusp gas in the parity class sits on the bound side at the physical field.
3. [numerical] In the full chain the gas exists only in a window above a sharp onset (below it, a gapped dressed vacuum with up to 0.1 virtual particles per site) and below a first-order jump to the blockade-packed phase (n ≈ 0.45): at r = 8, |w| = 1 the window is x ≈ 0.075 to 0.3, densities 0.05 to 0.3.
4. [numerical, derived] For C186's transverse-DM bond the uniform ΔN = 1 element vanishes exactly (antiunitary I·P·K), and the full-chain retuning of the valley degeneracy at V/Δ = 16 shifts J_x + J_y from −0.00525 to −0.00400.

The many-body tests are on Typhon: job 17604545 (the V/Δ = 8 points S1 and S2, L = 64 to 256 with parity gaps and controls, running) and job 17604499 (the 180-point screen including the transverse-DM bond, family G, and C181's off-symmetric point in the full chain, family H, with automatic confirmation to L = 256 and VUMPS FES, pending). Their pass criteria are registered below.

## Pre-registration (written before any C183 output was computed)

Seen beforehand: C177's local logs (F8 at L = 32, 48, 64 reads c_prof 2.0 to 2.2 and S̄ 1.17, 1.45, 1.47; its one-valley control C8 reads S̄ 1.01, 0.96; its generic G8 variants with λ = 0.05 and 0.1 read S̄ 0.59 and 0.51) and C177 PROOF §3 (the bound pair of C175 is absent in the full chain at V/Δ = 8 and present at V/Δ ≥ 12, |w| = 1). Nothing else of the full chain off C177's points.

Convention. h = J_x XX + J_y YY + D(XY − YX) + J_z ZZ + (h_z/2)(Z⊗1 + 1⊗Z) + J_xz(XZ + ZX) + K_xz(XZ − ZX) + (h_x/2)(X⊗1 + 1⊗X), Δ = J_x − J_y = 1, V = 4J_z, r = V/Δ. A particle is a down spin. The dressed one-particle band ε(k) of the full chain is measured by ring ED (C177 band.py), fitted to harmonics to order 4, and (J_x + J_y, D) are tuned so that its first harmonic is 2|w||t₂| cos(k − φ) with the requested |w|/|t₂| and φ (|t₂| = −b₂/2 from the measured cos 2k harmonic). h_z sets the density. The dressed density is n_d = n_mid − n_vac with n_vac the vacuum density of virtual pairs at the same r (ring ED far below the band).

### Stage T, the two-body map of the full chain (local)

At φ = π/2, r ∈ {4, 5, 6, 8, 10, 12, 16}, |w| ∈ {0.5, 1, 1.5, 2, 3}: full-chain two-magnon levels at K = π on rings L = 14 to 20, fitted by an effective two-body model with the exact band (C177 twob method, models A and B), then the pair binding at L = 400 and C175's dilute pseudospin anisotropy Δ_σ = (a_s − A_d)/|A_o| from C175's scat.py with the fitted V₂, V₃.
- Rule. A point is two-body favourable if the pair is unbound (E₂ − 2ε_min > −1e-3|t₂| at L = 400 in both models) and Δ_σ ∈ (−1, 1].
- Prediction T1. The (Δ/V)² corrections are repulsive in the unlike channel at every |w|, so each |w| has an unbinding ratio r_b(|w|) below which the pair is unbound, and r_b decreases with |w|: r_b(0.5) > r_b(1) (between 8 and 12, C177) > r_b(2). At |w| = 3 the pair stays bound for r ≥ 6.
- No prediction is registered for the sign of Δ_σ where the pair is unbound.

### Stage S, the many-body screen (Typhon, open chains of the full qubit chain)

Families (all with C166/C177's validated DMRG, fc183.py):
- A, on the J_x + J_y = 0 stratum: r ∈ {4, 6, 8, 10, 12, 16} × |w| ∈ {0.5, 1, 2, 3}, at three dressed densities n_d ≈ 0.08, 0.17, 0.28 (n_d ≈ 0.17 at r = 8 and 12, |w| = 1 is C177's and is skipped).
- B, off the stratum: φ ∈ {1.56, 1.53, 1.48, 1.40} at r ∈ {6, 10}, |w| ∈ {1, 2}, n_d ≈ 0.08 and 0.17.
- C, off the parity class: at r ∈ {6, 10}, |w| ∈ {1, 2}, φ = π/2, n_d ≈ 0.17, one of J_xz, K_xz, h_x equal to λ with λ/|t₂| ∈ {0.1, 0.3}.
- Density targeting: per family an L = 32, χ = 64 ladder in the dressed chemical potential, interpolated to the target n_d. A target inside a density jump is replaced by the first point above the jump and flagged.

Screen quantities at L = 32 and 64: S̄ (mean of the three middle cuts), c₂ = 6[S̄(64) − S̄(32)]/ln 2, c_prof(64) (C166 dmrg_c.profile_fits), the final truncation, n_d, the valley weights of the one-body row G(r).

Screen pass rule (promotion). c_prof(64) ≥ 1.4, c₂ ≥ 1.3 and truncation at L = 64 below 1e-6. The four parity-class points and the one generic (family C) point with the largest min(c₂, c_prof(64)) go to stage C whether or not they pass, and are labelled "passed" or "best failing".

Controls at the same sizes: the XX chain at n = 0.17 (c = 1) and two decoupled XX legs at n = 0.17 per leg as a d = 4 chain of L/2 rungs (c = 2), plus the full chain at φ = 1.30 (one valley, c = 1) at r = 6 and 10, |w| = 1 and 2, n_d ≈ 0.17. The screen is calibrated if the XX control reads c₂ within 0.25 of 1 and the two-leg control within 0.4 of 2.

Predictions for stage S.
- Q1. Every family B point with φ ≤ 1.53 fails the screen (one valley wins at |φ − π/2| ≥ 0.04 for any density reached).
- Q2. Family C at λ = 0.3|t₂| fails the screen at every point. At λ = 0.1|t₂| the reading is within 0.3 of the parity-class parent in c_prof(64).
- Q3. On the stratum at |w| ≥ 2 and r ≥ 10, where the pair and larger clusters bind, the points fail the screen or show a density jump (clustering). No prediction for r ≤ 6.

### Stage C, confirmation of the promoted points (same Typhon job, automatic)

For each promoted point: open chains at L = 128 and 256 (χ up to 256 and 320), the parity-sector ground states at L = 64 and 128 (parity-class points only), and damped VUMPS in the two-site cell to χ = 256 with the FES slope. A point is a c > 1 candidate only if all of these hold (C163 §4.2 as adapted by C177):
- P1. c₂(128 → 256) > 1.2, c_prof(128) and c_prof(256) above 1.2 with c_prof(256) ≥ c_prof(128) − 0.03, truncation below 1e-7 at both lengths.
- P2. FES slope (6 dS/d ln ξ over the last three χ, χ_max = 256) above 1.2.
- P3. |e_mid(256) − e₀| ≤ 2|e_mid(256) − e_mid(128)| + 1e-6, with e₀ the VUMPS energy per site and e_mid the mean bond energy over the central L/4 bonds.
- P4 (parity class). The parity gap Δ_p(L) = E_odd(L) − E_even(L) has L·|Δ_p| growing by less than 50 percent from L = 64 to 128.
- Controls at L = 128 → 256 must read within 0.2 of 1 (XX) and 0.3 of 2 (two legs), else the stage is reported as uncalibrated.

### Addendum (registered after seeing the stage T rows at r ≤ 6, before any stage S output existed)

Stage T at r = 4, 5, 6 shows unbound unlike pairs with Δ_σ inside (−1, 1] at |w| ≤ 1.5 (r ≤ 6) and |w| = 2 (r ≤ 5), and bound pairs at |w| = 3 (all r ≥ 4) and |w| = 2 at r = 6. So family D is added to the screen, aimed at the two-body favourable corner:
- D on the stratum: r = 5 at |w| ∈ {0.5, 1, 1.5} with n_d ≈ 0.06, 0.12, 0.17, 0.24, and r = 4 at |w| ∈ {0.5, 1} with n_d ≈ 0.06, 0.12.
- D off the stratum: r = 5, |w| = 1, φ ∈ {1.565, 1.55, 1.53}, n_d ≈ 0.06, 0.12.
- Selection: the D point with the largest min(c₂, c_prof(64)) also goes to stage C (so stage C has up to six points). Same pass rule.
- Prediction D1. If the dilute two-body criterion carries over to the many-body chain, the D points on the stratum at n_d ≤ 0.12 read c₂ and c_prof(64) above the C177-type c = 1 control and do not show a density jump. At L ≤ 64 and n_d ≈ 0.06 (two to four particles) the screen is not expected to resolve c, and a failure there is not counted against the two-body picture.

### Addendum E (registered after the coordinator's priority note and the full stage T map, before any stage S output existed; job 17604432 was cancelled while still pending and replaced by the job below)

Priority on the unbound corner at V/Δ = 6 to 12, on and slightly off J_x + J_y = 0 (off the line the two-pair vertex of C179 carries momentum).
- Family E: (r, |w|) ∈ {(6, 1), (8, 1), (8, 0.5), (12, 0.5), (6, 1.5)}, all on the unbound side in stage T, × φ ∈ {π/2, 1.565, 1.555, 1.545} (valley offsets about 0, 0.02, 0.06 and 0.10|t₂|), × n_d ≈ 0.08, 0.12, 0.17 (n_d ≈ 0.17 at r = 8, |w| = 1, φ = π/2 is C177's and is skipped).
- Forced confirmation: in each of (r, |w|) = (8, 1) and (6, 1), the on-stratum E point and the off-stratum E point with the largest min(c₂, c_prof(64)) go to stage C whatever their score, in addition to the six selected by the earlier rules. Same pass criteria P1 to P4.
- Prediction E1. If the two-body picture holds in the many-body chain, the unbound E points at n_d ≤ 0.12 on the stratum and at φ = 1.565 read c₂ and c_prof(64) above 1.3, and their parity gap L·Δ_p does not grow from 64 to 128. At φ = 1.545 (offset 0.10|t₂|) the dilute pseudospin field is comparable to the estimated saturation field at n_d ≈ 0.08, so I expect one valley there (c₂ and c_prof below 1.2) at n_d ≈ 0.08 and no prediction at higher density.
- The parity gap of the full chain is E_odd − E_even at fixed field. N is not conserved, only its parity, so this is the grand-canonical counterpart of E(N+1) + E(N−1) − 2E(N), equal to half of it when the field sits midway between the N ± 1 thresholds.

### Amendment to the density targeting (made while job 17604499 was still pending, before any stage S output)

The local ladders (local/p5, local/p8) showed that near the gas onset the dressed vacuum carries more virtual pairs than the far-below vacuum, so n_raw − n_vac reads 0.03 to 0.07 while the chain is still the gapped vacuum (at r = 5, |w| = 1, x = 0.45 the L = 64 chain has a flat S̄ = 0.44, no Friedel oscillation and |G(r)| decaying exponentially, local/Pr5_64.json). The targeting rule is therefore changed to count the density from the gas onset: the onset is the first ladder point whose S̄(32) exceeds that of the lowest x by 0.3 (or the lowest x if its S̄ is already above 0.6), n_d = n_raw − n_raw(point before the onset). A target whose bracketing upper ladder point is in the dense blockade-packed phase (S̄ lower by more than 0.3) is placed at the last gas point and flagged. The ladder grid is refined to x ∈ {−0.4, −0.2, −0.1, 0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.6, 2.4}. The forced (r, |w|) = (8, 1) confirmations of addendum E are dropped from the big job because stage F covers them. The screen pass rule and all criteria are unchanged.

### Correction and family G (registered after the coordinator's correction, while job 17604499 was still pending, before any of its output)

The coordinator (from C186, REPORT §9.146) corrected the premise of stage T: the two-body fits were taken at the far-below field (μ_eff = −12), and the effective V₂ depends on the field through the virtual-pair denominators, so at the gas onset (μ_eff ≈ −2.2) the parity-class pair is bound at V/Δ = 8 (C186: v ≈ 3(Δ/V)², resonance function −0.49). If this holds, the stage T map describes chains whose band lies about 10|t₂| above zero density, not the gas, and the "unbound corner" of stage T does not transfer to the gas. I check the field dependence myself (twob183.py at several μ_eff, below) before using it.
- Family E of addendum E is dropped from the big job (its premise is the far-below fit). Families A, B, C, D, K stay as registered (they are data on the full chain whatever the premise).
- Family G, the transverse-DM bond of C185 and C186, Γ(XZ − ZX) with Γ = Δ√v/4, which unbinds the pair at second order: r ∈ {12, 16}, v ∈ {0.3, 0.6} (Γ = 0.1369, 0.1936), |w| = 1.035, J_x + J_y retuned in the full chain so the two valleys are degenerate (offset 0) or split by 0.05|t₂| (v = 0.3 only), D tuned to the measured |w| (tune183.tune_dm, band read without a parity split since Γ breaks the Z₂ parity), densities 0.06, 0.12, 0.17 counted from the onset. At r = 16, v = 0.3 the full-chain retuning gives J_x + J_y = −0.00400 and D = 0.015212 against C186's second-order −0.00525 and 0.016174.
- Forced confirmation: the best G point at r = 16 and the best at r = 12 go to stage C (open chains 128 and 256, VUMPS FES). No parity sectors (the bond breaks parity). Criteria P1 to P3 as registered.
- Prediction G1. If C186's second-order unbinding survives in the full chain, the G points at n_d ≤ 0.12 with offset 0 read c₂ and c_prof(64) above 1.3 and above their offset-0.05 partners. If the ΔN = 1 vertex r Σ b†(n_{j+2} − n_{j−2}) that comes with v > 0 (C186) gaps or polarizes the gas, they read below 1.2. No prediction for which.

### Family H and two notes on family G (registered after two further coordinator notes, while job 17604499 was still pending at 12:45 EST, before any of its output)

- Family G keeps Γ_xz = −Γ_zx exactly with no h_x and no symmetric XZ + ZX. The uniform ΔN = 1 element is checked in tune_dm by ⟨X⟩ in the dressed vacuum (ring L = 10, full ED), and cancelled with h_x if nonzero. At r = 16, v = 0.3 it is 2e-16 (dm_uniform_check.out), against −0.0115 for a symmetric XZ + ZX of only 0.001 with or without the DM term, so the test is sensitive and no h_x is needed. [derived] It vanishes exactly on the ring because Θ′ = I·P·K (inversion, P = ΠZ, complex conjugation) is an antiunitary symmetry of the DM bond (every term is even: D picks up − from I and − from K, Γ(XZ − ZX) picks up − from I and − from P) and maps ⟨X_j⟩ to −⟨X_{−j}⟩. The G targets (0.06, 0.12, 0.17) avoid 1/5.
- Family H, C181's partially polarized point in the full chain: parity class, |w| = 1, φ ∈ {1.555, 1.545, 1.53}, r ∈ {12, 16, 24}, n = 3/16 counted from the onset (the physical μ), plus a one-valley control at φ = 1.30, r = 16 on the same chain. Forced confirmation (open chains 128 and 256, parity sectors 64 and 128, VUMPS): φ = 1.555 and 1.53 at r = 16 and 24. Reported in addition: the valley imbalance p from the ⟨XX + YY⟩ profile, and c₂(64 → 128) from the screen L = 64 and the confirmation L = 128.
- To make room, the A/B rule now sends two points to stage C instead of four (amended before any output).
- Prediction H1 (C181's effective-model reading carried to the full chain). At r = 24 and 16, φ = 1.555 and 1.53: c₂(64 → 128) above 1.6 and the one-valley control near 1. Prediction H2 (C179's two-pair vertex, detuned by q₄ ∝ |π/2 − φ|, with an amplitude growing as V/Δ falls): the reading drops at r = 12 relative to r = 24 at fixed φ, more at φ = 1.555 than at 1.53. No prediction for which holds.

### Stage F, the small V/Δ = 8 confirmation job (coordinator exception; registered before submission)

Points (couplings in typhon/f_points.json), r = 8, |w| = 1, stage T unbound with Δ_σ = 0.90:
- S1 on the stratum: J_x = 0.5, J_y = −0.5, D = 0.0290182, J_z = 2, h_z = −3.8811529 (gas density about 0.16 at L = 24, local/p8);
- S2 at φ = 1.555 (valley offset 0.062|t₂|): J_x = 0.5004616, J_y = −0.4995384, D = 0.0290146, J_z = 2, h_z = −3.8805730 (density about 0.17 at L = 24, local/p8b).

Runs: open chains at L = 64, 128, 256 (χ to 192, 256, 256), parity sectors at L = 64 and 128, controls (XX, two legs) at L = 64, 128, 256. Criteria: P1 (with the truncation bound 1e-7, and c₂ evaluated for 64 → 128 and 128 → 256), P4, the controls as in stage C, and n(k) from the one-body row.
- Prediction F1 (two-body picture). c₂(64 → 128) and c₂(128 → 256) both at or above 1.6 at S1 and S2, and L·Δ_p flat within 20 percent from 64 to 128.
- Prediction F2 (pairing or clustering re-emerging from three- and more-body effects, which stage T does not see). c₂ falls with L, below 1.4 at 128 → 256, and L·Δ_p grows.
- No prediction is registered for which of F1 and F2 occurs. Because Δ_σ = 0.90 is close to the isotropic end, the pseudospin backscattering is close to marginal, so slow logarithmic drifts in c₂ are expected in either case.

### Local pilot (registered after the full stage T map and the L = 24 density ladders, before any pilot open-chain output)

Parity gap Δ_p(L) = E_odd(L) − E_even(L) of the open chain at fixed h_z (fc.py --parity), L = 32 and 64, χ ≤ 160, on the stratum at n_d ≈ 0.12. Points: r = 5 with |w| = 1 and 2 (unbound, Δ_σ = 0.66 and 0.96), r = 12 with |w| = 2 (bound, two-body E_b = 0.076|t₂|), and r = 12 with |w| = 1 (bound on the resonance edge, E_b ≈ 0.001|t₂|, not expected to be resolved).
- Pairing prediction for r = 12, |w| = 2: L·Δ_p grows by more than 50 percent from L = 32 to 64, with Δ_p(64) of order E_b/2 or more.
- Two-component prediction for r = 5: L·Δ_p changes by less than 20 percent or falls, and the two-size c₂(32 → 64) is above the r = 12, |w| = 2 value.
- These sizes hold 4 to 8 particles, so this is a pilot for the Typhon stages, not a pass criterion.

## Typhon job

Job 17604499 (submitted 2026-11-09 after cancelling the still-pending 17604432, partition typhon, one node, 48 cores, time limit lowered from 72 h to 48 h while pending, after which the estimated start moved to 2026-11-09 13:45; pipeline amended twice while pending, as registered above), directory /scratch/lustre/jkudlerflam/c183/, script typhon/c183.slurm running pipeline183.py (phases 0 to 3 as 24 workers × 2 threads, phase 4 as 16 workers × 3 threads). Expected outputs:
- tunes.json, plan.json (density ladders and the chosen x per point), lad/*_32.json (ladders);
- scr/<tag>_32.json and scr/<tag>_64.json for the 180 screen points (families A, B, C, D, G, H, K; E was dropped), ctrl/xx_*.json and ctrl/leg2_*.json (32, 64, 128, 256);
- screen_summary.json and screen_summary.txt (stage S table, controls, the selected tags);
- conf/<tag>_128.json, conf/<tag>_256.json, conf/<tag>_{odd,even}_{64,128}.json, vc/<tag>_n2.json (stage C, up to ten points: two by the A/B rule, one generic C, one D, the best G at r = 16 and at r = 12, and four forced H points; parity sectors only for parity-class points);
- logs/ (per-task logs and c183_17604499.out).

To collect, from inside codex/work_C183:

    mkdir -p typhon_out && scp -r typhon:/scratch/lustre/jkudlerflam/c183/scr typhon:/scratch/lustre/jkudlerflam/c183/conf typhon:/scratch/lustre/jkudlerflam/c183/ctrl typhon:/scratch/lustre/jkudlerflam/c183/lad typhon:/scratch/lustre/jkudlerflam/c183/logs typhon:/scratch/lustre/jkudlerflam/c183/tunes.json typhon:/scratch/lustre/jkudlerflam/c183/plan.json typhon:/scratch/lustre/jkudlerflam/c183/screen_summary.json typhon:/scratch/lustre/jkudlerflam/c183/screen_summary.txt typhon_out/
    mkdir -p typhon_out/vc && scp typhon:/scratch/lustre/jkudlerflam/c183/vc/*.json typhon_out/vc/

then `python3 summ183.py typhon_out > typhon_out/summary_local.txt` (evaluates P1 to P4 and the controls). The vc/*.npz checkpoints need not be copied.

### Small job (stage F)

Job 17604545 (submitted 2026-11-09 as the coordinator's exception, partition typhon, one node, 16 cores, 12 h limit, estimated start the same day), directory /scratch/lustre/jkudlerflam/c183/f/, script typhon/c183f.slurm with task lists typhon/f_laneA.txt (the two L = 256 chains, 4 threads each) and typhon/f_laneB.txt (L = 128 and 64, parity sectors, controls, 4 slots × 2 threads). Expected outputs: oc/S1_{64,128,256}.json, oc/S2_{64,128,256}.json, oc/S1odd_{64,128}.json, oc/S1even_{64,128}.json, oc/S2odd_{64,128}.json, oc/S2even_{64,128}.json, ctrl/xx_{64,128,256}.json, ctrl/leg2_{64,128,256}.json, logs/.

To collect, from inside codex/work_C183:

    mkdir -p typhon_out_f && scp -r typhon:/scratch/lustre/jkudlerflam/c183/f/oc typhon:/scratch/lustre/jkudlerflam/c183/f/ctrl typhon:/scratch/lustre/jkudlerflam/c183/f/logs typhon_out_f/

then `python3 summ183f.py typhon_out_f` (c₂ for 64 → 128 and 128 → 256, c_prof, truncation, L·Δ_p, n(k) valley weights, the controls, and F1/F2).

## Results

### Stage T, the two-body map (twob183.py, twob_map.out, twob_map.json; PROOF §2)

Each cell: U (unbound side, A_d < 0) or B (bound side), E₂ − 2ε_min at L = 400 in |t₂| (model A; model B agrees to 1e-3 or better), and C175's Δ_σ with the fitted V₂, V₃.

| r = V/Δ | |w| = 0.5 | |w| = 1 | |w| = 1.5 | |w| = 2 | |w| = 3 |
| :-- | :-- | :-- | :-- | :-- | :-- |
| 4 | U, +0.0005, Δσ +0.36 | U, +0.0005, Δσ +0.57 | U, +0.0004, Δσ +0.74 | U, +0.0004, Δσ +0.84 | B, −0.0330, Δσ −1.06 |
| 5 | U, +0.0005, Δσ +0.43 | U, +0.0005, Δσ +0.66 | U, +0.0004, Δσ +0.84 | U, +0.0001, Δσ +0.96 | B, −0.0630, Δσ −1.23 |
| 6 | U, +0.0005, Δσ +0.51 | U, +0.0004, Δσ +0.75 | U, +0.0003, Δσ +0.93 | B, −0.0083, Δσ −1.08 | B, −0.0881, Δσ −1.41 |
| 8 | U, +0.0005, Δσ +0.64 | U, +0.0005, Δσ +0.90 | B, −0.0036, Δσ −1.09 | B, −0.0371, Δσ −1.26 | B, −0.1232, Δσ −1.74 |
| 10 | U, +0.0005, Δσ +0.75 | B (on resonance, A_d = 727), +0.0004, Δσ −1.00 | B, −0.0151, Δσ −1.20 | B, −0.0600, Δσ −1.39 | B, −0.1446, Δσ −1.99 |
| 12 | U, +0.0004, Δσ +0.83 | B, −0.0011, Δσ −1.08 | B, −0.0255, Δσ −1.28 | B, −0.0758, Δσ −1.49 | B, −0.1580, Δσ −2.19 |
| 16 | U, +0.0003, Δσ +0.94 | B, −0.0060, Δσ −1.17 | B, −0.0395, Δσ −1.38 | B, −0.0942, Δσ −1.60 | B, −0.1730, Δσ −2.43 |

Positive E₂ − 2ε_min of a few 1e-4 is the lowest scattering level at L = 400. The lowest two-body level over all K (bare band with the fitted V's) sits at K = π everywhere, and no like-valley bound state appears at any unbound point.

- [numerical, at μ_eff = −12 only, see the correction below] Prediction T1 holds. The pair unbinds below r_b(|w|) ≈ 10/|w| (r_b ≈ 16 to 20, 10, 7, 5.5 and below 4 at |w| = 0.5, 1, 1.5, 2, 3). In bond language the unbound side is D V²/Δ³ ≲ 2.5 (PROOF §2.4, the scaling form heuristic).
- [derived, PROOF §2.2] The bound pair and the easy-plane two-component gas are the two sides of one zero-energy resonance. Approaching it from the unbound side Δ_σ → +1 from below, from the bound side Δ_σ → −1 from below. The data follow this at every point.
- [numerical, at μ_eff = −12 only, with the dilute criterion derived by C175 and C180] Every unbound point has Δ_σ in (0.36, 0.96), inside C175's two-component window, and no bound state at any K. So at the two-body level the dilute full chain is a two-component (c = 2) gas on the whole unbound corner. Inputs not checked (PROOF §2.3): three-body binding, the extrapolation of the fitted interactions, the effective description at r = 4 and 5, and the dilute condition n|A| ≪ 1 (|A| is tens of sites near the resonance, so the controlled window is narrow there, and widest at small |w| and small r, e.g. |A| ≈ 2 to 4 at r = 4 and 5, |w| ≤ 1).
- [derived from C179 and C180, not computed here] On the stratum the momentum-2π two-pair vertex has dimension 2/K₊ in C180's convention, with K₊ → 1/2 in the dilute impenetrable gas, so it is irrelevant in the dilute limit and can only matter at the attractive, denser side (C179's K_ρ ≈ 3.2 to 3.6 was measured at r = ∞ with the pair bound).

### Stage T at the physical field: the correction holds (twob_mu.out, twob_mu_small.out)

The same fit at |w| = 1 with the field raised toward the band (μ_eff from −12 to −4; the one-particle band and the two-magnon levels are both measured at that field):

| r | μ_eff = −12 | −8 | −5 | −4 (fit rms 2 to 3e-2, less reliable) |
| :-- | :-- | :-- | :-- | :-- |
| 4 | V₂ −1.30, U, Δσ +0.57 | −1.50, U, +0.67 | −1.73, U, +0.81 | −2.21, B, −1.44 |
| 5 | −1.47, U, +0.66 | −1.64, U, +0.78 | −1.81, U, +0.92 | −2.35, B, −1.82 |
| 6 | −1.59, U, +0.75 | −1.72, U, +0.87 | −1.87, on resonance (A_d = 384), −1.006 | −2.40, B, −2.12 |
| 8 | −1.73, U, +0.90 | −1.83, on resonance (A_d = 641), −1.004 | −1.92, B (E_b = 6e-4), −1.11 | (−3.5: fit fails, rms 6e-2) |

- [numerical] V₂ grows more attractive as the field approaches the band, by about 0.03 to 0.05|t₂| per unit of μ_eff at r = 4 to 8, and the resonance crossing moves to smaller V/Δ: r_b(|w| = 1) ≈ 10 at μ_eff = −12, ≈ 8 at −8, ≈ 6 at −5, and below 4 at −4. The gas lives at the onset, μ_eff ≈ −2.2, beyond the range where the two-body extraction is clean (four-particle states come within reach and the fits degrade), but the trend leaves no unbound parity-class point there for r ≥ 4. This confirms the coordinator's correction (C186: bound at V/Δ = 8 at the onset) with an independent full-chain calculation.
- [numerical, correction to my own stage T claims] The stage T table and the "unbound corner r < 10/|w|" describe the chain far below the band, not the gas. Closer to the physical field the parity-class pair is bound (at μ_eff = −4 for r = 4 to 6, at −5 for r = 8), and the trend says it stays bound at the onset, so the two-component conclusion of stage T does not apply to the gas. The resonance structure of PROOF §2.2 and the three-body result (at the μ_eff = −12 interactions) stand as statements about those interactions.
- Consequence for the Typhon runs: families D and the stage F points were motivated by the far-field map. They remain data on the full chain but are no longer two-body favourable. Family E was dropped for this reason before the big job started.

### Three-body check (three183.py, three183.out; PROOF §4)

- [numerical] No trimer at any of the seven unbound points tested, (r, |w|) = (4, 1), (5, 1), (5, 0.5), (6, 1), (6, 1.5), (8, 1), (8, 0.5). The lowest three-body level at the aab and abb momenta is above 3ε_min and falls as 1/R² with the box (C175's code, validated to six digits against C175's bound trimer at r = ∞). So the cluster channel that C175 found at r = ∞ closes together with the pair channel, at the level of the fitted two-body interactions. Genuine three-body forces of the full chain are not included.

### Many-body onset and the local pilot (local/, PROOF §3)

- [numerical] Below the many-body onset the full chain is a gapped dressed vacuum whose density (virtual pairs) rises to 0.06 to 0.10 near the onset. The registered local pilot point r = 5, |w| = 1, x = 0.45 turned out to lie there: S̄ = 0.464 and 0.445 at L = 32 and 64, flat bulk S = 0.44, uniform density with no Friedel oscillation, exponentially decaying one-body correlations, truncation 1e-12 (local/Pr5_32.json, local/Pr5_64.json). So the pilot did not test the gas, and its registered parity-gap test was not carried out (stopped to free the CPU for the coordinator's r = 8 priority).
- [numerical, L = 24 ladders, local/p5, local/p8, local/p8b] The gas window is narrow and bounded above by a first-order jump. At r = 8, |w| = 1 the gas sets in at x ≈ 0.075 (φ = π/2) and 0.12 (φ = 1.555) and reaches n ≈ 0.3 by x = 0.3, then jumps to the packed phase (n ≈ 0.45) by x = 0.45. At r = 5 the jump lies between x = 0.6 and 1.0.
- These two findings changed the Typhon density targeting before that job started (amendment above).

### Report against the registered predictions so far

- T1: holds at μ_eff = −12 (r_b ≈ 10/|w|), but the premise (far-field interactions describe the gas) was wrong, see the correction.
- Local pilot: not carried out as registered. The r = 5 point lay below the onset (gapped vacuum), and the run was stopped. Nothing is claimed from it beyond the onset finding.
- D1, E1: E was dropped before any output; D runs in the big job.
- Q1 to Q3, F1/F2, G1, H1/H2, P1 to P4: pending on the Typhon jobs.

(Stage S, C and F results: to be filled in from the Typhon jobs.)

## What is left open

1. The many-body readings of the transverse-DM bond (family G, r = 12 and 16) and of C181's off-symmetric point (family H, r = 12, 16, 24) in the full chain, with confirmation to L = 256 and VUMPS FES (job 17604499), and the V/Δ = 8 points S1 and S2 to L = 256 with parity gaps (job 17604545).
2. A clean two-body (and four-body) extraction at the onset field itself, which needs larger rings or an open-chain two-excitation method, since ring ED at L ≤ 20 degrades above μ_eff ≈ −4.
3. The parity gap of the full chain is a fixed-field quantity (PROOF §3.2). Turning it into the symmetric E(N+1) + E(N−1) − 2E(N) needs a field scan at each L.
4. Families D and the stage F points were motivated by the far-field map. Their Typhon data remain informative about the full chain but no longer test a two-body-favourable regime.
