This run was done by a Claude Opus 5.5 agent (C177, 2026-11-09), not by Codex. Astra was not used (out of credits).

# C177 result: the full qubit chain at the Fibonacci cusp

Derivations are in PROOF.md and the log in PROGRESS.md.

## Pass criterion (registered before any full-chain production output was looked at)

Written at the start of the run, before the first full-chain DMRG or VUMPS output. It is C163 §4.2 as adapted in the spec, made numerical in the form C166 used.

Bonds. F_r for r = V/Δ in {8, 12, 16, 24}: h = J_x XX + J_y YY + D(XY − YX) + J_z ZZ + (h_z/2)(Z⊗1 + 1⊗Z) with Δ = J_x − J_y = 1, J_x + J_y = 0, J_z = V/4, D tuned so that the full chain's one-particle band has |w|/|t_2| = 1, and h_z tuned so that the bulk particle density is n ≈ 0.17. G: one generic point off the parity class and off J_x + J_y = 0 (small J_x + J_y, XZ + ZX and a transverse field h_x, all of order Δ²/V).

A bond is a c > 1 candidate only if all three hold.
- (P1) Open chains. With S̄(L) the mean entanglement over the three middle cuts, the two-size slope 6 [S̄(256) − S̄(128)]/ln 2 is above 1.2, the profile fit c_prof (C166 dmrg_c.profile_fits) is above 1.2 at L = 128 and at L = 256 with c_prof(256) ≥ c_prof(128) − 0.03, and the maximal discarded weight in the last sweep is below 1e-8 at both lengths.
- (P2) VUMPS. In the cell (n = 2 or 4) with the lowest variational energy at the largest χ, the FES slope (least squares 6 dS/d ln ξ over the last three χ, χ_max ≥ 256) is above 1.2.
- (P3) Energies. The open-chain bulk bond energy e_mid (mean of ⟨h_{j,j+1}⟩ over the central L/4 bonds) at L = 256 agrees with the lowest VUMPS e₀: |e_mid(256) − e₀| ≤ 2 |e_mid(256) − e_mid(128)| + 1e-6. If e_mid lies above e₀ by more, the open chain misses the VUMPS state. If it lies below by more than the χ-error of e₀, VUMPS is trapped.

Failing any of P1 to P3 at every bond tested means the full chain does not support the c = 2 candidate at the sizes reached. Passing at F_8 to F_24 but not at G means the c = 2 region is not shown to be open in B₂.

Registered expectations for the local stage (L ≤ 64), used to read the data and not as a criterion.
- If the effective model (C173 eq. 3.2) is a faithful description and is c = 2, then at band-tuned couplings and n ≈ 0.17 the full chain at V/Δ = 8 to 24 reproduces the effective-model S̄ at the same L within about 0.1 (C173: S̄(48) = 1.646, S̄(72) = 1.790, S̄(96) = 1.869), and two-size estimates between L = 32 and 64 lie between 1.6 and 2.3.
- If instead S̄ at L = 48 is below 1.4 for all V/Δ, or two-size estimates are below 1.2 at L ≤ 64 at all V/Δ while the effective model reads above 1.6, then the higher-order terms dropped in (3.2) change the low-energy physics at these V/Δ, and only the large-V/Δ runs of C173's job can test the candidate.

## Typhon job

Job 17603971 (submitted 2026-11-09, one node, 48 cores as 12 slots of 4 BLAS threads, 72 h limit), directory /scratch/lustre/jkudlerflam/c177/, tasks in typhon/tasks.txt, script typhon/c177.slurm. Expected outputs:
- oc/F8_{64,128,192,256}.json, oc/F12_{128,256}.json, oc/F16_{64,128,192,256}.json, oc/F24_{128,256}.json, oc/G8_{128,256}.json, oc/F8m23_128.json (F_8 at the lower dressed μ_eff = −2.3): open chains, χ up to 384 (F_8, F_16 at L = 256), 320 or 256 otherwise.
- vc/F8_n2.json (χ ≤ 384), vc/F8_n4.json, vc/F16_n2.json (χ ≤ 320), vc/F12_n2.json, vc/F24_n2.json, vc/G8_n2.json (χ ≤ 256): damped VUMPS with FES data and the order parameters.
- oc/C8_{128,256}.json: the single-valley c = 1 control (F_8 with J_x + J_y = 0.15|t_2|), added to the task list before the job started.
- logs/, and summary.txt written by summ177.py at the end of the job (section 4 evaluates P1 to P3).
- Not in the job: two parity-gap lines (F8p*, F16p* at L = 64 to 128) were appended to typhon/tasks.txt after the job had started and its task list had been read, so they will not run. The parity gap was done locally instead (below).

To collect, from inside codex/work_C177:

    mkdir -p typhon_out && scp -r typhon:/scratch/lustre/jkudlerflam/c177/oc typhon:/scratch/lustre/jkudlerflam/c177/vc typhon:/scratch/lustre/jkudlerflam/c177/logs typhon:/scratch/lustre/jkudlerflam/c177/summary.txt typhon_out/

then `python3 summ177.py typhon_out > typhon_out/summary_local.txt`. The vc/*.npz checkpoints are large and need not be copied (delete them from the copy or skip vc/ and copy vc/*.json by name if bandwidth matters).

## Results

All numbers below are from files in this folder or on Typhon (named). Units |t_2| of each bond unless stated.

### R1. The full chain maps onto the effective model with (Δ/V)² corrections [numerical, band.out, tune.json, PROOF 1]

The exact one-particle band (ring ED) keeps the k → π − k symmetry to all orders at J_x + J_y = 0 (antisymmetric harmonics below 1e-5), and the lattice symmetry behind it is proved (PROOF 1.3). |t_2| is renormalized by 1 − 5.5(Δ/V)² (0.926 at V/Δ = 8, 0.997 at 40), the NN hop by less than 0.3 percent, and extra symmetric harmonics stay below 2 percent of |t_2|. D and h_z were tuned on the exact band (|w|/|t_2| = 1, dressed μ_eff = −2.2). At V/Δ = 40 this reproduces C173's candidate bond (h_z to 3e-5, D to 0.3 percent). Raw densities at the tuned h_z are 0.17 to 0.18, including a virtual-pair background of 0.023 (V/Δ = 8) to 0.003 (V/Δ = 24).
Caveat found later (R2). The band was measured at μ_eff = −12. At the physical μ (number-filtered one-excitation states, twob_phys_*.out) the V/Δ = 8 band has |t_2|/t_2^LO = 0.98 to 1.00 and, with the production D, |w|/|t_2| = 0.93 (μ_eff = −2.4) to 0.95 (−2.6), not 1. So the production bonds sit on the symmetric stratum (that is exact) with |w|/|t_2| a few percent below 1 at V/Δ = 8 and closer to 1 at larger V/Δ. The 1 − 5.5(Δ/V)² renormalization is a statement at μ_eff = −12, where every energy denominator carries the extra −2μ.

### R2. The two-excitation bound state at K = π (coordinator's item i, revised after the C186 conflict) [numerical; twob_full.out, twob_fit.out, twob_phys_mu-2.4.out, twob_phys_mu-2.4_retune.out, twob_phys_mu-2.6.out; PROOF 3]

First pass (withdrawn as a statement about the physical binding). Ring ED at L = 14 to 20 matched to an N-conserving two-body model, with everything measured at μ_eff = −12 (particles far above the vacuum), then solved at L = 400. This gave v = V_2/|t_2| + 2 ≈ 0.27, 0.13, 0.08 at V/Δ = 8, 12, 16 (about 17 to 20 (Δ/V)²), unbound at V/Δ = 8 and bound at V/Δ ≥ 12. It is a correct description of the full chain at μ_eff = −12 (fit rms 2e-4|t_2|, one-excitation ⟨N⟩ = 1.01), but not at the physical chemical potential, for the reason C186 gives: the near-resonant momentum-zero ΔN = 2 vertex and the μ-dependence of every energy denominator make the N-conserving couplings depend on μ.

Direct test at the physical μ (coordinator's request). Same rings, μ_eff = −2.4 and −2.6 (0.15 and 0.35|t_2| below the band bottom, the dressed vacuum is the ring ground state with ⟨N⟩ equal to the virtual-pair background), states selected by momentum, parity and ⟨N⟩. Then B(L) = E_2(K = π) − E_vac − 2 min_k [E_1(k) − E_vac], with no model.
- The excitations are strongly dressed there: one-excitation ⟨N⟩ = 1.29 to 1.82 at V/Δ = 8 and 1.13 to 1.41 at V/Δ = 12 (against 1.01 at μ_eff = −12), two-excitation ⟨N⟩ = 2.24 to 2.53. This is C186's near-resonant ΔN = 2 mixing, seen directly.
- B(L)/|t_2| at V/Δ = 8 (μ_eff = −2.4, D retuned so |w|/|t_2| = 1 at this μ): +0.084, +0.025, −0.004, −0.018 at L = 14, 16, 18, 20. The second-order model (3.2) with the same band gives +0.036, −0.027, −0.058, −0.071. At V/Δ = 12 (μ_eff = −2.4): +0.058, +0.006, −0.019 at L = 14, 16, 18 against +0.042, −0.014, −0.040. With the production D (|w|/|t_2| = 0.93 to 0.95 at the physical μ) and at μ_eff = −2.6 the numbers are similar (V/Δ = 8, L = 20: −0.013 against −0.048).
- So the full chain is less attractive than (3.2) at the same L, by 0.02 to 0.05|t_2|, but its lowest K = π two-excitation level already lies below twice the one-excitation minimum at L = 18 to 20, as in the bound second-order model. An N-conserving fit at the physical μ is poor (rms 0.8 to 2.6e-2|t_2|, against 2e-4 at μ_eff = −12) and extrapolates to weak binding, E_2 − 2ε_min = −0.003 to −0.016|t_2| at V/Δ = 8 depending on μ and |w|, with v between −0.16 and +0.08 and g between 0.24 and 0.48 (D between −0.76 and −0.23).

Answer to the conflict. C186's number (v ≈ 3.5 (Δ/V)², an N-conserving part close to second order plus an explicit near-resonant ΔN = 2 vertex) is the one consistent with the full chain at the physical μ. My v ≈ 20 (Δ/V)² is an artefact of measuring at μ_eff = −12, where the ΔN = 2 and μ-dependent virtual processes are detuned by about 20|t_2| instead of about 0.5|t_2|, so it is absorbed differently into an effective V_2 (the coordinator's guess is right in substance). The physical two-excitation binding at V/Δ = 8 and 12 is not resolved at L ≤ 20, but the direct ring data at the physical μ point to a weak bound pair, not to an unbound one. My earlier statement "unbound at V/Δ = 8, D changes sign near V/Δ ≈ 10" is withdrawn. The pairing obstruction of C175 should be assumed present at V/Δ = 8 and 12, weakened relative to second order, and its quantitative size is not controlled by any N-conserving description (C186's warning stands).

### R3. Local open chains, L ≤ 96 (first look, not the criterion) [numerical, local/*.json, sfit.py, cmp_eff.py]

| bond | L | S̄ (3 middle cuts) | c_prof | trunc | n_mid (raw) |
| :-- | :-- | :-- | :-- | :-- | :-- |
| F_8 | 24 | 1.249 | 3.56 | 2e-11 | 0.177 |
| F_8 | 32 | 1.174 | 2.15 | 2e-10 | 0.162 |
| F_8 | 48 | 1.451 | 2.22 | 6e-9 | 0.184 |
| F_8 | 64 | 1.469 (Typhon 1.471) | 2.03 (2.04) | 2e-8 (7e-9) | 0.178 |
| F_8 | 96 | 1.639 | 2.16 | 1.4e-8 | 0.184 |
| F_8 | 128 (Typhon) | 1.763 | 2.42 | 1.0e-8 | 0.179 |
| F_12 | 128 (Typhon) | 1.791 | 2.20 | 7.5e-9 | 0.181 |
| F_16 | 48 | 1.440 | 2.19 | 2e-9 | 0.173 |
| F_16 | 64 (Typhon) | 1.525 | 2.14 | 4e-9 | 0.177 |
| C_8 (c = 1 control) | 32 | 1.007 | 0.98 | 7e-12 | 0.152 |
| C_8 (c = 1 control) | 48 | 0.961 | 1.11 | 6e-11 | 0.149 |

- The profile fits read c ≈ 2.0 to 2.2 at L = 48 and 64 for V/Δ = 8 and 16, and also with one or two fitted incommensurate oscillation wavevectors (2.1 to 2.3; the oscillations sit at q ≈ 2.7 ≈ π − 2k_a, the inter-valley interference, and at q ≈ 0.46 to 0.59 ≈ 2k_F per valley). The single-valley control reads 0.98 and 1.11 with the same estimator. Two-size estimates are noisy at small L (S̄ is not monotonic between L = 24 and 32), but between 48 and 96 they are stable: at the cuts C173 used (L/2 − 3, L/2, L/2 + 3) F_8 gives 6[1.6545 − 1.4100]/ln 2 = 2.12, the central-quarter average gives 2.27, and the profile fit at L = 96 gives 2.16 (2.29 to 2.31 with oscillation terms). The effective model reads 1.93 for 48 → 96 (C173).
- Against the effective model (C173, block cuts at L = 48: 1.632, 1.676, 1.632), the full chain at the same cuts is lower by about 0.25 (1.407, 1.417, 1.407 at V/Δ = 8). It is nearly V/Δ-independent between 8 and 16 (S̄ = 1.451 and 1.440). At L = 96 the gap to the effective model is 0.21 (1.654 against 1.866 at the same cuts), so the constant differs and the slope does not. So the expectation registered for "faithful description" (within 0.1) is not met for the constant, while the profile slopes agree. A boundary chemical potential of order |t_2| that the effective model omits (an edge particle blocks one bond, not two) can shift the constant [heuristic].
- No valley polarization. The open-chain profile of P_j = X_jX_{j+1} + Y_jY_{j+1} vanishes by symmetry, and its correlator ⟨P_{L/4}P_{L/4+r}⟩ at F_8 decays to below 5e-3 by r = 12 (L = 32 and 48), while the single-valley control has a plateau 0.23 (= P², L = 48). At L = 96 (F_8) it is −0.009, −0.003, −0.002, −0.001 at r = 12, 24, 32, 48. VUMPS at F_8, F_12, F_16 (χ ≤ 128) has ⟨P⟩ < 2e-4 in a two-site cell [numerical, Typhon vc logs]. The vector-chiral (single-valley) c = 1 state is excluded at these sizes.

### R4. Parity gap (coordinator's item ii) [numerical, local/par/*.json; sectors selected exactly with λ ∏Z, checked against ED to 1e-13]

Δ_P = E_odd − E_even on open chains at the tuned h_z (grand canonical in N, so the charge part enters with a sign that depends on L):

| L | F_8: Δ_P/|t_2| | L Δ_P/|t_2| | F_16: Δ_P/|t_2| | L Δ_P/|t_2| |
| :-- | :-- | :-- | :-- | :-- |
| 24 | +0.0082 | +0.20 | +0.0085 | +0.20 |
| 32 | −0.0021 | −0.07 | −0.0026 | −0.08 |
| 48 | +0.0028 | +0.13 | +0.0042 | +0.20 |

The ground-state parity alternates and L Δ_P stays bounded at L ≤ 48, so no pair gap is resolved at these sizes. This does not exclude a pair gap well below 0.2/L|t_2|, which is what C175's estimate (about 1e-4|t_2| at V/Δ = ∞) and R2 (weaker binding at V/Δ = 12 to 16, none at 8) would give.

### R5. VUMPS so far (job 17603971, χ ≤ 128 at the time of writing) [numerical, Typhon logs/vF*_n2.log]

F_8 (n = 2): S = 1.173, 1.364, 1.466, 1.540, 1.639, 1.704, 1.749 and ξ = 50, 87, 130, 172, 263, 349, 420 at χ = 16, 32, 48, 64, 96, 128, 160. Pairwise FES slopes 2.09, 1.56, 1.41, 1.38, 1.47, not yet settled. e₀(χ = 160) = −2.0057036, against the open-chain e_mid = −2.0057487, −2.0057576, −2.0056850 at L = 64, 96, 128 (the central-quarter average still oscillates by about 5e-5 with L). F_16: slopes 1.67, 1.78 at χ ≤ 64, but err 1.8e-3 at χ = 128 (not converged). P2 needs χ ≥ 256.

### R6. Openness in B₂: the generic point is gapped [numerical]

G_8 (V/Δ = 8 plus J_x + J_y = 0.02|t_2|, J_xz(XZ + ZX) = 0.25|t_2|, transverse field h_x = 0.25|t_2|, all below 0.13 Δ²/V): open chains at L = 128 and 256 both read S̄ = 0.371 with truncation 2e-20, and VUMPS gives S = 0.371, ξ = 16.8 at χ = 192 and 256, density locked at 0.2000, valley polarization ⟨P⟩ = 0.233 (oc/G8_{128,256}.json, vc/G8_n2.json on Typhon). This is a gapped state, c = 0.
Smaller parity-breaking terms at L = 48 (local): h_x = J_xz = 0.05 and 0.1|t_2| (with J_x + J_y = 0.02|t_2|) give S̄ = 0.59 and 0.51; h_x = 0.02|t_2| alone gives 0.79; J_xz = 0.02|t_2| alone gives 0.69; against 1.45 for F_8. Perturbations of 0.02|t_2| ≈ 1e-3 in absolute units therefore change the state at L = 48 by order one. The parity-preserving valley-offset term alone, J_x + J_y = 0.02|t_2|, does not (S̄ = 1.452, c_prof = 2.15, local/Gjp0.02_48.json), so the effect comes from the parity-breaking terms.
[heuristic] The parity-breaking terms create single particles at momentum 0, and the band minimum lies only k_a ≈ 0.25 away with a particle cost ε(0) − μ ≈ 0.2|t_2|, so they act almost resonantly on the dilute gas. Whether a floating (incommensurate) window survives for small enough amplitudes is not decided here. By the registered rule, the c = 2 region is not shown to be open in B₂.

### R7. Density sensitivity [numerical, oc/F8m23_128.json]

Lowering the dressed μ_eff from −2.2 to −2.3 at V/Δ = 8 drops the density from 0.18 to 0.042 (L = 128, S̄ = 0.21). n(μ) is very steep there, consistent with the attractive interactions and C175's clustering.

## Verdict (interim; P1 to P3 wait on job 17603971)

- The full qubit chain at V/Δ = 8 to 24 reproduces the qualitative reading of the effective model at L ≤ 64: profile fits c ≈ 2 (single-valley control ≈ 1), no valley polarization, no resolved parity gap. This is consistent with c = 2 and also with a pair gap or crossover below the scales reached [numerical].
- Open chains on Typhon at L = 128 read S̄ = 1.763 (V/Δ = 8) and 1.791 (V/Δ = 12) with c_prof = 2.42 and 2.20 at truncation 1.0e-8 and 7.5e-9. The two-size slope 64 → 128 at V/Δ = 8 is 2.5. Nothing in the open-chain data yet points to c = 1, and the L = 256 runs decide P1.
- C175's bound pair: after the C186 conflict (R2, revised), the direct ring data at the physical μ point to a weak bound pair at V/Δ = 8 and 12, less bound than at second order, and the N-conserving extrapolation is not controlled. So the pairing obstruction is not removed at small V/Δ. If it holds, the c ≈ 2 readings are a crossover over a pair correlation length longer than the sizes reached.
- The tested generic point off the parity class is gapped, and parity-breaking terms as small as 0.02|t_2| change the L = 48 state strongly (R6). By the registered rule the c = 2 region is not shown to be open in B₂, and the evidence points to the parity-breaking terms being dangerous, which would put the candidate (if c = 2 at all) on the parity-class stratum, i.e. at codimension ≥ 1, not open in B₂.
- The registered criterion P1 to P3 is not yet evaluated. Run summ177.py on the collected outputs.
