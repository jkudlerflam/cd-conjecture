This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

# C178 result: ring spectroscopy of the Fibonacci-cusp model (3.2) in the blockaded space

Derivations are in PROOF.md, the log in PROGRESS.md. Code: kring.py (momentum-resolved ED in the blockaded space, validated in validate.out and validate2.out), spec.py (first driver), spec2.py (flux-resolved, valley-labelled driver used for all results), analyze2.py, delta_R.py, casimir.py, parity.py, towers.py, xxz_fit.py, fit_runs.py. Data: local/*.jsonl and typhon_out/*.jsonl (Typhon job 17603981, partly collected). Analysis outputs: an_delta_R.out, an_casimir.out, an_parity.out, an_towers.out, an_analyze2.out.

The pre-registration below was written before any new spectrum was computed (the only spectra seen beforehand were C173's valley_ed.out and ring_levels.out, L ≤ 32). The results follow it.

## Pre-registration (written before looking)

Model. C173 PROOF eq. (3.2) on a ring of L sites, hard-core bosons with no two adjacent (the spin chain with periodic boundary conditions), units |t_2| = 1, w = |w| e^{iφ}, NNN hop −1, V_2 = −2, V_3 = −2/3, fixed N. The chemical potential only shifts levels at fixed N, so densities are set by N/L. The target density n ≈ 0.17 (μ = −2.2) is taken as n = 1/6 (L = 24, 30, 36, 42, 48, 54 with N = L/6), and the neighbouring densities n = 1/8 and n = 1/5. φ = π/2, 1.54, 1.50, |w| = 1.

Labels. M = N_a − N_b is the valley imbalance, measured from the JW-fermion momentum distribution (valley a near k_a = −s, valley b near k_b = π + s, s = arcsin(|w|/4) = 0.2527 at φ = π/2). At φ = π/2 the lattice symmetry S (inversion times the staggered sign on odd sites, L even) sends k → π − k, exchanges the valleys, and acts in the self-conjugate momentum sectors 2K ≡ Nπ.

What a c = 2 spectrum looks like (two-component Luttinger liquid, charge velocity v_c, pseudospin velocity v_s). In the dilute impenetrable picture of C173 PROOF 4 the pseudospin sector is an XXZ chain on the N particles, so v_s is expected to be much smaller than v_c, and the low part of the spectrum is pseudospin levels with energies of order v_s/L below the first charge phonon at 2π v_c/L.

- P1, pseudospin gap. Δ_ps(L) = energy of the lowest level with |M| ≈ 2 (one particle moved between valleys) above the ground state. Pass for c = 2: the ground state has M = 0 (balanced), and at fixed n the log-slope d ln Δ_ps/d ln L between the two largest even-N sizes lies in [−1.3, −0.7], so L·Δ_ps tends to a constant. Fail toward a gapped pseudospin (valley crystal or AF pseudospin, c = 1): the slope is above −0.5 and L·Δ_ps grows. Fail toward one valley (ferro pseudospin, c = 1): Δ_ps ≤ 0 (the polarized or partly polarized states drop below) at the largest size, or the slope is below −1.5 and heading through zero.
- P2, easy-plane versus easy-axis. In an easy-plane XXZ pseudospin chain the lowest |M| = 2 level (dimension 1/(4K_s)) lies below the lowest M = 0 level at the Néel momentum ΔP = π n (dimension K_s), so R = Δ(M = 0, ΔP = π n)/Δ_ps > 1, with R = 1 at the SU(2) point. In an easy-axis (gapped Néel, valley crystal) chain R < 1 and falls with L as the Néel level collapses onto the ground state. Pass for c = 2: R > 1 at the largest sizes and not decreasing with L. A fit of the low pseudospin levels to an N-site XXZ ring should give |J_z/J_⊥| < 1.
- P3, two velocities. The lowest M = 0 levels at ΔP = ±2π/L above the ground state should contain two single-boson levels with L·ΔE/2π → v_s and v_c, both stable (within 20 percent between the two largest sizes). For c = 1 there is one, and the next level at that momentum sits at three times the first. v_c should be close to the Girardeau value of the band curvature, v_c ≈ (π n/(1 − n))/m with 1/m = ε''(k_a) = 7.5 at |w| = 1, φ = π/2 (about 4.7 at n = 1/6), up to interaction corrections of order one.
- P4, Casimir. E_0(L) = e L − π(c v)_eff/(6L) at fixed n, with the boundary twist chosen to minimise E_0 (this removes the charge current zero mode). c = 2 predicts (c v)_eff = v_c + v_s, c = 1 predicts v_c. Because v_s is expected to be small and the incommensurate valley momenta leave a pseudospin zero mode of order v_s/L, this test separates the two only if v_s/v_c is above about 0.3. It is registered as secondary.
- Controls at the same sizes. c = 1: hard rods (blockade, hop w only, no NNN hop, no V_2, V_3). c = 2: two decoupled XX legs (w = 0, NNN hop that ignores the middle site, no blockade, no V_2, V_3). The controls must pass P3 and P4 with their known velocities (to 10 percent at L ≥ 36) for the method to count.

Verdict rule. c = 2 is supported at a point if P1 and P2 pass at the largest two sizes, P3 shows a stable second velocity, and P4 is not in conflict. A clean P1 failure (gap saturating, or a polarized ground state) is a c = 1 verdict at that point. A slope between −0.7 and −0.5, or R near 1, is inconclusive (an SU(2)-like pseudospin with a BKT-slow gap cannot be separated from a marginal c = 2 at these sizes, and this is stated as a limitation, not as a pass).

## Verdict

Inconclusive by the pre-registered rule, with a definite trend toward c = 1.

At every size reached (N ≤ 8 particles, L ≤ 40), the twist-relaxed ring spectrum at n = 1/6 to 1/4 is quantitatively that of a two-component liquid. It has two velocities, the combination levels of a two-boson tower, a Casimir term equal to their sum, and a pseudospin gap scaling as 1/L. But the one registered indicator that separates easy-plane from easy-axis, R = D_Néel/D_ps, falls monotonically with N at every density and |w|, and it does not slow down. In the twist-free series it goes 1.67, 1.51, 1.35 at n = 1/5 and 1.72, 1.42, 1.12 at n = 1/4, against R = 1 at the SU(2) point. At lower density (n = 1/12) it is already below 1 and falling, and in the strong-binding calibration (|w| = 2) it collapses to 0.02 or less, the charge-2 flux period of a pair liquid. Together with C175's inter-valley bound pair, the natural reading is a marginally relevant flow into a pseudospin-gapped pair liquid (c = 1), with a crossover scale beyond these rings that shrinks as the density falls. Rings of this size cannot exclude a flow that stops at an easy-plane fixed point.

## Method in one paragraph

Ring ED of (3.2) in the blockaded space, resolved in momentum and N, at up to 1.9 million states per sector locally and on Typhon (kring.py). Every level carries a valley label C1 = Σ_k cos k n_k (M ≈ C1/cos s). The flux θ through the ring is scanned. By the large gauge transformation, a level followed in θ over the real line visits the sectors m, m + N, ..., so minimising over θ and m ≡ κ mod N removes the charge current zero mode exactly. The residual label κ is the pseudospin chain momentum (PROOF 1, 2) [proved for the gauge identity, derived for the spin-charge reading]. Where the bare valley spacing makes the pseudospin chain twisted (Φ_s = (k_b − k_a)L mod 2π), |w| is retuned to 4 sin(πj/L) so that Φ_s = 0. These are the "twist-free" series (|w| = 1.0353 at n = 1/6, 1.2361 at n = 1/5, 1.5307 at n = 1/4).

## Key claims

1. [proved, numerical check] Validation. The union of the momentum sectors reproduces C173's full-basis ring ED to 1e-13 for L = 14 to 16, φ = π/2 and 1.3, flux 0 and π. n(k) agrees with C173's obdm to 1e-15, S commutes with H and squares to 1, and the flux and S images used to save sectors are exact (validate.out, validate2.out).
2. [proved] Gauge identity F(m, θ + 2π) = F(m + N, θ), so twist relaxation over κ = m mod N removes the charge zero mode for any interaction (PROOF 1c).
3. [derived, checked to 5 percent] Spin-charge factorisation. Energies at flux 0 and π obey E(0, π) − E(π, 0) = E(π, π) − E(0, 0) (−0.135 and −0.127 at L = 36, N = 6; +0.184 and +0.175 at L = 24, N = 4). The pseudospin q = π level alternates with N/2, like an antiferromagnetic or XY ring (PROOF 2; local/n6_pi2.out).
4. [numerical, an_casimir.out] Controls pass P3 and P4. Hard rods (c = 1) give (c v)_eff = 1.42 and 1.78 at n = 1/6 and 1/5 against v = 1.41 and 1.76. Two XX legs (c = 2) give 4.05, 4.78 and 5.81 at n = 1/6, 1/5, 1/4 against 2v = 4.00, 4.70 and 5.66.
5. [numerical, an_delta_R.out] P1 passes. The pseudospin gap D_ps (lowest |M| = 2 level) scales as 1/L in every φ = π/2 series. The log-slopes are −1.02 (n = 1/6, |w| = 1, L 24 → 36), −1.02 (n = 1/6, twist-free), −1.03 and −1.01 (n = 1/5 twist-free, L 20 → 30 → 40), −1.05 and −1.02 (n = 1/4 twist-free, L 16 → 24 → 32), and −1.02 (n = 1/12, L 24 → 48). L·D_ps = 3.72, 3.82, 4.46 and 4.93 at the largest sizes. The ground state is valley-balanced (M = 0) and the polarised states lie far above (D_4 ≈ 4 D_ps). C175 claim 9 reports a valley-flip energy with L·E = 2.5, 3.8, 4.2 growing at L = 12, 24, 36. That number is not twist-relaxed. At fixed flux the flipped state also pays a charge zero-mode energy of order v/L that varies with L, so the growth does not by itself indicate a gap. Twist-relaxed, L·D_ps = 3.75 and 3.72 at L = 24 and 36.
6. [numerical, an_delta_R.out] P2 fails its "not decreasing" clause everywhere.

   | series (φ = π/2) | Φ_s | R at N = 4, 6, 8 | Δ_R (XXZ ring with the same N, Φ_s) |
   | :-- | :-- | :-- | :-- |
   | n = 1/4, \|w\| = 1.5307 | 0 | 1.72, 1.42, 1.12 | 0.50, 0.72, 0.92 |
   | n = 1/5, \|w\| = 1.2361 | 0 | 1.67, 1.51, 1.35 | 0.53, 0.67, 0.78 |
   | n = 1/6, \|w\| = 1.0353 | 0 | 1.68, 1.57, pending | 0.52, 0.64 |
   | n = 1/6, \|w\| = 1 | 0.14π, 0.21π | 1.47, 1.22 | 0.61, 0.76 |
   | n = 1/12, \|w\| = 1 | 0.14π, 0.28π | 0.86 (N = 2), 0.64 (N = 4) | 1.29, 1.36 |
   | n = 1/6, \|w\| = 2 | 0 | 0.35, 0.003 | > 3 |
   | n = 1/12, \|w\| = 2 | 0 | 0.17 (N = 2), 0.019 (N = 4) | > 3 |

   R is above 1 (easy-plane-like) at n ≥ 1/6 and |w| ≈ 1 to 1.5, but it falls by 0.1 to 0.3 per step in N with no deceleration in ln N. At n = 1/4 a straight extension reaches R ≈ 0.9 at N = 10 [heuristic extrapolation]. At n = 1/12 and at |w| = 2 it is below 1 and falling. Runs with a large residual twist (Φ_s near π: n = 1/5 at |w| = 1, n = 1/8, n = 1/4 at |w| = 1, L = 32) give R values dominated by the twist and are not used.
7. [numerical, an_towers.out, an_casimir.out] P3 and P4 pass at the sizes reached. In the κ = 1 sector the two lowest M = 0 levels give v_1 and v_2. The combination level v_1 + v_2 appears twice at κ = 0, as a two-boson theory requires.
   - n = 1/5 twist-free. v_1 = 1.19, 0.99, 0.905 and v_2 = 1.96, 2.07, 2.05 at L = 20, 30, 40. At L = 40 the κ = 0 pair sits at 0.468 and 0.471, against 0.464 predicted. The Casimir fit over the three sizes (residual 2e-6) gives (c v)_eff = 2.83, against v_1 + v_2 = 2.96 at L = 40 and v_2 alone = 2.05.
   - n = 1/6 twist-free, L = 36. v_1 = 0.94, v_2 = 1.63, so v_1 + v_2 = 2.57, against (c v)_eff = 2.51.
   - n = 1/4 twist-free, L = 32. v_1 = 0.83, v_2 = 2.47, so v_1 + v_2 = 3.30, against (c v)_eff = 3.12.

   So the Casimir coefficient at L ≤ 40 needs both modes. The slow mode softens with N (v_1 = 1.15, 0.95, 0.83 at n = 1/4).
8. [numerical, an_parity.out] Parity gap. L·Δ_pb is flat at |w| = 1: 2.75 and 2.70 (n = 1/6, L = 24, 36), 2.77 and 2.68 (n = 1/5), 3.07 and 3.02 (n = 1/5 twist-free). No pairing gap is visible on these rings, so at n ≥ 1/6, |w| = 1 the pair gap is at most of order 2Δ_s ≲ Δ_pb(L = 36) ≈ 0.07 [numerical bound, not a measurement]. At |w| = 2 the parity gap is not a clean diagnostic either (2.73 → 2.16 at n = 1/6), while R and the flux period are (claim 6). The flux period δ_π = E_min(θ = π) − E_min(θ = 0) shows the charge-2 collapse at |w| = 2 (L·δ_π = 0.86 → 0.12 at n = 1/12) and not at |w| = 1, n = 1/6 (|L·δ_π| = 4.2, 4.2 at L = 24, 36).
9. [numerical, an_casimir.out] The charge sector is very soft. v_N = (L/π)∂²E/∂N² = 0.25 at n = 1/5 (twist-free, stable between L = 20 and 30) and 0.27 at n = 1/6, against v_J = 2.7. At |w| = 2, n = 1/6 the compressibility is already negative at L = 36 (clustering, as C175 found). Because the flux couples to a combination of the two modes, (v_J v_N)^{1/2} is not a mode velocity here (PROOF 4). But the ratio (v_J/v_N)^{1/2} ≈ 3.2 is far above the value 2 of two free components, which points the same way as C179's R4: the two-pair vertex has dimension 4/K_ρ, below 2 when K_ρ > 2 [heuristic link, C179's convention].
10. [numerical, local/phi.jsonl, an_analyze2.out] Off the symmetric point. At φ = 1.54 the M = −2 level drops toward the ground state (L·D_ps = 1.23 at L = 24, 0.81 at L = 30, n = 1/5), a valley field δ = 4|w| cos φ = 0.12 partly polarising the pseudospin. At φ = 1.50 (δ = 0.28) the relaxed ground state already has M = −2 at N = 4 and 6, heading to one valley. This matches C173's DMRG (c ≈ 2 readings persist to 1.54 and are gone by 1.50) if the pseudospin saturation field J(1 + Δ) is about 0.2 to 0.3 [heuristic, with J ≈ 0.15 to 0.19 from the XXZ fits].
11. [numerical, fit_runs.py] A nearest-neighbour XXZ ring does not fit the full set of pseudospin levels (rms 20 to 27 percent). The κ ≠ 0, N/2 levels sit at about half their XXZ values, because the slow mode mixes with the charge. Δ_R is therefore an effective label from one ratio, not a fitted anisotropy.

## Against the pre-registration

- P1 (D_ps ∝ 1/L, balanced ground state). Pass at every φ = π/2 series, including the pair-liquid calibration at |w| = 2, n = 1/6 (slope −1.02). So P1 does not discriminate at these sizes. That is a lesson about P1, stated after the fact.
- P2 (R > 1 and not decreasing). Fail. R > 1 at n ≥ 1/6 and |w| ≈ 1 to 1.5, but it decreases with N in every series.
- P3 (two stable velocities). Pass at n = 1/5 (v_2 stable to 1 percent, v_1 to 9 percent between L = 30 and 40), and seen at n = 1/6 and 1/4.
- P4 (Casimir consistent with v_1 + v_2). Pass to 3 to 6 percent at n = 1/6, 1/5 and 1/4. The controls pass to 1 to 3 percent.
- Verdict rule: c = 2 needs P1 and P2. P2 fails, so c = 2 is not supported. P1 passes, so this is not a clean c = 1 verdict either. Inconclusive, with the R drift as the main finding.

## What is left open

1. Whether R crosses 1. The pending Typhon outputs decide it at the next size. Registered now, before they are seen. If the flow continues, the n = 1/4 twist-free point at L = 40, N = 10 gives R < 1.0 (straight extension: 0.85 to 0.95), and the n = 1/6 twist-free point at L = 48, N = 8 gives R ≤ 1.50. If the flow stops (an easy-plane fixed point), R ≥ 1.1 at L = 40, N = 10 and R ≥ 1.55 at L = 48, N = 8.
2. The size of the pair gap at n ≈ 0.17, |w| = 1. The rings bound 2Δ_s only by about 0.07. An open-chain parity gap at L ≥ 200 (C175's suggestion) is the next test.
3. The identity of the slow mode (v_1 ≈ 0.9, softening with N) and why the nearest-neighbour XXZ description of the pseudospin fails at κ ≠ 0, N/2.
4. φ = 1.54 and 1.50 beyond L = 30, and the L = 48 runs at φ ≠ π/2. These were dropped from the job for time.

## Typhon job 17603981 (running at write-up)

One packed node (48 cores, partition typhon), started 2026-11-09 10:01 EST, 20 h limit, directory /scratch/lustre/jkudlerflam/c178/. Step .batch runs the original queue (typhon/c178.slurm with heavy.txt and medium.txt). The heavy L = 50, N = 10 task was killed, because it would not have finished. Step .1 (typhon/pool2.sh, srun --overlap inside the same job) runs medium.txt from the end. run_task.sh uses lock directories out/locks/<tag> so no task runs twice. The L = 72 tasks fail, because the code is limited to L ≤ 62, and are skipped. So are the L = 48 runs at φ ≠ π/2 and |w| = 0.5, and the rods at L = 48 and L = 40, N = 10.

Expected outputs: out/<tag>.jsonl (one JSON record per task, written at the end), out/<tag>.out, logs/done.txt. The most important pending tags are cusp_ppi2_w1.0353_L48_N8, cusp_ppi2_w1.0_L48_N8, cusp_ppi2_w1.5307_L40_N10, cusp_ppi2_w1.0_L48_N9 and cusp_ppi2_w2.0_L48_N9 (parity at L = 48), cusp_ppi2_w1.0_L56_N7, cusp_ppi2_w2.0_L48_N8, cusp_ppi2_w1.0_L40_N10, cusp_p1.54/1.50_w1.0_L40_N8, rods_ppi2_w1.0_L40_N8, cusp_ppi2_w1.0_L48_N6 and cusp_ppi2_w1.0_L42_N7. Collect with

    scp -r typhon:/scratch/lustre/jkudlerflam/c178/out typhon:/scratch/lustre/jkudlerflam/c178/logs "codex/work_C178/typhon_out/"

from the cd_conjecture directory, then rerun

    python3 delta_R.py local/v2_*.jsonl local/dil.jsonl local/tf.jsonl local/dN.jsonl local/phi.jsonl typhon_out/out/*.jsonl
    python3 casimir.py local/v2_*.jsonl local/dil.jsonl local/tf.jsonl local/dN.jsonl local/phi.jsonl typhon_out/out/*.jsonl
    python3 parity.py local/v2_*.jsonl local/dil.jsonl local/tf.jsonl local/dN.jsonl local/phi.jsonl typhon_out/out/*.jsonl

(local/n6_pi2.jsonl and local/legs.jsonl are in spec.py's older format and are read only by hand), then score item 1 above.
