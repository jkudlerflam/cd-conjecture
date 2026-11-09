This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

# C181 result: the partially polarized pair route off the symmetric line

Derivations are in PROOF.md and the log in PROGRESS.md. Copied unmodified: fibdmrg.py, analyze.py, run.py (C176), scat.py, nbody_box.py (C175).

## Pre-registered predictions, few-body part (written before running fewbody.py)

F1. The ab pair binding relative to ε_a + ε_b stays within 2 percent of its φ = π/2 value for φ in [1.48, π/2] at |w| = 1 and 2 (C175 phi_scan already shows this to 1.54).
F2. The aab trimer binding B3, relative to the pair at its minimum plus an a atom at the valley minimum, stays positive over φ in [1.50, π/2] at |w| = 2 and changes by less than 30 percent there. At |w| = 1 it stays of order 1e-3 or below.
F3. The valley offset δ(φ) is linear in π/2 − φ near π/2 with slope about 3.9 (|w| = 1) and 6.9 (|w| = 2), from C175's phi_scan.


## Pre-registered predictions and criteria, many-body part (written before any production output exists)

Registered 2026-11-09 before submitting the Typhon job, after the few-body runs and a design pilot (local/pilot_w2_L48.out, E(N) at φ = π/2 only, used to pick densities, and a 16-site smoke test).

Definitions (C176's, unchanged). S̄_W(L) is the mean of S(l) over cuts l in [L/2 − W, L/2 + W]. c_two(64 → 128) = 6[S̄_W(128) − S̄_W(64)]/ln 2 at fixed n = N/L, for W = 1, 6, 12. c_prof(L) is the smoothed-profile fit of C176's analyze.py. Truncation is the largest discarded weight of the final sweep. Parity gap Δ_P(L) = E(N+1) + E(N−1) − 2E(N) at even N, with the N ± 1 runs at the same final χ (384). Valley weights W₊ and W₋ are the fermion n(k) weights on cos k > 0 and cos k < 0 (C176's Wa_f, Wb_f), and the valley polarization is p = |W₊ − W₋|/(W₊ + W₋). Homogeneity: the density profile smoothed over 16 sites varies by less than 25 percent over the central half of the chain.

First-stage c = 2 criterion at a point (all must hold). c_two(64 → 128) > 1.5 for W = 6 and W = 12, c_prof(128) > 1.2, truncation below 1e-7 at L = 128 (otherwise the point is reported as unconverged), 0.05 < p < 0.9, L·Δ_P grows by less than 50 percent from L = 64 to 128, homogeneous. Calibration: the c = 1 control (|w| = 1, φ = 0, V₂ = V₃ = 0, n = 3/16) reads within 0.15 of 1 and the c = 2 control (two decoupled XX legs, n = 3/16) within 0.25 of 2 at the same sizes, otherwise every pass or fail is reported as uncalibrated. Passing this is a first stage only. A counterexample claim would further need L = 256, FES or ring spectroscopy, and the full qubit chain.

Predictions.
- M1, |w| = 1, n = 3/16, φ = 1.565, 1.555, 1.54. A partially polarized two-component liquid (PROOF 4): the first-stage criterion passes, with p increasing with π/2 − φ. At φ = 1.50 one valley: p > 0.9 and c_two(W = 6) < 1.3.
- M2, |w| = 1, n = 1/8, φ = 1.565 and 1.55. As M1, the criterion passes.
- M3, |w| = 2, n = 9/32, φ = π/2. Paired: c_two(W = 6) < 1.4, L·Δ_P grows by more than 50 percent from 64 to 128, p < 0.05.
- M4, |w| = 2 off π/2 (n = 9/32 and 3/8). No prediction for c is registered. The alternatives and their signatures are (a) partially polarized c = 2 (criterion passes, p continuous in φ), (b) first-order or phase separated (inhomogeneous profile at fixed N, or p jumping between neighbouring φ), (c) a trimer plateau (p near 1/3 over a range of φ, c_two below 1.4, a finite parity gap), (d) saturation (p > 0.9, c_two below 1.3). Prediction registered only for the far end: φ = 1.48 is (d).
- M5, equation of state at |w| = 2, L = 64. At φ = π/2 the energy per particle E(N)/N has its minimum at a density between 0.15 and 0.35 (a self-bound pair liquid, from the bound tetramer of C175). No prediction at φ = 1.535.
- Control V0 (|w| = 2, φ = 1.535, V₂ = V₃ = 0, n = 9/32). No pass or fail, it calibrates what an imbalanced two-valley liquid without the attractions reads at these sizes.
- M6 (added before any output, job pending). FFLO signature. At a point that passes the first-stage criterion, the pair correlation D(r) (separation 2, c frame, from the centre) is not a single wave at the pair momentum K_p(φ) of fewbody.py. Its Fourier weight sits at K_p ± q with q = π(n_a − n_b), the Luttinger (Yamanaka–Oshikawa–Affleck) value for two valley Fermi seas, and the single-particle |G_f(r)| decays as a power law. Pass if the dominant |D| Fourier peak is displaced from K_p by between 0.5 and 1.5 times π p n, with p from the n(k) weights. At M3's paired point the D peak sits at K_p within 2π/64.

## Pre-registered predictions, ring magnetization curves (written before running ring_mag.py beyond a 20-site timing test at φ = π/2)

Typhon job 17604449 has an estimated start several days away, so the φ dependence is first tested on rings with C178's kring.py (copied unmodified). For each φ, ring_mag.py finds the lowest level in every momentum sector at fixed N and the ground state's valley imbalance M from the fermion n(k) (cos k > 0 minus cos k < 0). Converting one b particle to a changes M by about 2.
- R1. |w| = 2, (L, N) = (28, 8), (32, 9), (36, 10), n ≈ 0.28. At φ = π/2 the even-N ground state is balanced (|M| < 1).
- R2. The test. Staircase (supports a continuous m(δ), the partially polarized route): at each of the two largest sizes there is a φ interval in which the ground state has an intermediate imbalance, 1 < |M| < N − 1.5, and the φ range of intermediate |M| does not shrink from L = 32 to 36. Jump (first order, no partially polarized phase): at both of the two largest sizes the ground state goes from |M| < 1 straight to |M| ≥ N − 1.5. No prediction is registered for which occurs at |w| = 2.
- R3. |w| = 1, (L, N) = (30, 5), (36, 6). Staircase, with the first intermediate |M| appearing by φ = 1.565 and saturation (|M| ≥ N − 1.5) by φ = 1.50.
- R4. Control, |w| = 2, V₂ = V₃ = 0, (32, 9). Staircase.
- Clarification written right after launch, before any output was read: for odd N the minimal imbalance is one unpaired particle, so there "balanced" means |M| < 2 and "intermediate" means 2 < |M| < N − 1.5.

## Registered reading of the local 64 -> 96 previews (written after the φ = 1.555 pair and the L = 64 runs at 1.53 and 1.50 were seen, before the L = 96 run at 1.53 and the controls finished)

The Typhon job has not started, so local open chains at L = 64 and 96 (n = 3/16, χ to 256) serve as a preview. The 1.555 pair is already seen and is reported as unregistered. For φ = 1.53: preview pass if c_two(64 → 96) > 1.5 for W = 6 and 12, truncation below 1e-6 at L = 96, both valley weights above 0.3 at both sizes, and the controls read within 0.15 of 1 (c = 1) and within 0.25 of 2 (c = 2) at the same sizes. Otherwise preview fail or uncalibrated. A preview pass is not the first-stage criterion, which stays on 64 -> 128 with truncation below 1e-7.

Amendment before any output (job still pending, 2026-11-09). The |w| = 2 main density was changed from 1/4 to 9/32, because n = 1/4 is commensurate with period 4 and could lock the total density (PROOF 4.3), which would confuse the test. Predictions M3, M4 and the V0 control now refer to n = 9/32. Nothing else changed.

## Verdict

The partially polarized route is not closed, and at |w| = 1 it is the most natural reading of the off-symmetric c ≈ 2 signal. It is not established. Off φ = π/2 the valley offset acts as a pseudospin field, the pairing term then oscillates, and an imbalanced two-valley liquid is c = 2 wherever it is homogeneous and away from composite plateaus and resonance tongues (the 1D FFLO mechanism). What can remove it is now located. In the dilute limit (n ≲ 0.01 at |w| = 2, a few times 0.01 at |w| = 1) mass imbalance and bound trimers make every pair–atom mixture first order, so there is no c = 2 window there, as C175 said. At |w| = 2 and the densities where the liquid is stable (self-bound near n ≈ 0.35), the off-symmetric open chain is a cluster liquid with short-ranged one-body and pair functions, and the L = 28 and 32 rings prefer trimer composition, so composite locking appears to pre-empt the two-component phase there at the sizes reached (the partial L = 36 ring is less clear). At |w| = 1, n = 3/16 the open chain does part of what the route needs: the valley imbalance grows continuously with the offset (p = 0, 0.12, 0.23 at φ = 1.555, 1.53, 1.50 for L = 64), with no jump to one valley and no parity gap at L = 64. But at the imbalanced point φ = 1.53 the entanglement grows like c ≈ 1 between L = 64 and 96, which is a negative indication. It is not decisive, because these chains hold only one to three excess particles, too few for the second mode's logarithm. The decisive runs (64 → 128 with N ± 1, both |w|, controls) are registered and queued on Typhon: the |w| = 1 core and the controls in backfill job 17604615, the rest in job 17604471. Neither had started at write-up.

## Results

Few-body (PROOF 2) [numerical, fewbody_w2_R100.out, fewbody_w1_R400.out, ad_scat_w2.out]
- F1 passes. E_b relative to ε_a + ε_b is constant to 1.5 percent from π/2 to 1.48 (0.12255 to 0.12432 at |w| = 2, 0.019527 at |w| = 1).
- F2 fails in part. B3 stays positive everywhere, but at |w| = 2 it falls by 42 percent from π/2 to 1.50 (0.0189 to 0.0109, 0.0093 at 1.48), more than the registered 30 percent. At |w| = 1 it is 6.0e-4 to 5.6e-4 (π/2 to 1.56).
- F3 passes. δ = 4ucϑ with c = √(1 − u²/16), slopes 3.873 and 6.928 [proved, PROOF 1.2].
- [derived] The valley-mass imbalance is η = δ/(16 − u²), at most 0.05 in the relevant range, so the offset is a clean pseudospin field (PROOF 1).
- [numerical] Atom–dimer lengths at |w| = 2 are of order the pair size, with a transmission part 2|a_o| ≈ 8 sites, so atoms and pairs exchange at threshold. At |w| = 1 the shallow trimer gives an atom–dimer length of about 95 sites.

Dilute phase diagram (PROOF 3) [derived]
- First species from the vacuum: pairs for δ < E_b − 2B3, aab trimers for E_b − 2B3 < δ < E_b + B3, a atoms above. At |w| = 2 the trimer window is φ in (1.5507, 1.5584), at |w| = 1 it is 5e-4 wide in φ.
- Mixtures of two dilute species phase separate (C164's mass-imbalance coefficient, pair–atom η = 0.38 and 0.45), so all transitions at small n are first order and the dilute c = 2 window is empty. This confirms C175 but confines the obstruction to n ≲ 0.01 (|w| = 2), and at |w| = 1 to densities below the inverse pair and trimer sizes (0.05 and 0.01).

Finite density (PROOF 4) [derived]
- Luttinger description. Two bosons with a general 2 × 2 Luttinger and velocity structure (Σ broken, drift terms allowed). Backscattering carries 2π(ρ_a − ρ_b), cut off by any imbalance. Composite plateaus at commensurate valley ratios (ρ_a = ρ_b pairs, ρ_a = 2ρ_b trimers) and total-density locking at n = 1/q are the only resonant structures. The magnetization curve has plateaus at m = 0, possibly m = n/3 and m = n, with c = 2 between unless a first-order jump intervenes.
- FFLO signature: pair correlations at K_pair ± π(ρ_a − ρ_b), parity gap closing as 1/L.

Ring magnetization curves (PROOF 5.1) [numerical, ringmag_*.out]
- R1 passes. R2's staircase holds formally at L = 28 and 32 (|w| = 2) and R4 holds for the control, but the control shows that the ring test does not discriminate at these sizes. The first conversion sits at δ ≈ 0.3 with and without the attractions. R3 fails as registered: at |w| = 1, L = 36 the first conversion is at φ ≈ 1.535 and the ring is not saturated at 1.50. No ring shows a jump from balanced to fully polarized.
- At |w| = 2, φ ≤ 1.52, the ground state has trimer composition (6, 3) at L = 32, N = 9 and the nearest composition (5, 3) at L = 28, N = 8. At L = 36, N = 10 (partial run, π/2 to 1.555) the first conversion already comes at δ ≈ 0.11, below E_b, so the ring gives no evidence of a pair gap of order E_b at n ≈ 0.28 and the conversion fields are dominated by finite-size shell effects.

Open chains, local previews (PROOF 5.2, 7) [numerical, local/dmrg/; not the registered Typhon test]
- |w| = 2, N/L = 9/32: self-bound liquid with bulk density 0.345 and depleted edges, one-body and pair functions short-ranged (pairs below 1e-3), n(k) a single broad hump, so the M4 alternatives (c) or a cluster liquid, not (a). No calibrated c.
- |w| = 1, n = 3/16: the valley imbalance grows continuously, W_a = 0.503, 0.557, 0.617 at φ = 1.555, 1.53, 1.50 (L = 64) and 0.537, 0.578 at 1.555, 1.53 (L = 96), with a homogeneous bulk, no jump to one valley, boson correlations at 0.02 to 0.08 out to r = 44 in both valleys, and no parity gap at L = 64 (L·Δ_P = −0.22 at φ = 1.555). But at the imbalanced point φ = 1.53 the entanglement grows like c ≈ 1: c_two(64 → 96) = 0.81, 0.99, 1.15 (W = 1, 6, 12), and 1.06 with the L-proportional window that reads 1.00 and 1.79 on the controls. The registered preview rule fails at φ = 1.53. At φ = 1.555 the same estimators read 1.5 to 1.7, but there the two sizes are on different sides of the first conversion. [heuristic] The chains hold only one to three excess particles, so a c = 2 partially polarized liquid would still read near 1 here; the reading is a negative indication at these sizes, not a closure.
- Controls at L = 64 and 96 (n = 3/16): c = 1 control c_two = 1.07, 1.01, 1.06 (W = 1, 6, 12), c_prof 0.93 and 0.95. c = 2 control c_two = 3.46, 1.54, 2.33, c_prof 1.87 and 1.99. The fixed-window two-size estimator is therefore not calibrated for c = 2 at 64 → 96, and by the registered rule the preview is uncalibrated. A window proportional to L (cuts with |l/L − 1/2| ≤ 1/4, relwin.py, chosen after the controls were seen, so not registered) reads 1.00 and 1.79 on the controls.

Qubit bond and vertices (PROOF 6) [derived]
- Off the stratum no vertex is resonant at generic (δ, n). The closest danger is the two-pair vertex O₄, detuned only by q₄ = 2u²ϑ/(8 − u²/2) (0.004 to 0.008 at |w| = 1 near φ = 1.54 to 1.555). It is cut off when q₄ξ₄ ≫ 1, which holds at V/Δ = 40 (factor 10 to 30 with C179's ξ₄) and fails at V/Δ = 8 for |w| = 1. The ΔN = 2 vertices (7.9e-4|t₂| and |t₂|/30) are off resonance except on tongues.

## Best candidate bond (to be tested, not established)

|w| = 1, φ = 1.53 to 1.555, n ≈ 3/16 to 0.22 bulk, at V/Δ = 40. At φ = 1.555, μ_eff ≈ −2.143 (local L = 64 E(N ± 1)), in h = J_x XX + J_y YY + D(XY − YX) + J_z ZZ + (h_z/2)(Z⊗1 + 1⊗Z):

    J_x = 0.500098723, J_y = −0.499901277, D = 0.006249220, J_z = 10, h_z = −19.97589125

(energy scale |t₂| = 0.0125; at φ = 1.53, J_x + J_y = (1/80) cos 1.53 = 5.10e-4 and D = 0.0062449). This is the parity class, off J_x + J_y = 0. A larger V/Δ improves the O₄ protection.

## Typhon backfill job 17604615 (coordinator's request, pending at submission)

The decisive |w| = 1 part of M1 to M6, split off so it can backfill: 16 cores, 32 GB, 11 h 55 min limit, 16 runs at once with one thread each, χ to 384 (so a run whose truncation stays above 1e-7 at L = 128 is reported as unconverged under the registered criterion). Runs: φ = 1.53 and 1.555 at n = 3/16, L = 64 (N = 11, 12, 13) and L = 128 (N = 23, 24, 25), and the c = 1 and c = 2 controls (ctl1, ctl2, n = 3/16) at L = 64 and 128. Configs in typhon_bf/ (task list typhon_bf/tasks.txt, built in the session, same run.py and fibdmrg.py). Scratch /scratch/lustre/jkudlerflam/c181/, outputs out_bf/bf_<tag>.jsonl and logs_bf/. Collect with

    mkdir -p typhon_out_bf && scp -r typhon:/scratch/lustre/jkudlerflam/c181/out_bf typhon:/scratch/lustre/jkudlerflam/c181/logs_bf typhon_out_bf/

then `python3 summarize.py typhon_out_bf/out_bf` (and `relwin.py` with DIR=typhon_out_bf/out_bf). The registered predictions M1 (at φ = 1.555; φ = 1.53 was not in the original M1 list and is read only against the first-stage criterion), M6 and the calibration rule apply unchanged.

## Typhon job 17604471 (pending, estimated start 2026-11-13)

Trimmed in place on scratch while pending (2026-11-09): the w1n316_p1.5550 runs and the ctl1 and ctl2 runs were removed from typhon/tasks.txt because the backfill job holds them. It now has 124 tasks.

One 48-core node, 36 h limit, 124 tasks (typhon/tasks.txt, from make_cfg.py minus the ten runs moved to the backfill job), 24 at a time with 2 threads each, C176's fibdmrg.py and run.py unmodified. Scratch /scratch/lustre/jkudlerflam/c181/. Expected outputs: out/<tag>.jsonl (one line per sweep, measurement on the last sweep of the N runs), out/ckpt/<tag>.pkl, logs/<tag>.log, logs/c181_<jobid>.out. Tags: w2n932_p{1.5708,1.5600,1.5530,1.5450,1.5350,1.5200,1.5000,1.4800}_L{64,128}_N{18,17,19 | 36,35,37}, w2n38_p{1.5708,1.5450,1.5200}_L{64,128}_N{24,23,25 | 48,47,49}, w1n316_p{1.5650,1.5400,1.5000}_L{64,128}_N{12,11,13 | 24,23,25}, w1n8_p{1.5650,1.5500}_L{64,128}_N{8,7,9 | 16,15,17}, ctlV0_w2_p1.5350 (n = 9/32) at L = 64 and 128, eos_w2_p{1.5708,1.5350}_L64_N{4..28 even}. Collect with

    mkdir -p typhon_out && scp -r typhon:/scratch/lustre/jkudlerflam/c181/out typhon:/scratch/lustre/jkudlerflam/c181/logs typhon_out/

then `python3 summarize.py typhon_out/out` gives c_two, c_prof, p, homogeneity, the pair momentum and the parity gaps, to be read against M1 to M6.

## What is left open

1. The registered test itself (Typhon). In particular M1 and M2 at |w| = 1 with N ± 1 at L = 128, and M3 and M4 at |w| = 2.
2. A calibrated c at |w| = 1 off π/2 at L ≥ 128 with truncation below 1e-7, plus FES or ring spectroscopy, and the full qubit chain at the bond above (C177's tools).
3. Whether the trimer plateau at |w| = 2 is a true plateau (c = 1) or a finite-size composition preference, and whether a c = 2 window opens between the pair and trimer plateaus.
4. The FFLO pair momentum shift, which needs L ≥ 128 to resolve (q ≈ 0.08 at φ = 1.53).
5. The size of ξ₄ in the imbalanced liquid, which decides openness in B₂ at moderate V/Δ.
