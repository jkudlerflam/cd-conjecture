This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

# C182 result: cluster bands and two-component pair liquids, one level up

Derivations are in PROOF.md and the log in PROGRESS.md.

## Pre-registered predictions and pass criteria (written before any C182 output existed)

Model: C173 PROOF eq. 3.2 in the c frame, units |t₂| = 1, ε(k) = 2|w| cos(k − φ) − 2 cos 2k, blockade, V₂ = −2, V₃ = −2/3. Valley a is the lower one off φ = π/2.

Q1. Pair band at φ = π/2. The opposite-valley (ab) pair band E₂(K) is symmetric about K = π (valley exchange maps K → 2π − K). A pair two-valley liquid at the symmetric point needs its curvature at K = π to change sign, so that the band has two minima π ± q₀.
- Registered expectation [heuristic]: 1/m_pair stays positive for all 0 < |w| < 4, because the free-pair value (8 − |w|²/2)/2 stays positive and binding has so far only reduced it by a modest factor (3.37 against 3.75 at |w| = 1, 2.25 against 3 at |w| = 2, collaborator's exact_pair.py). No prediction is made for |w| ≥ 3.
- Pass criterion for "two-minimum pair band": a bound pair band whose lowest two local minima in K are separated by at least 0.05 in K and lie within max(1e-3, 0.1 E_b) of each other, with each below the two-body continuum by at least 1e-4, stable between relative-box sizes R and 2R to 1e-6.

Q2. Like-valley (aa, bb) pairs. C175 reports a like-pair resonance near |w| ≈ 3.63 at φ = π/2. Registered expectation [heuristic]: aa and bb pairs are bound only at |w| ≳ 3.6, and wherever they are bound the ab pair is lower, so the aa/bb two-valley pair band is never the lowest two-body state for |w| ≤ 3.5.

Q3. Clusters at φ = π/2. Clusters with unequal valley content (aab, abb, pentamers) have two degenerate band minima by symmetry, balanced ones (ab, aabb) have one. Registered expectation [heuristic]: at every |w| the lowest energy per particle, E_N^min/N for N ≤ 5, sits on an even, balanced cluster (N = 2 at small |w|, N ≥ 4 at |w| ≳ 2), so the dilute onset species is never a two-valley cluster.

Q4. Off φ = π/2. The pair band minimum moves to K ≈ k_a + k_b and the pair stays bound relative to its own (ab) threshold with nearly the same E_b, while two a atoms cost 2ε_a. The onset species switches from ab pairs to a atoms at a valley offset δ* ≈ E_b. Registered expectation [heuristic]: no other cluster becomes the onset species anywhere in φ ∈ [1.3, π/2], |w| ∈ [0.3, 3.5], so the only near-degeneracies of onset species are the pair/atom line δ = δ*(|w|) (C181's route) and, if Q2 fails, the aa/ab line near |w| ≈ 3.6.

Q5. Many-body screen (effective-model DMRG, two sizes). A point is flagged for follow-up if c(L₁ → L₂) > 1.4 with both controls within 0.25 of 1 and 2 at the same sizes, and truncation below 1e-6. Registered expectation [heuristic]: the only flagged points lie within the pair/atom crossover band around δ*, at intermediate density, and none passes the converged stage.

Q6. Cluster ladder and collapse (registered before clusters.py and ladder6.py ran on Typhon). Registered expectation [heuristic]: at φ = π/2 the per-particle energy E_N^min/N − ε_min is lowest at N = 2 for |w| ≤ 1.5, and keeps falling with even N up to N = 6 for |w| ≥ 2, which would mean the dilute gas collapses (first-order onset into a dense liquid) there.

Q7. Typhon screen and stage-2 rule (registered before job submission, implemented in select_stage2.py).
- Screen: open chains at fixed N, L = 48 and 96, n ∈ {1/12, 1/6, 1/4, 1/3, 5/12}, |w| ∈ {0.5, 1, 1.5, 2, 2.5, 3, 3.5}, φ at δ/E_b ∈ {0, 0.5, 1, 2} (E_b the φ = π/2 pair binding) and φ = 1.45, 1.30; χ to 192 (L = 48) and 256 (L = 96); N ± 1 at L = 96 for the parity gap. Controls c1 and c2 (C176 definitions) at n = 1/4 and 1/6.
- Statistic: c_two(48 → 96) from the W = 6 window mean of S.
- Stage 2: among points with L = 96 truncation below 1e-6, the up to six highest with c_two > 1.4, filled up to six with 1.2 < c_two ≤ 1.4, are run at L = 96, 144, 192 with χ to 512, plus N ± 1 at L = 144 and 192, and both controls at n = 1/4.
- Amendment (added while the job was still pending, before any of its output existed, after the local runs of Q10 to Q12 showed droplets): stage 2 only considers screen points whose L = 96 mid-chain density is within 0.02 of N/L.
- Reading: a stage-2 point is a candidate only if c_two(96 → 192, W = 6) > 1.2 and not below c_two(48 → 96) by more than 0.3, truncation below 1e-7 at L = 192, controls within 0.15 (c1) and 0.25 (c2) at the same sizes, and L·Δ_P not growing by more than 50 percent from 96 to 192. A point failing any of these is reported as not a candidate.

Q8 (registered after submission of job 17604443, before any of its output or any local many-body output existed). Expected screen pattern [heuristic]:
- |w| = 0.5, δ ≤ 2E_b (E_b = 0.0017, pair size of order 60 sites): c_two(48 → 96) above 1.4 at n ≥ 1/6, as a crossover, because the pair coherence length exceeds L. If stage 2 selects these points, c_two(96 → 192) is expected to stay high too, and the parity gap L·Δ_P to stay flat at these sizes. Such a reading would not separate c = 2 from a pair gap below 1e-3, so it is not counted as a candidate unless it also passes the Q7 reading and an independent pair-gap estimate puts the gap above the finite-size spacing.
- |w| ≥ 2 at φ = π/2: low c_two (below 1.2) and inhomogeneous or self-bound behaviour at n ≤ 1/6 (collapse, §5 of PROOF), c ≈ 1 at n ≥ 1/4.
- φ = 1.45 and 1.30 (one valley or a large offset): c_two below 1.3 everywhere.

Q9 (registered before local batch b2 ran; written after the first local runs showed self-bound droplets at |w| = 2 and 3, local/mini). Local open chains L = 40 and 80, χ to 128, at densities above the droplet density.
- |w| = 3 (φ = π/2 and δ = E_b) and |w| = 2 (φ = π/2 and δ = E_b) at n = 0.3, and |w| = 3, π/2 at n = 0.4: c(40 → 80) below 1.3 (a paired or clustered liquid, c = 1).
- |w| = 0.5, π/2, n = 0.2: above 1.4 (pairing crossover, as in Q8). |w| = 0.5, φ = 1.45 and |w| = 1, δ = E_b, n = 0.2: below 1.3.
- Controls c1 and c2 at n = 0.2 within 0.25 of 1 and 2 at these sizes (otherwise reported uncalibrated).

Q10 (registered before local batch b3 ran, after local/mini showed c(32 → 64) = 2.5 to 2.6 at |w| = 2, δ = E_b, N/L = 1/6, inside a droplet of density 0.23 with valley ratio 0.77). Follow-up at |w| = 2, δ = E_b, uniform density n = 0.25 (above the droplet density): L = 40, 80, 120 with χ to 192, the parity gap at L = 80 and 120, and the neighbours δ = 0.75 E_b and 1.25 E_b at L = 40, 80.
- No value of c is predicted. Reading rule: the point is kept as a lead only if c(40 → 80) and c(80 → 120) (W = 1 and W = L₁/8) both exceed 1.4 at truncation below 1e-6, the profile is homogeneous (n_mid within 0.02 of N/L), and L·Δ_P does not grow by more than 50 percent from 80 to 120. Otherwise it is reported as a finite-size or droplet effect.

Q11 (registered before local batch b4 ran, after b2 and b3 showed that N/L = 0.25 to 0.3 at |w| = 2 is below the self-bound density, n_mid ≈ 0.35). Uniform liquids: |w| = 2 at φ = π/2 and δ = E_b, n = 0.4; |w| = 1 at φ = π/2 and δ = E_b and |w| = 0.5 at π/2, n = 0.35; L = 40 and 80. Prediction for |w| = 2: c(40 → 80) below 1.3 with n_mid within 0.02 of N/L. No prediction for |w| ≤ 1. Readings at n_mid more than 0.02 above N/L are reported as inhomogeneous and not used for c.

Q12 (registered before local/eos ran, after b4 showed a dense region, n_mid ≈ 0.48, at |w| = 0.5, N/L = 0.35). Equation of state on open chains, L = 40, φ = π/2, N = 1 to 20, |w| = 0.5, 1, 2, 3, χ to 96: μ(N) = E(N) − E(N − 1). Prediction: μ(N) decreases over a range of low N at |w| = 2 and 3 (self-bound liquid, first-order onset), and is increasing up to n ≈ 0.3 at |w| = 1 (as C176's L = 48 scan). No prediction at |w| = 0.5. A decreasing stretch of μ(N) (beyond the even-odd alternation, so compared on even N) is read as phase coexistence.

Criterion for a c > 1 candidate (adapted from C163 PROOF §4.2, as the C181 to C185 spec states it): converged open-chain fits above 1.2 that do not fall with L, at truncation below 1e-7; FES or ring-spectroscopy confirmation; a check on the full qubit chain; controls at c = 1 and c = 2 at the same sizes; the relevant gap (parity gap, pseudospin gap) measured against L.

## Results

### Verdict

No c > 1 opportunity was found on the cluster-band side, and none of the local many-body diagnostics produced a lead that survives the registered homogeneity and size conditions. One level up the two-valley question does arise, but only for odd clusters at φ = π/2, and those are never the onset species. Their unlike channel is also reactive (two unlike trimers rearrange into three pairs). The pair band itself has a single minimum at K = π for every 0 < |w| < 4, with a positive exact inverse mass. Pair–atom coexistence at δ = E_b is pre-empted by bound trimers, and at |w| ≳ 2 by ever larger clusters, so the dilute gas is self-bound there. The uniform liquids above the coexistence region read c ≈ 0.85 to 0.95 at |w| = 2 and 3, both at φ = π/2 and at δ = E_b. The only region with high readings is the one already under test (|w| ≲ 1 near φ = π/2 at intermediate density, the BCS-crossover regime of C173, C175, C176), which this run extends to n ≈ 0.35 at |w| = 1 without being able to separate c = 2 from a small pair gap. Best candidate bond: none. The Typhon screen (job 17604443, below) is the systematic phase-diagram check and was still queued when this was written.

### Key claims

1. [proved, §1 of PROOF] At φ = π/2 cluster bands satisfy E_N(K) = E_N(Nπ − K). Every odd cluster band (trimers, pentamers) therefore has two degenerate minima, and balanced even clusters sit at the symmetric momentum.
2. [proved algebra, numerical evaluation; PROOF 2.2, pairmass_exact.out] The ab pair inverse mass at K = π, from the collaborator's exact wavefunction, is positive and monotone on (0, 4), from 4 down to 0.643 (3.374 at |w| = 1, 2.251 at 2, 1.296 at 3, 0.931 at 3.5, 0.695 at 3.9). The banded numerics agree (0.93138, 0.69469 at 3.5, 3.9). [numerical, sym_pairs.out] The full K scan finds a single bound minimum at K = π for all 14 values of |w| from 0.3 to 3.9. So the pair band never has two minima at the symmetric point. Q1 passed as registered.
3. [numerical, sym_pairs.out] Like-valley (aa, bb) pairs are unbound for |w| ≤ 3.6 and bound at 3.7 (by 0.0022), 3.8 (0.0125) and 3.9 (0.0323). They then form a two-minimum pair band, but it lies above the ab pair by 0.13, 0.10 and 0.06. Q2 passed.
4. [proved, elementary; PROOF 3] Pre-emption lemma: two clusters with equal per-particle energy that bind to each other are pre-empted by their bound state. A dilute two-species liquid at the onset therefore needs two clusters on the lower envelope that do not bind. [derived] With the valley offset acting as a field, e(m, n) = [nδ − B(m, n)]/(m + n), and pure a atoms win only for δ > sup B(m, n)/n.
5. [derived; PROOF 3 to 5] Pair–atom coexistence at δ = E_b is pre-empted by the aab trimer whenever B₃ > 0 (true at every |w| checked), which opens a trimer window E_b − 2B₃ < δ < E_b + B₃. At |w| = 3 the ladder keeps falling through N = 5 at δ = 0.75 E_b.
6. [derived, PROOF 4] Trimers are the onset species at φ = π/2 only if B₃ > E_b/2. Under the same inequality reversed, two unlike trimers (total momentum π) can rearrange into three ab pairs with energy release E_b − 2B₃ > 0. [numerical] B₃ = 0.0006 (|w| = 1, C175), 0.019 (|w| = 2, C175, C181), 0.068 (|w| = 3, R = 40 and 60 equal) against E_b/2 = 0.0098, 0.061, 0.098. Q3 passed for N ≤ 3. At |w| ≥ 2 the onset species is not the pair but larger balanced clusters (next claim), which is also consistent with Q3's "even, balanced".
7. [numerical, upper bounds; PROOF 5, local/ladder_w3_summary.out] Per-particle energies at |w| = 3, φ = π/2: −0.0975, −0.0877, −0.1107, −0.1087, −0.1194 for N = 2 to 6 (the hexamer at R = 14 is already below the tetramer). Q6 is borne out at |w| = 3. At |w| = 2 the local hexamer box (R = 14) is too small, and the Typhon ladder has R = 16 and 18.
8. [numerical, local/eos/eos.out; PROOF 6.1] Equation of state at φ = π/2, L = 40. At |w| = 1, μ rises monotonically up to n = 0.5, so the C173/C176 candidate region is a uniform phase (Q12 passed). At |w| = 3 the energy per particle falls through N = 8 and μ falls, a self-bound liquid (Q12 passed). At |w| = 2 clusters of 4 to 6 particles are self-bound (droplets of mid-chain density about 0.35). At |w| = 0.5, μ is flat to 0.003 for 0.2 ≲ n ≲ 0.3 and the N/L = 0.35 ground state is a droplet of density 0.47, so small |w| has a coexistence region between a dilute pair gas and a dense phase. [heuristic] The legs are the ferromagnetic Heisenberg chain at w = 0, the edge of phase separation, which is why small attractions tip the model into clustering.
9. [numerical, local/b2, local/b4; PROOF 6.2] Homogeneous liquids above the coexistence region read one mode: c(40 → 80) = 0.85 to 0.96 (W = 1, 5) at |w| = 3 (n = 0.3, 0.4; φ = π/2 and δ = E_b) and |w| = 2 (n = 0.4; φ = π/2 and δ = E_b), with truncation below 1e-10 and balanced valleys. The controls read 1.10/1.01 (c1) and 2.17/1.80 (c2). Q9 and Q11 passed where they could be tested.
10. [numerical; PROOF 6.3] High two-size readings (2.4 to 3.9) at |w| ≤ 1 and inside the |w| = 2 droplet regime come with mid-chain densities 0.03 to 0.12 above N/L that drift with L. The |w| = 2, δ = E_b lead from local/mini fails Q10 (n_mid = 0.328, 0.313, 0.297 at L = 40, 80, 120; the parity gap is negative at L = 80 and positive at 120, as for droplet rearrangements). Q9's predictions at |w| = 0.5 and at |w| = 1, δ = E_b could not be scored, because those chains were inhomogeneous.

### Against the registration

- Q1, Q2: passed. Q3: passed for N ≤ 3 at φ = π/2, and the onset species at |w| ≥ 2 is a larger balanced cluster. Q4: the pair/atom line is pre-empted by trimers (claim 5), so "the only near-degeneracy is the pair/atom line" is replaced by the pair/trimer and trimer/atom lines, both subject to the lemma. The complete off-symmetric ladder is in the Typhon job. Q5 to Q8: see "Typhon results" at the end of this file. Q9: passed at |w| = 3, inhomogeneous elsewhere (not scored). Q10: not a lead by the registered rule. Q11: passed at |w| = 2; |w| ≤ 1 inhomogeneous. Q12: passed at |w| = 1 and 3. The |w| = 2 μ₂ is nearly flat rather than decreasing, and at |w| = 0.5 a coexistence region appeared that was not predicted.

### Typhon job

Job 17604443 (partition typhon, one 48-core node, 320 GB, 3 days), submitted 2026-11-09, completed 2026-11-10 02:30 EST (results in the last section of this file). Scratch /scratch/lustre/jkudlerflam/c182/.
- Stage 1: 48 cluster-ladder tasks (out/cl_w{1,1.5,2,2.5,3,3.5}_d{0,0.25,0.5,0.75,1,1.25,1.5,2}.json), 8 hexamer tasks (out/l6_w{2.0,2.5,3.0,3.5}_R{16,18}.json), and the DMRG screen (out/s_*.jsonl and controls out/c1_*.jsonl, out/c2_*.jsonl, 848 runs), logs in logs/.
- Stage 2: typhon/stage2_ranking.txt, typhon/tasks2.txt, out2/*.jsonl, summary1.txt, summary2.txt.
- Collect with

      mkdir -p typhon_out && scp -r typhon:/scratch/lustre/jkudlerflam/c182/out typhon:/scratch/lustre/jkudlerflam/c182/out2 typhon:/scratch/lustre/jkudlerflam/c182/typhon typhon:/scratch/lustre/jkudlerflam/c182/logs typhon:/scratch/lustre/jkudlerflam/c182/summary1.txt typhon:/scratch/lustre/jkudlerflam/c182/summary2.txt typhon_out/

  then run `python3 summarize.py clusters typhon_out` and `python3 summarize.py screen typhon_out` in this folder. The checkpoints in out/ckpt are large and can be left out by copying out/*.jsonl and out/*.json instead of out/.

### What is left open

1. The Typhon screen and the stage-2 converged runs. These complete the cluster ladder over |w| = 1 to 3.5 and δ/E_b = 0 to 2, decide whether the off-symmetric trimer window ends in a non-binding pair of species anywhere, and map c across (|w|, φ, n).
2. The near-symmetric |w| ≲ 1 region at intermediate density (n up to 0.35 at |w| = 1). It is uniform, but at the sizes reached here a pair gap cannot be told from c = 2 there. That is the subject of C176, C178 and C181.
3. The small-|w| coexistence region: its boundaries in (|w|, n) and whether the dense phase near n ≈ 0.47 is a commensurate (gapped) state.
4. The full qubit chain was not run here. With no candidate, there is nothing to transfer.

### Files

pairband.py, nbody.py, clusters.py, ladder6.py, pairmass_exact.py, sym_pairs.py, make_tasks.py, select_stage2.py, summarize.py, mini_an.py, eos_an.py, typhon/; C176's fibdmrg.py, run.py and analyze.py were copied unmodified. Outputs: sym_pairs.out, pairmass_exact.out, local/ (cl_w3_*, l6_small.jsonl, ladder_w3_summary.out, mini/, b2/, b3/, b4/, eos/, each with summary.out or eos.out).

## Typhon results (collected 2026-11-10 by the coordinator, Claude Fable 5)

The job ran both stages to completion. Outputs are in typhon_out/, with the tables typhon_out/clusters_summary.out (`python3 summarize.py clusters typhon_out`), typhon_out/screen_summary.out (`python3 summarize.py screen typhon_out`) and typhon_out/stage2_reading.out (`python3 -E -s stage2_read.py typhon_out`). The verdict above changes in one respect. The DMRG screen has four points that pass the registered stage-2 reading, so Q5 is not borne out as written. None of them is a new mechanism, and none is a converged candidate (see the assessment at the end).

A caveat on the job's own summary2.txt. It was made with C176's analyze.py, which files the N + 1 runs under the same series and size as the main runs, so its L = 144 and 192 rows are the N + 1 states. It is not a valid stage-2 reading. stage2_read.py selects the main-N runs by tag and applies Q7 as registered.

### Cluster ladder (Q3, Q4, Q6)

- [numerical] At φ = π/2 the onset species among N ≤ 5 is the pair at |w| = 1 (e₂ = −0.00976, e₃ = −0.00649, e₄ = −0.00879) and 1.5 (e₂ = −0.03179, e₄ = −0.03076), the tetramer at |w| = 2 (e₂ = −0.06127, e₄ = −0.06345, e₅ = −0.05691), 2.5 (e₄ = −0.09550) and 3 (e₄ = −0.11073), and the pentamer at 3.5 (e₅ = −0.09520). Here e_N = E_N/N − ε_min.
- [numerical] The hexamer per particle at R = 16 and 18 is −0.06228, −0.06334 at |w| = 2, −0.10155, −0.10178 at 2.5, −0.11954, −0.11959 at 3, and −0.10128, −0.10130 at 3.5. It lies below the tetramer and pentamer for |w| ≥ 2.5. At |w| = 2 it is level with the tetramer (−0.06345) and still moving with R. Q6 is borne out for |w| ≤ 1.5 and |w| ≥ 2.5 and is marginal at |w| = 2.
- [numerical] Off the symmetric line the onset species against δ/E_b = 0, 0.25, 0.5, 0.75, 1 is N = 2, 2, 2, 2, 1 at |w| = 1, then 2, 2, 2, 2, 3 at 1.5, then 4, 4, 4, 5, 5 at 2, then 4, 4, 5, 5, 5 at 2.5, then 4, 5, 5, 5, 5 at 3 (and still 5 at 1.25), and 5 throughout at 3.5. For larger offsets the onset is the single atom. Two rows (|w| = 1.5 at δ = 2E_b and 2.5 at 1.5E_b) print N = 2 with e₂ = 0 to five digits, which is a tie with the atom threshold. So the trimer is the onset species at one point only (|w| = 1.5, δ = E_b, inside the trimer window of claim 5). For |w| ≥ 2 the lowest cluster among N ≤ 5 is the tetramer or the pentamer up to δ = E_b, and these tasks (N ≤ 5 only) do not decide where that ladder stops.
- Flag. Three rows at |w| = 3.5 (δ/E_b = 1.25, 1.5 and 2) are identical, with φ = 1.54226 and δ printed as 0. The bisection in clusters.phi_for_delta appears to have converged to the φ where the second valley minimum disappears, so offsets that large do not exist at |w| = 3.5 and these three rows are one point at that edge.

### DMRG screen (Q5, Q8)

- [numerical] 210 screen points. 101 are homogeneous at L = 96 (|n_mid − N/L| ≤ 0.02). 98 have c_two(48 → 96, W = 6) > 1.4 with truncation below 1e-6, and 85 of those are inhomogeneous (self-bound or phase-separated, n_mid far from N/L). The homogeneous points fall as 7 below 0.8, 81 between 0.8 and 1.2, none between 1.2 and 1.4, and 13 above 1.4.
- [numerical] The 13 homogeneous flagged points are s_w2.5_p1.45_n1o6 (7.031, with the c = 1 window reading 7.220, an artefact), s_w1_x2_n1o12 (1.723), s_w1.5_p1.45_n1o6 (1.576), s_w1_p1.45_n1o6 (1.574), s_w1.5_p1.3_n1o4 (1.550), s_w0.5_{x0, x0.5, x1, x2}_n1o12 (1.531 each), s_w1_{x0, x0.5, x1}_n1o12 (1.486, 1.485, 1.485) and s_w3_p1.45_n1o4 (1.444, with c_prof = 0.33).
- [numerical] Screen controls, c_two(48 → 96) at W = 1 and 6: c1 reads 1.047, 0.979 at n = 1/4 and 1.117, 1.016 at n = 1/6, and c2 reads 2.121, 1.913 and 2.218, 1.808.
- Q8, second bullet (|w| ≥ 2 at φ = π/2). Every density up to n = 1/4 or 1/3 is inhomogeneous (n_mid 0.15 to 0.35 against N/L from 1/12 to 1/3), so "self-bound" holds and the c_two values there (−2.5 to 6.5) carry no meaning. At n = 5/12 the reading is 0.94 at every |w|.
- Q8, third bullet (φ = 1.45 and 1.30, c_two below 1.3 everywhere) is wrong. Five homogeneous off-line points read above 1.3 at the screen (1.55, 1.58, 1.57, 7.03 and 1.44, listed above).

### Stage 2 (Q7), from the main-N runs

All stage-2 runs finished their schedules (23 sweeps, χ = 512). Controls at n = 1/4, c_two(96 → 192) at W = 1, 6, 12: c1 = 1.037, 0.979, 1.017 and c2 = 2.056, 1.923, 1.978 (truncation 8.7e-9), so the size pair is calibrated. n_mid − N/L is about +0.009 at L = 192 for every point.

| point | c_two(48 → 96) | c_two(96 → 144) | c_two(144 → 192) | c_two(96 → 192) | c_prof(192) | W_min/W_max | L·Δ_P at 96, 144, 192 | trunc(192) | Q7 |
|---|---|---|---|---|---|---|---|---|---|
| s_w0.5_x0_n1o12 | 1.531 | 1.625 | 1.712 | 1.661 | 2.253 | 1.000 | 1.514, 1.413, 1.356 | 3.6e-9 | passes |
| s_w1.5_p1.3_n1o4 | 1.550 | 2.289 | 2.591 | 2.414 | 3.719 | 0.169 | −1.013, 0.602, 0.571 | 1.2e-8 | passes |
| s_w1.5_p1.45_n1o6 | 1.576 | 2.297 | 0.976 | 1.748 | 2.159 | 0.203 | 0.893, 1.061, 1.279 | 1.7e-8 | passes |
| s_w1_p1.45_n1o6 | 1.574 | 1.101 | 2.221 | 1.566 | 2.188 | 0.306 | 1.481, 0.439, 0.089 | 3.0e-8 | passes |
| s_w1_x2_n1o12 | 1.723 | 1.765 | 1.721 | 1.747 | 2.159 | 0.797 | −0.814, 1.043, 1.252 | 4.4e-9 | fails (gap) |
| s_w2.5_p1.45_n1o6 | 7.031 | 1.076 | 0.975 | 1.034 | 7.031 | 0.079 | −0.596, 0.914, 0.918 | 8.1e-11 | fails |

All c_two entries are at W = 6. The step |S̄₆(χ = 512) − S̄₆(χ = 384)| at L = 192 is 0.0010, 0.0115, 0.0051, 0.0055, 0.0011 and 0.0000 in the order of the table. For s_w1.5_p1.3_n1o4 the gap changes sign between L = 96 and 144, and the rule was applied to its absolute value.

### Against the registration, and assessment

- Q5 ("none passes the converged stage") is not borne out. Four points pass the Q7 reading as registered. Q7's own rule therefore names them candidates of this screen.
- s_w0.5_x0_n1o12 (|w| = 0.5, the symmetric line, n = 1/12) is the crossover that Q8 anticipated. It passes Q7, but Q8 also requires an independent pair-gap estimate above the finite-size spacing, and none exists, so by the registered rule it is not counted. Its three size pairs agree (1.63, 1.71, 1.66), its valleys are balanced and L·Δ_P falls slowly, which is what a pair gap far below the level spacing looks like at these sizes. It is the same question as the main candidate at |w| = 1, at weaker binding.
- The three off-line points (φ = 1.45 and 1.30, |w| = 1 and 1.5) are valley-imbalanced (W_min/W_max = 0.17 to 0.31), their two size pairs disagree by more than 1 (0.98 against 2.30, 2.22 against 1.10), and c_prof(192) is 2.2 to 3.7. These are not converged readings of a central charge. They sit in the partially polarized two-valley regime that C181 and C187 already study closer to the line (REPORT §9.156, §9.157, §9.173), so they are further leads in that family and not a new class.
- None of this changes the verdict on the cluster-band side. The |w| ≥ 2 ladder beyond N = 6 is still open.

### Left open after the Typhon job

1. An independent pair-gap estimate at |w| = 0.5, n = 1/12 (the Q8 condition), and L = 256 or 384 for that point.
2. The off-line points at larger sizes with N ± 1 and N ± 2 partners, to see whether the erratic size dependence is a shell effect of the imbalanced valleys.
3. The ladder beyond N = 6 at |w| = 2, where the hexamer and tetramer are level.

### Files added

stage2_read.py; typhon_out/ (summary1.txt, summary2.txt, stage2_ranking.txt, tasks2.txt, clusters_summary.out, screen_summary.out, stage2_reading.out, out/ with the cluster and hexamer JSON files, the controls and the six stage-2 screen series, out2/ with the stage-2 runs). The 840 raw first-stage screen files (121 MB) are kept out of git by typhon_out/out/.gitignore. Their numbers are in summary1.txt and screen_summary.out, and the files remain on Typhon under /scratch/lustre/jkudlerflam/c182/out/.
