This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

# C185 result: degenerate points where second-order dynamics can win

Derivations are in PROOF.md and the log in PROGRESS.md.

Scripts: classify.py, fib_first.py, sw_dm.py (imports C179's sw_exact.py), pseudo_v.py and trimer_v.py (import C175's scat.py and nbody_box.py), ring_kappa.py, degen_opt.py with degen_test.py, eff/blockk.py (extends C173's block_dmrg.py), eff/ copies of C176's fibdmrg.py, run.py, analyze.py, full/ copy of C177's fc.py (with a --dmy option) and its dependencies. Outputs: classify.out, fib_first.out, sw_dm.out, pseudo_v.out (pseudo_v_partial.out is an aborted first run), trimer_v.out, ring_kappa_L20.out, ring_kappa_L24.out, eff/pilot.log, eff/pilot2.log, and the validation files full/test_fc.out, eff/test_eff.out, eff/test_blockk.out.

## Verdict

No closure and no confirmed candidate, but a new route. The list of points where second-order dynamics can win is short and complete at the physics level, and one first-order-invisible direction at the Fibonacci point produces exactly the repulsive inter-valley channel the spec asked for. Its price is a ΔN = 1 term of the same size, so whether c = 2 survives is a many-body question that is set up on Typhon and not answered here.

- Classification [proved for classical and frustration-free bonds, C143 for local integrals of motion, unbiased numerical search pending]. A nearest-neighbour qubit bond has an exponentially degenerate ground space only on the Fibonacci orbit (J ZZ + b Z at |b| = 2J, J > 0, the rank-one projector onto |aa⟩) and the XZ orbit. Every other frustration-free kernel gives L + 1 or at most 2 states, every other classical bond at most 2. The polynomial points (ferromagnetic multiplet, the Bravyi–Gosset set, product projectors |ab⟩, KTM) cannot host the mechanism, because their degenerate space has zero entropy density and the low-energy extensive physics is H₀'s own single-minimum magnon band (PROOF 2.5).
- XZ [C144, C157]. Closed. The second-order B-chain interaction is a positive-semidefinite kernel, attractive in every channel.
- Fibonacci, complete tuned locus [proved]. Modulo rotations and scale, the first-order-invisible normal space is four-dimensional: the pairing Δ ∈ ℂ and the transverse DM vector ℓ = λ_y + iλ_x, with bond terms λ_y(XZ − ZX) + λ_x(YZ − ZY). In the parity class only Δ survives, so the second-order model is rigid (C173's eq. 3.2) and the pair is bound for every bond in that class (PROOF 3.2).
- The transverse-DM cusp [derived, checked against the full chain]. λ(XZ − ZX) adds V₂′ = +8λ²/V (repulsive), a real NN hop 4λ²/V and nothing at distance 3. In units |t₂| the model is V₂ = −2 + v with v = 16λ²/Δ², plus a ΔN = 1 conditioned creation κ Σ σ⁻_s(n_{s+2} − n_{s−2}) + h.c. with κ = 4λ/Δ = √v. The uniform single-particle creation cancels exactly.
- Two-body and three-body problems of its U(1) part [numerical with the proved boundary]. For |λ|/Δ > u/8 the pair is unbound, the dilute pseudospin is easy-plane (Δ_σ between 0.27 and 0.95 for every v > u²/4 tested at u = 0.5 to 2), and there is no trimer. This is the first setting in the Fibonacci corner where the dilute two-valley gas is two-component by C175's criterion.
- The catch [derived, numerical]. κ = √v > u/2 wherever the pair is unbound, so the ΔN = 1 term is as large as the hops. It breaks U(1) and the valley symmetry. Ring ED at κ = 1 shows the valley-balance point moved to φ* ≈ 1.64 and a ground state with n ≈ 0.2 and strong number dressing. Local open chains at L = 36 (χ ≤ 64) are too small to read c: the C173 calibration point itself gives a profile fit of 7.6 there. The test is in the Typhon job.

## Best candidate bond (not established)

V/Δ = 40, u = |w|/|t₂| = 1, λ = 0.3 (v = 1.44, κ = 1.2), in h = J_x XX + J_y YY + D(XY − YX) + J_z ZZ + (h_z/2)(Z⊗1 + 1⊗Z) + λ(X⊗Z − Z⊗X):
J_z = 10, J_x − J_y = 1, J_x + J_y = −4λ²/V − 2D₀ sin δφ = −0.009 − 0.0125 sin δφ, D = D₀ cos δφ with D₀ = 0.00625, h_z = −19.97625 − 4λ²/V = −19.98525, λ = 0.3, and the valley offset tuned by δφ ≈ 0.07 to 0.1 (ring ED at κ = 1 gives 0.065 to 0.075). Energy scale |t₂| = 0.0125. It is off the parity class and off J_x + J_y = 0. At V/Δ = 8 the corresponding full-chain points are in the Typhon job (C177's tuned D and h_z plus the same shifts).

## Key claims

1. [proved] Classical bonds: exponential ground degeneracy only at J > 0, |b| = 2J (PROOF 2.1, classify.out).
2. [proved for |aa⟩, numerical otherwise] Frustration-free bonds: exponential only for the rank-one kernel complement |aa⟩; every other rank-one kernel gives exactly L + 1, ranks 2 and 3 give at most 2 (PROOF 2.2, classify.out, to L = 14 with 40 random kernels per class).
3. [derived] Polynomially degenerate points cannot host second-order-dominated extensive physics (PROOF 2.5).
4. [proved, fib_first.out] At the Fibonacci point every bond term projects onto w, Ω, μ; the invisible directions are Δ and the transverse DM ℓ (PROOF 3.1).
5. [derived] Parity-class rigidity: the second-order model on the parity-class tuned locus is uniquely C173's, so the pair is bound there for every bond (PROOF 3.2, with C175 and the contrib theorem).
6. [derived, numerical sw_dm.out to 1e-4] Transverse DM at second order: V₂′ = 2(λ_L − λ_R)²/V in general, +8λ²/V for the DM, hop 4λ²/V, V₃′ = 0, uniform ΔN = 1 creation zero, conditioned creation 2Δλ/V (PROOF 3.3, 3.4).
7. [numerical] Above v = u²/4 the U(1) part is unbound, two-component in the dilute limit, and has no trimer (pseudo_v.out, trimer_v.out).
8. [derived] v = κ²: the repulsion and the ΔN = 1 term come together (PROOF 3.4).
9. [numerical, small sizes] κ = 1 moves the valley balance to φ* ≈ 1.64 and raises the density at fixed μ from 0.08 to about 0.2 (ring_kappa_L20.out, ring_kappa_L24.out). L = 36 open chains (eff/bk_pilot.jsonl, eff/bk_pilot2.jsonl) give profile fits of 0.92 at φ = π/2 and 1.63 at φ = 1.64 for κ = 1, but 7.6 at the C173 calibration point (S_mid = 1.479), so these fits carry no information at this size. They are recorded, not used.
10. [C144, C157] XZ is closed by the sign of a positive kernel.

## What is left open

1. Whether the full model (3.5) at κ ≈ 1 and φ ≈ φ* is a two-component liquid. Tests: the blockk runs of jobs 17604568 and 17604455 (L = 48, 96), then a blockade-exact non-U(1) DMRG at L ≥ 128 with n(k), the valley polarization and FES; the parity gap is not defined without U(1), so the pseudospin gap has to be read from the valley-imbalance correlator.
2. The full qubit chain with the DM term at V/Δ = 8 to 16 (jobs 17604568 and 17604455, φ scan at L = 64 only), and then C163's criterion.
3. Whether κ can be reduced at fixed v. At second order v = κ² is forced by the single source Δ for the NNN hop. A third-order source, or a second dominant bond, is not available in B₂ at this order.
4. The unbiased degeneracy search (R3) on Typhon, which would close the completeness of the extensive-point list numerically outside the classes covered by proofs.
5. Overlap: the transverse DM is a term C186 (general bond near the Fibonacci point) should also find. The route belongs equally to C186's question, and C181's partially polarized route meets it at φ ≠ π/2.

## Typhon jobs

Two jobs, split at the coordinator's request so that the highest-value runs can backfill. They share the scratch directory /scratch/lustre/jkudlerflam/c185/ and its code (eff/, full/, typhon/), and write to separate output directories. No run appears in both.

Backfill job 17604568 (partition typhon, 12 cores, 12 h limit, submitted 2026-11-09, PENDING at submission). Script typhon/c185bf.slurm, tasks typhon/bf_block_tasks.txt and typhon/bf_full_tasks.txt.
- blockk (bf/outbk/bk_*.jsonl): κ = 1, v = 1, μ = −2.1 at φ = 1.625, 1.640, 1.655, and the κ = 0, v = 1 reference at φ = π/2, each at L = 48 and 96 (8 runs).
- full chain, V/Δ = 8, L = 64 (bf/oc/*.json): λ = 0.3 at φ shift 0.06 (F8dm_p06_64.json) and the λ = 0 reference (F8ref_64.json).
- logs in bf/logs/.

Collect from inside codex/work_C185 with

    mkdir -p typhon_out_bf && scp -r typhon:/scratch/lustre/jkudlerflam/c185/bf/outbk typhon:/scratch/lustre/jkudlerflam/c185/bf/oc typhon:/scratch/lustre/jkudlerflam/c185/bf/logs typhon_out_bf/

Main job 17604455 (partition typhon, one 48-core node, 2-day limit, submitted 2026-11-09, PENDING, estimated start shown as 2026-11-13 earlier and 2026-11-09T19:40 at the last check). After the split it holds:
- 10 fibdmrg runs of the U(1) part (out/*.jsonl): dmv1 at L = 64, 128 with N ± 1, dmv2 and mgs at L = 64, 128;
- 4 blockk runs (outbk/bk_*.jsonl): κ = 1 at φ = π/2 and the C173 calibration point κ = 0, v = 0, each at L = 48, 96;
- 4 full-chain runs (oc/*.json): λ = 0.3 at φ shift 0, 0.10, 0.14 (L = 64), and the λ = 0 reference at L = 128;
- 22 degeneracy-search tasks (degen/*.jsonl), logs in logs/, and summary_eff.txt at the end.

Collect from inside codex/work_C185 with

    mkdir -p typhon_out && scp -r typhon:/scratch/lustre/jkudlerflam/c185/out typhon:/scratch/lustre/jkudlerflam/c185/outbk typhon:/scratch/lustre/jkudlerflam/c185/oc typhon:/scratch/lustre/jkudlerflam/c185/degen typhon:/scratch/lustre/jkudlerflam/c185/logs typhon:/scratch/lustre/jkudlerflam/c185/summary_eff.txt typhon_out/

out/ckpt holds large checkpoints and can be deleted from the copy. R6′ and R7 are scored on the union of both jobs.

## Report against the registration

- R1: passed (classify.out). R2: passed, with rank 2 giving 2 (structured) or 0 (random), and rank 3 at most 1.
- R3: pending (Typhon job 17604455). Only the calibration (degen_test.py) has run: W = 0 at Fibonacci, XZ and the ferromagnetic L + 1 count.
- R4: passed. Every predicted coefficient agrees to 0.5 percent at V/λ = 40 (Δ = 0) and to 1e-4 at V = 320, with deviations falling as 1/V². The registration did not predict the ΔN = 1 term, which the extraction found (2Δλ/V conditioned creation, uniform part zero).
- R5: nothing run, as registered.
- R6, R6′, R7, R8: pending (Typhon job). The local pilots are too small to score them.

## Pre-registered predictions and pass criteria (written before any computation, 2026-11-09)

R1 (classical bonds). Among classical bonds of B_2 (diagonal in one uniform product basis, so J_zz ZZ + b Z after a rotation), the ground space is extensive (dimension growing exponentially in L) only at J_zz > 0, |b| = 2 J_zz, the Fibonacci point. Checked by exact transfer-matrix counting for L up to 16.

R2 (frustration-free bonds, classified by the kernel S of the bond). Rank 1 (dim S = 3, S = psi^perp): the open-chain ground space grows exponentially only when psi is a product |aa> of two equal states. Every other psi, entangled or a product |ab> with a not parallel to b, gives dim G_L = L + 1. Rank 2 and rank 3: dim G_L stays bounded by L + 1 for all L, with no exponential growth. Pass: dimensions computed exactly (rank of the intersection) for L = 2 to 14 on structured and random kernels.

R3 (unbiased completeness check over B_2). Minimizing the normalized width of the lowest K(L) = round(g^L) states over the 8 parameters of B_2, jointly at L = 8 and 10, for several growth rates g, lands only on the Fibonacci orbit (g near the golden ratio) and the XZ orbit (g = sqrt 2), or fails to reach small width. A new extensive point is declared if the final normalized width is below 1e-3 at both L, at a bond further than 0.05 (in normalized coupling distance, after optimizing over rotations) from both known orbits. No new point is predicted.

R4 (first-order dynamics at the Fibonacci point). Projected onto the blockaded space P, every bond perturbation reduces to a nearest-neighbour hop w, a PXP creation amplitude Omega and a chemical potential mu. The first-order-invisible directions are the pairing (XX - YY, XY + YX), the antisymmetric off-diagonal terms XZ - ZX and YZ - ZY (DM along y and x), and the combinations of transverse field with XZ + ZX or YZ + ZY whose PXP amplitude cancels. Prediction for the second-order effect of lambda (XZ - ZX) on its own: a distance-2 interaction +8 lambda^2/V (repulsive), a real nearest-neighbour hop of size 4 lambda^2/V, no distance-3 interaction and no change of the NNN hop. Pass: exact Schrieffer-Wolff extraction from the full open chain (L <= 12) agrees within 3% at V/lambda >= 40, with deviations falling like 1/V.

R5 (XZ point). No new numerics planned. C144 and C157 cover it.

No expectation is registered for the valley-pseudospin sign of any unbound gas found.

### Addendum registration R6 to R8 (written after the two-body results of pseudo_v.out and before any many-body run of the modified model)

Context. Section 3 of PROOF finds that the transverse DM term λ(X⊗Z − Z⊗X) is invisible at first order at the Fibonacci point and adds a repulsive V₂′ = +8λ²/V, i.e. V₂ = −2 + v with v = 16λ²/Δ² in units |t₂|. pseudo_v.out shows that for v above the proved binding boundary u²/4 the dilute pseudospin criterion of C175 reads two-component.

R6 (effective model, U(1) part, C176's blockade-exact DMRG unmodified, φ = π/2, V₃ = −2/3, n = 3/16 unless stated). Points: |w| = 1 with v = 1 (dmv1) and v = 2 (dmv2), |w| = 2 with v = 3 (dm2v3), |w| = 1, v = 1, n = 1/8 (dmv1n8); control V₂ = V₃ = 0 at |w| = 1 (mgs, the Mishra–Greschner–Santos limit, two-component by C175). First-stage criterion, as C176's: truncation below 1e-7 at L = 128; c_two(64 → 128) > 1.2 for windows W = 1, 6, 12; c_prof(128) > 1.2; min(W_a, W_b)/max ≥ 0.8; C_P(L/4) below 0.05 P_pol² and falling. Calibration: the mgs control must read c_two(64 → 128) within 0.25 of 2. Prediction: dmv1, dmv2 and dm2v3 pass, with c_two(64 → 128) ≥ 1.6 at W = 6. Parity gap at dmv1 (N ± 1 at L = 64, 128): L Δ_P flat or falling within 20 percent (gapless pseudospin). A rise of L Δ_P by more than 50 percent is read as a pair gap, against the two-body result.

R7 (full qubit chain, V/Δ = 8, C177's fc.py with the DM option, tuned by C177's band fit at λ = 0 and the leading-order shifts J_x + J_y = −4λ²/V, h_z → h_z − 4λ²/V). λ = 0.3 (v ≈ 1.4 at leading order) at L = 64 and 128, and the λ = 0 reference at the same L. No pass is expected at these sizes, because the ΔN = 1 conditioned-creation terms are of order λ/Δ relative to |t₂| (PROOF 3.3) and the valley offset is tuned only at leading order. A reading of S̄(128) − S̄(64) more than 0.05 above the λ = 0 reference is noted as support, not as a pass.

R8 (degeneracy search, R3 above) is part of the same job.

### Amendment to R6 to R8 (written after sw_dm.out confirmed the coefficients, before any many-body output)

sw_dm.out shows a ΔN = 1 term (2Δλ/V) Σ_s σ⁻_s (n_{s+2} − n_{s−2}) + h.c., i.e. κ = 4λ/Δ in units |t₂|, tied to the repulsion by v = κ² at leading order. It is of order one where unbinding needs it (κ > u/2). So the U(1)-only runs (dmv1, dmv2) test only part of the model, and I added blockk.py, the full effective model with the κ term, 3-site blocks with C173's penalty U = 40 (validated against direct ED, test_blockk.out). The fibdmrg list is trimmed to dmv1 (L = 64, 128, with N ± 1), dmv2 and mgs at L = 64, 128 (no L = 256).

R6′ (blockk, open chains L = 48 and 96, two-size estimator c₂ = 6[S_mid(96) − S_mid(48)]/ln 2 at block cuts, grand canonical μ). Points: κ = 1, v = 1 at μ = −2.2, −2.1, −2.0; κ = 0, v = 1 at μ = −2.1; κ = √2, v = 2 at μ = −2.1; and κ = 0, v = 0 at μ = −2.2, which must reproduce C173's penalty-model reading (S̄ = 1.646 at L = 48 and 1.869 at L = 96, c₂ ≈ 1.9) within 0.1 as a code calibration. Registered expectation, honestly uncertain: κ = 0, v = 1 reads c₂ ≥ 1.6 (two-component U(1) part). For κ = 1 I register no direction. A κ = 1 reading of c₂ ≥ 1.6 at two of the three μ, with final truncation below 1e-6, is called support for the transverse-DM route at the effective-model level (not a pass of the C163 criterion, which needs L ≥ 128 and controls). A reading of c₂ ≤ 1.3 at all three μ while κ = 0, v = 1 reads ≥ 1.6 is called a failure caused by the ΔN = 1 term.

### Second amendment (job 17604455 still pending, no many-body output of the modified model seen)

Ring ED of (3.5) (ring_kappa.py, ring_kappa_L20.out, ring_kappa_L24.out) shows that κ = 1 moves the valley-balance point to φ* ≈ 1.635 to 1.645 (W_a/W_b crosses 1 between φ = 1.631 and 1.661 at L = 24, between 1.621 and 1.671 at L = 20). So the blockk runs were changed before the job started to κ = 1, v = 1, μ = −2.1 at φ = π/2, 1.625, 1.640, 1.655 (L = 48, 96), with the κ = 0 controls unchanged at φ = π/2 and κ = √2 dropped. R6′ now reads: support for the route if c₂ ≥ 1.6 at one of φ = 1.625, 1.640, 1.655 with final truncation below 1e-6 and a lower reading at φ = π/2; failure caused by the ΔN = 1 term if c₂ ≤ 1.3 at all four φ while κ = 0, v = 1 reads ≥ 1.6. The full-chain stream is now a φ scan at L = 64 (λ = 0.3, V/Δ = 8, φ shifted by 0, 0.06, 0.10, 0.14 through J_x + J_y and D at fixed |w|), plus the λ = 0 references at L = 64 and 128. R7 is unchanged in spirit: no pass is expected at L = 64, and a reading of S_mid above the λ = 0 reference by more than 0.1 at some φ shift is noted as support only.
