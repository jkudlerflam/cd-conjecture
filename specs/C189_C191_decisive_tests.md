# C189 to C191: the decisive tests of the Fibonacci-cusp candidate (2026-11-09)

Coordinator: Claude Opus 5.5. Context: REPORT §9.137 to §9.159 and STATUS item 1.

## Where things stand

C173's second-order effective model of the parity-class cusp (C173 PROOF eq. 3.2) reads c ≈ 1.9 on open chains from L = 128 to 256 on the symmetric line φ = π/2, |w| = 1, at n = 1/8, 5/32 and 3/16 (§9.158). It also reads that off the line at φ = 1.50 to 1.565 (§9.156, §9.157). The controls are calibrated, the parity gap L·Δ_P is flat (about 2.4), and single-particle and pair correlations both decay as powers.

The two-body pair of one particle from each valley is bound (C175, §9.139, exact boundary §9.140). So the c = 1 alternative is a pseudospin (relative-valley) gap far below the finite-size spacing, exponentially small as in a weakly attractive 1D two-component gas. C175 estimates about 1e-4 at n ≈ 0.17. Open chains at L ≤ 256 cannot see that.

The real qubit chain is barely tested at the candidate. The collaborator's full-qubit run at the φ = 1.53 parity-family bond (w = e^{1.53i}/80, Δ = 1, V = 40, h_z = −19.97589125) reads c = 2.15 from L = 32 to 64 but 0.74 to 0.93 from 64 to 128. Those L = 128 states are not fully optimized, and the bulk density drifted to about 0.22 (§9.155, contrib/d2_dmrg_priority_20261109/, contrib/d2_tenpy_20261109/).

## Common rules for C189 to C191

- Work in codex/work_C18x/ (your number below). Keep RESULT.md (main file, with the registration), PROOF.md (derivations), PROGRESS.md (log).
- Register every pass/fail criterion in RESULT.md before you look at any production output that decides it.
- Physics standard (RG, perturbation theory, calibrated numerics), but show every step and run consistency checks: limits, controls with known c, matching to known results.
- Typhon: `ssh -o LogLevel=ERROR typhon` (retry with ~30 s waits, the bastion rate-limits). SLURM only, never compute on the login node. Partition typhon, CPU. Scratch /scratch/lustre/jkudlerflam/c18x/<yours>/. Submitting is always OK. Never touch other jobs, and never release held jobs. Name jobs c189, c190, c191 (the name must start with c1).
- Local machine: keep local load modest (at most about 8 cores). Large runs go to Typhon.
- Don't pip install anything. numpy and scipy are available.
- Don't commit or push. The coordinator does that. No probabilities in any file. No absolute local paths in files you write.
- No Fable subagents. Astra (gpt_solve.sh) is out of credits, so don't call it.
- Prose: commas, not em-dashes, few colons, no italics.
- Report at the end: what you established, with numbers, what failed, and what's left. Put Typhon job ids and collection commands in RESULT.md.

## C189: is the pseudospin backscattering relevant? K_s from the data, plus infinite-chain entanglement scaling

Goal: decide between the c = 2 liquid and a hidden exponentially small pseudospin gap in the effective model, on the symmetric line, |w| = 1, n = 3/16 (with 5/32 and 1/8 as checks).

1. Extract the Luttinger parameters K_c and K_s from existing converged data, without new runs first.
   - Data: codex/work_C176/typhon_out/out/*.jsonl, at L = 128, 192 and 256, χ up to 1024. The last line of each file holds 'meas', with density, the boson and fermion one-body matrices, C_nn, C_P and Dpair (separations 2 and 3). analyze.py shows the format.
   - C181's off-line data (codex/work_C181/typhon_out*/) give the valley imbalance p against the offset δ(φ), which is linear in π/2 − φ with slope about 3.9 at |w| = 1 (C181 RESULT F3). That gives the pseudospin susceptibility.
   - Several independent estimators:
     - the density structure factor at small q, which gives K_c;
     - the exponents of the single-particle function (about −2.0) and the inter-valley pair function (about −0.78), with the operator content of each derived for this two-valley blockade model (state which frame, c or d, and how the Jordan–Wigner string enters);
     - the pseudospin susceptibility against v_s;
     - bipartite fluctuations of the valley imbalance if you can build them from the one-body matrix.
   - Derive the relations you use in PROOF.md. Check each estimator on the c = 2 control (two XX legs, K_s = 1 exactly) and on the hard-rod c = 1 control, using the same files.
   - The decisive question is K_s against 1, and whether the estimates drift with L in the way a marginal flow would (logarithmic corrections, K_s → 1).
2. Infinite-chain test. Write an iDMRG or VUMPS for C173's effective model, blockade-exact and with U(1) number conservation at fixed density, or at fixed μ with the density monitored. Reuse C176's fibdmrg.py conventions and validate against it on finite chains and against C176's bulk energy density. Run S against log ξ over a χ sequence at n = 3/16 and 1/8, plus both controls. Separate the correlation lengths in the charge and pseudospin sectors if you can, using the transfer-matrix spectrum by quantum number.
   - A c = 2 liquid gives S = (c/6) log ξ with c ≈ 2, and ξ growing without bound.
   - A gapped pseudospin sector shows a pseudospin correlation length that saturates while the charge one keeps growing, and the slope tends to 1/6.
   - Register the criterion before running. Run on Typhon (one node).
3. Also prepare codex/work_C189/desktop_run/ for Jonah's desktop: a self-contained package (code, configs, run_desktop.sh using about 8 cores, summarize.py, README.txt and PROMPT.txt in the style of codex/work_C187q/desktop_run3/). It should run the χ sequence that Typhon can't reach in time, or a second density.

## C190: a converged full-qubit test at the candidate bond

Goal: a converged real-chain answer at the parity-family candidate, so that the full-qubit reading at L = 128 and 256 can be compared with the effective model at the same bulk density.

- Bonds:
  - The collaborator's φ = 1.53 parity-family bond and their rational simplification, 10 ZZ + (XX − YY)/2 + (XX + YY)/4000 + (XY − YX)/160 with uniform field −199759 Z/10000.
  - The symmetric-line bond (φ = π/2) at the same V/Δ = 40.
  - Optionally V/Δ = 80, where corrections to the effective model are smaller.
  - Derive the mapping (C173 PROOF, C177 RESULT): the effective w, t₂ = −Δ²/(2V), V₂, V₃, and how the qubit magnetization maps to the effective density. Tune h_z so that the bulk mid-chain density matches the effective model's bulk density at n = 3/16 and at n = 1/8 (C176's n_mid values, about 0.196 and 0.133 at L = 256).
- Tools: C177's dmrg2.py (Davidson DMRG with an exact fermion-parity selector, matches ED to 1e-11) or C183's code. Read the collaborator's contrib/d2_dmrg_priority_20261109/RESULT.md and SIMPLE_PARITY_CANDIDATE.md for their convergence problems (local residuals 1e-5 to 1e-6, energy changes comparable to the parity splitting at L = 128) and fix them:
  - more sweeps;
  - χ to at least 512 at L = 128 and 768 at L = 256;
  - subspace expansion or noise early, then off;
  - both parity sectors;
  - several initializations;
  - energy and entropy convergence against χ reported;
  - a residual target registered in advance.
- Measure c_two over windows W = 1, 6, 12, c_prof, the bulk density profile, the parity splitting, and the valley weights, with the XX and two-leg XX controls at the same sizes and χ. Compare directly with C176's effective-model S̄ and L·Δ_P at the matched density.
- Register in RESULT.md before production. The verdict should say whether the full chain tracks the effective model's c ≈ 1.9 at L = 128 → 256, or falls toward 1, and whether the collaborator's L = 128 downturn reproduces once converged.

## C191: the pseudospin RG from the scattering data

Goal: an independent analytic prediction of whether the pseudospin backscattering is relevant at the candidate densities, to compare with C189's numbers.

- Inputs:
  - C175 (codex/work_C175/RESULT.md, PROOF.md, smatrix.out, scan_w.out): the two-valley scattering data and its pseudospin criterion, c = 2 if A_d − |A_o| < a_s < A_d + |A_o|.
  - The exact binding boundary in REPORT §9.140, which comes from the collaborator.
  - C179's REFEREE.md and PROOF.md (codex/work_C179/): the resonant operators on the symmetric line, the two-pair vertex O₄, and the K_ρ estimate 3.2 to 3.6 in its convention.
  - C181 PROOF (the FFLO structure off the line).
- Task:
  - Build the low-energy theory of the two-valley gas at finite density: two Fermi seas at the valley points, the forward-scattering couplings (intra- and inter-valley), the inter-valley backscattering g₁ (momentum transfer between valleys), and any umklapp or pair-hopping terms allowed at π/2.
  - Derive K_c, K_s and the initial g₁ from the scattering lengths (dilute limit, with the density corrections you can control) and from the bare couplings V₂ = −2 + v, V₃ = −g.
  - Then solve the KT flow of (K_s, g₁) and say on which side of the separatrix the candidate densities start.
  - Account for the fact that pseudospin is not SU(2) here, since intra- and inter-valley interactions differ.
  - Give the gap scale if it flows to strong coupling, and the crossover length.
  - Check the construction on known limits: free two-leg XX (decoupled, c = 2), the Gaudin–Yang attractive gas (gapped for any attraction), and the hard-rod chain.
- Register your prediction for C189's K_s and its drift in RESULT.md before C189's numbers exist. You may read C189's RESULT.md only after your prediction is written. Keep the derivation in PROOF.md, with every step.
