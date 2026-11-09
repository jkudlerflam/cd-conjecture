# C175 to C180: a swarm on the Fibonacci-cusp candidate for c = 2 at d = 2

## The candidate

C173 (REPORT §9.137, `codex/work_C173/`) found a corner of bond space where a nearest-neighbour qubit chain may have c = 2. If it holds, c ≤ 1 at d = 2 is false.

- Setting. Take J ≈ J_z ZZ with the field near its blockade value B ≈ 2J_z ẑ. Tune the first-order transverse hop w = (J_x + J_y)/2 + iD down to order Δ²/V, with Δ = J_x − J_y and V = 4J_z. The pairing Δ then acts only virtually across the Rydberg blockade.
- Effective model (C173 PROOF 3.2, units |t₂| = Δ²/(2V) = 1). Hard-core particles with no two adjacent, a NN hop w = |w|e^{iφ}, a NNN hop −1, interactions −2 n_j n_{j+2} and −(2/3) n_j n_{j+3}, and a chemical potential μ. It conserves N at this order. Particles never pass, so it is sign-free.
  - The band 2|w| cos(k − φ) − 2 cos 2k has two valleys.
  - At φ = π/2 an exact lattice symmetry (inversion × staggered rotation, present when J_x + J_y = 0) makes them degenerate with equal masses, and k_a + k_b = π.
  - The gauge c_j → i^j c_j maps the φ = π/2 model to a frustrated J1–J2 chain, with real NN hop |w| and NNN hop +1, plus the blockade.
- Reading. Effective-model DMRG at φ = π/2, |w| = 1, μ = −2.2 (n ≈ 0.17) gives a mid-chain entropy of 1.646, 1.790 and 1.869 at L = 48, 72 and 96, so c ≈ 1.93 for 48 → 96.
  - The controls read 0.90 for c = 1 and 1.88 for c = 2.
  - The runs are not converged: the truncation is 4e-6 at L = 96, and the estimate moves between 1.8 and 2.1 with the choice of cut window.
  - Ring ED gives a valley-balanced ground state, with the fully polarized states 0.07 to 0.10 higher at N = 4.
- Candidate bond, in the convention h = J_x XX + J_y YY + D(XY − YX) + J_z ZZ + (h_z/2)(Z⊗1 + 1⊗Z): J_x = 0.5, J_y = −0.5, D = 0.00625, J_z = 10, h_z = −19.97625. This is V/Δ = 40, with an energy scale of 0.0125.
- What decides it. Is the valley pseudospin easy-plane (c = 2), or does it choose one valley or order (vector-chiral, valley crystal, c = 1)? Does the full qubit chain agree with the effective model? Is the region open in B₂?
- Already queued. Typhon job 17603887 (C173's, pending) runs C173's block_dmrg.py to L = 192 at χ = 192 for φ = π/2, 1.54 and 1.52, both controls at L = 144, and the full chain at V/Δ = 40 (L = 96, 144) and 80 (L = 96). Do not duplicate those runs.

Six agents work on this at once. Your role is the section with your number, and you do only that. Other agents' folders are in progress, so you may read them but do not edit them.

## Shared reading for every agent

- `specs/_qubit_c_context.md` and STATUS.md
- REPORT.md §9.126 to §9.137
- `codex/work_C173/` in full: RESULT.md, PROOF.md (sections 3, 4 and 5 especially) and the scripts (blockade.py, block_dmrg.py, dmrg_d.py, full_dmrg.py, ring_ed.py, valley_ed.py, twobody.py)
- `codex/work_C162/RESULT.md` and `codex/work_C164/RESULT.md`, on two-valley gases and their pseudospin sign
- `codex/work_C163/PROOF.md` §4.2, the pre-fixed rules for a c > 1 candidate

## Shared rules for every agent

- Physics standard. RG, bosonization, perturbation theory and calibrated numerics count, and "not rigorous" is not a gap. Show all steps and include consistency checks (limits, special cases, matching to known results).
- Tag every claim proved, derived, numerical, heuristic or assumed. No numerical claim without the file it came from.
- Register any numerical prediction or pass criterion in RESULT.md before you look at the output, and report honestly against it.
- Your folder is `codex/work_C1xx/`, with your number. RESULT.md opens with a line saying which model did the work. Derivations go in PROOF.md and the log in PROGRESS.md.
- Write only inside your folder. Do not edit STATUS.md, REPORT.md, specs or other agents' folders. Do not commit or push. Never stash, reset or switch branches, do not pip install, and do not launch agents.
- Astra is out of credits, so do the work yourself. Never call codex directly, and never use paid APIs.
- Local CPU is shared by six agents. Set OMP_NUM_THREADS, OPENBLAS_NUM_THREADS and MKL_NUM_THREADS to 2 for local runs, and send anything longer than about 20 minutes to Typhon.
  - Use at most one packed 48-core node job per agent: `ssh -o LogLevel=ERROR typhon`, SLURM only, never compute on the login node, scratch /scratch/lustre/jkudlerflam/c1xx/.
  - The bastion rate-limits, so keep ssh calls few and retry with waits.
  - Record in RESULT.md the job id, the expected outputs and the exact scp -r collection command, and do not block on the job.
- Prose: commas not em-dashes, few colons and semicolons, no italics for emphasis, no filler intensifiers.
- Final report: the verdict, key claims with tags, and what is left open.

## C175: the pseudospin sign, analytically

Decide whether the dilute two-valley gas at φ = π/2 is two-component (easy-plane pseudospin, c = 2) or picks one valley or orders (c = 1), from the microscopic model (3.2).
- Solve the two-body problem exactly on the lattice, in the blockaded space with V₂ and V₃, for one particle in each valley and for two in the same valley. Large rings are cheap, and C173's twobody.py is a start.
- Extract the low-energy inter- and intra-valley scattering data (1D even-channel scattering length, or the T-matrix at the valley momenta), including the valley-exchange amplitude. Check for a zero-energy resonance.
- Build the pseudospin chain couplings J_⊥ and J_z to leading order in n, in the impenetrable regime (C173 PROOF 4, C162), with the three-body input if two-body does not fix the sign. State whether |J_z| < J_⊥ (easy-plane), and the sign of J_z.
- Map the answer against |w| and n, and give the predicted c = 2 window in (φ, μ) at small n.
- Cross-check against exact ring ED of (3.2) at N = 3 to 6. The pseudospin levels and their 1/L scaling should match your J_⊥ and J_z.
- If the dilute limit is inconclusive, carry the calculation to the density C173 used (n ≈ 0.17) by any controlled route you can find, for example a two-component Bethe-ansatz-like or bosonization estimate with the computed scattering data.

## C176: a blockade-exact DMRG and converged open chains of the effective model

The C173 runs used a penalty and were unconverged. Write a DMRG for (3.2) that is exact in the blockaded space and conserves N. One route is a local basis with the constraint in the MPO and a U(1) charge, or the Fibonacci-constrained local Hilbert space. Validate it against C173's ED and ring ED at L ≤ 24.
- Run φ = π/2, |w| = 1, μ = −2.2 open chains at L = 48, 64, 96, 128, 192 and 256, with truncation below 1e-8 where reachable, plus the c = 1 and c = 2 controls at the same sizes.
- Report the two-size estimates with several cut windows, the profile fits, S(x) at fixed x against L, and the energy convergence in χ.
- Also run φ = 1.54, 1.52 and 1.50, and two densities, μ = −2.0 and −2.4.
- Measure the competing orders:
  - the valley polarization (the order parameter of the k → π − k symmetry, from n(k) or from the vector chirality of the J1–J2 picture);
  - its correlations;
  - the valley-density-wave correlations.
- The pre-fixed first-stage criterion is C173 PROOF 5.4: two-size and profile c above 1.2 at L = 128 with truncation below 1e-7, a two-valley momentum distribution, and no long-range valley polarization. Register it before you run.

## C177: the full qubit chain

The decisive question is whether the actual nearest-neighbour qubit bond, not the effective model, reads c = 2. C173's job tests V/Δ = 40 and 80 at L ≤ 144 with a plain DMRG that may stall.
- Work at smaller V/Δ (for example 8, 12, 16, 24), where the scales are less separated and the full chain is tractable.
  - Tune h_z to the density n ≈ 0.17. Read C173's typhon/cfg/F*.json for the mapping from μ_eff to h_z.
  - Use J_x + J_y = 0, and D chosen so that |w|/|t₂| = 1 at each V/Δ.
  - Tie the reading to the effective model as V/Δ grows.
- Use C163/C166's validated dmrg.py and vumps.py (`codex/work_C166/`), with a parity-conserving or energy-shifted variant if it helps. Open chains to L = 256, and VUMPS FES at χ up to 256 or more.
- Compare the open-chain bulk energy with the VUMPS e₀, as in C166, so that a cat state or a trapped state is caught.
- Also take one generic point off the parity class and off J_x + J_y = 0 (add small J_xz-type and field-tilt terms of order Δ²/V), to test openness in B₂.
- Pre-register the criterion of C163 §4.2: open-chain fits above 1.2 that do not fall from L = 128 to 256 at truncation below 1e-8, a FES slope above 1.2 at χ ≥ 256, and energy agreement.

## C178: ring spectroscopy in the blockaded space

Exact diagonalization of (3.2) on rings, using the Fibonacci-constrained basis with momentum and N resolved. Push L as far as memory allows, about 36 to 44 sites.
- Identify the low-energy spectrum at φ = π/2, μ = −2.2 and its neighbours. Ask whether it is a two-component Luttinger liquid, with two velocities and a compactification lattice for two bosons, or one boson plus a gapped or slowly closing pseudospin sector.
- Track the valley-pseudospin gap (the lowest excitation that changes N_a − N_b, or the antisymmetric sector of the k → π − k symmetry) against 1/L. Fit the ground-state energy to e L − π(v₁ + v₂)/(6L) with independently measured velocities.
- Compare with the c = 1 and c = 2 controls at the same sizes.
- Do this at φ = π/2 and at φ = 1.54 and 1.50, and at two densities.
- Pre-register what a c = 2 spectrum looks like (towers, a 1/L pseudospin gap, the Casimir coefficient) before you look.

## C179: referee and hidden-symmetry check of C173

Referee C173 PROOF sections 3 to 5 independently.
- Re-derive the second-order Schrieffer–Wolff model: the signs, t₂, V₂, V₃, the diagonal shifts, and the order of every omitted term, including the ΔN = 2 processes and the parity-odd terms off the parity class.
- Check the claimed lattice symmetry at J_x + J_y = 0.
- Check the openness argument of PROOF 5.4. Every inter-valley vertex must carry nonzero momentum off resonance, so list the vertices and their momenta. At φ = π/2, k_a + k_b = π, so a vertex creating two A–B pairs carries 2π ≡ 0. Find the lowest order at which it appears, its scaling dimension in the two-component liquid, and whether it can lock a mode, as C153's resonant conversion vertex locks the density at w = 0.
- Do C163 §4.2's hidden-symmetry check. Is the candidate bond equivalent, by an exact site unitary, a duality or a known map, to a free, integrable or decoupled model (for example the staggered-charge U(1) class of C153 at D = 0, via particle-hole on one sublattice)? Does the D term break that structure at the order that matters?
- Write REFEREE.md with each claim marked correct, gap (fixable), or wrong, with corrected statements where needed.

## C180: field theory and literature

- Bosonize the two-valley gas at φ = π/2 and nearby: two bosons with a 2 × 2 Luttinger matrix.
  - List every operator allowed by the lattice symmetries (translation, the k → π − k symmetry at φ = π/2, N conservation in the parity class and its breaking to N mod 2 or nothing off it), with its momentum.
  - Give the scaling dimensions as functions of the Luttinger matrix.
  - State the exact conditions under which c = 2 is stable on an open set, and what the competing c = 1 phases are (vector chiral, valley crystal, pair or bound-state liquids).
- Use the frustrated J1–J2 chain near saturation as the calibration. Find in the literature what is known about the two-component versus vector-chiral competition there, including Rydberg and constrained versions and models with NNN hopping generated at second order. Give citations you have checked (arXiv numbers), not recalled ones.
- Check whether a c = 2 phase of a translation-invariant nearest-neighbour spin-1/2 chain is already claimed or excluded anywhere. Look for "two-component Luttinger liquid" in nearest-neighbour chains, Rydberg chains with pairing, the XYZ chain in a field, and constrained PXP-type models.
- Say where the candidate sits relative to Latorre and Sierra's conjecture (arXiv:2403.17242) at d = 2.
