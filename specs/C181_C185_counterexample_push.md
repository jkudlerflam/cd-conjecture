# C181 to C185: pursue every route to a c > 1 nearest-neighbour qubit chain

## Where things stand

The goal is a translation-invariant nearest-neighbour qubit chain whose low-energy theory is a unitary CFT with c > 1. That would break Latorre and Sierra's c ≤ d − 1 at d = 2. A fine-tuned point counts. A multicritical point at any finite codimension in bond space B₂ is a counterexample if its low-energy theory has c > 1.

**The Fibonacci cusp** (REPORT §9.137 to §9.139, `codex/work_C173`, `C175`, `C179`, `C180`). This is strong ZZ with the field near its blockade value. The transverse hop is tuned to the size of the generated next-nearest hop t₂ = −Δ²/(2V). The effective model is C173 PROOF eq. 3.2: hard-core particles with no two adjacent, NN hop w = |w|e^{iφ}, NNN hop −1, interactions −2 n_j n_{j+2} and −(2/3) n_j n_{j+3}, and chemical potential μ, in units |t₂| = 1. Its band has two valleys, degenerate at φ = π/2 (the J_x + J_y = 0 line).
- At φ = π/2 one particle from each valley binds into a pair at total momentum π. Binding is 0.0195 at |w| = 1 and 0.1226 at |w| = 2. Trimers bind too. So the dilute gas is a pair liquid with c = 1 (C175, checked independently).
- At n ≈ 0.17 the pair gap may be tiny, and C173's c ≈ 2 reading at L ≤ 96 is probably a crossover, but this is not excluded.
- In the full qubit chain on J_x + J_y = 0, a two-pair vertex at momentum 2π appears at order Δ⁴/V³. With K_ρ ≈ 3.2 to 3.6 from rings it is relevant and locks the total density (C179).
- Off J_x + J_y = 0 both resonant operators carry momentum.
- C179 suggests an untested heuristic route. A valley offset acts as a Zeeman field on the pairs, as in the attractive Gaudin–Yang gas, where the partially polarized phase is a two-component liquid with c = 2.

**Running now** (read their folders, do not edit them):
- C176, a blockade-exact DMRG of the effective model;
- C177, the full qubit chain at V/Δ = 8 to 24;
- C178, ring spectroscopy in the blockaded space;
- C179, finishing.

Typhon also has C173's and C175's effective-model runs.

**Other routes already examined:**
- C169 found no exact map from a known c > 1 model, and argued Z₅ parafermions and Potts × boson (c = 9/5) need at least two tunings at amplitudes of order the exchange.
- C165 closed the staggered Ising mode inside a liquid near the hidden Heisenberg point.
- C162 and C164 closed symmetric two-valley soliton gases.
- C171 listed routes R3 (composite second bosons) and R4 that remain open.
- STATUS notes that isolated multicritical points have not been inventoried.

Six agents work in parallel (C181 to C186). Your role is the section with your number, and you do only that.

## Shared reading for every agent

- `specs/_qubit_c_context.md` and STATUS.md
- REPORT.md §9.126 to §9.139
- `specs/C175_C180_fibonacci_swarm.md`
- RESULT.md of C173, C175, C179 and C180, plus whatever your section names

## Shared rules for every agent

- **Standard.** Use the physics standard. RG, bosonization, perturbation theory and calibrated numerics count, and "not rigorous" is not a gap. Show all steps and include consistency checks (limits, special cases, matching to known results).
- **Tags.** Tag every claim proved, derived, numerical, heuristic or assumed. No numerical claim without the file it came from.
- **Registration.** Register predictions and pass criteria in RESULT.md before looking at output, and report honestly against them.
- **Criterion for a c > 1 candidate.** Adapt C163 PROOF §4.2:
  - converged open-chain fits above 1.2 that do not fall with L, at truncation below 1e-7;
  - FES or ring-spectroscopy confirmation;
  - for an effective model, a check on the full qubit chain;
  - controls at c = 1 and c = 2 at the same sizes.
  - Watch for long crossovers from small gaps. Measure the relevant gap (parity gap, pseudospin gap) against L, not only the entropy.
- **Folder.** Your folder is `codex/work_C18x/`, with your number. RESULT.md opens with a line saying which model did the work. Derivations go in PROOF.md and the log in PROGRESS.md.
- **Writing only in your folder.** Do not edit STATUS.md, REPORT.md, specs or other agents' folders. Do not commit or push. Never stash, reset or switch branches, do not pip install, and do not launch agents.
- **No outside compute.** Astra is out of credits, so do the work yourself. Never call codex directly, and never use paid APIs.
- **Local CPU.** It is shared by about nine agents. Set OMP_NUM_THREADS, OPENBLAS_NUM_THREADS and MKL_NUM_THREADS to 2, and send anything longer than about 20 minutes to Typhon.
- **Typhon.**
  - At most one packed 48-core node job per agent, through `ssh -o LogLevel=ERROR typhon`, SLURM only, never computing on the login node.
  - Scratch is /scratch/lustre/jkudlerflam/c18x/.
  - Keep ssh calls few and retry with waits.
  - Record the job id, the expected outputs and an exact scp -r collection command in RESULT.md. List explicit remote paths, because brace expansion fails. Do not block on the job.
- **Reuse code.** C176's blockade-exact DMRG, once it exists in `codex/work_C176/`, is the right tool for the effective model, so check there before writing your own. C166's dmrg.py and vumps.py and C163's tools are validated for the full chain.
- **Prose.** Commas not em-dashes, few colons and semicolons, no italics for emphasis, no filler intensifiers. No probabilities in any file.
- **Final report.** The verdict, the best candidate bond if any (explicit couplings in the convention h = J_x XX + J_y YY + D(XY − YX) + J_z ZZ + (h_z/2)(Z⊗1 + 1⊗Z), or a general 4×4 bond), key claims with tags, and what is left open.

## C181: the partially polarized pair route off the symmetric line

Test C179's heuristic. Off φ = π/2 the valley offset δ acts as a field on the ab pairs.
- **The analogy.** In the attractive Gaudin–Yang gas, paired (c = 1), partially polarized (two-component, c = 2) and fully polarized phases occur as the field grows.
- **The difference here.** Trimers bind (C175), and the valleys have unequal masses off π/2 (C164's immiscibility input). Either could remove the mixed phase.

Tasks:
- Derive the dilute phase diagram in (δ, μ):
  - the pair–atom scattering length and the trimer binding as functions of φ;
  - whether a stable mixture of pairs and excess atoms exists, and over what window;
  - its Luttinger description.
- Use the exact two- and three-body problems on rings or in the infinite-line scattering formalism (C175's scat.py, nbody_box.py and analytic_res.py are a start).
- Then test numerically on the effective model:
  - φ from 1.48 to 1.57;
  - several μ, at |w| = 1 and |w| = 2 (compact pairs make the phases cleaner);
  - c from two sizes;
  - the parity gap;
  - n(k) (the atom peaks);
  - pair correlations.
- If a c = 2 window exists, give the best explicit qubit bond off J_x + J_y = 0 and check the C179 vertex list there.

## C182: one level up, cluster bands and two-component pair liquids

The ab pairs, and aa or bb pairs if bound, have their own bands. C180 noted that a pair band with two minima repeats the two-valley question one level up.
- Compute the exact dispersion of the bound pair against total momentum K, and of trimers, at |w| from 0.3 to 3.5 and φ from 1.3 to π/2.
- Find where a cluster band has two degenerate or near-degenerate minima, or where pairs and single atoms coexist at the same chemical potential.
- For each such region, run C175's analysis one level up:
  - cluster–cluster scattering;
  - the pseudospin criterion;
  - bound states of clusters.
- Map the effective model's full phase diagram in (|w|, φ, μ) for every c > 1 opportunity:
  - two-component atoms;
  - pair–atom mixtures;
  - pair two-valley liquids;
  - coexisting pair and atom Fermi seas.
- Use cheap diagnostics first, then converged DMRG at the best points.
- Read C181's folder as it develops so you don't duplicate the pair–atom route. Your target is the cluster-band and phase-diagram side.

## C183: the full qubit chain near the cusp, targeted

The effective model may miss what decides c, such as C179's two-pair vertex, and higher-order terms at moderate V/Δ. Search the full qubit chain directly.
- Work at V/Δ = 4 to 16, where it is tractable.
- Cover the neighbourhood of the cusp in all its directions:
  - |w|/|t₂| from 0.5 to 3;
  - φ from 1.4 to π/2;
  - h_z across the density range n = 0.05 to 0.35;
  - J_x + J_y off zero;
  - generic symmetry-breaking terms of order Δ²/V (XZ, ZX, field tilts) that leave the parity class.
- Use a cheap two-size entropy screen with controls, then confirm the best points with converged open chains to L = 256, VUMPS FES and the parity gap.
- Moderate V/Δ is where the generated terms are not small, which is a chance for something the effective model does not contain.
- Do not duplicate C177's points (read its folder and its registered list). C177 tests the C173 bond family at V/Δ = 8 to 24. You search the surrounding space.

## C184: multicritical points with c > 1

A fine-tuned point counts. Inventory the ways a nearest-neighbour qubit bond can reach a multicritical point with c > 1 at finite codimension, and try to build one.
- **Candidates:**
  - Ising × Gaussian (c = 3/2) where an Ising surface crosses a gapless liquid or a Gaussian transition surface;
  - Ising × Ising × boson;
  - tricritical Ising × Ising (c = 6/5);
  - three-state Potts × boson (c = 9/5);
  - Z₅ parafermions (c = 8/7);
  - SU(2)₂ or other WZW multicritical points;
  - the c = 2 point where two Gaussian transition surfaces meet;
  - Lifshitz-type points, ruled out where z ≠ 1.
- **For each:**
  - What does the crossing generically flow to? Coupled at the crossing, decoupled, or first order. Do the scaling dimensions of the coupling operators make the decoupled product stable or unstable?
  - Is there a nearest-neighbour bond family where both transitions exist and can be tuned to cross? C155 and C163 numerics and the strata of C154 PROOF 7 tell you where Ising and Gaussian surfaces are.
- **Then the most promising crossing.** Locate two transition surfaces in a two- or three-parameter slice of B₂, tune to the crossing, and measure c there with converged DMRG and ring spectroscopy.
  - At a crossing, open-chain c can read c₁ + c₂ only if both sectors are critical at the same point. Use two-size estimates with careful tuning, or a finite-size crossing analysis.
- **Read first:** C169's RESULT and PROOF on parafermions and Potts, C154 PROOF 7 on strata and the closure lemma (does the closure lemma exclude multicritical c > 1? check its reasoning critically), and C174's REDTEAM.md.

## C185: other degenerate points where second-order dynamics can win

The Fibonacci cusp exists because a dominant nearest-neighbour term has an extensive degenerate ground space, inside which the first-order dynamics can be tuned to the size of the second-order dynamics. Find every such point.
- **Classify.** Find nearest-neighbour qubit bonds H₀ whose ground space on the chain is extensively degenerate. Examples:
  - diagonal bonds a ZZ + b Z at their critical field (Fibonacci);
  - rank-1 and rank-2 projector bonds;
  - frustration-free bonds with exponentially many ground states;
  - classical bonds with degenerate manifolds after a local rotation;
  - the ferromagnetic multiplet, which is only L + 1 dimensional, so say why it differs.
- **For each,** derive the constrained low-energy space and the first-order dynamics allowed in it. Find where the first-order part can be tuned to vanish or be small, and derive the second-order effective model (Schrieffer–Wolff, as in C173 §3, checked against the full chain at L ≤ 12).
- **Then ask whether the effective model has:**
  - two valleys;
  - a second species;
  - a known c > 1 phase (four-Fermi-point liquids, two-leg ladders with independent legs, a decoupled spin sector).
  - Look specifically for cases without the inter-valley binding that kills the Fibonacci case at φ = π/2. A repulsive inter-valley channel is what you want.
- Test any promising effective model with small DMRG or ring ED, and give the qubit bond.

## C186: unbind the pair by changing the microscopic bond

**New shared input, 2026-11-09 (commit `3c7dda4`).** The collaborator has
now carried this route further in
[`contrib/counterexample_push_20261109/FACILITATED_ROUTE.md`](../contrib/counterexample_push_20261109/FACILITATED_ROUTE.md).
Read it before repeating the reduction. A genuine transverse DM term
`(lambda/2)(XZ-ZX)` generates both the desired repulsion and an
unavoidable leading density-assisted creation term. The complete model,
literal 4-by-4 bonds, independent operator checks, dressed atom/pair
resolvents and all-number ring spectra through L=24 are supplied.
After balancing the two dressed valley energies, the calculated
order-r² direction weakens binding; the apparent finite-ring crossing
has strong size drift. No infinite-volume unbinding or c=2 is claimed.

For this family the next unresolved step is **dressed scattering with
particle-number mixing**, including larger-cluster channels and
nonperturbative valley retuning. A fixed-N two-body calculation or the
replacement `v=4(lambda/Delta)^2` alone omits a term of the same leading
scale and does not test the actual microscopic deformation. The
exploratory small-u path in the linked calculation additionally needs
chemical-potential retuning and a remainder bound on the binding-energy
scale. The general-bond classification below remains useful for finding
other independent interaction knobs.

A collaborator (GitHub: salpace) proved the binding boundary exactly (`contrib/cusp_binding_certificate_20261109/`, `contrib/cusp_scattering_20261109/`, `contrib/cusp_validation_20261109/`, read them first).
- In the symmetric effective model take NN hop u, NNN hop +1, V₂ = −2 + v and V₃ = −g, with 0 < u < 4 and 0 < g < 2. A bound opposite-valley pair exists exactly when v < g u²/[2(2 − g)].
- At the Schrieffer–Wolff values (v = 0, g = 2/3) the pair is bound. With V₂ = V₃ = 0 the gas is two-component (C175, Mishra–Greschner–Santos).

So the c = 2 route at the cusp needs a nearest-neighbour qubit bond whose generated interactions land on the unbound side, with v > g u²/[2(2 − g)] or g reduced, while keeping the two valleys.
- Derive the second-order effective model for the most general nearest-neighbour bond near the Fibonacci point: a general 4×4 bond = dominant J_z ZZ + field at its blockade value + all small terms at first and second order (XX, YY, XY, YX, XZ, ZX, YZ, ZY, X and Y fields, and the deviation of h_z). Find which terms enter V₂, V₃, the NNN hop and the NN hop, and whether any choice moves (v, g) across the binding boundary at fixed |w|/|t₂| and φ = π/2 (or near it).
  - Example. Terms like XZ + ZX hop a particle conditioned on its neighbour, and at second order they shift the interactions at distances 2 and 3.
- Check the symmetry (the valley degeneracy needs the C173 PROOF 3.3 symmetry or a substitute), the number-changing terms, and the two-pair vertex of C179 for any modified bond.
- If a modified bond unbinds the pair, test it. Run the two-body problem (the collaborator's exact_pair.py or C175's code), then C175's pseudospin criterion, then effective-model DMRG with the parity gap, then the full qubit chain at moderate V/Δ. Give the explicit qubit bond.
- If no nearest-neighbour bond can unbind it, prove that, and say what the obstruction is.
