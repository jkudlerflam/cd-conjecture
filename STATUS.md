# Where the c–d conjecture stands

Last updated 2026-10-27. Read this page first. It is the current state of the project in one place. [REPORT.md](REPORT.md) is the chronological record behind it, and where the two disagree, this page is the one to trust. Every claim carries a tag saying how it was established (see "Tags" at the end).

## The question

Latorre and Sierra ([arXiv:2403.17242](https://arxiv.org/abs/2403.17242)) consider H = Σ_i h_{i,i+1} on (ℂ^d)^⊗L, with one fixed Hermitian h on ℂ^d ⊗ ℂ^d, periodic boundary conditions and one-site translation invariance. They conjecture that if the low-energy spectrum is described by unitary CFT(s), the total central charge obeys c ≤ d − 1. The SU(d) swap chain (Uimin–Lai–Sutherland, "ULS") flows to SU(d)₁ with c = d − 1, so the claim is c_max(d) = d − 1. We also ask the weaker question of whether any bound c ≤ f(d) exists.

The working standard is theoretical physics (standard RG, controlled perturbation theory, numerics calibrated with at least two estimators), not mathematical rigor.

## Summary

1. **No counterexample at small d.** Nothing above d − 1 at d = 2, in symmetric families at d = 3 to 7, or at the best d = 4 target. One d = 3 hint is being checked.
2. **No general proof.** Every lattice quantity we know that bounds c grows like (a length scale) × log d, never like d. Proofs exist in restricted classes.
3. **A candidate construction against it at large d.** If it holds, c is unbounded at a fixed local dimension d₀, so there is no bound c ≤ f(d) and c–d fails for every d ≥ d₀ (direct sums propagate a violation upward). It is not yet a counterexample, because one ingredient, a fixed-alphabet computation track, is missing and d₀ is not yet fixed.

## The construction, step by step

The idea fits in a paragraph. Sites spontaneously group into long cells. Inside each cell a local "computation" prepares a degenerate code space, which acts as one big effective site of dimension N ≫ d. Neighbouring cells interact only through the single bond at their junction, so perturbation theory in that bond gives an effective nearest-neighbour chain of big sites with a logical Hamiltonian we choose. We choose a stable CFT with large c. The difficulty is doing all of this with one repeated h on a fixed ℂ^{d₀}.

1. **Spontaneous cells.** An explicit d = 4 chain in which a "head" symbol hops freely inside a segment, so a segment of length b costs the kinetic energy −2t cos(π/b). Tuning one chemical potential then selects any cell length b ≥ 5, with couplings independent of L and a gap for each b. [checked] `codex/work_C7/C7_result.md`.
2. **A code in each cell, with ports that stay put.** The cell runs a history-state computation (Feynman–Kitaev type) that prepares the code. If the computation never acts on the cell's two end sites (the "ports"), every port operator commutes with the circuit, so the code and the compression of port operators onto it are exact, and the effective coupling of two cells converges to its target as 1/J_c in the clock coupling. [checked] `contrib/programmable_cell/`.
3. **Coupling the cells.** There are two routes.
   - *Route A, second order (the main line).* Choose a code on which every traceless port operator compresses to zero, so first order vanishes, and engineer the cell's excited spectrum so the second-order bond is the logical bond we want. The G1 cell gives aI + bP_N (a swap) for every N and any port dimension, confirmed by full diagonalization [checked, rerun from scratch] (`contrib/unbounded_c/G1/`). It works with identical cells that have distinct left and right ports [Astra] (`codex/work_C8/`). Any Hermitian logical bond can be made through two qubit ports in a cell of dimension 4⌈(6N³ + 13N)/4⌉ [Astra] (`codex/work_C9/RESULT.md`). Route A makes d₀ = 8k_c, independent of N, where k_c is the alphabet of the computation track (step 5) [outline, unchecked].
   - *Route B, first order.* The ports compress to scaled anticommuting Majoranas, αγ_a on the left and αiγ_aP on the right, so the inter-cell bond is free Majorana hopping with c = r/2 for r channels [checked]. An exactly degenerate cell that does this exists [checked] (`contrib/programmable_cell/`). Route B is capped by bond rank at c ≲ q²/4 for port dimension q, so it cannot give unbounded c. It can beat d₀ − 1, but only for q > 16k_c, which means d₀ ≈ 64k_c². Its flavour-symmetric version is unstable for m ≥ 7 channels (step 4).
4. **The logical target and its stability.** The target is an n-flavour chain of logical complex fermions with distinct, incommensurate Fermi momenta, giving c = n. Detuning the momenta removes every marginal non-Gaussian channel, so the target needs no exact symmetry beyond an on-site Z₂ that acts as logical fermion parity. Sidon-set momenta give margins δ₂ ≥ π/(6n) and δ₄ ≥ π/(576n²). A small iDMRG toy agrees: the symmetric chain dimerizes and the detuned one does not. [Astra, checked at the physics standard] `contrib/ir_control/`. Symmetric targets do not work. SU(N)₁ (the ULS-type target G1 makes) and the SO(2n) Majorana point sit on separatrices and run away under small anisotropy, and an exact SU(N) cannot be protected through a small port for N > q at any order [Astra]. The symmetric route-B family turns marginally relevant at m = 7, and keeping U(1)^m flavour charges with distinct fillings repairs it [Astra] (`codex/work_C11/`). The infrared analysis was done for the encoded chain in general. That it applies to route A's second-order realization unchanged has not been checked separately.
5. **The computation track, which fixes d₀.** An explicit translation-invariant nearest-neighbour track runs any hard-coded gate program, with proved bounds on illegal sectors and on the selector margin. A small example has k_c = 24, q = 2, d₀ = 192 [Astra] (`codex/work_C10/RESULT.md`). But k_c grows with the cell size and the program length, so as written d₀ grows faster than c and gives no counterexample. **This is the missing ingredient.** It needs a fixed-alphabet compiler of the Gottesman–Irani or Kohler–Piddock–Bausch–Cubitt type, or a much cheaper cell. An alphabet count based on the smallest known translation-invariant history machine (Bausch–Cubitt–Ozols, 42 dimensions) gives d₀ ≈ 130 for route A with the C7 selector tensored on, and d₀ ≈ 45 if the machine's own clock can select the cell length [heuristic; the key lemma, that one static Bell-check term is the only term that must touch the port, is proved by Astra and checked] (`contrib/minimal_d0/`). These counts are for frozen cells. A minor point at the same step: at finite clock coupling the history states are ground states only up to O(1/J_c), which is acceptable if that error only renormalizes the logical couplings.
6. **Frozen versus mobile cells.** If the cell markers do not move, translation is spontaneously broken. Even on rings with b | L the ground state is b-fold degenerate (momenta 2πk/b), and on other ring lengths there are defects, with at least L ground states on large prime rings. That is a gapless theory with c = n on top of a broken translation, not a CFT in the strict sense the conjecture intends. The fix is to let the cells move. A mobile selector lets cell lengths fluctuate between b and b + 1 as a free c = 1 "length fluid" on every ring [Astra numerics: c ≈ 1.004 from entropy, 0.9999997 from the Casimir energy]. The logical payload rides on it, as in spin–charge separation, which should give a unique ground state and c = n + 1 on every ring. Two pieces of this are written. An "equalizer" makes the first-order logical bond independent of the lengths of the two adjacent cells, at local dimension 96q + 48 (`codex/work_C8/helix_equalizer.md`). With it, the complete first-order spectrum of length fluid plus any payload is known exactly: the two sectors decouple, with their momenta glued by one global constraint (`codex/work_C8/joint_fluid.md`; every eigenvalue matches direct diagonalization to 10⁻¹⁴, and the Casimir energy gives c = 2.0000 for an XX payload) [Astra; check scripts rerun by the coordinator]. What is not yet shown is that this survives the corrections at finite coupling (leakage, residual length–payload terms, commensurate locking of the length fluid). That is being worked on now (`codex/work_C13/`, and microscopic numerics in `contrib/mobile_cells/`). A separate scope question remains, whether the conjecture could be read as requiring conformal behaviour only on compatible rings, but the mobile version would make it moot.
7. **In progress.** Finite-coupling numerics on both routes (`contrib/routeA_numerics/`, `contrib/encoded_chain_numerics/`), the length fluid plus payload (`contrib/mobile_cells/`), the smallest d₀ (`contrib/minimal_d0/`), an adversarial review (`contrib/redteam_construction/`), and the refined statement that should replace c–d (`contrib/refined_conjecture/`). These folders are drafts until their READMEs are filled in.

## How the construction sits in parameter space

At fixed d₀ the bonds h, modulo shifts and overall scale, form a compact space, so any sequence h_n with c_n → ∞ must accumulate at some h_∞. By the c-theorem, c is locally bounded near any point whose infrared is a CFT with finite c (a nearby flow passes close to that CFT before leaving it), and a gap is stable, so h_∞ must be gapless but not conformal. The construction is consistent with this [heuristic]. It is not a continuous family. Each n gives an open region of bonds (within those keeping the on-site Z₂) where c is constant, and moving from n to n + 1 changes the cell length, which is a first-order level crossing. In a fixed-alphabet version the natural form is one line h(δ), with δ selecting the cell length, along which c is a staircase whose steps shrink and accumulate at the threshold where cells become infinitely long. There the logical energy scales vanish and the dilute cell walls have z = 2. So unbounded c at fixed d is equivalent to c failing to be locally bounded near some gapless non-conformal point, and the construction says that happens at a commensurate–incommensurate threshold whose "particles" carry unboundedly complex internal states.

## Small d

- **d = 2.** c ≤ 1 is proved on the free-fermion slice, and holds given standard input on the U(1) and XYZ slices (codimensions 5, 5 and 6 of the 8-parameter family) [Astra] (`codex/work_C6_d2/`). The generic residual (a single Z₂, only antiunitary symmetry, or none) is open analytically. Numerically nothing exceeds c = 1. The last two flagged points, 220 and 1282, are incommensurate floating spirals with c = 1.00 and 1.02 from finite-entanglement scaling, calibrated on Ising, XX, XXZ and two decoupled XX chains [checked] (`contrib/d2_followup_220_1282/`).
- **d = 3 to 7.** Targeted symmetric families at d = 3 and 4 (C4-D3), and about 900 further points in spin-orbital, spin-3/2, SO(5)-vector and Sp(4) families, show nothing above d − 1 (a G₂ family at d = 7 was inconclusive). Every apparent excess was traced to rotor towers, dimer doublets, ferromagnetic multiplets or incommensurate sizes [numerical] (`contrib/counterexample_search/`). At d = 3, a qutrit chain with a three-fluid mixture of atoms and bound pairs gives entropy fits c ≈ 3.05 on rings up to L = 21, but its gaps also fit z = 2 [Astra, unchecked]. An independent check is running (`codex/work_C12/`, `contrib/qutrit_candidate/`).
- **d = 4.** Pure Sp(4)_k cannot be the IR of an Sp(4) chain on ℂ⁴ at any level, since it contains SU(2)_{10k}, whose even level cannot carry the spin-3/2 anomaly [checked]. The best target, a tricritical Ising singlet on top of the c = 3 spin-3/2 phase (c = 16/5), needs two tunings, and a map of the whole two-parameter SU(2) family at L = 8 and 12 shows no sign of it [Astra field theory, checked; numerical] (`contrib/tci_target_d4/`).

## Bounds that hold

- Quasi-free chains of any period: c ≤ log₂ d [proved] (REPORT §4, `codex/work_C1_astra/`).
- SU(d)-symmetric fundamental chains (h = a + bP), toroidal IR theories whose bosons come from on-site U(1)s, and WZW IR theories with level at most the embedding index [proved, given the IR identification] (REPORT §3.4, §4).
- Conditional bounds [checked]: c ≤ 1 + (6/π) K(d) log d if c-many conformal levels appear on rings of length at most K(d); c < d² − 1 if the IR stress tensor is built from on-site symmetry currents (`codex/work_C7/`). Also c ≤ ⌊r(h)/2⌋, with r(h) the operator Schmidt rank of the bond, for Gaussian modes resolved by commuting on-site charges [Astra].
- What cannot work: entanglement, thermodynamic and Lieb–Robinson arguments give at best c ≲ ℓ* log d, with ℓ* the length at which conformal behaviour sets in, which is too weak for large d (REPORT §3.1 to §3.2, `tools/`).

## Superseded

- REPORT's "Bottom line", §4 and §5 date from 2026-10-21, before the construction existed. They are kept as a record.
- The README's open problems aim at a proof. They are still open, but most of the current work goes into the construction.
- Approaches that failed: the ground-state Koo–Saleur ratio, single-scale Lin–McGreevy fixed points, Green's-function zero counting, Sp(4)₂ and the d = 4 tricritical point, exact product codes, symmetric logical targets (SU(N)₁, SO(2n)), and the conjecture c ≤ r(h), which a G1 cell declared as one big site appears to violate.

## Where help is most useful

1. A fixed-alphabet computation track, or a much cheaper programmable cell. This fixes d₀ and turns the construction into an explicit counterexample.
2. An independent check that the second-order logical chain of route A flows to c = n at finite coupling, including whether leakage out of the code only renormalizes the logical couplings.
3. Mobile cells: show that the length fluid plus payload stays gapless with c = n + 1 and a unique ground state at finite coupling.
4. An independent rederivation of G1, the C9 synthesis, and the joint Majorana compression.
5. The generic d = 2 residual, and the d = 3 qutrit hint.

## Glossary

- **Astra / Sol.** OpenAI models run closed-book through the Codex CLI: gpt-6-astra for almost everything, gpt-5.6-sol for two early comparison runs.
- **C1 to C12.** Problem statements given to Astra, in `specs/`, with outputs in `codex/C*.out` and work in `codex/work_C*/`. C1 free-fermion bound. C2 lattice Virasoro (Koo–Saleur). C3 Lin–McGreevy fixed points at d = 2. C4 a five-run campaign with sub-runs P1 and P2 (proof attempts) and D1 to D3 (searches for counterexamples). C5 sub-runs F1 and F2 (interacting generalizations of C1). C6 the d = 2 classification. C7 whether any bound c ≤ f(d) exists. C8 a long-running autonomous run (runs 1 to 6). C9 cell realization. C10 the computation track. C11 the infrared fate of route B. C12 a simple counterexample at moderate d.
- **F1 in `contrib/tci_target_d4/`** is a different, local spec (the glued tricritical point), not C5's F1. **A1** is a spec inside `contrib/counterexample_search/specs/`.
- **G1.** The second-order exchange cell of route A (`contrib/unbounded_c/G1/`).
- **Route A / route B.** Second-order and first-order coupling of cells (step 3).
- **Port.** One of a cell's two end sites, through which it touches its neighbours.
- **q, k_c, d₀, J_c.** Port dimension, alphabet of the computation track, total local dimension of the construction, clock coupling.
- **N, n, m, r.** Code dimension, number of complex fermion flavours (c = n), number of flavours in route B, number of Majorana channels (c = r/2).
- **K_FF(q).** The finite witness length of C8 run 3: for frustration-free selection, a zero-energy state at length K_FF(q) implies one at every length.
- **c_half, c_S, c_Cas, c_Δ.** Central-charge estimators: half-chain entanglement, ring entanglement fit, Casimir energy, and the Lin–McGreevy entropy combination (REPORT §3.7).
- **Verdict (b)** in `contrib/ir_control/` means standard but not rigorous; type (c) would be an obstruction.
- **ULS, LSM.** Uimin–Lai–Sutherland (the SU(d) swap chain), Lieb–Schultz–Mattis.

## Tags

- [proved] has a proof in the repository.
- [checked] means the coordinating Claude session reran or rederived the result itself. No human has refereed anything yet, and the stricter definition in [CONTRIBUTING.md](CONTRIBUTING.md) (checked by someone other than the author) applies to new contributions.
- [Astra] means produced by a closed-book Astra run and not independently rechecked.
- [numerical] means numerical, with sizes and estimators given in the linked folder.
- [read] or [recalled] are literature claims, from reading the source or from memory.

Everything here was produced with AI assistance (Claude Opus 5.5 coordinating, with gpt-6-astra runs) and has not been refereed.
