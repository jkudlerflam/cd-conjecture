# Where the c–d conjecture stands

Last updated 2026-10-28 UTC. Read this page first. It is the current state of the project in one place. [REPORT.md](REPORT.md) is the chronological record behind it, and where the two disagree, this page is the one to trust. Every claim carries a tag saying how it was established (see "Tags" at the end).

## The question

Latorre and Sierra ([arXiv:2403.17242](https://arxiv.org/abs/2403.17242)) consider H = Σ_i h_{i,i+1} on (ℂ^d)^⊗L, with one fixed Hermitian h on ℂ^d ⊗ ℂ^d, periodic boundary conditions and one-site translation invariance. They conjecture that if the low-energy spectrum is described by unitary CFT(s), the total central charge obeys c ≤ d − 1. The SU(d) swap chain (Uimin–Lai–Sutherland, "ULS") flows to SU(d)₁ with c = d − 1, so the claim is c_max(d) = d − 1. We also ask the weaker question of whether any bound c ≤ f(d) exists.

The working standard is theoretical physics (standard RG, controlled perturbation theory, numerics calibrated with at least two estimators), not mathematical rigor.

## A d=24 conditional counterexample candidate (2026-10-28)

[proved local ingredients; conditional infrared conclusion] The
[one-vacancy construction](contrib/counterexample_campaign/COUNTEREXAMPLE_D24.md)
is an explicit one-site invariant nearest-neighbour family on 24 states
with predicted c=24>d-1=23. The earlier statement that a counterexample
was verified is withdrawn after fresh cross-review. The finite bond,
frozen-sector bounds, band certificate and conditional vertex analysis
pass, but finite-J matching and a uniform dressed-defect comparison
are not yet established for the full microscopic chain. No fatal error
or actual refutation has been found. The current status is a credible
analytic candidate, not an established counterexample.

[proved and checked] The campaign replays C33's n=23 finite-g band
certificate, proves all-length/illegal-sector exclusion on compatible
material rings at nonzero port coupling, and checks the literal 576-by-576
bond, exact charge-sea projection and cyclic gluing. A local phase-model
representation handles the gapless vacancy mode. Physical translation
excludes dangerous fractional marker harmonics; bounded occupied-number
rounding handles arbitrary physical ring lengths. Tagged-defect numerics
through L=20 support the separate continuum per-defect comparison.

[read; proved application; standard RG input] Published local
Schrieffer-Wolff bounds provide an exact finite-order change of basis with
a volume-independent quasi-local remainder. A block-diagonal stability
argument gives a massive auxiliary threshold above the reference ground.
Taking g small, then lambda small controls the frozen local reduction;
the additional J/t hierarchy is consistent but still needs the uniform
vacancy matching and defect comparison.
This avoids requiring C40's stronger exact all-order kernel for the
physics inference. The
[source audit](contrib/counterexample_campaign/route_a/M2_REVIEW.md) and
[independent source review](contrib/counterexample_campaign/infrared/SOURCE_REVIEW.md)
state precisely what is proved and what remains conventional RG input.
The separate collaborator [C41 report](codex/work_C41/RESULT.md) agrees
with the harmonic and finite-ring rules and gives conditional acceptance,
retaining local matching and positive dressed defect costs as assumptions.

[fresh audit; corrected verdict] The
[adversarial audit](contrib/counterexample_audit_20261028/README.md)
reconstructs the literal small-cell effective Hamiltonian, regenerates the
full n=23 certificate, checks the finite-ring quotient and dangerous
quartic/marker branches, and strengthens quantum-interface and Gaussian
control. The second-order estimates do not by themselves establish the
finite-J common-bulk-subtracted defect comparison. Cross-review agrees
with the separate referee's conditional classification. New workers are
targeting that missing comparison directly. The finite results remain
valid; no microscopic c=24 measurement is claimed.

[proved partial reductions; numerical stress tests] The latest
[matching attack](contrib/counterexample_audit_20261028/matching_attack/README.md)
disproves an overly strong two-sided remainder bound that counted an
arbitrarily long wrong cell as one defect; the growing correction favors
selection and does not refute the candidate. It excludes a whole-ring giant
cell directly. A [quantum-selector test](contrib/counterexample_audit_20261028/selector_numerics/README.md)
finds no undercut in 264 exact diagonalizations at b=2,3 through L=20,
including noncommuting b=3 walkers. The
[new matching route](contrib/counterexample_audit_20261028/uniform_matching/README.md)
derives exact cancellations but leaves the required spatially weighted
infrared estimate open. The overall classification remains conditional.

## Lead claim: unbounded c at fixed local dimension (how to reproduce)

The current fixed-d candidate predicts c = n + 1 for every n at local dimension d₀ = 4464, which would exclude every bound c ≤ f(d). This is conditional on the stated Wilsonian input S, not an established refutation. The C44/C45 referee finds no fatal error and accepts the reduction to S; C46 explicitly leaves the full infrared step open. No human has checked it yet. To reproduce, read and rerun in this order:
1. `codex/work_C39/` (the 46-letter frustration-free core and its all-n lacunary schedule; `checks.py`, `certify.py n`), refereed in `contrib/referee_C39/`.
2. `codex/work_C42/` (explicit finite-D, finite-g bounds for every n).
3. `codex/work_C44/` (the mobile wrapper at d₀ = 4464, audited for every n; `bounds_checks.py`, `consistency_checks.py`), built on `codex/work_C8/` and `codex/work_C32/`.
4. `codex/work_C45/` (the infrared step and the coupling schedule; `checks.py`), refereed with C44 in `contrib/referee_C44_C45/`.
5. Supporting evidence at small n: `contrib/routeA_demonstrator/`, `contrib/mobile_cells/`, `contrib/routeA_n2_chains/`, `contrib/routeA_n34_chains/`, `contrib/margins_audit/`; the infrared literature in `contrib/lit_ir_step/`.
Verification done since the referee report (`contrib/verify_C44/`, REPORT §9.28): an independent rebuild of C44's bare bond, second-order source formula (11), transport overlap (9) and vanishing first-order port coupling agrees with C44 at n = 1, 2 and extends to n = 3. The second-order projected n = 1 chain at the scheduled point gives calibrated c = 2 from finite-entanglement scaling (1.97 and 2.00) and from interval entanglement (2.03 and 2.00), with a c = 1 control. This tests the Gaussian-stability half of S at n = 1, not the elimination of the massive sectors. Two cautions from it: C44's n = 2 traceless length–payload numbers are below double-precision roundoff, and at C44's diagnostic ratios the attractive length interaction would phase-separate the length fluid, so the s/δ → 0 ordering of the schedule is essential. Still open: a human read of C39 to C45. Where input S stands after the round of proof attempts C47 to C54 (REPORT §9.26): the massive-sector elimination now has a resummed, volume-uniform matching theorem on the enlarged cell model (C48, building on the collaborator's operator-valued theorem), and the step from a complete multichannel Gaussian limit to c = n + 1 is proved, with several velocities and moving zero modes (C51). Exact decoupling is ruled out for product-cell encodings with invertible logical hopping (C52), so S has to be proved perturbatively. No mechanism breaking S inside the coupling schedule was found (C54). What remains is the all-scale constructive RG for the matched retarded action (C50 isolates the missing comparison estimate and counterterm map), the temperature-sensitive scalar bound that yields the T² coefficient (C51), and the transfer to the physical fixed-length ensemble. Round two (REPORT §9.27) closed the ensemble transfer in thermodynamic form (C58: the physical fixed-length ring and the matched logical action have the same free energy up to activated errors), so S now reduces to two estimates about the logical fermion action alone: a normalized renormalized-tree bound for retarded endpoints (Lemma N of C55, the same estimate C56 needs for the multichannel comparison) and the scalar thermal bound P (C57). The self-contained write-up of both mobile forms (the d = 24 vacancy chain and the d₀ = 4464 wrapper chain) is `contrib/routeA_final/README.md`.

## Smallest candidate so far (2026-10-28, updated)

[proved finite ingredients; conditional infrared conclusion; not yet refereed] Astra C80 reduced the alphabet further to 12 letters with a certified n = 12 member (24 simple Fermi points, positive pair-sum margins): predicted c = 12 > 11 for the frozen chain at d = 12 (REPORT §9.43). For a mobile version, C82 repaired the vacancy lift with two coloured vacancies; on the 15-letter core this gives a certified conditional candidate c = 17 > 16 at d = 17, and C83 applied it to the 12-letter bond: two vacancy colours, d = 14, with a certified n = 13 member, giving the mobile conditional candidate c = 14 > 13 (REPORT §9.44). An independent referee (REPORT §9.45) found no fatal gap in C80, C82 or C83; the remaining serious items are the inherited infrared inputs. Smallest mobile candidate now: d = 13 with conditional c = 14 > 12 (C88, REPORT §9.51; eleven letters plus two vacancy colours, using C85's certified n = 13 member; mobile ground-state selection conditional on an extended dressed-defect comparison; not yet refereed). Since then C84 reached eleven letters with a certified frozen n = 11 member, predicted c = 11 > 10 at d = 11 (REPORT §9.46, not yet refereed). Earlier steps follow.

### d = 15 (C75)

[proved finite ingredients; conditional infrared conclusion; not yet refereed] Astra C75 reduced the route-A alphabet from 23 to 15 letters (`codex/work_C75/`, REPORT §9.41) with a certified member at n = 15: 30 simple Fermi points, distinct sea centres, exact U(1), positive pair-sum margins. The predicted infrared is c = 15 > d − 1 = 14 for the frozen chain, or c = 16 > 15 at d = 16 with one vacancy. It is conditional on the same steps as the d = 24 candidate. An independent referee (`contrib/referee_C75/`, REPORT §9.42) found no fatal gap and confirmed the finite claims; the frozen d = 15 chain stands on the same footing as the 23-letter frozen chain, while the d = 16 vacancy lift has a serious open gap (a domain wall can form next to a vacancy without penalty). A survey of simpler mechanisms at d = 3 to 12 (C76, §9.39) and a hard-core multi-species analysis (C77, §9.38) found nothing beyond c = d − 1.

## Summary

1. **No counterexample at small d.** Nothing above d − 1 at d = 2, in symmetric families at d = 3 to 7, or at the best d = 4 target. The d = 3 hint failed C12's own follow-up; the independent check is still pending in the recorded status. New [qubit proof reductions](contrib/d2_gap_reduction/README.md) close symmetric ferromagnetic rank-one exchange with arbitrary uniform fields (c ≤ 1/2), add a rigorous interacting gapped region, and exhaust additive on-site spin symmetries, leaving one staggered-field family beyond the known controlled classes.
2. **No general proof.** Every lattice quantity we know that bounds c grows like (a length scale) × log d, never like d. Proofs exist in restricted classes.
3. **A leading conditional candidate at d=24.** The one-vacancy construction predicts c=24 using the explicit 23-state chain and its n=23 certificate. Its finite ingredients and conditional infrared argument have been checked; uniform vacancy matching and dressed-sector stability remain open. The older sparse species-graph proposal remains a separate route toward an unbounded ratio c/d; its complete controller still needs composition and review.
4. **The stronger fixed-d candidate now has a finite spatial construction.** C39/C42 and the C44 wrapper supply the construction at d₀=4464; C45 closes a coupling schedule conditional on S. Its remaining central gap is microscopic infrared control, not the older missing spatial compiler. C46 proves partial thermal and finite-order results without establishing the full CFT limit. The sparse first-order family with growing d remains a separate construction.

## Direct attack on the remaining matching steps (2026-10-28)

[proved local reduction; infrared step open] The
[direct-attack campaign](contrib/direct_attack_20261028/README.md) gives an
explicit ultraviolet change of basis with a volume-independent remainder
carrying one spatial moment. A translation-covariant version improves the
cutoff power from 18 to 10. For s=J/t and a Gaussian cutoff delta=s^a,
0<a<1/10, its weighted remainder is o(J). This resolves a
technical ultraviolet estimate in the vacancy matching route, but not the
inhomogeneous low-frequency energy comparison or the CFT inference.

[proved finite-cell theorem; independently audited] The same campaign's
[massive-matching theorem](contrib/direct_attack_20261028/infrared/FINITE_CELL_MATCHING.md)
uses finite spectral embedding and a convergent determinant expansion to
bound the logical action correction linearly in leakage, with time and space
moments and no inverse logical-bandwidth factor. The hypotheses include a
fixed finite cell, a locally block-diagonal reference, and strong coefficient
norms. It does not establish the final logical RG/CFT limit.

[proved abstract extension; independently AI-audited] The new
[operator-valued thermal theorem](contrib/direct_attack_20261028/infrared/INFINITE_CELL_PROGRESS.md)
allows an infinite excited cell space with a finite code, positive gap and
finite heat trace. Its full thermal influence functional has a correction
linear in leakage, with spatial/time moments, source and scalar control,
and exponential convergence to the zero-temperature kernels. Constants are
independent of auxiliary cutoffs and the logical bandwidth. C46 supplies a
candidate heat bound on its specified legal rooted space; the complete
microscopic sector and physical-ensemble application remains conditional.
The arbitrary-degree retarded multichannel RG and central-charge inference
are still open. This strengthens the analytic case for pursuing C44/C45;
it is not a completed counterexample.

[proved ensemble comparisons; infrared conclusion open] The
[physical-ensemble campaign](contrib/thermal_ensemble_20261028/README.md)
extends pressure transfer to arbitrary length-dependent interacting bonds
by pinning a marker and controlling the root seam. A length-fugacity
pressure-zero theorem and a tilted C46 thermal comparison control the
all-length reference at fixed physical length. Its new
[positive completion](contrib/thermal_ensemble_20261028/ensemble/ILLEGAL_SECTORS.md)
adds helix pins only on illegal selector pairs, preserving the complete
legal Hamiltonian at d0=4464 while supplying the full all-string cell
gap and heat trace. An independently checked two-charge fermion
dictionary and [explicit parameter reduction](contrib/thermal_ensemble_20261028/MOBILE_MATCHING_APPLICATION.md)
instantiate massive matching and the physical-pressure zero for this
completed bond. The interacting logical RG/current or thermal-pressure
conclusion remains open; no central-charge violation is asserted.

[proved pressure bound; conformal interpretation conditional] The parallel
[frozen-chain campaign](contrib/decisive_attack_20261028/README.md) bounds
the full frozen-defect thermal contribution exponentially, strengthening
the d=23 bulk route. Its selected-cell RG/CFT limit remains open, and the
pressure result does not remove forced defects or translation degeneracy
on incompatible ring lengths.

## Parallel sparse first-order continuation

**Reversible single-head compilation, 2026-10-28:**
[digit_compiler](contrib/digit_compiler/README.md) supplies all signed
arithmetic operations (4,749 active controls, 28,494 digit rules) and serial
GET/SET/APPEND/LEN (97 controls, 207 rules). A proved spatial compiler runs
these finite tables with one reversible physical head and a finite history
alphabet. Actual lowered arithmetic takes 648,266 and 3,055,626 microsteps
in two checked runs, with exact complete cleanup. Explicit full-space
onsite/two-site permutation gates pass 63 autonomous compute/reset cycles.
Real-generator traces also validate 8,165 memory events and arithmetic for
p=2,3,5. The full interpreter's stack, heap aliases, allocation, instruction
dispatch and initialization still need a joined table. These component
results do not give a final high-channel dimension or a new CFT conclusion.

**Digit compilation and explicit stability budget, 2026-10-28:**
[digit_compilation](contrib/digit_compilation/README.md) removes the finite
parameter-reader-table gap: 96 collision entries and eight controls,
with exact output and autonomous cleanup. Eleven arithmetic kernels now
have 4,428 width-independent adjacent-cell rules. The actual p=2, 21-mode
circuit generator runs eight arithmetic opcodes through these rules:
1,175,426 local calls reproduce all 287,798 gates. This remains a hybrid
execution: memory, dispatch, operand placement and other arithmetic are
not thereby compiled. A separate proved Gaussian matrix inequality gives
an explicit relative stiffness threshold delta<5/13 protecting all
reference-dimension-three vertices, conditional on actual filling and
momentum margins. Full microscopic composition, the final alphabet count,
and local massive-sector elimination remain open.

**Autonomous continuation, 2026-10-28 UTC:**
[autonomous_controller](contrib/autonomous_controller/README.md) supplies one
3,200-entry reversible collision table (80 neutral symbols, 40 controls)
that starts with blank dynamic tracks, generates a binary-digit program, executes noncommuting
quantum gates, exposes both endpoints, and resets. All 125 lengths 4–128
pass complete cycles. Actual C8-clock refinement checks initialization and
every intermediate vertex; a full cycle fits for n>=9. The initialized
demonstration core is explicitly specified, but its program is not the
many-fermion sparse encoding. A frozen 3,115-instruction word machine now
generates that encoding exactly at p=2,3 in the checked abstract streams;
its arithmetic, arrays and journal still need local digit-level compilation.
Explicit neutral parameter queries can supply p from the already
p-dependent matrices, removing length-to-prime decoding. No final
high-channel constant C, violating d, or new infrared conclusion is claimed.

**Latest implementation work, 2026-10-28 UTC:**
[controller_instantiation](contrib/controller_instantiation/README.md) replaces
several algorithms with executable finite data. The encoding circuit now
lowers to 17 two-site gate families; an explicit preloaded dispatcher has
88 neutral symbols and 150 controls, and a separate blank-to-length counter
has 12 symbols and 16 controls. A complete d=240 one-mode fixture bond is
saved and independently checked, including its physical hopping projection.
An audit found and repaired an ambiguity in the neutral input-format check.
The full autonomous program generator, accepted program sector and endpoint
readout are still uncomposed; these component counts are not the full core's
constant C. No new violating example or thermodynamic determination is claimed.

[proved finite ingredients; conditional infrared conclusion] Let the species graph be the point-line incidence graph over F_p: q=2(p²+p+1) species and N=(p+1)(p²+p+1) logical fermions. Girth six forbids every flavour-changing quartic. Explicit endpoint maps realize all N hopping channels. A sparse reversible preparation, exact indexed algebraic rotations, a single-record moving head, and the existing fifteen-letter full-sweep clock give an initialized core of dimension at most Cq(1+log₂q). All-time FAIL penalties select exactly 2^N inputs; cyclic READY readout avoids an invalid absorbing halt. The complete mobile wrapper costs 96d_core+48 and preserves the species torus and marker number.

[physics-level inference] Choose k_a=π/4+aπ/[12(N+1)] and length density in [5/12,7/12]. Payload mismatches are at least π/[6(N+1)], mixed mismatches at least π/6, and all other non-Gaussian resonances have free dimension at least three. At each fixed p and cell length, choose a sufficiently small **fixed nonzero** physical coupling, including all visibility factors and small gaps. Standard local integration and RG then give N+1 compact bosons, hence c>d−1 for sufficiently large p. The [infrared audit](contrib/research_continuation/IR_AUDIT.md) states the input and the distinction from a global band Schrieffer–Wolff argument.

The new [research continuation](contrib/research_continuation/README.md) records the proofs, algorithms and reproducible finite checks. The constant C is not numerically optimized, and there is no expanded full-machine transition table or thermodynamic simulation of a violating chain. Three AI workers cross-reviewed the construction; no human review is claimed. A parity-compatible repair of C9 and a separate-left/right-port realization also strengthen route A without completing its spatial compiler.

## Route A (fixed d) and the d = 2 program, latest (coordinator, 2026-10-27 evening)

- **All-n margins (unrefereed).** Astra C39 builds a 46-letter frustration-free core whose bands have distinct centres and positive pair-sum margins for every n (a lacunary schedule, with exact rational certificates at n = 1, 2, 3); with the mobile wrapper this gives a fixed d₀ = 4464 candidate for unbounded c, conditional on the infrared step. An independent referee confirms the finite all-n claim (`contrib/referee_C39/`); what remains is an audit of the mobile wrapper on this core and the thermodynamic elimination of massive sectors. A one-vacancy mobile completion (collaborator's proposal, `contrib/counterexample_campaign/doping/`) would make the mobile versions much cheaper (d = 24 for the 23-letter core), and is under test.

- **The route-A mechanism works end to end in one explicit bond.** A ten-letter nearest-neighbour Hamiltonian forms cells spontaneously, its ports are silent at first order, cells couple only at second order, and the full microscopic chain measures c = 1.00 to 1.01 by finite-entanglement scaling at three couplings, with matching ring spectra and z = 1 (`contrib/routeA_demonstrator/`; frozen cells, one flavour). Going from one flavour to n is the open problem, with candidate cells in `contrib/routeA_spatial/` and `contrib/routeA_compile/`. The lead n-flavour candidate is now written out completely as one bond on ℂ²³ (`contrib/explicit_candidate/README.md`, checked against the code); its bands provably cross for every n (C32, which also gives a mobile version at d₀ = 2736), and distinct band centres with non-resonant margins are now certified for every n from 1 to 23 (C33, C37), so the frozen chain's finite-dimensional input at d = 23 is complete. The frozen chain is not a counterexample to the conjecture as stated, because of spontaneous translation breaking (red team), so the work now targets the mobile fixed-d₀ version, which needs the same certificates for every n.

- Route A's finite cell is now uniform: C19 gives n logical fermion flavours with exact parity, independent hopping phases and chemical potentials from second-order exchange, and only n + 1 port-accessible energies [Astra, checks rerun]. What is missing is a translation-invariant constant-alphabet controller. The cheapest frustration-free port interface provably fails (scalar source resolvents, C18 and C19 §9), so the interface must be non-frustration-free; C8 run 7 supplies a 15-letter initialized clock. Four Astra runs (C21 to C23, C8 run 8) and a microscopic demonstrator are on this now.
- Route A numerics: a symmetric target (SU(3)₁) drifts at larger coupling, and an XX target with Fermi points at ±k_F is gapped by a generated pairing term, so cells must keep an exact U(1) or displace each flavour's Fermi-sea centre (`contrib/routeA_numerics/`).
- The mobile-cell infrared claim (c = n + 1 on an open set) is now confirmed by a blind closed-book run (C20) as well as C13.
- d = 2: the two numerical hunts found nothing above c = 1 (multicritical junctions are gapped or non-conformal; the atoms-plus-pairs two-fluid mechanism is excluded with a U(1) and fails without one because a dense canted state intervenes; `contrib/d2_multicritical/`, `contrib/d2_two_fluid/`). The first analytic round (C15 to C17) found no general proof and no counterexample. The [2026-10-28 continuation](contrib/d2_gap_reduction/README.md) now proves an arbitrary-field ferromagnetic rank-one subclass using the Bravyi–Gosset–Liu–Wong gap theorem. Exact additive on-site charges reduce exhaustively to known XXZ/commuting cases or a(XX−YY)+uZZ+bZ. After a=1, b>0, its unresolved conformal candidates lie in −1<u<sqrt(1+b²/4), u≠0; the upper boundary is nonconformal and both outer regions are gapped by explicit proofs. The neutral-sector possibility for the full symmetry-free family remains open. See REPORT §9.17.
- The species-graph candidate above: Astra's review (C24) finds no fatal error and accepts the selection rule, ports, circuit and conditional Gaussian stability, but not yet an established refutation, because the complete controller and the local elimination for the mobile microscopic model are unverified. An independent Claude referee (`contrib/review_species_graph/`) reaches the same two serious gaps, finds no fatal error, and notes that with the components already explicit the first violation would occur only at d of order 10³⁵ or larger.
- Route A is now assembled end to end (`contrib/routeA_assembly/`): nothing found kills it, the scale hierarchy closes at each fixed n, and the decisive gap is a fixed-alphabet port interface that gives logical-operator-valued second-order exchange while keeping the code frustration-free (required by the mobile wrapper). Mobile route A costs d₀ above about 2000. Since then, fixed-alphabet cells with exact codes, exact parity and growing second-order hopping rank have been built (16 and 17 states, C21 and C22); the one missing ingredient is independent flavour fields and hopping phases (or an exact U(1)) at fixed alphabet, so the Fermi points are not pinned to 0 and π. Astra C28 is on it. See REPORT §9.19. A Claude agent has since built a fixed-alphabet interface with exact fermion number (mobile carriers through fermionic Bell-dimer ports, frustration-free, exact code, no pairing at any order) and a central-switch cell whose effective chain has c = n at n = 1 to 4, with displaced Fermi seas shown at n = 2; open are an equal-arm selector, displaced seas for larger n and the scales (REPORT §9.20, `contrib/routeA_compile/`). A second agent built an explicit computation-free chain with only 23 letters, selector included, with exact U(1) and parity, exact code, zero first order, and an effective chain whose bands all cross with distinct centres through n = 3 (n = 4 checked on the hopping matrix), verified by full ED of two and three cells (`contrib/routeA_spatial/`). It is a frozen-cell construction: making it mobile (equal payloads at two cell lengths) and the exponentially shrinking band scales are the open points.

## Earlier construction and the separate fixed-d program, step by step

The steps below describe the earlier route-A program and its Majorana alternative. Their missing fixed-alphabet spectral compiler remains an issue for that stronger program; the new sparse first-order core above uses a different exact preparation and readout.

The idea fits in a paragraph. Sites spontaneously group into long cells. Inside each cell a local "computation" prepares a degenerate code space, which acts as one big effective site of dimension N ≫ d. Neighbouring cells interact only through the single bond at their junction, so perturbation theory in that bond gives an effective nearest-neighbour chain of big sites with a logical Hamiltonian we choose. We choose a stable CFT with large c. The difficulty is doing all of this with one repeated h on a fixed ℂ^{d₀}.

1. **Spontaneous cells.** An explicit d = 4 chain in which a "head" symbol hops freely inside a segment, so a segment of length b costs the kinetic energy −2t cos(π/b). Tuning one chemical potential then selects any cell length b ≥ 5, with couplings independent of L and a gap for each b. [checked] `codex/work_C7/C7_result.md`.
2. **A code in each cell, with ports that stay put.** The cell runs a history-state computation (Feynman–Kitaev type) that prepares the code. If the computation never acts on the cell's two end sites (the "ports"), every port operator commutes with the circuit, so the code and the compression of port operators onto it are exact, and the effective coupling of two cells converges to its target as 1/J_c in the clock coupling. [checked] `contrib/programmable_cell/`.
3. **Coupling the cells.** There are two routes.
   - *Route A, second order (the fixed-d program).* Choose a code on which every traceless port operator compresses to zero, so first order vanishes, and engineer the cell's excited spectrum so the second-order bond is the logical bond we want. The G1 cell gives aI + bP_N (a swap) for every N and any port dimension, confirmed by full diagonalization [checked, rerun from scratch] (`contrib/unbounded_c/G1/`). It works with identical cells that have distinct left and right ports [Astra] (`codex/work_C8/`). Any Hermitian logical bond can be made through two qubit ports in a cell of dimension 4⌈(6N³ + 13N)/4⌉ [Astra] (`codex/work_C9/RESULT.md`). Route A makes d₀ = 8k_c, independent of N, where k_c is the alphabet of the computation track (step 5) [outline, unchecked].
   - *Route B, first order.* The ports compress to scaled anticommuting Majoranas, αγ_a on the left and αiγ_aP on the right, so the inter-cell bond is free Majorana hopping with c = r/2 for r channels [checked]. An exactly degenerate cell that does this exists [checked] (`contrib/programmable_cell/`). Route B is capped by bond rank at c ≲ q²/4 for port dimension q, so it cannot give unbounded c at fixed q. It can beat d₀ − 1, but only for q > 16k_c, which means d₀ ≈ 64k_c². Its flavour-symmetric version is unstable for m ≥ 7 channels (step 4).
4. **The logical target and its stability.** The target is an n-flavour chain of logical complex fermions with distinct, incommensurate Fermi momenta, giving c = n. Detuning the momenta removes every marginal non-Gaussian channel, so the target needs no exact symmetry beyond an on-site Z₂ that acts as logical fermion parity. Sidon-set momenta give margins δ₂ ≥ π/(6n) and δ₄ ≥ π/(576n²). A small iDMRG toy agrees: the symmetric chain dimerizes and the detuned one does not. [Astra, checked at the physics standard] `contrib/ir_control/`. With a code of dimension N = 2ⁿ this gives c = n, so route A still gives unbounded c. Symmetric targets do not work. SU(N)₁ (the ULS-type target G1 makes) and the SO(2n) Majorana point sit on separatrices and run away under small anisotropy, and an exact SU(N) cannot be protected through a small port for N > q at any order [Astra]. The symmetric route-B family turns marginally relevant at m = 7, and keeping U(1)^m flavour charges with distinct fillings repairs it [Astra] (`codex/work_C11/`). The infrared analysis was done for the encoded chain in general. That it applies to route A's second-order realization unchanged has not been checked separately.
5. **The computation track, which fixes d₀.** An explicit translation-invariant nearest-neighbour track runs any hard-coded gate program, with proved bounds on illegal sectors and on the selector margin. A small example has k_c = 24, q = 2, d₀ = 192 [Astra] (`codex/work_C10/RESULT.md`). But k_c grows with the cell size and the program length, so as written d₀ grows faster than c and gives no counterexample. **This remains missing for the route-A spectral construction.** It needs a fixed-alphabet compiler of the Gottesman–Irani or Kohler–Piddock–Bausch–Cubitt type, or a much cheaper cell. An alphabet count based on the smallest known translation-invariant history machine (Bausch–Cubitt–Ozols, 42 dimensions) gives d₀ ≈ 130 for route A with the C7 selector tensored on, and d₀ ≈ 45 if the machine's own clock can select the cell length [heuristic; the key lemma, that one static Bell-check term is the only term that must touch the port, is proved by Astra and checked] (`contrib/minimal_d0/`). These counts are for frozen cells. A minor point at the same step: at finite clock coupling the history states are ground states only up to O(1/J_c), which is acceptable if that error only renormalizes the logical couplings.
6. **Frozen versus mobile cells.** If the cell markers do not move, translation is spontaneously broken. Even on rings with b | L the ground state is b-fold degenerate (momenta 2πk/b), and on other ring lengths there are defects, with at least L ground states on large prime rings. That is a gapless theory with c = n on top of a broken translation, not a CFT in the strict sense the conjecture intends. The fix is to let the cells move. A mobile selector lets cell lengths fluctuate between b and b + 1 as a free c = 1 "length fluid" on every ring [Astra numerics: c ≈ 1.004 from entropy, 0.9999997 from the Casimir energy]. The logical payload rides on it, as in spin–charge separation, which should give a unique ground state and c = n + 1 on every ring. Two pieces of this are written. An "equalizer" makes the first-order logical bond independent of the lengths of the two adjacent cells, at local dimension 96q + 48 (`codex/work_C8/helix_equalizer.md`). With it, the complete first-order spectrum of length fluid plus any payload is known exactly: the two sectors decouple, with their momenta glued by one global constraint (`codex/work_C8/joint_fluid.md`; every eigenvalue matches direct diagonalization to 10⁻¹⁴, and the Casimir energy gives c = 2.0000 for an XX payload) [Astra; check scripts rerun by the coordinator]. At finite coupling, a bosonization analysis of every correction (residual length–payload terms, payload-dependent marker hopping, three-cell terms, leakage) finds no relevant mixed resonance: with exact marker conservation, exact payload parity, displaced and jointly non-resonant Fermi points (including the length fluid's) and a positive stiffness matrix, the mobile chain is an (n + 1)-boson phase, with an explicit open region at weak coupling [Astra, soluble controls rerun by the coordinator] (`codex/work_C13/RESULT.md`). A microscopic toy (d = 6, each cell carrying a spin payload) has a unique ground state on every ring tested including primes, z = 1 and c = 2 from three estimators, where the frozen version is degenerate [numerical: ring entropy, finite-entanglement scaling and Casimir energy with measured velocities] (`contrib/mobile_cells/`). At irrational densities the zero modes drift with L, so the spectrum converges to CFT towers in L-dependent charge sectors, as for the XXZ chain in a generic field. The sparse first-order continuation above now supplies an existential high-channel core. A fully expanded numerical alphabet and generic interacting numerics at sizes that meet the two-estimator standard are still absent; the conformal inference uses the stated physics inputs. A separate scope question remains, whether the conjecture could be read as requiring conformal behaviour only on compatible rings, but the mobile version would make it moot.
7. **Adversarial review.** A red-team agent reimplemented the decisive finite-dimensional checks and found no fatal flaw in them. It ranks the problems as frozen-cell ring degeneracy (step 6), the unstable SU(N)₁ target (dropped in favour of step 4's target), the unwritten computation track (step 5), and the absence of a rigorous Schrieffer–Wolff reduction (not needed at the physics standard). It also found that route B's cell, left untuned, puts every flavour's Fermi points at δk and π − δk, which leaves inter-flavour pair tunnelling resonant. Independent hopping phases per flavour remove that, as step 4's prescription requires (`contrib/redteam_construction/`).
8. **In progress.** Finite-coupling numerics on both routes (`contrib/routeA_numerics/`, `contrib/encoded_chain_numerics/`; route B with r = 2 already gives c = 1.00 ± 0.03 at λ/Δ up to 0.8), the mobile-cell toy (`contrib/mobile_cells/`), the smallest d₀ (`contrib/minimal_d0/`), and the refined statement that should replace c–d (`contrib/refined_conjecture/`). These folders are drafts until their READMEs are filled in.

## How the construction sits in parameter space

At fixed d₀ the bonds h, modulo shifts and overall scale, form a compact space, so any sequence h_n with c_n → ∞ must accumulate at some h_∞. By the c-theorem, c is locally bounded near any point whose infrared is a CFT with finite c (a nearby flow passes close to that CFT before leaving it), provided the crossover scale stays uniformly controlled nearby, and a gap is stable, so h_∞ must be gapless but not conformal, or the crossover scale must diverge there. The construction is consistent with this [heuristic]. It is not a continuous family. Each n gives an open region of bonds (within those keeping the on-site Z₂) where c is constant, and moving from n to n + 1 changes the cell length, which is a first-order level crossing. In a fixed-alphabet version the natural form is one line h(δ), with δ selecting the cell length, along which c is a staircase whose steps shrink and accumulate at the threshold where cells become infinitely long. There the logical energy scales vanish and the dilute cell walls have z = 2. So unbounded c at fixed d is equivalent to c failing to be locally bounded near some gapless non-conformal point, and the construction says that happens at a commensurate–incommensurate threshold whose "particles" carry unboundedly complex internal states.

## Small d

- **d = 2.** c ≤ 1 is proved on the free-fermion slice, and holds given standard input on the U(1) and XYZ slices (codimensions 5, 5 and 6 of the 8-parameter family) [Astra] (`codex/work_C6_d2/`). The generic residual (a single Z₂, only antiunitary symmetry, or none) is open analytically. Numerically nothing exceeds c = 1. The last two flagged points, 220 and 1282, are incommensurate floating spirals with c = 1.00 and 1.02 from finite-entanglement scaling, calibrated on Ising, XX, XXZ and two decoupled XX chains [checked] (`contrib/d2_followup_220_1282/`).
- **d = 3 to 7.** Targeted symmetric families at d = 3 and 4 (C4-D3), and about 900 further points in spin-orbital, spin-3/2, SO(5)-vector and Sp(4) families, show nothing above d − 1 (a G₂ family at d = 7 was inconclusive). Every apparent excess was traced to rotor towers, dimer doublets, ferromagnetic multiplets or incommensurate sizes [numerical] (`contrib/counterexample_search/`). At d = 3, a qutrit chain designed as a three-fluid mixture of atoms and bound pairs gave entropy fits c ≈ 3.05 on rings up to L = 21. Astra's own follow-up finds that the bound-pair minima the mechanism needs disappear in a saddle-node at δ ≈ 0.570, below the δ = 0.75 and 0.9 points that gave the large estimates, and that those estimates fail a size-scaling check, so it is not a counterexample [Astra] (`codex/work_C12/RESULT.md`). An independent check is still running (`contrib/qutrit_candidate/`).
- **d = 4.** Pure Sp(4)_k cannot be the IR of an Sp(4) chain on ℂ⁴ at any level, since it contains SU(2)_{10k}, whose even level cannot carry the spin-3/2 anomaly [checked]. The best target, a tricritical Ising singlet on top of the c = 3 spin-3/2 phase (c = 16/5), needs two tunings, and a map of the whole two-parameter SU(2) family at L = 8 and 12 shows no sign of it [Astra field theory, checked; numerical] (`contrib/tci_target_d4/`).

## Bounds and proposed refinements

- Quasi-free chains of any period: c ≤ log₂ d [proved] (REPORT §4, `codex/work_C1_astra/`).
- SU(d)-symmetric fundamental chains (h = a + bP), toroidal IR theories whose bosons come from on-site U(1)s, and WZW IR theories with level at most the embedding index [proved, given the IR identification] (REPORT §3.4, §4).
- Conditional bounds [checked]: c ≤ 1 + (6/π) K(d) log d if c-many conformal levels appear on rings of length at most K(d); c < d² − 1 if the IR stress tensor is built from on-site symmetry currents (`codex/work_C7/`). Also c ≤ ⌊r(h)/2⌋, with r(h) the operator Schmidt rank of the bond, for Gaussian modes resolved by commuting on-site charges [Astra].
- A proposed replacement for c–d, the Casimir bound Σ_a c_a v_a ≤ 2π‖h‖ (`contrib/refined_conjecture/`), is false: the integrable SU(N) chains with the rank-k symmetric representation on each site have Σ c v/(2π‖h‖) = k(N² − 1)/(N(N + k)H_k), above 1 for k = 2 and N ≥ 7 (d = 28; for example SU(8) with k = 2 has d = 36, c = 12.6, v = π/4, ‖h‖ = 3/2 and ratio 1.05), and growing like N/(2 log N) for k = N [Astra; velocity and bond spectrum checked by the coordinator; confirmed directly at N = 7 and 8 from finite-size Bethe roots on rings up to L = 120, `contrib/casimir_survey/`] (`codex/work_C14/`). It does hold for all quadratic fermion chains (A ≤ π‖h‖). Whether some bound A ≤ C‖h‖ log d survives is open; every chain computed so far obeys A ≤ π‖h‖ log₂ d, with equality only for the Heisenberg chain [numerical] (`contrib/casimir_survey/`).
- What cannot work: entanglement, thermodynamic and Lieb–Robinson arguments give at best c ≲ ℓ* log d, with ℓ* the length at which conformal behaviour sets in, which is too weak for large d (REPORT §3.1 to §3.2, `tools/`).

Additional material from "Prove c d conjecture" is in [contrib/chat_followup/](contrib/chat_followup/README.md): the species-torus selection rule, a direct-sum massive-auxiliary embedding, exact extensions of controlled classes, and a detailed qubit catalog. The chat's proposed c = r², d = O(r log r) family remains conditional and is distinct from unbounded c at fixed d. Its compiler and general infrared transfer are not established by the imported notes; the newer sparse-graph continuation above supplies a different finite-core construction and explicitly states its infrared inputs.

## Superseded

- REPORT's "Bottom line", §4 and §5 date from 2026-10-21, before the construction existed. They are kept as a record.
- The README's open problems aim at a proof. They are still open, but most of the current work goes into the construction.
- Approaches that failed: the ground-state Koo–Saleur ratio, single-scale Lin–McGreevy fixed points, Green's-function zero counting, Sp(4)₂ and the d = 4 tricritical point, exact product codes, symmetric logical targets (SU(N)₁, SO(2n)), and the conjecture c ≤ r(h), which a G1 cell declared as one big site appears to violate.

## Where help is most useful

1. Independently audit the new sparse first-order compiler, initialized-core composition and local infrared argument. Expand its finite alphabet to obtain a concrete first violating d. For the separate fixed-d goal, realize route A's excited-spectrum construction spatially.
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
