# The c–d conjecture

Latorre and Sierra ([arXiv:2403.17242](https://arxiv.org/abs/2403.17242), J. Stat. Mech. (2024) 113103) conjecture that a critical, translation-invariant, nearest-neighbour Hermitian chain on (ℂ^d)^⊗L flows to unitary CFT(s) with total central charge c ≤ d − 1. The bound is saturated by the SU(d) Uimin–Lai–Sutherland chain h = P (the swap), whose IR is SU(d)₁ with c = d − 1, so the conjecture is equivalent to c_max(d) = d − 1.

This repository is a shared notebook for proving or refuting it. Everything we know goes here, and everyone with access pushes to and pulls from it as the work develops. Maintained by Jonah Kudler-Flam.

## How to use this repository

### 1. Getting oriented

Read [STATUS.md](STATUS.md) first. It lists the main results and what is open, on one short page. The earlier long version, with the main construction written out step by step and a glossary of the run labels (C1 to C12, G1, route A and B, and so on), is kept in [STATUS_DETAIL.md](STATUS_DETAIL.md). [REPORT.md](REPORT.md) is the chronological record behind it, organized as follows.

- §1 states the conjecture precisely.
- §2 surveys the literature.
- §3 analyzes the proof routes.
- §4 and §5 give the verdict and the next calculations as of 2026-10-21. Both are superseded by STATUS.md.
- §6 records provenance.
- §7 and §8 give the results of the first model campaigns.
- §9 is the log of the campaign of 2026-10-27: the cell construction, the small-d searches and the infrared analysis.

[lit.md](lit.md) is the literature map. Every claim in this repository carries a tag saying how it was established, so check the tag before relying on it.

### 2. Getting a copy and staying in sync

```bash
git clone https://github.com/jkudlerflam/cd-conjecture.git
cd cd-conjecture
```

Before you start any work, pull with rebase so your work sits on top of everyone else's:

```bash
git pull --rebase
```

At meaningful checkpoints (roughly every 30 minutes during longer sessions)
and before finishing, commit in small pieces, pull again, and push:

```bash
git add <files>
git commit -m "Short description of what changed"
git pull --rebase
git push
```

If the push is rejected because someone else pushed first, run `git pull --rebase` again, resolve any conflicts, and push. Everyone with write access works directly on `main`. If you don't have write access, fork the repository and open a pull request.

Every project chat follows this routine under [AGENTS.md](AGENTS.md). Preserve
other chats' edits and verify that the branch is synchronized before finishing.

### 3. Adding a result

Each new piece of work gets its own folder, `contrib/<short-name>/`. Copy [contrib/_template/README.md](contrib/_template/README.md) into it and fill it in.

- State the claim with every hypothesis.
- Give its status: proved, numerical, heuristic, or refuted.
- Say who checked it and how, including any AI assistance and the model used.
- Say how to reproduce any numerics.

Put the code and data in the same folder. If the result changes the overall picture, also update "Where things stand" below and add a short section to REPORT.md.

### 4. Discussing

Use [GitHub issues](https://github.com/jkudlerflam/cd-conjecture/issues) for open problems, proposed attacks, and claims that need an independent check. Keep one question per issue and link the files or commits it concerns. The open problems listed below are good starting points.

### 5. Running the existing code

The code is plain Python and needs numpy and scipy unless a folder says otherwise. Each `codex/work_*/` folder holds the scripts, data and write-up of one run, and its report or README says how the scripts were run and how long they took. `codex/work_C4_D2/` and `codex/work_C4_D3/` also contain SLURM scripts for larger scans. To rerun a whole calculation from scratch, use the matching problem statement in `specs/`.

### 6. Conventions

These are spelled out in [CONTRIBUTING.md](CONTRIBUTING.md). In brief, tag every claim, name which central charge you mean (Virasoro c, or c_eff for non-compact or non-unitary limits), and keep files under about 50 MB. Never commit passwords, tokens or other private data.

## Where things stand

**Qubit update (2026-11-23):** The earlier c=3/2 interpretation of the
simple chirality chain was withdrawn on November 19. Its proposed c
approximately two region now has both infinite-chain and open-chain
entropy evidence, with the cleanest new point at delta=3/16,h=1.783.
The latest correlator estimates support the two-valley interpretation
there, but do not settle long-distance locking or establish a microscopic
CFT. The older large-blockade and phase-1.53 candidates remain conditional.
The [current review](contrib/qubit_c2_status_20261123/README.md) compares
the evidence, convergence limits and coefficient simplicity, and
reproduces the central-window fits from the committed tables.

**Stability round completed (2026-11-18):** The modified d=9 family now
has finite-order local interaction bounds independent of chain length
and exact nonzero-port ground-sector selection on compatible rings.
Qubit work computes a connected finite-density quartet response and
finds that saved carrier shifts disagree in sign with bare-band estimates;
the continuum matching needed to settle c=2 remains missing. The original
d=10 mobile model gains a uniform bound on one half-occupancy response
channel, with full thermodynamic selection still open. All four GPT-6 Astra
Max chats completed. No new microscopic c is established and the
coefficient hierarchies are not improved. The
[results and remaining gaps](contrib/research_coordination_20261116/STABILITY_RESULTS_20261118.md)
also distinguish the collaborator's separate, simpler c=3/2 numerical
qubit candidate. Long production remains stopped.

**Focused round completed (2026-11-18):** A new moment-sign proof gives
an arbitrary-n leading quadratic band family at physical d=9, using a
modified biased clock and L1 readout. The original d=10 mobile model gains
integer-certified finite sector comparisons and an exact conserved-charge
field improvement. Qubit work computes finite interacting quartet mixing
and proves crystal preemption in the tested number-conserving three-valley
leading model. No new microscopic c or unbounded-c CFT is established, and
the small-d coefficient hierarchy remains severe. All four GPT-6 Astra Max
chats have completed; [results and scope](contrib/research_coordination_20261116/FOCUSED_RESULTS_20261118.md).

**Ten-state family located (2026-11-18):** The reported d=10 candidate
family is [contrib/markerless_d9_d10/](contrib/markerless_d9_d10/README.md),
dated 2026-11-12 and published today. It removes the marker letters of the
d=11 and d=13 bonds: the predictions are c=n on nine-state frozen cells
and c=n+1 after adding one uncoloured vacancy, at d=10. The legal sector
is exactly invariant, eliminating finite-K leakage. Logical band data
are certified at n=9,10,11,12,14,16,20, not for all n. The infrared step
remains conditional, and the corrected mobile bond still needs a
quantitative J/t budget and multi-cell sector checks. This is now the
priority fixed-d family. An AI referee report is in the folder.

**Earlier completed round (2026-11-18):** The separate single-rail nine-state three-flavour
bond has six certified simple Fermi roots, but its generated interactions
destabilize the tested three-mode Gaussian state in the coupled leading
RG flow. The eleven-state loop gives full pair-sum separation for a new
analytic n=1,...,4 family and an explicit remaining all-n nonvanishing
criterion, without adding states. The qubit audit permits quartet and
mixed locking, corrects C180's drift covariance, and identifies the
momentum-pi pair exponent as a useful future diagnostic. No microscopic
central charge or unbounded-c theorem follows. See the
[completed-round summary](contrib/research_coordination_20261116/ROUND5_RESULTS_20261118.md).

**Earlier source audit (2026-11-18, before publication above):** Direct
GitHub branch, issue, comment and history checks found no d=10 arbitrary-c
source then; that source search is now superseded.
New work supplies a literal nine-state two-flavour resonant seed with a
favorable leading RG flow, and transfers clock, zero-port selector and
simple-band arguments to eleven states for arbitrary finite n. The first
does not establish arbitrary flavours; the second retains exact Cooper
resonances in its tree and an open all-n loop-splitting problem. Neither
establishes a microscopic CFT or unbounded c. The new tasks address those
gaps, plus qubit mode-locking criteria without restarting DMRG. See the
[source check and follow-ups](contrib/research_coordination_20261116/GITHUB_D10_AUDIT_20261118.md).

**Smaller-d priority and source audit (2026-11-17, historical):** The
reported d=10 unbounded-c construction had not yet been located; the d=10
models identified then were one- or two-flavour demonstrations. The
markerless contribution above supersedes this source-search conclusion.
Existing smaller candidates include d=11 frozen and d=13 mobile bonds
with certified finite members, and a d=12 mobile proposal with additional
heuristic dynamics. These are distinct from the recent d=20/22 all-n
program. The new round prioritizes <=10-state mechanisms and transfer
to the smallest concrete alphabet while preserving all CFT assumptions.
See the [source ledger and active assignments](contrib/research_coordination_20261116/SMALLEST_D_FOLLOWUP_20261117.md).

**V=4 size-scaling follow-up (2026-11-17):** The
[bulk-phase and entropy study](contrib/v4_bulk_diagnosis_20261117/README.md)
reaches L=144 in both parities and reproduces the L=96 state independently.
The latest interrupted chi=320 entropy slope is 2.57 for L=96 to 144,
versus 2.82 for L=48 to 96; size/window drift and incomplete accuracy
checks prevent identifying a thermodynamic c. Density/bond and secondary
valley-channel diagnostics do not establish persistent bulk order.
The user requested a battery-saving stop; all jobs are stopped and
checkpoints are preserved.

**Mobile length-selection test (2026-11-17):** A
[tail-sensitive comparison](contrib/mobile_length_selection_20261117/README.md)
reduces legal vacuum length words to explicit scalar vacancy problems at
J/t<=1/80. For every length, the second-order clock correction favors
the intended word in the weak marker-changing test. The scalar marker
response does not cancel; uniform selection and retuning remain unresolved.
Bounded checks retain the frozen selector and need no new penalty or
smaller J/t. The full construction's coefficient hierarchy remains.

**Two-term qubit diagnostic (2026-11-17):** The
[line `2XY+gZZ`](contrib/two_term_qubit_phase_20261117/README.md)
has an exact obstruction to the standard boost third charge, including
all range-two repairs, and a controlled large-|g| bulk-gap exclusion
with a non-explicit threshold. Near the flat point, a slow hopping is
generated at first order; this is not the tuned cancellation line.
Sector-resolved rings through L=20 and an audit of existing L=128 data
favor c=1 at g=1/2,1; stronger-coupling controls favor density order.
No all-g central-charge bound or new counterexample is established.

**Simpler qubit candidates (2026-11-17):** A
[46-run calibrated campaign](contrib/simple_qubit_candidates_20261117/README.md)
closes the missing local chi check for the old V=4 point at L=96.
Its common-basis form `(9/16)XY+(7/16)YX+ZZ-(7/4)Z` has four
monomials and coefficient spread 4 (six monomials and spread 28 in
the original basis). It is the simplest retained follow-up lead;
strong density modulation and unresolved size scaling prevent a CFT
claim. V=8 reaches L=128,chi=256 but still fails the chi-step gate;
small parameter probes show substantial density/sector changes.
No robust c>1 phase or new counterexample is established.

**Mobile vacuum-clock control (2026-11-17):** A
[coherent-gradient comparison](contrib/mobile_vacuum_clock_20261117/README.md)
proves a positive compressed clock-pin cost and bounds vacuum-clock
excitation weight with the terminal selector and vacancy dynamics retained.
Its sufficient J/t<=1/80 allowance is independent of the charged mass D
and cell length, in the stated tilted ensemble and its canonical
minimizing blocks. A parameter audit distinguishes genuine simplifications
from proof-induced hierarchies. Dressed dynamical matching, wrong-length
and illegal-sector selection, and the CFT remain open.

**Minimal qubit terms (2026-11-17):** The
[four-face classification](contrib/qubit_minimal_terms_20261117/README.md)
of `a(XX-YY)+b(XY-YX)+cZZ+hZ` excludes central charge above one when
pairing or interaction is absent, using the XXZ and free classifications.
The staggered-field and zero-field interacting faces retain explicit
unresolved regions. The latter is an asymmetric Ashkin–Teller model,
with a controlled C=1 neighbourhood of its XXZ limits. A common rotation
reduces six Pauli monomials to four but leaves all four independent
couplings; it cannot remove the reference candidate's invariant hierarchy.
The two-term subfamily `2a XY+cZZ` is a simple residual test, not a newly
established counterexample.

**First d=22 mobile auxiliary estimate (2026-11-16):** The
[mobile matching contribution](contrib/small_d_mobile_matching_20261116/README.md)
derives the literal contact-weighted representation and proves a positive
charged-rail pin cost with vacancy dynamics retained, in a positive
length-tilt ensemble and its canonical minimizing blocks. It also bounds
charged admixture at nonzero port and handles physical cyclic gluing with
one vanishing global energy allowance. Vacuum-clock matching, arbitrary
prescribed densities, wrong-length selection, and the interacting CFT
remain open; a separate canonical scalar-impurity lemma is not full matching.

**Exact qubit resonance hierarchy (2026-11-16):** The
[all-period parent analysis](contrib/qubit_resonance_hierarchy_20261116/README.md)
constructs exact nearest-neighbour qubit parents for every resonance
period m. On compatible rings their ground space has dimension L+1,
its singular limit is a span of single-residue Dicke multiplets, and
their gap has an explicit O(L^-2) upper bound. Exact local calculations
and larger complete-space checks find eighth-order motion at m=5,6.
These exact parents do not supply m positive-velocity conformal modes;
the general detuned phase problem and uniform local bound remain open.

**Resonant qubit detunings (2026-11-16):** The
[many-body and scattering analysis](contrib/qubit_resonant_liquids_20261116/README.md)
shows that the exceptional leading models sit at ferromagnetic density
thresholds. Their natural pair-amplitude detuning also creates particles
and retains a positive parent structure; treating it as repulsion alone
is incorrect. A separate hopping detuning mixes residue sectors and
passes the stated few-particle tests. An explicit fixed microscopic test
point is supplied, but its interacting infrared phase remains unverified.

**Coordinated follow-up (2026-11-16):** An
[independent effective-model check](contrib/research_coordination_20261116/README.md)
shows that the new three-/four-valley qubit branches, at their specified
sixth-order threshold, are positive configuration-graph Laplacians with
ferromagnetic single-residue sectors and exactly L+1 ground states on
compatible rings. Their band minima alone do not supply three or four
conformal modes. Follow-up work studies admissible detunings, the general
resonance hierarchy, and mobile matching for the d=22 family. The original
large-chi qubit convergence test remains blocked on production resources.

**Qubit cancellation classification (2026-11-16):** The
[new perturbative classification](contrib/qubit_cancellation_classification_20261116/README.md)
identifies both dilute and dense classical cancellation branches, gives
their complete fourth-order operators, and finds a sixth-order obstruction
to scalar cancellation on the dense branch. Two exceptional dilute paths
produce distance-four and distance-three hopping; the earlier parent paths
instead extend the blockade through distance eight. Exact symbolic and
full-sector finite-ring checks are supplied. Neither a two-valley argument
nor parent tangency proves a finite bound on total interacting c; universal
local boundedness remains open.

**Frozen d=20 defect stability (2026-11-16):** An
[exact cutting comparison](contrib/small_d_defect_stability_20261116/README.md)
now selects the all-length-n ground sectors at strictly nonzero port
coupling on rings divisible by n+1. It gives a positive cost per wrong
cell and per forbidden bond, uniformly in volume and wrong-cell length,
with finite grammar penalty and the required positive readout. No
chemical-potential retuning is needed. The interacting CFT identification,
incompatible-ring qualification, and d=22 mobile matching and defect
questions remain; this is not an unbounded-c counterexample.

**Smaller fixed-d candidate (2026-11-16):** The
[new 22-state mobile family](contrib/smaller_d_unbounded_20261116/README.md)
has an explicit bond, an exact all-length cell selector and an arbitrary-n
leading-band argument, with proposed c=n+1. Its frozen precursor has d=20.
The infrared conclusion remains conditional on mobile ground-sector stability,
vacancy matching and multichannel conformal control; it has not had
independent review. The note distinguishes this unbounded-c target from
the smaller finite-c candidates and records the
[compactness argument and its missing hypothesis](contrib/smaller_d_unbounded_20261116/COMPACTNESS.md).

**Small-d continuation (2026-11-16):** The
[stacking and vacancy note](contrib/small_d_lifting_20261116/README.md)
proposes attainable c=2,3,4 at d=2,3,4 while leaving their maximal values
open. Two onsite copies give d=4, c=4 exactly if the qubit c=2 seed is
confirmed. An explicit d=3 vacancy lift targets c=3, with additional
Gaussian stability and finite-coupling matching hypotheses. Exact closure
proofs, an entropy identity, and finite-matrix checks are supplied; no new
thermodynamic central charge is measured.

**Proposed revised qubit conjecture (2026-11-16):** For a fixed Hermitian
nearest-neighbour qubit bond repeated with **one-site translation symmetry**,
the [new working conjecture](contrib/d2_updated_conjecture_20261116/README.md)
is total unitary conformal central charge c≤2, conjecturally sharp but now
explicitly assessed as a tentative upper-bound guess. This is
not a proved upper bound or a promotion of the numerical c≈2 candidate to
an established counterexample. The note gives the exact microscopic-U(1)
reduction to c≤1, checks the two-valley mechanism, and tests repeated
blockade tuning, including two exact nonconformal frustration-free families.
Additional composite or neutral critical sectors remain the principal gap
in a general argument.

**Independent continuation (2026-11-09):** [new checked contributions](contrib/counterexample_push_20261109/README.md)
give an exact binding boundary and pair wavefunction for the Fibonacci
qubit candidate, identify its leading omitted number-changing term, and
improve the d=13 legal projected pin bound from J/4096 to J/135.
The qubit result agrees with C175's new two-body obstruction; its
finite-density c remains open. The [route audit](contrib/route_audit_20261109/README.md)
ranks d=13 as the strongest mature construction and keeps its
model-specific matching and relative defect hypotheses explicit.
The same continuation supplies a [new transverse-DM qubit family](contrib/counterexample_push_20261109/FACILITATED_ROUTE.md)
whose complete effective Hamiltonian weakens binding perturbatively;
the apparent small-ring crossing is strongly size dependent and is not
a counterexample claim.

**d=2 update (2026-11-10):** More accurate full-spin reruns recover c≈2
at the phase-1.53 point through L=256. The earlier unconverged L=128
downturn does not reproduce. Six tuned bonds give ground-branch entropy
estimates 1.94–2.03 from 128 to 256 sites, with calibrated controls.
The L=256 entropy still changes too much with bond dimension, so no
registered pass or established d=2 counterexample is declared. L=512
runs are recorded as in progress; infinite-chain microscopic convergence
remains open. See [current status](STATUS.md) and the [ground-sector
selection audit](contrib/d2_status_audit_20261110/README.md). The transverse-DM
points previously tested read near c=1.

**Conditional large-d constructions (scope clarified 2026-11-17).** Explicit
chains predict c=14 at d=13 and unbounded c at fixed d=4464, subject to
their stated matching, sector-selection and infrared inputs. The finite
certificates do not by themselves establish those microscopic CFTs.
The proposed replacement, Σ_a c_a v_a ≤ π M log₂ d, survives the recorded
tests and is proved in the stated classes, but remains open in general.
See [STATUS.md](STATUS.md) and each construction's hypothesis ledger.

A d=24 conditional counterexample candidate is recorded with independent AI reviews and a [fresh adversarial audit](contrib/counterexample_audit_20261028/README.md). The earlier verified-counterexample wording is withdrawn: finite-J matching and a uniform dressed-defect comparison remain open. The working standard is theoretical physics (standard RG, controlled perturbation theory, calibrated numerics), not mathematical rigor.

- **Small-alphabet candidates with finite certificates (2026-10-29): d = 11 frozen with predicted c = 11, and d = 13 mobile with predicted c = 14.** Explicit route-A bonds with certified finite members, refereed, and red-teamed for d = 13; both rest on the same inherited infrared inputs as the d = 24 candidate. Full bonds and status in [contrib/small_d_note/](contrib/small_d_note/note.pdf) and STATUS.md. The logical vacancy model gives calibrated c = 2, 2.9 ± 0.1 and 3.98 ± 0.05 at n = 1, 2, 3 (REPORT §9.34, §9.55, §9.67). The smaller d=12 mobile proposal is separately recorded with additional heuristic dynamics.

- **Unbounded c at fixed d0 = 4464 (the C39/C44/C45 wrapper).** c = n + 1 for every n, conditional on the Wilsonian input S. Several pieces of S are now proved and refereed (massive matching C48, ensemble transfer C58, c from a Gaussian limit C51); the open core is the constructive RG (Lemma N and lemma TP). The full microscopic demonstrator flows to its matched XXZ Luttinger liquid numerically (REPORT §9.58). Expert summary: [contrib/expert_note/note.pdf](contrib/expert_note/note.pdf).

- **d = 2.** No c > 1 found in about 103,000 sampled qubit chains of the residual strata, on the staggered branch to L = 256, or inside floating phases (REPORT §9.37, §9.40, §9.62); the analytic covering argument is not complete.

- **d=24, predicted c=24: a one-vacancy candidate.** Add one mobile vacancy to the explicit 23-state chain, with its original bond acting on occupied pairs. The local phase-model representation, certified 23-band reference, nonzero-coupling sector bounds and physical momentum rules support 23 payload bosons plus one charge boson, conditional on finite-coupling matching and defect stability. Published local perturbation bounds control the original cell machinery. The conclusion uses conventional Gaussian phase stability; no microscopic c=24 measurement or rigorous CFT-limit theorem is claimed. Read [the construction and independent reviews](contrib/counterexample_campaign/COUNTEREXAMPLE_D24.md).

- **Independent vacancy continuation.** A proved density window removes exact half-filling tuning from the proposed d=24 construction. New reviews check the local perturbation argument and C39's all-n finite input, including a new exact four-mode certificate. A full microscopic d=3 control has entropy and independently calibrated energy scaling consistent with c=2; this tests the vacancy mechanism, not the c=24 chain. See [the proofs, reviews and numerical data](contrib/vacancy_continuation/README.md).

- **Small d: no counterexample.** Nothing above d − 1 at d = 2, in symmetric families at d = 3 to 7, or at the best d = 4 target. The d = 3 hint failed its own follow-up (C12); an independent check remains pending in the recorded status.
- **New d = 2 proof reductions.** Arbitrary uniform fields in the symmetric ferromagnetic rank-one class give at most c = 1/2, using a recent rigorous gap theorem and the exact Ising solution. An exhaustive classification of additive on-site spin charges leaves one additional two-parameter staggered-field family; exact gap arguments exclude both outer regions. The full qubit conjecture remains open. See [the proofs and remaining region](contrib/d2_gap_reduction/README.md).
- **No general proof.** Every lattice quantity we know that bounds c grows like (a length scale) × log d. Proofs exist for quasi-free chains (c ≤ log₂ d), SU(d)-symmetric bonds, and some classes of IR theories, and there are conditional bounds.
- **Large d: a candidate physics-level refutation.** The new sparse species-graph family has predicted c=(p+1)(p²+p+1)+1 and microscopic d=O(p² log p), so c/d diverges. Exact finite constructions now supply the ports, reversible circuit, O(q log q) controller and initialized positive core; the existing mobile wrapper applies. The conformal conclusion uses stated standard massive-sector integration and Gaussian RG inputs. It is not a rigorous thermodynamic theorem or a measured numerical violation. Read [the complete candidate argument](contrib/research_continuation/PROPOSED_COUNTEREXAMPLE.md) and [the infrared audit](contrib/research_continuation/IR_AUDIT.md).
- **Unbounded c at one fixed d remains conditional.** C39/C42 and the C44 wrapper now supply a finite spatial construction at d₀=4464. C45 gives a coupling schedule conditional on a stated Wilsonian input S; C46 proves partial rigorous results but leaves the full infrared step open. Read [the referee](contrib/referee_C44_C45/README.md) and [C46](codex/work_C46/RESULT.md). The sparse first-order family with growing d is a separate route.

[STATUS.md](STATUS.md) has the list of where help is most useful, and [STATUS_DETAIL.md](STATUS_DETAIL.md) the details.

The [direct-attack continuation](contrib/direct_attack_20261028/README.md) supplies ultraviolet matching bounds, an exact exclusion of a hidden free-fermion qubit mechanism, and an independently AI-audited thermal massive-matching theorem allowing infinite excited cell spaces. The retarded multichannel RG and central-charge conclusion remain open.

The [physical-ensemble continuation](contrib/thermal_ensemble_20261028/README.md)
proves interacting physical-to-rooted pressure transfer, an all-length
pressure-zero theorem and a length-tilted thermal reference comparison.
Its [positive completion and matching application](contrib/thermal_ensemble_20261028/MOBILE_MATCHING_APPLICATION.md)
repair the mobile bond on illegal words, retain d0=4464 and the complete
legal dynamics, and supply the full-space massive and ensemble inputs.
The parallel [frozen-chain continuation](contrib/decisive_attack_20261028/README.md)
controls the frozen-defect pressure while retaining the incompatible-ring
qualification. The logical infrared conclusion remains open in both.

**Reversible digit compiler (2026-10-28):**
[digit_compiler](contrib/digit_compiler/README.md) gives fixed tables for all
generator arithmetic and selected-array operations, then compiles them to
one physical reversible head. Actual arithmetic runs restore every work and
history cell. Full-space local permutation gates also pass autonomous
compute/reset cycles. Whole-interpreter storage and dispatch, initialized
Hamiltonian composition, and the final violating dimension remain open.

**Digit compilation and infrared progress (2026-10-28):**
[digit_compilation](contrib/digit_compilation/README.md) supplies an explicit
96-entry parameter-reader table, eleven finite local arithmetic kernels,
and checked execution of eight arithmetic opcodes inside the complete
21-mode circuit generator. All 287,798 abstract gates match. A separate
matrix theorem gives the explicit Gaussian stability budget delta<5/13,
including density-current couplings and unequal velocities. The remaining
memory/controller composition and microscopic low-energy reduction are
still required; no new violating dimension or central-charge measurement
is claimed.

**Microscopic implementation progress (2026-10-28):**
[executable circuit and controller work](contrib/controller_instantiation/README.md)
now supplies a full sparse-encoding gate generator, explicit reversible
dispatcher and length-counter tables, a repaired input-format check, and
a saved d=240 one-mode test Hamiltonian. The autonomous program generator,
verified program initialization and endpoint readout still need to be
composed into the high-channel core. The saved fixture is not a c-d violation.

The [autonomous-controller continuation](contrib/autonomous_controller/README.md)
now composes blank program generation, coherent execution, endpoint readout
and complete reset in one explicit local machine for a demonstration program.
Its initialized clock is also instantiated. A separate fixed 3,115-instruction
word program reproduces the full sparse-circuit generator; digit-level local
compilation and high-channel composition remain. Two explicit parameter-query
gates remove the need to decode the chosen prime from cell length. None of
these demonstration dimensions is a newly claimed violating d.

Additional derivations from the chat "Prove c d conjecture" are collected in
[contrib/chat_followup/](contrib/chat_followup/README.md): a species-torus
selection rule for the encoded target, a massive-auxiliary construction,
exact extensions of controlled classes, and a detailed qubit catalog. The
imported proposal remains conditional; the [new continuation](contrib/research_continuation/README.md) gives a different sparse-graph construction. Also see [C14's progress
log](codex/work_C14/PROGRESS.md) for its reported physics-level counterexample
to the proposed velocity-weighted Casimir bound; that does not refute c-d.

## Open problems

These were written when the main line was a proof, and they are still open. For the tasks that matter most for the construction, see STATUS.md.

1. Find a locality-based budget for neutral sectors: starting from a fixed two-site bond on ℂ^d ⊗ ℂ^d, bound all emergent stress-tensor sectors by d − 1.
2. Find an improvement-independent version of the velocity-free stress response F_h(q) = ∫ S_h(q,ω)/ω² dω → C|q|/(24π), and a microscopic bound on it.
3. Prove a nearest-neighbour level ceiling k ≤ ind(R) for chains with an on-site symmetry flowing to WZW models, or find a chain that violates it (for example Sp(4)₂ on ℂ⁴, or SO(N)₃ on an SO(N) vector).
4. Run larger sizes on the two numerical loose ends: qubit chains with entropy slopes 1.35 to 1.7, and small-angle Sp(4) (REPORT.md §7).
5. Establish Bethe-ansatz completeness for fundamental gl(d) chains, with control of massless strings. That would give c ≤ d − 1 for this integrable class.

## Layout

| Path | Contents |
|---|---|
| `STATUS.md` | The current state on one page; start here |
| `REPORT.md` | The running report, with every claim tagged and provenance in §6 |
| `lit.md` | Literature map, claims tagged [read] or [recalled] |
| `routes_draft.md` | Short derivations behind REPORT.md §3 |
| `specs/` | Self-contained problem statements given to the model runs (`_background.md` is shared) |
| `codex/*.out` | The answers those runs returned |
| `codex/work_*/` | Code, data and full write-ups of each run |
| `codex/opus_check_C3/` | Independent check scripts |
| `contrib/` | New results from collaborators, one folder each |
| `tools/` | Literature digests of techniques that might apply (entanglement bootstrap, McGreevy's notes) |

## Provenance

The material up to October 2026 was produced with AI assistance and has not been refereed. Claude (Anthropic, Opus 5.5) planned the work, wrote REPORT.md and checked results, so "I" in REPORT.md refers to Claude. Closed-book derivations and numerics were run on GPT through the Codex CLI, almost all on gpt-6-astra. Two early runs on gpt-5.6-sol are kept for comparison. Contributions should say in the same way who or what produced them.
