# The c–d conjecture

Latorre and Sierra ([arXiv:2403.17242](https://arxiv.org/abs/2403.17242), J. Stat. Mech. (2024) 113103) conjecture that a critical, translation-invariant, nearest-neighbour Hermitian chain on (ℂ^d)^⊗L flows to unitary CFT(s) with total central charge c ≤ d − 1. The bound is saturated by the SU(d) Uimin–Lai–Sutherland chain h = P (the swap), whose IR is SU(d)₁ with c = d − 1, so the conjecture is equivalent to c_max(d) = d − 1.

This repository is a shared notebook for proving or refuting it. Everything we know goes here, and everyone with access pushes to and pulls from it as the work develops. Maintained by Jonah Kudler-Flam.

## How to use this repository

### 1. Getting oriented

Read [STATUS.md](STATUS.md) first. It is the current state of the project on one page, with the main construction written out step by step, what is checked, what is open, and a glossary of the run labels (C1 to C12, G1, route A and B, and so on). [REPORT.md](REPORT.md) is the chronological record behind it, organized as follows.

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

When you are done, commit in small pieces and push:

```bash
git add <files>
git commit -m "Short description of what changed"
git push
```

If the push is rejected because someone else pushed first, run `git pull --rebase` again, resolve any conflicts, and push. Everyone with write access works directly on `main`. If you don't have write access, fork the repository and open a pull request.

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

**Conclusion (2026-10-29), at the theoretical-physics standard: the c–d conjecture is false.** Explicit chains have c = 14 at d = 13 and c unbounded at fixed d = 4464, conditional only on standard Wilsonian RG input whose main pieces are now proved; the natural replacement, Σ_a c_a v_a ≤ π M log₂ d, survives every test. See the conclusion at the top of [STATUS.md](STATUS.md).

A d=24 conditional counterexample candidate is recorded with independent AI reviews and a [fresh adversarial audit](contrib/counterexample_audit_20261028/README.md). The earlier verified-counterexample wording is withdrawn: finite-J matching and a uniform dressed-defect comparison remain open. The working standard is theoretical physics (standard RG, controlled perturbation theory, calibrated numerics), not mathematical rigor.

- **Smallest candidates (2026-10-29): d = 11 frozen with predicted c = 11, and d = 13 mobile with predicted c = 14.** Explicit route-A bonds with certified finite members, refereed, and red-teamed for d = 13; both rest on the same inherited infrared inputs as the d = 24 candidate. Full bonds and status in [contrib/small_d_note/](contrib/small_d_note/note.pdf) and STATUS.md. The logical vacancy model gives calibrated c = 2, 2.9 ± 0.1 and 3.98 ± 0.05 at n = 1, 2, 3 (REPORT §9.34, §9.55, §9.67).

- **Unbounded c at fixed d0 = 4464 (the C39/C44/C45 wrapper).** c = n + 1 for every n, conditional on the Wilsonian input S. Several pieces of S are now proved and refereed (massive matching C48, ensemble transfer C58, c from a Gaussian limit C51); the open core is the constructive RG (Lemma N and lemma TP). The full microscopic demonstrator flows to its matched XXZ Luttinger liquid numerically (REPORT §9.58). Expert summary: [contrib/expert_note/note.pdf](contrib/expert_note/note.pdf).

- **d = 2.** No c > 1 found in about 103,000 sampled qubit chains of the residual strata, on the staggered branch to L = 256, or inside floating phases (REPORT §9.37, §9.40, §9.62); the analytic covering argument is not complete.

- **d=24, predicted c=24: a one-vacancy candidate.** Add one mobile vacancy to the explicit 23-state chain, with its original bond acting on occupied pairs. The local phase-model representation, certified 23-band reference, nonzero-coupling sector bounds and physical momentum rules support 23 payload bosons plus one charge boson, conditional on finite-coupling matching and defect stability. Published local perturbation bounds control the original cell machinery. The conclusion uses conventional Gaussian phase stability; no microscopic c=24 measurement or rigorous CFT-limit theorem is claimed. Read [the construction and independent reviews](contrib/counterexample_campaign/COUNTEREXAMPLE_D24.md).

- **Independent vacancy continuation.** A proved density window removes exact half-filling tuning from the proposed d=24 construction. New reviews check the local perturbation argument and C39's all-n finite input, including a new exact four-mode certificate. A full microscopic d=3 control has entropy and independently calibrated energy scaling consistent with c=2; this tests the vacancy mechanism, not the c=24 chain. See [the proofs, reviews and numerical data](contrib/vacancy_continuation/README.md).

- **Small d: no counterexample.** Nothing above d − 1 at d = 2, in symmetric families at d = 3 to 7, or at the best d = 4 target. The d = 3 hint failed its own follow-up (C12); an independent check remains pending in the recorded status.
- **New d = 2 proof reductions.** Arbitrary uniform fields in the symmetric ferromagnetic rank-one class give at most c = 1/2, using a recent rigorous gap theorem and the exact Ising solution. An exhaustive classification of additive on-site spin charges leaves one additional two-parameter staggered-field family; exact gap arguments exclude both outer regions. The full qubit conjecture remains open. See [the proofs and remaining region](contrib/d2_gap_reduction/README.md).
- **No general proof.** Every lattice quantity we know that bounds c grows like (a length scale) × log d. Proofs exist for quasi-free chains (c ≤ log₂ d), SU(d)-symmetric bonds, and some classes of IR theories, and there are conditional bounds.
- **Large d: a candidate physics-level refutation.** The new sparse species-graph family has predicted c=(p+1)(p²+p+1)+1 and microscopic d=O(p² log p), so c/d diverges. Exact finite constructions now supply the ports, reversible circuit, O(q log q) controller and initialized positive core; the existing mobile wrapper applies. The conformal conclusion uses stated standard massive-sector integration and Gaussian RG inputs. It is not a rigorous thermodynamic theorem or a measured numerical violation. Read [the complete candidate argument](contrib/research_continuation/PROPOSED_COUNTEREXAMPLE.md) and [the infrared audit](contrib/research_continuation/IR_AUDIT.md).
- **Unbounded c at one fixed d remains conditional.** C39/C42 and the C44 wrapper now supply a finite spatial construction at d₀=4464. C45 gives a coupling schedule conditional on a stated Wilsonian input S; C46 proves partial rigorous results but leaves the full infrared step open. Read [the referee](contrib/referee_C44_C45/README.md) and [C46](codex/work_C46/RESULT.md). The sparse first-order family with growing d is a separate route.

[STATUS.md](STATUS.md) has the details and the list of where help is most useful.

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
