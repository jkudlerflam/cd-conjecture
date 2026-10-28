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

A new candidate refutation is recorded below and needs independent review. The working standard is theoretical physics (standard RG, controlled perturbation theory, calibrated numerics), not mathematical rigor. In brief:

- **Small d: no counterexample.** Nothing above d − 1 at d = 2, in symmetric families at d = 3 to 7, or at the best d = 4 target. The d = 3 hint failed its own follow-up (C12); an independent check remains pending in the recorded status.
- **No general proof.** Every lattice quantity we know that bounds c grows like (a length scale) × log d. Proofs exist for quasi-free chains (c ≤ log₂ d), SU(d)-symmetric bonds, and some classes of IR theories, and there are conditional bounds.
- **Large d: a candidate physics-level refutation.** The new sparse species-graph family has predicted c=(p+1)(p²+p+1)+1 and microscopic d=O(p² log p), so c/d diverges. Exact finite constructions now supply the ports, reversible circuit, O(q log q) controller and initialized positive core; the existing mobile wrapper applies. The conformal conclusion uses stated standard massive-sector integration and Gaussian RG inputs. It is not a rigorous thermodynamic theorem or a measured numerical violation. Read [the complete candidate argument](contrib/research_continuation/PROPOSED_COUNTEREXAMPLE.md) and [the infrared audit](contrib/research_continuation/IR_AUDIT.md).
- **Unbounded c at one fixed d is a separate, unfinished goal.** Route A still requires a spatial implementation of its engineered excited spectrum. The new first-order family has growing d and does not prove that no bound c≤f(d) exists.

[STATUS.md](STATUS.md) has the details and the list of where help is most useful.

**Microscopic implementation progress (2026-10-28):**
[executable circuit and controller work](contrib/controller_instantiation/README.md)
now supplies a full sparse-encoding gate generator, explicit reversible
dispatcher and length-counter tables, a repaired input-format check, and
a saved d=240 one-mode test Hamiltonian. The autonomous program generator,
verified program initialization and endpoint readout still need to be
composed into the high-channel core. The saved fixture is not a c-d violation.

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
