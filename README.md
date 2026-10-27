# The c–d conjecture

Latorre and Sierra ([arXiv:2403.17242](https://arxiv.org/abs/2403.17242), J. Stat. Mech. (2024) 113103) conjecture that a critical, translation-invariant, nearest-neighbour Hermitian chain on (ℂ^d)^⊗L flows to unitary CFT(s) with total central charge c ≤ d − 1. The bound is saturated by the SU(d) Uimin–Lai–Sutherland chain h = P (the swap), whose IR is SU(d)₁ with c = d − 1, so the conjecture is equivalent to c_max(d) = d − 1.

This repository is a shared notebook for proving or refuting it. Everything we know goes here, and everyone with access pushes to and pulls from it as the work develops. Maintained by Jonah Kudler-Flam.

## How to use this repository

### 1. Getting oriented

Read "Where things stand" below, then [REPORT.md](REPORT.md). The report is the running account of what has been tried, organized as follows.

- §1 states the conjecture precisely.
- §2 surveys the literature.
- §3 analyzes the proof routes.
- §4 gives the current verdict.
- §5 lists the next calculations.
- §6 records provenance.
- §7 and §8 give the results of the model campaigns.

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

The conjecture is open. We know of no proof and no counterexample within its scope, and it appears to be true.

- Entropic arguments cannot prove it. A monotone lattice c-function bounded by K log d at the lattice scale would give c ≤ K log d, which SU(d)₁ violates at large d. Thermodynamic, Lieb–Robinson and sum-rule bounds only control c times a power of the conformal window.
- It is proved in restricted classes:
  - quasi-free nearest-neighbour chains of any period, where c ≤ log₂ d;
  - SU(d)-symmetric bonds;
  - Gaussian IR theories whose modes are resolved by commuting on-site charges, where c equals the rank of their susceptibility matrix and is at most d − 1;
  - WZW IR theories obeying a level ceiling.
- Interacting invariants bounded by d − 1 exist (the affine dimension of a parton occupation code, the susceptibility rank), but they equal c only when the IR has no neutral sector. Every route so far stops at the same missing lemma, a locality-based bound on neutral (Majorana, parafermion, coset) stress-tensor sectors.
- Numerical and analytic searches found nothing above d − 1. They covered all nearest-neighbour qubit chains, Sp(4), spin-3/2 SU(2), spin-orbital and spin-1 families, and integrable and non-compact candidates. The non-compact candidates turned out to be non-Hermitian.

## Open problems

1. Find a locality-based budget for neutral sectors: starting from a fixed two-site bond on ℂ^d ⊗ ℂ^d, bound all emergent stress-tensor sectors by d − 1.
2. Find an improvement-independent version of the velocity-free stress response F_h(q) = ∫ S_h(q,ω)/ω² dω → C|q|/(24π), and a microscopic bound on it.
3. Prove a nearest-neighbour level ceiling k ≤ ind(R) for chains with an on-site symmetry flowing to WZW models, or find a chain that violates it (for example Sp(4)₂ on ℂ⁴, or SO(N)₃ on an SO(N) vector).
4. Run larger sizes on the two numerical loose ends: qubit chains with entropy slopes 1.35 to 1.7, and small-angle Sp(4) (REPORT.md §7).
5. Establish Bethe-ansatz completeness for fundamental gl(d) chains, with control of massless strings. That would give c ≤ d − 1 for this integrable class.

## Layout

| Path | Contents |
|---|---|
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
