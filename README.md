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

Put the code and data in the same folder. Add a short dated paragraph at the top of [UPDATES.md](UPDATES.md) and a short section to REPORT.md. Edit the summary in "Where things stand" below only if the overall picture changes, and keep it short.

### 4. Discussing

Use [GitHub issues](https://github.com/jkudlerflam/cd-conjecture/issues) for open problems, proposed attacks, and claims that need an independent check. Keep one question per issue and link the files or commits it concerns. The open problems listed below are good starting points.

### 5. Running the existing code

The code is plain Python and needs numpy and scipy unless a folder says otherwise. Each `codex/work_*/` folder holds the scripts, data and write-up of one run, and its report or README says how the scripts were run and how long they took. `codex/work_C4_D2/` and `codex/work_C4_D3/` also contain SLURM scripts for larger scans. To rerun a whole calculation from scratch, use the matching problem statement in `specs/`.

### 6. Conventions

These are spelled out in [CONTRIBUTING.md](CONTRIBUTING.md). In brief, tag every claim, name which central charge you mean (Virasoro c, or c_eff for non-compact or non-unitary limits), and keep files under about 50 MB. Never commit passwords, tokens or other private data.

## Where things stand

Last updated 2026-11-25. No counterexample to c ≤ d − 1 is established, and no general proof exists. There are explicit chains that predict a violation, at d = 2 from numerics and at d ≥ 9 conditionally on standard infrared inputs, and there is a proposed replacement bound that is proved in several classes.

The working standard is theoretical physics (standard RG, controlled perturbation theory, calibrated numerics), not mathematical rigor. Nothing here has had independent human review, and a conditional infrared prediction is not counted as an established counterexample. [STATUS.md](STATUS.md) has more detail and the list of where help is most useful, [REPORT.md](REPORT.md) is the chronological record, and [UPDATES.md](UPDATES.md) keeps the dated announcements of each research round that used to sit in this section.

### d = 2: a candidate with c ≈ 2

At d = 2 the conjecture says c ≤ 1. One family of qubit chains reads c ≈ 2. Strong ZZ repulsion near its blockade field leaves a dilute gas of particles whose direct hop is tuned down to the size of the second-order next-nearest hop. The particles then move on a zigzag ladder with two band valleys, and the proposed infrared theory is a two-component Luttinger liquid. The family has two well-studied points.

| Bond | Evidence | Where |
|---|---|---|
| Large-blockade bond, and Pace's neighbour with hopping phase 1.53 (written out in STATUS.md) | Open chains read c = 2.00 from L = 128 to 256, converged in bond dimension, and 1.95 to 1.98 from 256 to 512 at χ = 2048. The effective model reads 1.99 to 2.00 on the infinite chain. | [STATUS.md](STATUS.md), REPORT §9.178 to §9.205 |
| Four-term chain Σ[((1+δ)/2) XY + ((1−δ)/2) YX + ZZ] − h Σ Z at δ = 3/16, h = 1.783 | Infinite chain 1.97 ± 0.01, open chains L = 96 to 384 1.96 ± 0.08, inside a narrow wedge below a chirality-ordering line. It lies in the same family, at V/Δ = 4. The c = 3/2 reading reported on November 18 was withdrawn the next day. | [contrib](contrib/d2_xy_yx_chirality_imps_20261118/README.md), [review](contrib/qubit_c2_status_20261123/README.md), REPORT §9.209b to §9.215 |

At both points the pseudospin Luttinger parameter exceeds one (K_s ≈ 1.4 to 1.6 and ≈ 1.8 to 2.0), which makes the backscattering between valleys irrelevant. What stands between this and a counterexample is mostly long-distance stability. The chain conserves only parity, and a number-changing term can lock the charge mode on scales beyond the reach of L ≤ 512. Its coupling has been estimated but not computed. The operators that carry the measured exponents are also not yet identified (REPORT §9.215 and the [list of missing measurements](contrib/qubit_ir_cft_team_20261123/TODO.md)).

Away from this corner no c > 1 has been found, in about 103,000 sampled qubit chains or in targeted searches, and c ≤ 1 is derived at generic points of the symmetric classes, with one time-reversal sheet still open.

The [November 25 quartet audit](contrib/quartet_typhon_20261125/README.md)
finds a substantial filtering bias in the historical x4≈2.05±0.2 estimate;
it is not a calibrated charge-four dimension. Operator matching and the
locking question remain unresolved (REPORT §9.219).

### Unbounded c at a fixed small d (conditional)

These are explicit cell constructions. Gapped cells of length n carry logical particles whose band has 2n simple Fermi points, so the predicted infrared theory has c growing with n while d stays fixed. In each case the finite band data are certified. The step from the microscopic chain to the multicomponent Luttinger liquid rests on the same stated infrared inputs, namely matching to the logical model, selection of the right ground sector, and Gaussian RG stability.

| d | Predicted c | Certified | Where |
|---|---|---|---|
| 9 (frozen cells) | n | Band data at n = 9, 10, 11, 12, 14, 16, 20. A modified d = 9 bond has an all-n band proof, local interaction bounds independent of chain length, and exact ground-sector selection on compatible rings. | [markerless_d9_d10](contrib/markerless_d9_d10/README.md), [all-n](contrib/markerless_all_n_20261118/README.md), [uniform bounds](contrib/markerless_uniform_matching_20261118/README.md) |
| 10 (one mobile vacancy) | n + 1 | Same n as d = 9, plus finite mobile sector comparisons. The mobile selector and its J/t budget are open. | [markerless_d9_d10](contrib/markerless_d9_d10/README.md), [mobile control](contrib/markerless_mobile_control_20261118/README.md) |
| 11 (frozen) and 13 (mobile) | 11 and 14 | Members n = 11 and n = 13. d = 13 was red-teamed. | [small_d_note](contrib/small_d_note/note.pdf) |
| 4464 | n + 1 | Every n. Massive matching, ensemble transfer and c from a Gaussian limit are proved and refereed. The open core is a constructive RG. | [expert note](contrib/expert_note/note.pdf), [referee](contrib/referee_C44_C45/README.md) |

Numerics on the logical vacancy model give c = 2.0, 2.9 ± 0.1 and 3.98 ± 0.05 at n = 1, 2, 3, and the projected d = 4464 chain gives c = 3.00 ± 0.02 at n = 2. Intermediate versions (a d = 12 mobile proposal, d = 20 and 22 all-n families, a d = 24 one-vacancy chain) are kept in `contrib/` and described in UPDATES.md.

### Growing d with c/d → ∞ (conditional)

A sparse species-graph family has predicted c = (p+1)(p²+p+1) + 1 at d = O(p² log p), so c/d diverges. The finite pieces (ports, reversible circuit, controller, initialized core) are constructed explicitly, and the conformal conclusion uses standard massive-sector integration and Gaussian RG. See [the candidate argument](contrib/research_continuation/PROPOSED_COUNTEREXAMPLE.md) and [its infrared audit](contrib/research_continuation/IR_AUDIT.md).

### Proofs, and a replacement bound

c ≤ d − 1 is proved for quasi-free chains (c ≤ log₂ d), SU(d)-symmetric bonds, and some classes of infrared theories. At d = 2 it is proved in the symmetric ferromagnetic rank-one class (c ≤ 1/2) and in the outer regions of the one remaining staggered-field family ([proofs](contrib/d2_gap_reduction/README.md)).

Since the cell constructions suggest that c alone is not bounded at fixed d, we also propose a bound with the velocities left in, A = Σ_a c_a v_a ≤ π M log₂ d, with M the gauge-optimized bond half-spread. It holds in every solvable family and is proved for quasi-free chains (A ≤ 2M, sharp), free-fermion and classical chains, and WZW sectors coupled to on-site symmetries. The general case reduces to showing C(d) = sup_L L(E_{2L} − 2E_L)/M is finite.

- At d = 2, C(2) is finite on the whole U(1) class (XXZ plus DM along z plus field), with an explicit constant of 6.91 on −1 < J_z/J ≤ −1/√2 next to the ferromagnet (REPORT §9.207, §9.210).
- At general d, (L e₀ − E_L)/M → 0 uniformly over all bonds, but with no rate at reachable sizes. A ≤ 2M holds for two-species exchange chains, and a crossover-free bound holds for the critical ferromagnetic Potts chains at Q = 3, 4 (REPORT §9.207).
- For the cell constructions above, the sum of Fermi speeds is at most C d⁴ times the bond norm, so free cells cannot make A unbounded. With interactions this is not controlled (REPORT §9.210, §9.211).

A finite C(d) for general chains is open.

### Other small d

At d = 3 to 7 nothing above d − 1 has been found in symmetric families. If the qubit c = 2 is confirmed, two on-site copies give c = 4 at d = 4, and a vacancy lift targets c = 3 at d = 3 ([note](contrib/small_d_lifting_20261116/README.md)).

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
| `UPDATES.md` | Dated announcements of each research round, newest first |
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
