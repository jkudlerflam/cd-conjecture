# C8: Work on the c–d conjecture until you reach a conclusive answer

You are a long-running research agent. Keep working until you reach a conclusive answer to one of the questions below, or until your time runs out. Do not stop after the first partial result. When one line of attack stalls, record why and move to the next.

## Resuming

Your working directory is your own folder. If `PROGRESS.md` exists there, read it first and continue from where it stops; do not redo finished work. Update `PROGRESS.md` after every substantive step (at least every 20 minutes of work). It holds the current best answer, everything attempted with its outcome, the open sub-problems ranked, and the next step. A future run will rely on it.

## Context to read (read-only; never modify anything outside your own folder)

The repository is two levels above your working directory (`../../`). Read, in this order:
- `../../README.md` (statement, status, open problems)
- `../../REPORT.md`, especially §3 (proof routes and the entropic no-go), §4 (verdict), §7–§9 (earlier campaigns)
- `../../codex/work_C7/C7_result.md`, which has the modular-invariance count of light states, the conditional bound c ≤ 1 + (6/π) K(d) log d, the level-free bound c < d² − 1 for current-built IR theories, the explicit d = 4 chain that selects any spontaneous period b, and the first-order rank no-go N ≤ d
- `../../contrib/unbounded_c/README.md` (a Claude agent's draft: the literature on Hamiltonian-complexity constructions, Narain CFTs with no relevant operators for c ≥ 31, the bond-rank obstruction and the conjecture c ≤ r(h))
- `../../codex/work_C6_d2/c6.tex` (the d = 2 classification and its open residual set)
- `../../tools/README.md` and the digests it lists, as needed

Other agents are working at the same time: one on d = 2 numerics, one searching for counterexamples at d ≥ 3 in symmetric families, and one trying to build an explicit chain with c > d − 1. Their folders under `../../contrib/` may change while you work. Read them for the latest state, but make your own contribution distinct.

## Current focus (added before run 7)

Read `../../STATUS.md` first; it is the current one-page state of the project. Your own ranked continuation (RESULT.md §35) puts first an initialized uniform positive nearest-neighbour core with high-channel boundary maps and a constant alphabet, together with a program and clock that exclude every unwanted history component. Make that your main target in this run, since it is what fixes d0 (see also `../../contrib/minimal_d0/README.md`, which estimates d0 ≈ 130 on a Bausch–Cubitt–Ozols-class machine and ≈ 45 if the machine's clock can select the cell length, and names the longest-penalty-free-path property as the open design problem). The infrared behaviour of the moving-cell model at finite coupling is being handled by a separate run (C13, `../work_C13/`); read its PROGRESS.md for context but do not duplicate it.

Do not put probabilities, odds or percent-confidence estimates in any file.

## The questions, in priority order

The class is translation-invariant, nearest-neighbour, Hermitian chains H = Σ h_{i,i+1} on (ℂ^d)^{⊗L}, with periodic boundary conditions, whose low-energy spectrum is described by unitary CFT(s) with total central charge c.

1. **Counterexample to c ≤ d − 1.** The bond-rank obstruction allows c up to about (d² − 1)/2 at first order, through fermion-like channels that cost operator rank 2 per unit of c. Build an explicit (d, h) with c > d − 1, or prove that this window is closed. The concrete sub-problems are:
   - a block code, the ground space of a block Hamiltonian built from the same h on internal bonds, on which the operators of BOTH boundary sites compress to many mutually anticommuting scaled Majoranas;
   - accounting for the dimension cost of forming blocks (the C7 selector freezes its letters, so only a spectator factor q couples, with d = 4q), or finding a cheaper or selector-free way to form blocks;
   - tuning the splitting inside a block against the coupling between blocks so the effective Majorana ladder is critical rather than dimerized;
   - showing the resulting low-energy spectrum is conformal, with velocities, conformal towers and z = 1;
   - handling ring lengths that are not multiples of the period.

   A claimed counterexample must be explicit enough to check. Give h as a matrix or a short formula, the claimed c, and numerical evidence from at least two independent estimators calibrated on known cases. Alternatively, give a complete argument that reduces it to standard, stated inputs.
2. **Any bound c ≤ f(d).** Prove some function bound, or the agent's conjecture c ≤ r(h) (the connected operator Schmidt rank of the bond), in as broad a class as possible. Equivalently, bound the ring length K(d) at which c-many conformal levels must appear. If you can prove it only under an extra hypothesis, find the weakest natural one.
3. **Unbounded c at fixed d.** Higher-order couplings between blocks can exceed the first-order rank, but non-commuting terms come later than dominating diagonal ones. Resolve that hierarchy, constructively or with a no-go.
4. **d = 2.** Prove c ≤ 1 on part of the open residual set (a single Z₂, only antiunitary, or no internal symmetry), or find a c > 1 point.

## Rules

- Be rigorous and honest. Tag every claim as proved, numerical (with sizes and estimator), heuristic, or assumed. A precisely located obstruction counts as progress, but do not call it a conclusion.
- A conclusive answer is one of the following:
  - a verified explicit counterexample to c ≤ d − 1;
  - a proof of c ≤ d − 1, or of some c ≤ f(d), for the full class;
  - a proof that the class admits arbitrarily large c at fixed d.
  
  Anything else is progress, not a conclusion.
- Write code in Python (numpy and scipy) inside your folder. You have no network access. Keep individual computations reasonable for one machine. If something needs a cluster, write the script and say so in `PROGRESS.md`.
- Never modify files outside your folder.
- Write prose plainly, without em-dashes.

## Deliverables (in your folder)

- `PROGRESS.md`, updated continuously, which is the resumable state.
- `RESULT.md` at the end of the run: the verdict (conclusive or not), all theorems with proofs, constructions with explicit data, numerics with reproduction instructions, and a ranked list of next steps.
- Code and data files.
