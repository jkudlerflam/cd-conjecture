# C18: A history machine whose own clock selects the cell length

You are a long-running research agent. Your working directory is your own folder. If `PROGRESS.md` exists there, read it first and continue from where it stops. Update `PROGRESS.md` after every substantive step, and write `RESULT.md` at the end.

## Background

The repository is two levels above your folder (`../../`). Read `../../STATUS.md` first (the c–d conjecture and the cell construction, steps 1 to 8), then `../../contrib/minimal_d0/README.md` in full (alphabet accounting; this run attacks its open direction (i)), `../../codex/work_C7/C7_result.md` (the C7 selector: a head hopping freely in a segment of b − 1 sites costs −2t cos(π/b), and convexity in b lets one chemical potential select any b ≥ 5), `../../codex/work_C10/RESULT.md` (an explicit translation-invariant track for a hard-coded program, and its §9 obstruction), and `../../contrib/programmable_cell/README.md`. Astra C8 (`../work_C8/`) and C19 (`../work_C19/`) work on related parts of the computation track; read their PROGRESS.md and do not duplicate them.

## The problem

minimal_d0 estimates d0 ≈ 45 for unbounded c if the history machine's own clock selects the cell length, against ≈ 130 when C7's walker is tensored on. Running the history clock in hopping form, −J Σ (U|t+1⟩⟨t| + h.c.), gives a cell whose legal history has K(b) steps the energy −2J cos(π/(K(b) + 2)), so C7's convexity argument could select b with no extra walker letters. The obstruction: in hopping form a penalty at one node barely raises a long path, so every locally legal configuration must lie on a path no longer than the legal one, and internal lengths (number of modes, counter sizes) must be rigid or derived from b, or the cell picks whichever split has the longest path.

## Tasks

1. Design a translation-invariant nearest-neighbour history machine on marker-bounded segments with the longest-path property: in every segment of length b, the legal history is the unique longest penalty-free path of the transition graph, with every competing component either strictly shorter or penalized by an amount that survives in hopping form. Use monotone counters started at zero, domain-wall letters for a unique head, oblivious control flow, and any other device you need. Prove the property, including all illegal and partially legal sectors.
2. Make it compute what the cell needs: for route A (qubit ports, second-order exchange, `../../contrib/unbounded_c/G1/` and `../work_C9/RESULT.md`) a cell program that depends on the target size only through quantities derived from b, with the port idle throughout (`../../contrib/programmable_cell/` §2).
3. Check C7's convexity selection with K(b) in place of b, including the spacing of competing cell energies and the gap.
4. Count the alphabet exactly and give d0. Test the legality and selection claims numerically on small segments.

## Rules

The standard is theoretical physics, but this part is combinatorial and should be proved where possible. Tag every claim as proved, numerical (with sizes), heuristic, or assumed. A precisely located obstruction counts as progress. Do not put probabilities, odds or percent-confidence estimates in any file. Python (numpy and scipy) only, inside your folder; no network; never modify files outside your folder; keep each computation to minutes. Write prose plainly, without em-dashes. Deliverables: `PROGRESS.md`, `RESULT.md` (construction with explicit transition rules and alphabet, proofs, numerical checks with reproduction instructions, d0, ranked next steps), code and data.
