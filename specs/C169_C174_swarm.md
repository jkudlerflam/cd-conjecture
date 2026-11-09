# C169 to C174: a swarm on the critical path to c ≤ 1 at d = 2

## The critical path

Target: c ≤ 1 for every translation-invariant nearest-neighbour qubit chain, at the theoretical-physics standard.

- Derived: every symmetric class, see STATUS.md "Is c ≤ 1 at d = 2?".
- Still assumed:
  - (A_W) on the time-reversal strata T and K_I (C163). Agent C167 is attacking it.
  - (H)/(A6) on the symmetry-free stratum N (C154 PROOF 7). Agent C168 is attacking it.
  - A_S near the hidden Heisenberg point in N (C165). A Typhon screen is testing it.
- A derivation also needs the load-bearing steps to survive refereeing. C159 (Lemmas A, B, Proposition C), C163 (bound-state and uniform-tilt exclusion) and C165 §4 (the SU(2)_1 identification) have not been refereed.

Six agents run alongside C167 and C168, each with one role. Your role is the section with your number. Do only that.

## Shared reading for every agent

- `specs/_qubit_c_context.md`
- STATUS.md
- REPORT.md §9.110 to §9.126
- `codex/work_C154/RESULT.md`

After that, read what your section names.

## Shared rules for every agent

- Physics standard. RG, bosonization, perturbation theory and calibrated numerics count, and "not rigorous" is not a gap. Show all steps and include consistency checks (limits, special cases, matching known results).
- Tag every claim proved, derived, numerical, heuristic or assumed.
- Your folder is `codex/work_C1xx/`, using your number. RESULT.md opens with a line saying which model did the work. Derivations go in PROOF.md and the log in PROGRESS.md.
- Write only inside your folder. Do not edit STATUS.md, REPORT.md, specs or other agents' folders. Do not commit or push. Never stash, reset or switch branches, do not pip install, and do not launch agents.
- Astra is out of credits, so do the work yourself. Never call codex directly, and never use paid APIs.
- Numerics stay small and local (ED, DMRG at L ≤ 128, χ ≤ 128) unless your section says otherwise. If a large run is essential, use at most one packed 48-core Typhon job: `ssh -o LogLevel=ERROR typhon`, SLURM only, never on the login node, scratch /scratch/lustre/jkudlerflam/c1xx/. Record the job id and do not block. The queue is shared, so avoid Typhon if you can.
- Prose: commas not em-dashes, few colons and semicolons, no italics for emphasis, no filler intensifiers.
- Final report: the verdict, key claims with tags, and what is left open.

## C169: constructive counterexample hunter

Try hard to build a d = 2 nearest-neighbour bond with c > 1, by structure rather than sampling. Sampling is done (C155, about 150,000 evaluations), so do not repeat it.

Look for exact maps that send a known c > 1 model to a one-site nearest-neighbour qubit bond, or show that none exist:
- products of single-site unitaries, possibly site-dependent with a period, as in C154's twisted strata;
- Kramers–Wannier, Jordan–Wigner and cluster-entangler dualities restricted to those that keep nearest-neighbour range;
- blocking.

Then consider deformations of a c = 2 or c = 3/2 model with a two-site cell (for example two XX chains, or Ising × XX) that become nearest neighbour on a one-site cell. Ask whether the obstruction is range, translation by one site, or dynamics.

Also read `codex/work_C165/RESULT.md` (the classical coexistence that fails quantum mechanically), `codex/work_C157` (XZ, the one bond with a second decoupled Majorana sector) and `codex/work_C156` (literature). Output either an explicit candidate tested by C163 PROOF §4.2's pre-fixed rules, or a precise statement of why each construction fails.

## C170: referee of C165 §4 (the SU(2)_1 identification)

Referee `codex/work_C165/PROOF.md`, especially §1.2 (exact site rotations to the XXZ chain), §2 (the Landau quartic) and §4 (spin 1/2: staggered and spiral orders as one SU(2)_1 primary, quartic of dimension ≥ 8, to all orders in weak harmonics).

Check every step independently, re-deriving where you can. Check the claim "to all orders in weak harmonics" in particular: does any harmonic generate a relevant operator at momentum 0 or π that splits the multiplet? Check the scope statements. Small numerical checks are fine.

Write REFEREE.md with each claim marked correct, gap (fixable), or wrong, and give the corrected statement where needed.

## C171: referee of C159 and C163 (the T and K_I chain)

Referee:
- `codex/work_C159/PROOF.md`: Lemma A (half-floating liquids under an antiunitary that commutes with translation), Lemma B (the Kramers LSM anomaly on T) and Proposition C (any c > 1 candidate filling a hypersurface or an open set of T or K contains an attractive Luttinger sector);
- `codex/work_C163/PROOF.md`: the exclusion of bound-state and uniform-tilt second sectors in the attractive window, and the statement of (A_W).

Check each step independently, re-deriving where you can, and check that the reduction to (A_W) is complete: every route to c > 1 on T and K_I ends in (A_W) or in a derived exclusion. `codex/work_C161/REFEREE.md` is a model of the format. Write REFEREE.md with each claim marked correct, gap (fixable), or wrong.

## C172: closing A_S analytically

A_S (C165 RESULT.md): in N, with harmonics comparable to the exchange near the hidden Heisenberg point (parent Δ near 1, fields up to the Néel plateau edge), no lattice-scale excitation independent of the liquid's boson orders at momentum π while the liquid survives. C165 derived it to all orders in weak harmonics. The Typhon screen (job 17601089) tests amplitudes 0.2 to 0.6 numerically, so do not redo that screen.

Find an argument that works at order-one amplitude. Some leads:
- a continuity or non-crossing argument for the SU(2)_1 conformal towers along paths of growing amplitude;
- an exact identity of the rotated bond that keeps a hidden SU(2) or U(1) along a submanifold;
- a bound on the gap of the momentum-π sector in terms of the liquid's own vertex operators.

Use small DMRG or ED only to check a specific derived statement, for example conformal-tower spectra at L ≤ 24. Read `codex/work_C165/` in full.

## C173: two floating bosons at codimension 0, a coupled-wire route

C154's kinematic count allows two independent floating bosons on an open set of N. Attack this with a specific picture, different from C168's general approach.

A nearest-neighbour qubit chain is a zigzag ladder (even and odd sites as legs) with only rung couplings, so any intra-leg dynamics is generated at second order. Two floating bosons would need both legs, or two collective channels, to become independent Luttinger liquids.

Ask:
- whether the Wilsonian generation of intra-leg terms from rung terms can ever make two channels independent and gapless at codimension 0;
- or whether the rung coupling, being first order, always locks one combination or gaps it.

Consider the regimes of strong rungs, weak rungs and the classical spiral. Compare with known c = 2 zigzag ladders, which need next-nearest-neighbour couplings. State the result as a derived obstruction with its regime, or as a mechanism that evades it (a candidate, to be tested with C163 PROOF §4.2's rules).

Read `codex/work_C154/PROOF.md` (the codimension count) and `codex/work_C168/` if it exists (it is in progress, so read it, do not edit it).

## C174: red team of the conditional theorem

Attack the overall logic, not individual calculations:
- C154 PROOF 7: the stratification of B_2, the closure lemma ("c is constant on a critical component and c at a conformal frontier point is at least c of the component"), its reliance on C141's weak c-theorem, and whether the controlled set really bounds every component.
- Whether the list of strata is complete.
- Whether the controlled regions' boundaries are covered.
- Whether "c ≤ 1 on every symmetric class" plus (H) truly implies c ≤ 1 on all of B_2. Look for non-sequiturs, missing cases (non-conformal gapless points, Lifshitz points, first-order lines, spontaneously broken translation with a larger cell) and any illegal interchange of limits.

Read `codex/work_C154/` and `codex/work_C141` in full, and the RESULT.md of C152, C153, C159, C160, C162 and C163. Write REDTEAM.md, giving each hole you find its severity and the narrowest fix.
