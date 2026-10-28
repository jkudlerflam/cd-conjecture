# C24: Independent review of the species-graph candidate refutation of c ≤ d − 1

Your working directory is your own folder. If `PROGRESS.md` exists there, continue from it. Update it as you go and write `RESULT.md` at the end.

## What to review

The repository is two levels above your folder (`../../`). A collaborator's Codex workers have proposed a physics-level counterexample to the c–d conjecture. Read `../../STATUS.md` (its "Latest continuation" section), then `../../contrib/research_continuation/PROPOSED_COUNTEREXAMPLE.md`, `IR_AUDIT.md`, `SCALE_HIERARCHY.md` and `README.md` in that folder, and the supporting folders `../../contrib/torus_bottleneck/`, `../../contrib/exact_route/` and `../../contrib/graded_route_a/`. The claim: with the point-line incidence graph of the projective plane over F_p as a species graph (q = 2(p² + p + 1) species, N = (p + 1)(p² + p + 1) logical fermions, girth six), an explicit first-order encoded-cell chain has c = N + 1 at local dimension d ≤ 192C(p² + p + 1)[1 + log₂(2(p² + p + 1))] + 48, so c/d → ∞.

## Your job

Act as an expert referee at the theoretical-physics standard. Do not take any step on trust.

1. **The selection rule.** Re-derive that girth six implies every neutral monomial with at most four fermion fields preserves each flavour, and check what happens at six fields: dimensions at the Gaussian point and at the interacting point with the stated stiffness, and whether some six-field or mixed length–payload operator can be relevant or marginal and non-oscillating.
2. **The exact finite ingredients.** Check the port isometry (W†F_L W = εf†, the right-port analogue with parity), the charge-preserving circuit, the initialized core and its FAIL/READY selection of exactly 2^N inputs, and the claimed core dimension C q (1 + log₂ q). Rerun or rewrite the supplied check scripts, and write independent checks for the most load-bearing identities.
3. **The infrared argument.** Audit IR_AUDIT.md: the local massive-sector elimination, the momentum choices k_a = π/4 + aπ/[12(N + 1)] and length density in [5/12, 7/12], the claimed mismatch bounds, the free dimension of every other non-Gaussian resonance, the positivity of the stiffness matrix, and the order of limits (fixed p, fixed cell length, small fixed coupling). In particular check that the scale hierarchy in SCALE_HIERARCHY.md is consistent and that the coupling can be fixed and nonzero for each p while all the needed margins hold.
4. **The pieces it inherits** (the mobile wrapper of `../work_C8/`, the analysis of `../work_C13/` and `../work_C20/`): check that the interfaces actually match (charges, grading, wrapper cost 96 d_core + 48).
5. **Verdict.** Classify every gap you find as fatal, serious, fixable, or cosmetic, at the physics standard. State plainly whether a careful theoretical physicist should accept this as a refutation of c ≤ d − 1, and if not, exactly what would make it acceptable.

## Rules

Tag every claim as proved, numerical (with sizes), heuristic, or assumed. Do not put probabilities, odds or percent-confidence estimates in any file. Python (numpy and scipy) only, inside your folder; no network; never modify files outside your folder; keep each computation to minutes. Write prose plainly, without em-dashes.
