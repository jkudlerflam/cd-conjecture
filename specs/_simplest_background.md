## Background for the search for the simplest counterexample (C75 to C77)

Goal: the simplest one-site translation-invariant nearest-neighbour Hermitian bond h on C^d with c > d − 1, at the theoretical-physics standard (controlled perturbation theory, standard RG, calibrated numerics where possible). "Simplest" means first the smallest local dimension d, then the fewest and most transparent bond terms and the mildest coupling hierarchy.

What exists (read `../../STATUS.md` first):
- `../../contrib/explicit_candidate/README.md` and `bond23.py`: one 23-letter bond whose member with cell length n + 1 has c = n for every n (frozen cells), with C33's finite certificate at n = 23 (`../work_C33/`).
- `../../contrib/counterexample_campaign/COUNTEREXAMPLE_D24.md`: the collaborator's one-vacancy lift, d = 24, predicted c = 24 (23 payload modes plus a charge boson).
- `../work_C39/`, `../work_C44/`: a 46-letter core with all-n margins and the mobile wrapper, d0 = 4464, c = n + 1 for all n.
- Demonstrators: `../../contrib/routeA_demonstrator/` (10 letters, c = 1 measured), `../../contrib/routeA_demonstrator_n2/` (13 letters, c = 2 measured).
- `../../contrib/minimal_d0/README.md` (alphabet accounting), `../work_C30/` (switch-cell selector), `../work_C22/` (17-letter cell), `../work_C8/` (constant-alphabet clock, mobile cells).
- Bounds: free-fermion chains obey the bound (C1). At d = 2 nothing exceeds c = 1 (`../work_C63/`, `../work_C67/`, `../work_C68/`). Integrable SU(N)_1 and SU(N|M) chains saturate c = d − 1.
The key structural fact: in the route-A family the alphabet is fixed while the number of flavours n grows with the cell length, so a fixed alphabet of size a yields a violation as soon as some member with c ≥ a is under control. The cost is a hierarchy of scales that grows with n.

Rules: tag claims proved, numerical (with sizes), heuristic, assumed. No probabilities, odds or percent-confidence in any file. Python (numpy, scipy) only, inside your folder; no network commands (a browsing tool may be used for primary literature); never modify files outside your folder. Plain prose without em-dashes. Keep PROGRESS.md updated and write RESULT.md (and PROOF.md for derivations).
