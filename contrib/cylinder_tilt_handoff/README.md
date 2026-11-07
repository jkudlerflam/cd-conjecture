# Real-tilt c-theorem and the velocity-weighted bound (session handoff)

## Claim

A collection of results from one session, each tagged in `HANDOFF.md`. The main ones for the repository: (i) in a Lorentz-invariant QFT the ratio bound R_L <= K + kappa M beta/L is a consequence of the real-tilt c-theorem, since R_L is exactly a crossed-channel Renyi-2 entropy; (ii) for chains with a positive spatial transfer operator the same identity holds on the lattice, and R_L and E_{2L} - 2E_L are nonincreasing in L; (iii) Lieb-Robinson plus tilt positivity cannot replace Euclidean rotation invariance, because the reachable rapidity is capped at v_s/v_LR < 1; (iv) monotonicity of the velocity-weighted Casimir function c_A(L) cannot give A <= f(d) M for route-A chains, whose first vacuum-compatible ring length grows at fixed d, so the bound must hold at every length, which is the reconnection/ratio bound itself; (v) an explicit convex model shows that thermal-channel convexity plus lattice bounds do not bound A, and a d = 4 quasi-free chain has a fast species with arbitrarily long crossover.

## Status

Mixed: items are tagged proved, numerical, heuristic or assumed in `HANDOFF.md`. The positivity of the quantum transfer matrix for reflection-positive chains (ii) is assumed and needs checking.

## How it was established and checked

- Author(s): Claude (Anthropic) in a session with Jonah Kudler-Flam, 2026-11-07.
- AI assistance: the continuum duality and the claim Theta(q) = -L E_0(L) convex in L^2 come from a parallel session; here they were rechecked (Hessian at finite boost, strip-analyticity argument, sign via J(tau)).
- Checked by: nobody else yet. The numerics were rerun when this folder was assembled.

## Argument or method

See `HANDOFF.md`, sections 1 and 2.

## Reproducing

All scripts use numpy and scipy, run from this folder, and take under a minute: `python3 -B ceff_free_np.py`, `tilt_check.py`, `xx_tilt.py`, `lattice_cA.py`, `avoided_cA.py`. Their outputs are quoted in `HANDOFF.md`.

## Relation to the rest of the repository

Addresses the open velocity-weighted bound (STATUS.md item 5, REPORT.md sections 9.68 to 9.86). It does not change "Where things stand": it identifies the ratio bound as the right target from the continuum side and rules out two shortcuts (monotonicity with a short anchor, Lieb-Robinson in place of rotation invariance). Spec `../../specs/C125_crossed_c_theorem.md` poses the transfer-matrix question to Astra.
