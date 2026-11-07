# C146: The last three inequalities on the ferromagnetic strip

Your folder is `../work_C146/`. Paths below are repository-relative, so prefix them with `../../` from your folder. Read `specs/_Cd_target.md`, then all of `codex/work_C145/` (RESULT.md, PROOF.md, scripts: your starting point), `codex/work_C142/` PROOF.md Part II, `codex/work_C138/` (RESULT.md, PROOF.md) and the RESULT.md and REFEREE.md files of `codex/work_C136/`. Rules of `_bound_background.md` apply (tag claims proved, numerical, heuristic or assumed; no probabilities in files; Python only in your folder; browsing tool only for literature; plain prose without em-dashes; PROGRESS.md, RESULT.md, PROOF.md).

## Where things stand

On the attractive side −1 < Δ = −cos η ≤ 0 of the XXZ ring, C145 Theorem A.2 gives L(E_{2L} − 2E_L) ≤ 686 M at every field, DM angle and length, conditional on three statements. Everything else is proved. In C145's notation (Section 1 of its PROOF):

(G) On the shell S = {x ∈ [−Q, Q] : 2 ≤ λ(x) < 4C_D(η)}, λ = ηL r_Q, with C_D(η) = (4/π)(log(1/η) + π/2) from C142 Theorem II.6: σ_L ≥ (3/4) r_Q and η·Var_S(r_Q/σ_L) ≤ C_GS = 1. Numerical: min σ/r = 0.995 and η·Var ≤ 0.0018 (36 cases, L to 16384, 1 + Δ to 1e-6). The minimal energy-level statement behind it is |Σ over cells in S of [ε_Q(x_j) − L ∫_{X_j} ε_Q σ_L]| ≤ C/L (C145 PROOF line near "minimal energy-level statement"). It is implied by C138's (D) with C_D ≤ 1/2 on {λ ≥ 2} (comparability half).

(W1') v_Q = ε_Q'/r_Q is nondecreasing on [0, Q] for every η and Q. Gives (W1) with C_W = 4π. Numerical: 55 of 56 cases, the exception a 3.4e-6 mesh-level dip.

(W3') η Z_Q(Q)² ≥ c_Z tanh(2Q) for every η and Q, with c_Z = 3/4, Z_Q = (I − T_Q)^(−1) 1. Gives (W3) with C_e = π/c_Z. Numerical minimum 0.808, approached in the dilute window as a function of Q/η alone (a limit function of Q/η with infimum near 0.8). C145 notes the heuristic Wiener–Hopf square-root edge, Z_Q(Q) of order (min(Q,1)/η)^(1/2).

## Tasks

1. (W3'). Prove the edge lower bound, or any replacement with an explicit c_Z > 0, uniformly in η and Q. In the η → 0 window the restricted equation tends to a limit equation in the variable Q/η, so a natural route is: (a) identify that limit kernel and prove the bound there by Wiener–Hopf or by an explicit subsolution; (b) control the approach to the limit with an explicit error; (c) cover Q of order 1 and larger separately (where C145's reflection idea, or a subsolution built from r_0 and Z_∞ = 1/(1 − q), may suffice). The positivity of (I − T_Q)^(−1) and the derivative identity (I − T_Q) Z' = −Z(Q) J_0 of C145 Section B are available.
2. (W1'). Prove monotonicity of v_Q on [0, Q], or replace it by the bound Var(v_Q) ≤ C_W η that it is used for, with any explicit C_W. A total-variation bound may be easier than monotonicity. The edge-cancellation structure of C145 (r'/r and Z'/Z carry the same edge-singular part) and C136 (5.1) are the starting points.
3. (G). Prove the shell estimate, or the minimal energy-level statement behind it, with explicit constants. C145 and C142 show that first-order charging and the naive bootstrap fail there by a log. Ideas: a second-order bound on the shell that uses (W2) (proved) and (W3) for the sparse side, or a comparison of the shell's Bethe roots with an explicit dilute-gas solution in the scaling variable Q/η, since on the shell the density is low and the problem is close to the dilute limit.
4. Assemble the strip theorem with explicit constants for whatever you prove, and say precisely which hypothesis is left if one is not proved, with its numerical margin.
5. Check every new inequality numerically (existing solvers in work_C145). Keep numerics on the critical path.
