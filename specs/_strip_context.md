# Shared context for C149 to C151: the ferromagnetic strip

Paths are repository-relative. Read `specs/_Cd_target.md`, then `codex/work_C146/` and `codex/work_C145/` (RESULT.md, PROOF.md), `codex/work_C142/` PROOF.md Part II, `codex/work_C138/` (RESULT.md, PROOF.md, Theorems A and B), `codex/work_C147/REFEREE.md`, and the RESULT.md and REFEREE.md files of `codex/work_C136/`. Copy Bethe solvers from `codex/work_C145/` or `codex/work_C142/` into your own folder.

The target. For the XXZ ring H = J Σ (SˣSˣ + SʸSʸ + Δ SᶻSᶻ) at flux Lπ, in every magnetization sector, with Δ = −cos η and η → 0 (the attractive side next to the ferromagnet), prove

    L [F_{2L}(2n) − 2 F_L(n)] ≤ C J       for all η in (0, η₀], L ≥ 2, 0 ≤ n ≤ L/2,

with explicit C and η₀ > 0. Together with C129's reductions (field cancels in matched sectors, DM is a twist costing at most π²M/4, M ≥ J/2), this closes the U(1) class with C = 10781 everywhere. It is enough to bound the matched defect directly. A sector theorem |F_L(n) − L e_bulk| ≤ K/L is one way, not the only way.

What is known. (i) C138 Theorem A, unconditional: |F_L(n) − L e_bulk| ≤ (5π/4) η J, so L·D ≤ C ηL J, fine for ηL bounded; the C131 cover lemma gives L·D ≤ 4L²(1 + Δ)J·const similarly. (ii) C142 Theorem II.6: ηL |σ_L − r_Q| ≤ (4/π)(log(1/η) + π/2) everywhere on the occupied interval. (iii) With (W2), (W3') proved and (W1') replaced by a proved variation bound (C145, C146), the only missing input of the current proof is the shell comparability σ_L ≥ (3/4) r_Q on S = {2 ≤ ηL r_Q < 4C_D(η)}. The log enters because first-order charging of the shell costs about π λ_max/L with λ_max = 4C_D(η), and second-order charging needs regular spacing of the actual roots, which a bootstrap from (ii) cannot supply. (iv) Numerically ηL sup |σ − r| ≤ 0.32, min σ/r on the shell is 0.995, and the true L·D is of order η J: there is a factor of order 1/η of slack against the target.

Notation (C136, C138, C142): q = 1 − 2η/π, T ≥ 0 even unimodal with ∫T = q, T(0) = cot(η)/π, a₁ the bare density, σ_L = a₁ + (1/L) Σ_k T(· − x_k) the counting density, z_L the counting function with z_L(x_j) = I_j/L, r_Q = (I − T_Q)^(−1) a₁, Z_Q = (I − T_Q)^(−1) 1, λ = ηL r_Q the number of roots per kernel width.

Rules of `specs/_bound_background.md` apply (tag claims proved, numerical, heuristic or assumed; no probabilities in files; Python only in your folder; plain prose without em-dashes; PROGRESS.md, RESULT.md, PROOF.md). Numerics only on the critical path.
