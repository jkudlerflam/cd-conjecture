# C149: The ferromagnetic strip by a discrete comparison principle

Your folder is `../work_C149/`. Read `specs/_strip_context.md` first; it states the target and what is known.

## The idea

On the attractive side the Bethe map F(x)_j = L p(x_j) − Σ_k θ(x_j − x_k) − 2π I_j has ∂F_j/∂x_k = θ'(x_j − x_k) = −2π T(x_j − x_k) ≤ 0 for k ≠ j, and diagonal 2π L σ_L(x_j) − 2π T(0) with row sums 2π L a₁(x_j) > 0. So the Jacobian is a strictly diagonally dominant Z-matrix, a nonsingular M-matrix, everywhere, and F is an off-diagonally antitone, diagonally isotone map (an M-function in the sense of Rheinboldt, J. Math. Anal. Appl. 1970; Ortega–Rheinboldt Ch. 13.5). Then F(y) ≥ 0 = F(x) componentwise implies y ≥ x, and F(y) ≤ 0 implies y ≤ x. The current proof fails because second-order quadrature needs regular spacing of the true roots, which cannot be bootstrapped. With a comparison principle, regularity is needed only for explicit trial configurations, which can be designed to be regular.

## Tasks

1. State and prove the comparison principle for the ground-state Bethe system at flux Lπ in every sector, with the exact hypotheses (ordering of quantum numbers, real roots, the domain), citing or reproving Rheinboldt's theorem.
2. Build explicit sub- and supersolutions y^± near the ground state, for example the quantiles y_j of corrected densities r_Q ± ψ (z̃(y_j) = I_j/L with z̃' = r_Q ± ψ), and compute the residuals F(y^±). The residual is a quadrature error of the kernel sum on a configuration whose spacing you control, plus a linear term in ψ. Choose ψ (through (I − T_Q)ψ = bound on the quadrature error, or a better localized choice) so the residuals have the required signs, with explicit control in the sparse shell and tails. Watch the amplification by (I − T_Q)^(−1): it is of order 1/η for uniform sources but smaller for sources localized near the edge (Z_Q(Q) is of order (min(Q,1)/η)^(1/2), C146).
3. Convert the root sandwich y⁻ ≤ x ≤ y⁺ into what the energy estimate needs: either the shell comparability σ_L ≥ (3/4) r_Q (then C146's assembly gives C_strip = 773), or directly a bound on the energy error, or directly on the matched defect L·D.
4. Assemble the strip theorem with explicit η₀ and C, and check every inequality numerically down to 1 + Δ = 1e-6.
