# C150: The ferromagnetic strip by comparing L with 2L directly

Your folder is `../work_C150/`. Read `specs/_strip_context.md` first; it states the target and what is known.

## The idea

The matched defect never needs the continuum. If each root of the (L, n) system is doubled, the kernel sum doubles and the counting function is exactly z_L. The (2L, 2n) roots then sit, to leading order, on the refined grid z_L(x) = I'/(2L). So D = F_{2L}(2n) − 2F_L(n) splits into (a) a midpoint-rule difference of one function at steps h and h/2, in which the first-order (O(1)) terms cancel identically and the leading term is the edge second-order term, the Casimir energy, and (b) a self-consistency correction from the change of counting function z_{2L} − z_L, itself a difference of two kernel sums. In the sparse shell both systems have comparable structure, so the errors that defeat the discrete-to-continuum comparison may cancel in the difference.

## Tasks

1. Write the exact decomposition of D along these lines at flux Lπ, every sector, with every term explicit, including the precise relation between the quantum numbers of (2L, 2n) and (L, n).
2. Bound the self-consistency correction. Natural tools: the linearization of the 2L Bethe equations around the refined-grid configuration (its Jacobian is a diagonally dominant M-matrix on the attractive side, so its inverse is positive and controlled), convexity of the Yang–Yang action (C136), and the counting-function monotonicity. Show how the 1/η amplification of (I − T)^(−1) acts on a source that is itself a difference of two quadratures.
3. Bound term (a) uniformly in η, including the sparse shell and tails, using that it is a difference of midpoint rules for the same function (Richardson structure) rather than a quadrature error against the continuum.
4. Assemble L·D ≤ C J on the strip with explicit η₀ and C, and check every term numerically down to 1 + Δ = 1e-6, L up to a few thousand, all fillings.
