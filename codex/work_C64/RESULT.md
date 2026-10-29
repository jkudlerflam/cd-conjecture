# C64 result

[assumed, unresolved target] The global cut-coverage construction remains open. I did not prove C60 (17) or exhibit a legal nonvanishing tree violating it. This run must not be cited as completing Lemma N.

[proved, local spatial realization] [PROOF.md Sections 2-5](PROOF.md) establish an exact weighted-line substitution and embed it in a quartic-to-quadratic GN/BBF configuration. Route the parent moment before transposing the child remainder. The child then differentiates and relocates the entire weighted covariance. For the case with one child derivative and two parent moments, the two terms are

    r integral_0^1 [2 Z_u g_s(Z_u)+Z_u^2 g_s'(Z_u)] du.

Both have F_child=1 and F_parent=2. When the derivative consumes a parent moment, its child path is concatenated with a parent path. The resulting gain is

    2^[-(t-s)-2(s-h)].

The same calculation works for child order c<=2 and line weight b<=2. Every Leibniz term has counts c and b and the same scale exponent. Its local coefficient mass is at most seven in one coordinate.

[proved, local repair] Keeping the weight at its old argument while relocating only the covariance changes the expression by

    (Z_0^b-Z^b) g_s(Z_0).

This is generally nonzero. Distributional transposition proves the whole-weight rule before interpolation. In the displayed BBF configuration it also preserves the anchor Jacobian, the moment bound and the determinant representation. This repairs that local shortcut, not the complete routing.

[assumed, remaining step] A general quartic child can retain several ports in an ancestor Taylor field. The child operator can then act on that exterior field as well as on covariance factors. C61's normal form gives algebraic counts, but a simultaneous spatial construction with bounded loads and total coefficient mass C^(sum_i p_i+n) is still missing. The restricted substitution proved here assumes and explicitly verifies the absence of those additional coordinate dependencies.

[numerical, finite exact-arithmetic controls] `python3 -B checks.py` passes 117 weighted Leibniz identities, 390 anchored Taylor identities, 195 operator-error identities, 396 two-cut incidence cases and 24 BBF determinant checks. These checks do not test arbitrary GN routing or interacting RG.

[record] PROGRESS.md, PROOF.md, STATE.md, checks.py and checks.json are local to work_C64. No files outside this folder were modified. No subagents, large computations or independent review were used. No source, thermal or central-charge conclusion is claimed.
