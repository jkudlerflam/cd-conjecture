# Independent audit of the ultraviolet normal form

**Status: checked algebra and locality counting; no closure of vacancy matching.**

I independently reviewed [UV_NORMAL_FORM.md](../vacancy/UV_NORMAL_FORM.md), including the sign of the filter, the finite-depth support growth, and the weighted remainder estimate. I found no error in the stated conditional theorem.

The distributional identity is consistent: with the author's Fourier convention, \(F'=-i\delta_0+ig\) and \([\tau_u(v),H_0]=i\partial_u\tau_u(v)\), integration by parts gives \([S,H_0]=-V+W\), including the displayed endpoint signs.

A layer of disjoint radius-\(R\) gates enlarges the support interval by at most \(O(R)\). There are \(O(R)\) layers, so the causal interval has length \(D=O(R^2)\). Each layer contributes at most \(O(D/R+1)\) relevant gates, yielding \(O(R^2)\) gates in total. Each generator has norm \(O(\delta^{-1})\). The second derivative is therefore bounded by \(O(R^4\delta^{-2})\) per original anchored term. The interaction norm introduces the separate overlap and diameter factors \(D^2=O(R^4)\), giving \(O(s^2R^8\delta^{-2})=O(s^2\delta^{-18})\). None of these steps requires a small global operator norm or a finite-volume gap.

Infinite local Hilbert spaces do not obstruct this argument when every local interaction used is bounded. The proof uses only bounded-operator norms and the dimension-independent finite-range Lieb–Robinson estimate. An application containing an unbounded on-site chemical potential or number term would require a separate check; the displayed bounded phase hopping and contact projectors do not have that problem.

The Gaussian filter retains all sufficiently low **energy transfers**, including transfers between high-energy states. It is not a low-energy projector and cannot be replaced by a vacuum expectation. The source note already states this correctly.

One application qualification is essential. A common geometric coloring and layer order makes transformed remainders identical when the original anchored coefficient arrays agree on their entire \(O(R^2)\) causal neighborhoods. It does **not** automatically make translated or rephased cell bulks coincide. A residue coloring or greedy coloring need not be covariant under the translation used to identify those bulks. A per-interface ground-energy comparison must therefore specify a common anchored identification, or supply an extra comparison between the corresponding circuits. The source note explicitly leaves bulk cancellation and low-frequency ground-energy control open, so this does not invalidate its theorem.

Reviewed 2026-10-28 by a separate GPT-6-based Codex agent. This was an independent algebraic and locality audit, not a new numerical test or a human review.
