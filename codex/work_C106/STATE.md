# C106 local state

2026-10-29. Task completed as a negative result for the uniform thermal lemma, with conditional counting theorems and the unresolved operator lemma identified.

[Proved] The nearest-neighbor two-orbital chain with B = -i t sigma_z and X = epsilon sigma_x-a epsilon 1 has d = 4, two c = 1 sectors, v = 2t sqrt(1-a^(-2)) sqrt[1-epsilon^2(a^2-1)/(4t^2)], and gauge-optimized M = sqrt(4t^2+epsilon^2). Equal marginals of the extremal Gaussian bond eigenstates prove the gauge optimization.

[Proved] Tune a(epsilon) so v/M = r exactly, using x=a^2-1 and epsilon^2 x^2 + [r^2(4t^2+epsilon^2)-4t^2]x + r^2(4t^2+epsilon^2)=0. The small positive root approaches r^2/(1-r^2). Then s/T -> 2pi/(3v) at fixed epsilon in the infrared, but s/T -> pi/(3t) when epsilon << T << t. The ratio to the infrared entropy prediction approaches r.

[Proved] For any fixed matching fraction alpha in (0,1), choose r < alpha. The full low-temperature matching endpoint is O(epsilon), although d and v/M are fixed. This rules out every strictly positive universal lower bound M F_alpha(d,v/M) for that endpoint. It does not rule out a velocity-dependent matching fraction or all other entropy inequalities. A/M = 2r, so there is no counterexample to the velocity-weighted conjecture.

[Proved conditional results] Velocity-window entropy matching at T_w = kappa M(w/M)^p gives the desired integrable tail for 1 < p < 2 and a logarithmic tail for p = 1. At p = 2, a stronger entropy-increment assumption on disjoint temperature intervals gives A <= 24 M log d/(alpha pi kappa). The universal versions of these specific matching assumptions are themselves refuted by the two-band family. See PROOF.md Sections 3 and 4.7.

[Open] The general interacting N(u) bound is not proved. The operator route needs dim W_u >= eta N(u) and a surjective microscopic scaling map with k <= K(d)+beta(d) log(V/u), 2 beta(d) log d < 1. Neither requirement is proved. Nonzero overlap is weaker than independent access, and stress tensors alone count sectors rather than central charge.

[Proved file map] RESULT.md is the concise handoff. PROOF.md has the complete derivations. checks.py is the only Python program and CHECKS.md its final passing output. PROGRESS.md records the work. Project STATE.md remains unchanged because writes outside work_C106 were prohibited.
