# C145 progress log

Run by a Claude Opus agent (2026-11-07), not by Codex. J=1 throughout.

## Step 0
Read the spec, _Cd_target.md, C142 (RESULT, PROOF, PROGRESS), C138 (RESULT, PROOF), C136 (RESULT, REFEREE, PROOF Sections 1-5), C134 (RESULT, REFEREE). Copied bethe_ferro.py, continuum.py, check_attractive.py from work_C142.

## Step 1, analysis (route (c) first)
- Theorem B cost split: core term with M = 4C_D(eta) is bounded (ratio C_D1(eta)/C_D(eta) is bounded), the log sits only in the first-order tail cost 4aM(1+2/(kappa b))/eta, i.e. on the shell {M_0 <= eta L r < 4C_D(eta)}. The shell holds about M/eta roots, not few, so it cannot be charged first order. Route (c) therefore reduces to second-order control on the shell only.
- (W2) PROVED by a reflection supersolution: r_0 - r_Q <= r_0(2Q-x) + r_0(2Q+x) on [-Q,Q], so r_Q >= (1 - 4 e^{-2b d}) r_0 at distance d from the Fermi point. kappa = 1-4e^{-2} = 0.459 at d=1, uniform in eta and Q.
- Referee of C138 Theorem B: every line re-derived, no error found (details go to PROOF Part C).
- (W1), (W3) both hinge on an edge cancellation between r_Q'/r_Q and Z_Q'/Z_Q at the Fermi point (log-derivative of omega = r_Q/Z_Q). Not proved here.

## Step 2, numerics for (W2)
check_W2.py (30 cases, eps 0.5 to 1e-6, Q 0.3 to 6, N=900): the reflection bound and r_Q >= (1-4e^{-2bd}) r_0 hold everywhere except at Q >= 4 for eps <= 1e-4, where the apparent violation (up to 0.021) comes with max r_Q/r_0 > 1, which contradicts the proved r_Q <= r_0. check_W2_refine.py: under mesh refinement N = 900, 1800, 2700 these violations shrink as 0.0105, 0.0025, 0.0007 (eps=1e-6, Q=4) and 0.0020, 0.00024, 0.00005 (eps=1e-4, Q=6), so they are mesh artifacts. min r_Q/r_0 at d >= 1 is 0.762 (proved 0.459).

## Step 3, Part A written, shell numerics
PROOF Part A: Theorem A.2 replaces (D), (D1) by the shell hypothesis (G) (sigma >= 3r/4 and eta Var_S(r/sigma) <= C_GS on S = {2 <= eta L r < 4C_D(eta)}); the core piece on {eta L r >= 4C_D(eta)} is proved log-free, eta Var_K1(r/sigma) <= C_G1 <= 13.19 (10.11 numerically).
check_shell.py, check_shell_big.py: eps 1e-2, 1e-4, 1e-5 with L to 8192 (and 1e-6 at L=16384 running). On K_0 = {eta L r >= 2}: min sigma/r >= 0.995 (the 0.995 is at n = L/2-2, eps=1e-5, where the continuum is mesh-limited; elsewhere >= 1.0000), eta Var_K0(r/sigma) <= 0.0018. So (G) holds with comparability margin 0.995 against the needed 3/4 and C_GS <= 0.002 in the sampled cases. At eps=1e-6 the set K_0 is empty for L <= 8192 (eta L r(0) < 2).
check_W2 done (Step 2). check_W13.py running.

## Step 4, final
check_W13.py done: v_Q monotone in 55 of 56 cases (one mesh-level dip), min eta Z(Q)^2/tanh(2Q) = 0.808. constants_C145.py: strip constant 686 conditional on (G), (W1'), (W3') (C_GS = 1, c_Z = 3/4), 661.6 with sampled inputs; C138's Theorem B constant with proved (W2) is 128.99. PROOF.md Parts A to E and RESULT.md written. Outcome: log not removed; remaining inequality isolated as the shell bound (G) plus (W1'), (W3'); (W2) proved; C138 Theorem B refereed, correct.
