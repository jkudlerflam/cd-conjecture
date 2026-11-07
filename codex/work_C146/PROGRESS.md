# C146 progress log

Run by a Claude Opus agent (2026-11-07), not by Codex. J=1 throughout. Notation of C145 PROOF Section 1.

## Step 0
Read specs/C146_shell_and_edge.md, specs/_Cd_target.md, specs/_bound_background.md, all of work_C145 (RESULT, PROOF, scripts), work_C142 PROOF Part II, work_C138 RESULT and PROOF, work_C136 RESULT and REFEREE. Copied continuum.py, bethe_ferro.py, check_attractive.py from work_C145.

## Step 1, analysis
- (W3'). Idea: W(Q) = int_{-Q}^{Q} Z_Q satisfies W'(Q) = 2 Z_Q(Q)^2 (domain-derivative identity), and Z_Q(Q) is nondecreasing in Q (translation plus positivity), so W is convex with W(0) = 0 and Z_Q(Q)^2 >= W(Q)/(2Q). A variational lower bound on W = <1,(I-T_Q)^(-1)1> with the semicircle trial function and the symbol bound 1 - T^(w) <= eta(|w| + 2/pi) gives eta Z_Q(Q)^2 >= pi Q/(4(1 + 16Q/(3 pi^2))). This is the Wiener-Hopf square-root edge in its natural lower-bound form, uniform in eta and Q.
- (W1). Freezing v(x0) in u' - v(x0) r' gives an exact first-order linear identity v' = F_1 - (r(Q)/r)[(v(Q)-v)G + (v(Q)+v)H], with G, H the edge Green's functions and F_1 = a[(R a1')^2 - (R a1)(R a1'')]/(R a1)^2. Only Var(v) <= const (not const * eta) is used in Theorem A.2. The edge terms integrate to a multiple of a(Z(Q) - 1) <= a(B-1) <= pi^2/2, F_1 is bounded below by -4a everywhere and is nonnegative up to an exponentially small error away from the edge (reflection bound of C145 Theorem B.1 applied to derivatives). Gives Var(v_Q) <= V_W with an explicit V_W.
- (G). Not closed; see later steps.

## Step 2, numerics
- check_W3.py (161 cases, 1+Delta 0.5 to 1e-6, Q 0.003 to 4, N = 900, about 15 minutes): Z_Q(Q) monotone in Q, dW/dQ = 2Z^2 to 1.5 percent, eta Z^2 / bound >= 1.578, min eta Z^2/tanh 2Q = 0.808.
- check_W1.py (30 cases): identity of Lemma 2.1 to 1e-3 (1.3e-2 mesh-limited at 1e-6, Q >= 2), int_0^Q (G+H) = Z(Q)-1, F_1 >= 0 in all resolved cases, Var(v)/a <= 3.0, Theorem 2 bound loose by 1.02 to 215.
- constants_C146.py: c_Z = 0.373, C_e = 8.43, I(b) <= 5.998, V_W <= 92.87, K = 154.0, C_strip = 772.5 (conditional on (G)); loglog table for Proposition 3.1.

## Step 3, (G)
- Not proved. Proposition 3.1: derivative half from comparability half with a log log(1/eta) loss through the shell length. Residual (G_c): sigma_L >= (3/4) r_Q on the shell. The log in the first-order charge sits in the sparse cells' first moments; a local bound turns it into log log, removing it needs spacing regularity of sparse roots.

## Step 4, final
PROOF.md and RESULT.md written.
