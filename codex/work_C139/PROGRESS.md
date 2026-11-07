# C139 progress log

Run by a Claude Opus agent (not Astra). J=1, Delta=cos gamma, 0<gamma<=pi/2 on the repulsive side.

## Done

1. Read spec, _Cd_target.md, C134 (RESULT, PROOF, REFEREE), C136 (RESULT, PROOF), C129 RESULT, C138 RESULT and PROGRESS, TOOLBOX sections 4 and 7.
2. Copied bethe_scan.py and check_c129.py from work_C136 (solver handles 0<Delta<1 and Delta=1 already).

## Key early observation

The C134 cumulative root bound can be rederived without the kernel sign and without monotonicity of the counting function: at each root the cumulative error G jumps between 1/(2L)+H and -1/(2L)+H exactly, G decreases between roots because rho>0, and H=-K*G has sup norm at most q sup|G|. Only the ordering of the roots is used.

## Next

Numerics on the sign of sigma_L (globally and at roots) for Delta in (0,1); then the identification of the ground branch for Delta>0.

## Update 2 (findings so far)

3. sigma_sign.py: sigma_L>0 on the whole line and at roots in every tested case, Delta in {0.1,0.5,0.9,0.99}, L to 1024, four fillings. Its e^{-2|x|} tail coefficient sin(g)-sin(2g)(1/L)sum e^{2x_k} is positive but tends to zero at half filling like L^{-q}.
4. Literature: Kozlowski 1508.05741 (intro, Prop 2.1) says Yang-Yang uniqueness holds only for -1<Delta<=0; the YY functional is not convex for Delta>0 (still coercive, so a minimizer exists). YY did show the ground state comes from a real solution continuous in Delta. So the spec's "Yang-Yang convexity for all -1<Delta<1" is not available.
5. Key new tool: the hole (dual) kernel J=(I+K)^{-1}K has Jhat=[sinh((B-C)w)/sinh(Bw)][1/(2cosh(Cw))], B=(pi-g)/2, C=g/2, a product of two transforms of positive even unimodal functions, so J=k1*r0>=0, even, unimodal, ||J||_1=(pi-2g)/(2(pi-g))<1/2. jkernel.py confirms numerically. Consequences (proved): r0<=r_Q<=(1+4/pi)r0, 1/(1+q)<=Z_Q<=1, |eps_Q'|<=ab r_Q, eps_Q'>0 on (0,Q], Var(v_Q)<=2ab+40ab/pi, e_bulk''<=ab(1+q)^2/2. The continuum side is better than C136 (kappa=1).
6. Plan: Theorem A (Delta in [0,1-delta], all fillings) with sign-free cumulative bound (B/(2L)), core trapezoid where sigma_L>=r/2, first-order ends. Conditional only on identification (I). Theorem C: under global monotonicity (M) all constants are uniform up to Delta=1.

## Update 3 (final)

7. continuum_rep.py and check_continuum.py: all dual-kernel continuum lemmas pass (18 Q values, Delta 0.1 to 0.99); underresolved large-Q cases at large b are skipped explicitly.
8. check_finite.py (572 cases, Delta 0.1 to 0.999, L to 2048, doubled 4096), ed_check.py (all sectors L<=12, max error 6.2e-14), heis_scan.py (half filling to L=4096, Delta to 1). Summary in finite_summary.json.
9. constants.py: K_rep(gamma) explicit, sup over delta of delta*sup K_rep=1027, analytic majorant 1086/delta.
10. PROOF.md and RESULT.md written: Theorem A (conditional on identification (I)), Theorem C (conditional on monotone counting function (M), uniform through Heisenberg), sign audit of 21 items.

Open: (M), and (I) at odd L. No referee pass yet.
