# C139 result

This run was done by a Claude Opus agent, not by Astra. 2026-11-07. Proofs are in [PROOF.md](PROOF.md), the log is [PROGRESS.md](PROGRESS.md). No independent referee pass has been made.

## Bottom line

[Proved, conditional on the identification hypothesis (I)] For 0<delta<=1, every bond of the U(1) class with J_z/J in [0,1-delta], every field, DM angle and L>=2,

    L(E_(2L)-2E_L) <= (5400/delta) M,

from the sector theorem |F_L(n)-J L e_bulk(Delta,n/L)| <= (1086/delta) J/L at every 0<=n<=L/2 and both parities. The constant grows only like 1/(1-Delta), much more slowly than C134 and C136 on the attractive side.

[Proved, conditional on (M)] If the ground-state counting function is monotone on its Fermi interval, every constant is uniform in Delta in [0,1). The bound L(E_(2L)-2E_L) <= (2215+pi^2/4) M then holds on all of J_z/J in [0,1], Heisenberg included by finite-matrix continuity.

[Assumed, unresolved] Hypothesis (I) is not proved at odd L, and (M) is not proved. Both are true in every computed case.

## Task 1: the sign audit

[Proved] PROOF Section 1 lists 21 places where C134 and C136 use K<=0, with a replacement for each. The essential discovery is a hole kernel J=(I+K)^(-1)K. For 0<Delta<1 it factorizes in Fourier space into two transforms of positive, even, unimodal functions, so J=k1*r0>=0 and ||J||_1<1/2. Duality, s_Q=r0+J*(1_out s_Q), turns every positivity argument of C136 into a positive Neumann series on the complement of the Fermi interval. The continuum side comes out better than on the attractive side: r0<=r_Q<=(1+4/pi)r0 (kappa=1), 1/(1+q)<=Z_Q<=1, Var(Z_Q)<=q_J, |epsilon_Q'|<=ab r_Q, epsilon_Q'>=0, Var(v_Q)<=(2+40/pi)ab, 0<=e_bulk''<=ab(1+q)^2/2. The restricted inverse has ||(I+K_Q)^(-1)||<=3 uniformly in Q and Delta.

[Proved] On the finite side, the sign enters only through global monotonicity of z_L. A new cumulative bound ||C_mu-C_nu||<=B/(2L) uses only the ordering of the roots, at every filling. The quadrature then runs on the core, where sigma_L>=r/2 is guaranteed, with first-order bounds at the ends.

[Literature, correction to the spec] Yang-Yang convexity does not hold for Delta>0 (Kozlowski 1508.05741, Section 2). Kozlowski's introduction says Yang and Yang proved the ground state comes from a real solution continuous in Delta, with uniqueness open for Delta>0, for even length. The functional remains coercive, so an ordered solution exists. Its Hessian is 2pi L diag(sigma_L(x_j)) plus a strictly positive Gram matrix (Bochner), so sigma_L>=0 at the roots gives the identification at both parities through C134's continuation argument.

## Task 3: Heisenberg

[Proved, tracked] The only non-uniform estimate as Delta->1 is the cumulative root bound. Its factor B=pi/(2gamma) comes from the fixed-point inequality ||D||<=1/(2L)+q||D||. All continuum constants are uniform, and so is the true inverse norm.

[Proved, conditional on (M)] With a monotone counting function, the sawtooth has sup 1/(2L), and the uniform restricted inverse gives ||G||<=27/(4L). This yields the uniform constant above and decides the Heisenberg endpoint. An O(1/L) bound is not obstructed by the marginal operator.

[Assumed] The single remaining estimate is (M): sigma_L=a1-(1/L)sum_k K(.-x_k)>=0 on [-Q_L,Q_L] for the ground-state roots, at 0<Delta<1. A weaker form, a negative mass of at most c/L, suffices for uniformity given (I).

## Numerics

[Numerical] ED matches Bethe for every sector at L<=12 and Delta in {0.1,0.5,0.9,0.99,0.999}, with maximum error 6.2e-14. There were 572 finite-L checks of every intermediate inequality, at Delta in {0.1,0.5,0.9,0.99,0.999}, L up to 2048 (doubled 4096), with fillings 1/L to 1/2. L||D|| lies in [0.50,0.522] against the proved B/2, which reaches 17.6. L|F-Le_bulk|/K_rep<=9e-4. sigma_L>0 holds on R in every case. The Casimir value -pi v/6 is reproduced. Delta=1/2 reproduces -3L/8 (odd L) and the twisted identity (even L).

[Numerical] Half-filled defect L(F_(2L)-2F_L)/M_0:

| L | 0.5 | 0.9 | 0.99 | 0.999 | 1 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 3.255 | 3.864 | 3.987 | 3.999 | 4.000 |
| 8 | 2.100 | 2.471 | 2.550 | 2.557 | 2.558 |
| 1024 | 2.04053 | 2.38514 | 2.46065 | 2.46810 | 2.46892 |
| 4096 | 2.04053 | 2.38481 | 2.46007 | 2.46752 | 2.46834 |
| pi v/2 | 2.04052 | 2.38459 | 2.45917 | 2.46658 | 2.46740 |

The defect decreases in L at every Delta. Its supremum is 4, at L=2 and Delta=1.

## Files

Solvers are bethe_scan.py and check_c129.py (copied from C136) and continuum_rep.py, the restricted dressing. The checks are jkernel.py, check_continuum.py, check_finite.py with summarize_finite.py, ed_check.py, heis_scan.py, sigma_sign.py and constants.py. Their outputs are the matching .json files. Typhon was not used.
