This run was done by a Claude Opus agent (2026-11-07), not by Codex.

# C150 result: the ferromagnetic strip by comparing L with 2L directly

Proofs are in [PROOF.md](PROOF.md), the log in [PROGRESS.md](PROGRESS.md). No independent referee pass. J = 1, eta = pi - gamma, Delta = -cos(eta), a = pi sin(eta), D = F_2L(2n) - 2F_L(n).

## Bottom line

[Assumed, not closed] The strip is not closed unconditionally. On this route it is reduced to one first-order, integrated statement about the L roots on the shell:

    (H_I)   int over S' of dy/(eta L sigma_L(y)) <= C_I,     S' = {2 <= eta L r_Q < 2 C_D(eta)}.

[Proved, conditional on (H_I)] For 0 < eta <= eta0 = 1/2, every L >= 2 and 0 <= n <= L/2,

    L |F_2L(2n) - 2F_L(n)| <= 112.2 + 4.87 C_I.

(H_I) follows from sigma_L >= c_* r_Q on S' with any c_* > 0 (C_I = 3.18/c_*). C145's (G_c) (c_* = 3/4) gives L|D| <= 133, two-sided. No derivative bound, no variation of r/sigma and no log log are needed, unlike C145's (G) and C146 Proposition 3.1.

[Proved, unconditional] For 0 < eta <= 1/2, L|D| <= 85.8 log(1/eta) + 180, and for every Delta in (-1, 0], |D| <= pi sin(eta) n/L.

## The decomposition and the key lemma

[Proved] Quantum numbers: t'_(2j-1) = t_j - 1/(4L), t'_(2j) = t_j + 1/(4L), I'_(2j-1) + I'_(2j) = 4I_j. On doubled configurations S_2L(x, x) = 4S_L(x) - n Th(0), E_2L(x, x) = 2E_L(x), and the Yang-Yang gradient at the doubled L ground state is +-pi on each pair. Trial yt: the refined grid Y(t_j -+ 1/(4L)) (Y = z_L^(-1)) on a central block of cells, doubled roots (x_j, x_j) outside. Then D = A + B, A = E(yt) - 2E_L(x) (a Richardson difference of midpoint rules in the counting coordinate, zero outside the block) and B = E(y*) - E(yt).

[Proved] Key lemma (homotopy). For every trial yt, |E(y*) - E(yt)| <= 2a sum_k |rho_k|, rho_k = z_2L[yt](yt_k) - t'_k. Proof: follow the minimizers of S_2L - 4pi L(1-s) rho.y; the Jacobian J(y) is a symmetric M-matrix with row sums a1 at every configuration, and |e'| <= 2a a1, so |J^(-1) e'| <= 2a. The 1/eta amplification of (I - T)^(-1) is absorbed exactly; there is no linearization error.

[Proved] The block residual is a second difference of theta(x - Y(t)) (a T-convolution of the pair sawtooth). Per-cell Euler-Maclaurin with the kernel v, ||v|| = 9/512, gives the master bound (Theorem 4.4)

    L|D| <= a[N_sp + (H+1)/4] + (a/eta)(9/512)(2 C_3 H + C_3e) I_1,

with N_sp the number of doubled cells, I_1 = int over the block of dy/lambda_L, and local constants H = 2.58, C_3 = 12.8, C_3e = 22.0 at eta = 1/2 (1.81, 10, 10 as eta -> 0). N_sp <= 30.87/eta + 2 is proved by C136's cumulative bound and C138's tail lemma. I_1 is controlled on the deep core by C142 Theorem II.6, and on the shell it is (H_I).

Proved versus numerical: everything above is proved; the suprema defining tau_2, alpha_i and H are of explicit one-variable functions evaluated on fine grids (constants_C150.py), increasing in eta.

## Numerics

doubling.py, 216 cases: 1 + Delta in {0.1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6}, L in {8, 32, 128, 512, 2048} paired with 2L up to 4096, eight fillings from L/64 to L/2 (summary_C150.json).

| 1+Delta | max L D | L D / (pi v/4) at half filling, largest L | max (4.1) with actual N_sp, I_1 | max eta N_sp | max I_1 | min lambda_L/lambda_Q on S' |
| --- | --- | --- | --- | --- | --- | --- |
| 0.1 | 0.200 | 1.0000 | 16.6 | 4.51 | 0.82 | 1.0002 |
| 1e-2 | 0.058 | 0.99997 | 15.5 | 4.25 | 0.95 | 0.99995 |
| 1e-3 | 0.0178 | 0.9997 | 15.4 | 4.20 | 0.97 | 0.9993 |
| 1e-4 | 0.0056 | 0.997 | 15.3 | 4.13 | 0.89 | 0.9998 |
| 1e-5 | 0.0017 | 0.968 | 14.4 | 4.58 | 0 (no block) | no shell sampled |
| 1e-6 | 0.0004 | 0.703 | 4.6 | 1.45 | 0 (no block) | no shell sampled |

Proved ingredients checked numerically in every case: D = A + B to rounding; max|J^(-1)e'|/(2a) <= 0.584 (about 0.51 at small eta, 0.995 at Delta = 0); pair half-width / delta <= 0.99999.

Consistency checks. [Proved] Delta = 0: rho = 0 and B = 0 exactly, and Euler-Maclaurin gives L D -> (pi/4) sin(pi nu), the CFT value; [Numerical] L D = 0.790 at L = 16. [Numerical] CFT: L D/(pi v/4) = 0.99997 at 1 + Delta = 1e-2, L = 2048. [Proved and Numerical] eta -> 0: D -> eps 2n(L-n)/((2L-1)(L-1)) (log-gas sum rule), matched to 1.4e-3 relative at eta L < 0.2. At small eta the Casimir energy sits in B, not in A (L A_ref = 0.0002, L B_ref = 0.058 at 1e-2, L = 1024).

## What remains

(H_I), with numerical margin: sampled lambda_L/lambda_Q >= 0.9993 on S' (any c_* > 0 suffices), and I_1 <= 0.97 against the 6.45 used for C = 133. The shell is only reached at 1 + Delta >= 1e-4 for L <= 2048; C145 sampled it at 1e-6 up to L = 16384 (min sigma/r = 0.995). The bounds of the route hold for any configuration with counting function z_L, so the shell density must come from the L Bethe equations themselves; Theorem II.6 supplies it only with the log of Corollary 5.3. The sparse count (97 of the 133) is crude by about a factor 7 against the sampled eta N_sp.

Files: doubling.py, summarize.py, constants_C150.py, sym_check.py, doubling_A.json, doubling_B.json, summary_C150.json, constants_C150.json, log_*.txt; bethe_ferro.py and continuum.py copied from work_C145. All runs local, each scan under fifteen minutes; Typhon was not used.
