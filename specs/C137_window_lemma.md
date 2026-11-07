# C137: The window lemma, a universal telescoping certificate at d = 2

Your folder is `../work_C137/`. Paths below are repository-relative, so prefix them with `../../` from your folder. Read `specs/_Cd_target.md` first (lessons 1 to 5, and the sine-square-deformation route set aside there), then `specs/C135_lti_rate.md` and all of `codex/work_C135/` (RESULT.md, PROOF.md section 5 and the XX certificate, NUMERICS.md, restricted_dual.py and its outputs), `contrib/marginal_toolbox_map/MAP.md`, `codex/work_C105/RESULT.md` (weighted XX inequality, tent weights, the open-interval obstruction g*), `contrib/cd_free_route/PROOF.md` sections 1 to 3 and `specs/_bound_background.md` with all updates. Rules of `_bound_background.md` apply (tag claims proved, numerical, heuristic or assumed; no probabilities in files; Python only in your folder, numpy, scipy and cvxpy if present, do not install; browsing tool only for literature; plain prose without em-dashes; PROGRESS.md, RESULT.md, PROOF.md). Your sandbox cannot reach Typhon; write any large computation as a SLURM-ready script and say so.

## The lemma

C135 reduced the target, at fixed d, to a telescoping certificate: if for every bond h on C^d (x) C^d and every n there is a Hermitian X on n - 1 sites with

    h_{12} (x) 1 + X (x) 1 - 1 (x) X >= [e_0(h) - K M(h)/n^2] 1,          (T)

then summing translates around any ring of length L >= n removes X and gives E_L >= e_0 L - K M L/n^2, hence (taking n = L) E_L >= e_0 L - K M/L and A <= (6/pi) K M. C135 proved (T) for every quadratic qubit bond with K = pi^2/2, and gave an explicit XX certificate with K = pi in which X is a weighted sum of translated bonds built from a sine window. It also found numerically that restricting X to weighted translated bonds loses nothing for XX through n = 8 and for Heisenberg through n = 6 (1.8e-6 at n = 8), and about 0.002 for critical TFIM at n = 8.

With X = sum_j w_j h_{j,j+1} plus telescoping single-site terms, the left side of (T) is an open chain with tapered couplings, sum_j c_j h~^(j)_{j,j+1}, with sum_j c_j = 1 and each bond in its own gauge representative. The window lemma is

    (W) for every bond h at d = 2 and every n there are weights c_j and gauge representatives such that
        lambda_min( sum_{j=1}^{n-1} c_j h~^(j)_{j,j+1} ) >= e_0(h) - K M(h)/n^2,

with K independent of h and n. Make the admissible class precise first (signs of c_j, which single-site terms telescope, whether the weights may depend on h).

## Why it might be provable, and the risk

A tapered open chain is a statement about one explicit finite Hamiltonian, not about crossed spectra or infinite-volume states. Several known mechanisms give 1/n^2 rates from smooth weights. (i) Sine-square deformation: for free fermions (Katsura) and in CFT the sine-square-deformed open chain has the periodic vacuum as its ground state, with the Casimir term removed. (ii) Gosset and Mozgunov improved the Knabe finite-size criterion for frustration-free chains from 1/n to order 1/n^2 by using quadratic weights. (iii) IMS-type localization with a smooth partition of unity has error of order the squared gradient of the cutoff, which is 1/n^2 for windows of width n. (iv) The momentum-box and covariance smoothings of cd_free_route and C135 are windows of this kind. The risk is lesson 2: a fixed envelope (sine-square deformation) is a local energy functional and may see near-critical massive modes at the order of the anomaly. Here the weights may depend on h, which may evade it. Decide this.

## Tasks

1. Formalize the window class and prove it is a restriction of the LTI dual (so (W) implies (T) and the target). Check the C135 XX certificate is in it.
2. Optimal windows numerically. For fixed h and n, maximizing lambda_min over the weights is a concave maximization (lambda_min is concave in c). Compute optimal windows by exact diagonalization for n up to about 16 (DMRG beyond if you can write it in the sandbox), for Heisenberg, XX, XXZ at Delta = 0.5, 2 and -0.9, -0.99, critical TFIM and TFIM at g = 1 + c/n with c in {-2, 0.5, 5}, the chiral bond of C135, and 30 random bonds on the M = 1 gauge sphere. Report n^2 (e_0 - lambda_min)/M for the optimal window, for the sine-square window and for the C135 sine window, with e_0 bracketed rigorously (LTI lower bounds from C135, Petz or MPS upper bounds). Describe the optimal windows.
3. Lesson 2. On TFIM at g = 1 + c/n and on any near-critical massive family you can construct at d = 2, decide numerically and, if possible, analytically whether the optimal-window constant stays bounded while the fixed sine-square constant grows. Explain the mechanism either way.
4. Proof attempts at d = 2, in order of promise as you judge it after reading: Gosset–Mozgunov-type weighted sum-of-squares (check what they need beyond frustration freeness), many-body IMS localization on the lattice, a deformation argument from the free certificate (C135) along a path of bonds with control of lambda_min, and any argument that uses that h-dependent weights can track the scale where massive modes appear. Prove (W) on the largest class you can (a stratum of `specs/C131_d2_cover.md`, frustration-free bonds, the U(1) class, a neighbourhood of the quadratic class) and state the exact obstruction otherwise.
5. State what is proved, for which class, and the single remaining lemma.
