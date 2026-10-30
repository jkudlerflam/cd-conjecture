# Ring-aware bootstrap bounds on the Casimir coefficient

## Claim

Setting. A one-site translation-invariant nearest-neighbour chain H = sum_i h_{i,i+1} on a ring of L sites with local dimension d and a real bond h, with Casimir coefficient A = sum_a c_a v_a defined by E_0(L) = L e_0 - pi A/(6L) + o(1/L), and M the gauge-optimized bond half-spread. Bonds tested: Heisenberg S.S (d = 2, M = 1/2, A = pi/2), XXZ at Delta = 1/2 (d = 2, M = 1/2, A = 3 sqrt(3)/4), the SU(3) swap (d = 3, M = 1, A = 4 pi/3), and the seeded real random bonds of C117 at d = 2 (M = 0.5916, gapped) and d = 3 (M = 2.1047, ED suggests a very small A).

A translation-invariant moment bootstrap for the ring ground energy was built that contains ring-wide objects: two-block marginals at every separation r of the ring, and positivity of <O_q^dag O_q> and <O_q^dag [H, O_q]> for the momentum modes O_q = sum_j e^{iqj} O_j at every q = 2 pi m/L, together with stationarity and positivity of a contiguous k-site marginal. Its optimum ell(L) is a lower bound on E_0(L). The findings are the following.

1. [numerical] With one-site operators in the momentum modes (s = 1), the bootstrap does not resolve the Casimir term. The bound per site converges by L of about 20 to a constant ell_inf that sits 0.012 to 0.022 below e_0 for the critical d = 2 chains (0.07 below for SU(3)), and its own 1/L coefficient is small and of the wrong sign (between -0.4 and 0). The certified quantity A_cert(L) = (6L^2/pi)(e_0 - ell(L)/L) grows like L^2 (to 98 at L = 48 for Heisenberg, against A = 1.571).

2. [numerical] Adding the energy-density modes h_q (two-site operators, s = 2) changes the picture qualitatively. The bound then tracks the finite-size spectrum: it agrees with E_0(L)/L to 1e-5 per site up to L = 10 (1e-6 at k = 7), the gap E_0(L) - ell(L) is 2 percent of the Casimir energy pi A/(6L) at L = 12 and 20 percent at L = 20, and the bound's own 1/L coefficient is A_relax = 1.46 (Heisenberg, true 1.571) and 1.20 (XXZ, true 1.299), 0.014 +- 0.006 for the gapped random bond, and 1.10 for SU(3) at k = 5 (true 4.19). The single operator h_q (with the identity) already carries essentially all of this gain; the two-site one-site-products and higher momenta add almost nothing.

3. [numerical] Even with h_q, a per-site offset e_0 - ell_inf of about 5e-4 to 8e-4 (critical d = 2 chains), 1.5e-4 (gapped random d = 2) and 0.037 (SU(3), k = 5) survives at every fixed operator support. Hence A_cert(L) still grows like (6/pi)(e_0 - ell_inf) L^2. At fixed support, no bootstrap of this kind certifies a finite upper bound on A.

4. [heuristic] The offset cannot vanish at fixed support for a structural reason (periodization argument below): any feasible point of the infinite-line program with decaying correlations periodizes to an almost feasible point of every ring program, so ell(L)/L is pinned near the infinite-line relaxation value, which lies strictly below e_0 when the supports are finite. Certifying A requires the offset to fall faster than 1/L^2, which in practice means supports growing in proportion to L and exponential cost.

5. [numerical, consistent with but not evidence for a bound] Every bootstrap 1/L coefficient found satisfies A_relax <= pi M log_2 d (ratios 0.93 Heisenberg, 0.76 XXZ, 0.22 SU(3), about 0 for the random bonds). This is compatible with C = pi and f(d) = log_2 d, but the bootstraps supply no certificate of it, and they say nothing about the constant beyond what the exact data already say.

## Status

Numerical, with one heuristic structural argument. The exact ring ground state was checked to satisfy every constraint (to 1e-13) in eight representative cases covering the s = 1, s = 2 and h-only program types; the d = 3 variant with the one-site-plus-h operator set was not checked this way (the check ran out of memory). The SDP values are interior-point primal values with a precision of about 1e-6 per site (Clarabel, reduced-accuracy termination in most runs); they are not rationally certified. What would upgrade it: a dual certificate extraction for the s = 2 program at a few L, and a proof (or disproof) of the periodization bound in point 4.

## How it was established and checked

- Author: Claude (Opus 5.5) as a numerics agent in the c-d campaign, for Jonah Kudler-Flam. No Codex delegation was used; the code is self-written.
- All SDPs, ED and analysis ran as SLURM jobs on Typhon, with cvxpy 1.9.2 and the Clarabel and SCS solvers that were already installed in the user environment. Nothing was installed.
- Checks performed. (i) Feasibility of the exact ground state: `ed.py --mode feasible` assigns the ED ground-space expectation values to every string variable and evaluates every equality and every PSD block; worst equality residual 4e-13 and smallest PSD eigenvalue -4e-15 over eight cases (Heisenberg L = 8, 10 at k = 5, s = 1 and L = 12 at k = 6, s = 2 with the h-only operator set; XXZ L = 10 at k = 7, s = 2; random d = 2 at L = 10; random d = 3 at L = 7; SU(3) at L = 6, 9). The energy computed from the string variables reproduces E_0(L)/L to 1e-15. (ii) Monotonicity in the constraint set holds within 4e-6 per site (for example k = 4, 6, 7 at s = 2). (iii) At small L where the program covers the whole ring the bound reproduces E_0(L) (L = 6, 8 for d = 2 and L = 6 for SU(3)). (iv) ED reproduces the known A: 1.2993 for XXZ against 1.2990 and 4.214 for SU(3) against 4.189 from fits; Heisenberg gives 1.58 at L = 24 because of its logarithmic corrections. (v) Early SCS runs at s = 2 disagreed with each other by up to 1e-4 per site; all d = 2, s = 2 numbers quoted here come from the regularized Clarabel reruns, which agree with monotonicity and with ED to 1e-6. The d = 3, s = 2 numbers are SCS values and carry errors of order 1e-4 per site, which does not affect the conclusions drawn from them (offsets of 0.04 to 0.07).
- Not checked: no independent re-implementation, no human review.

## Method and results

### The program

The variables are the expectation values x_c of real operator strings P_{mu_1} (x) ... (x) P_{mu_n}, one variable per class under ring translation (sites mod L), with P_mu a real orthogonal basis of d x d matrices (identity, diagonal Gell-Mann, symmetric and antisymmetric off-diagonal). Strings with an odd number of antisymmetric factors are dropped, which is valid because for a real bond the averaged ground-space projector is a real, translation-invariant state that satisfies every ground-state constraint. In this parametrization translation invariance, marginal compatibility and the identification of separation r with L - r are automatic, and no consistency equalities are needed. Constraints, each satisfied by the exact ground space:

- rho_k >= 0, the contiguous k-site marginal (for k > L/2 it wraps around the ring);
- G_loc = <O_a^T [H, O_b]> >= 0 and Hermitian for O on m_loc consecutive sites (local stationarity);
- P_r >= 0 for the two-block marginal on {0..s-1} u {r-1..r+s} at every separation r whose sites do not fit in a k-arc (`--pair-s` sets the block size independently of s);
- Gamma(q) = <O_q^dag O_q>/L >= 0 and G(q) = <O_q^dag [H, O_q]>/L >= 0 at q = 2 pi m/L, m = 0..L/2, for O in an operator set on s sites; the Hermiticity of G(q) at every q, equivalently g_r = g_{-r}^T for the real-space blocks, is imposed as an equality (bilocal stationarity).

Operator sets: `full` (all d^{2s} matrix units on s sites), `onesite_h` (one-site matrix units plus the bond h), `h` (identity and h only). The objective is <h_01>. Linear dependencies among the equalities are removed by pivoted QR, and Clarabel is run with static regularization 1e-7 (needed because the s = 2 programs have nearly flat directions); SCS is a fallback.

### ED baseline

| bond | d | M | e_0 used | A (exact) | A from ED |
| --- | ---: | ---: | ---: | ---: | ---: |
| Heisenberg | 2 | 0.5000 | -0.4431472 (exact) | 1.5708 | 1.58 (L = 24, log corrections) |
| XXZ, Delta = 1/2 | 2 | 0.5000 | -0.3750000 (fit) | 1.2990 | 1.2993 |
| SU(3) swap | 3 | 1.0000 | -0.703196 (fit, L <= 15) | 4.1888 | 4.214 |
| random d = 2 | 2 | 0.5916 | -0.769658 (L = 22) | 0 (gapped, two-fold quasi-degenerate) | 0 |
| random d = 3 | 3 | 2.1047 | -1.594502 (fit) | unknown | 0.04 (rough) |

### Which ingredients matter (L = 16, gap measured in units of the Casimir energy pi A/(6L))

| program | Heisenberg ell/L | gap/Casimir | XXZ ell/L | gap/Casimir |
| --- | ---: | ---: | ---: | ---: |
| exact E_0(L)/L | -0.4463935 | 0 | -0.3776727 | 0 |
| local only, k = 4 | -0.5000000 | 16.7 | -0.4215352 | 16.5 |
| local only, k = 5 or 6 (equals E_0(6)/6) | -0.4671293 | 6.45 | -0.3946723 | 6.40 |
| s = 1, k = 4 | -0.4652592 | 5.87 | -0.3934166 | 5.93 |
| s = 1, k = 6 | -0.4574747 | 3.45 | -0.3872248 | 3.60 |
| s = 1, k = 6, no pair PSD | -0.4631994 | 5.23 | | |
| s = 2, k = 4 | -0.4467540 | 0.112 | -0.3779506 | 0.105 |
| s = 2, k = 6 | -0.4467577 | 0.113 | -0.3779548 | 0.106 |
| s = 2, k = 6, momenta m = 0, 1, 2 only | -0.4467576 | 0.113 | -0.3779565 | 0.107 |
| s = 2, k = 6, operator set {1, h} | -0.4467853 | 0.122 | -0.3779808 | 0.116 |
| s = 2, k = 6, one-site units and h | -0.4467853 | 0.122 | -0.3779808 | 0.116 |
| s = 2, k = 6, no pair PSD | -0.4472066 | 0.253 | -0.3783895 | 0.270 |
| s = 2, k = 7 | -0.4467000 | 0.095 | -0.3779242 | 0.095 |

Reading: the local programs see a shorter ring, as C117 proved. Ring-wide two-block constraints with one-site modes help only modestly. The energy-density mode h_q at the lowest momenta, with positivity of the two-block marginals, is what makes the bound track the finite-size energy. The size of the contiguous marginal is nearly irrelevant once h_q is present.

### L dependence

s = 1 (Clarabel and SCS, precision about 2e-5 per site), two-parameter fits ell(L)/L = ell_inf - (pi A_relax/6)/L^2 over L in [20, 48]:

| bond | k | ell_inf | e_0 - ell_inf | A_relax | A_cert(48) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Heisenberg | 4 | -0.465458 | 0.0223 | -0.05 +- 0.01 | 98 |
| Heisenberg | 5 | -0.462632 | 0.0195 | -0.15 +- 0.04 | 85 |
| Heisenberg | 6 | -0.457889 | 0.0147 | -0.4 +- 0.1 | 64 |
| XXZ | 4 | -0.393526 | 0.0185 | -0.01 +- 0.01 | 82 |
| XXZ | 6 | -0.387404 | 0.0124 | -0.15 +- 0.02 | 42 (L = 42) |
| random d = 2 | 4 | -0.771187 | 0.0015 | 0.000 +- 0.001 | 6.7 |
| SU(3) | 4 | -0.77155 (L <= 27) | 0.068 | consistent with 0 | 95 (L = 27) |

The small negative A_relax values reflect the bound approaching its plateau from above at small L (the ring constraints bite harder on short rings), not a Casimir term.

s = 2 (regularized Clarabel, precision about 1e-6 per site), per L:

| L | Heis. k = 7 ell/L | E_0/L | gap/Casimir | A_cert | XXZ k = 7 ell/L | E_0/L | gap/Casimir | A_cert |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 10 | -0.4515436 | -0.4515446 | 0.000 | 1.60 | -0.3819025 | -0.3819033 | 0.000 | 1.32 |
| 12 | -0.4490531 | -0.4489492 | 0.018 | 1.62 | -0.3798524 | -0.3797727 | 0.017 | 1.33 |
| 14 | -0.4476206 | -0.4473964 | 0.053 | 1.67 | -0.3786779 | -0.3784970 | 0.052 | 1.38 |
| 16 | -0.4467000 | -0.4463935 | 0.095 | 1.74 | -0.3779242 | -0.3776727 | 0.095 | 1.43 |
| 18 | -0.4460740 | -0.4457083 | 0.144 | 1.81 | -0.3774134 | -0.3771092 | 0.145 | 1.49 |
| 20 | -0.4456388 | -0.4452193 | 0.204 | 1.90 | | | | |

With the operator set {1, h} at k = 6 the series continues to L = 24: gap/Casimir = 0.43 (Heisenberg) and 0.41 (XXZ), A_cert = 2.25 and 1.83. (Bound values a few times 1e-7 above E_0 at L = 10 are within solver precision.) Two-parameter fits:

| bond | program | L range | ell_inf | e_0 - ell_inf | A_relax | A_relax/A | A_relax/(pi M log_2 d) |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| Heisenberg | s = 2, k = 7 | 14-20 | -0.443726 | 5.8e-4 | 1.457 +- 0.005 | 0.93 | 0.93 |
| Heisenberg | s = 2, k = 6, {1, h} | 16-24 | -0.443918 | 7.7e-4 | 1.401 +- 0.004 | 0.89 | 0.89 |
| XXZ | s = 2, k = 7 | 14-18 | -0.375475 | 4.8e-4 | 1.199 +- 0.002 | 0.92 | 0.76 |
| XXZ | s = 2, k = 6, {1, h} | 16-24 | -0.375609 | 6.1e-4 | 1.159 +- 0.003 | 0.89 | 0.74 |
| random d = 2 | s = 2, k = 4 | 14-24 | -0.769809 | 1.5e-4 | 0.014 +- 0.006 | | about 0 |
| SU(3) | s = 2, k = 5, one-site units and h, pairs of size 1 (SCS) | 15-21 | -0.740200 | 0.037 | 1.10 +- 0.03 | 0.26 | 0.22 |

The quoted errors are fit errors only; corrections of order 1/L^4 and the Heisenberg logarithms are not included, so the A_relax values are good to perhaps 5 percent. The offset shrinks from k = 6 to k = 7 by a factor of about 0.75 to 0.78, roughly as 1/k^2 (heuristic, two points).

For d = 3 the k = 4 marginal is the bottleneck. At k = 4 the energy-density mode adds nothing for SU(3) (bounds identical to s = 1 to 1e-7), and the random d = 3 bond gains under 1e-4 per site. At k = 5 (SU(3), SCS) the bound per site rises with L (-0.74558, -0.74458, -0.74274, -0.74200, -0.74148 at L = 9 to 21), so a Casimir-like term with the right sign appears, but the offset of 0.037 per site dominates.

### Constraints that grow with L (ring wraparound, k > L/2, s = 1)

| bond | (L, k) | gap/Casimir | A_cert |
| --- | --- | ---: | ---: |
| Heisenberg | (8, 5), (10, 6), (12, 7), (14, 8) | 0.24, 0.26, 0.62, 0.53 | 1.99, 2.02, 2.58, 2.42 |
| Heisenberg | (10, 7), (12, 8) | 0.01, 0.17 | 1.61, 1.86 |
| XXZ | (8, 5), (10, 6), (12, 7), (14, 8) | 0.31, 0.48, 0.74, 0.73 | 1.73, 1.94, 2.28, 2.26 |
| SU(3) | (9, 5), no pair PSD, SCS | 0.57 | 6.70 |

A contiguous marginal of more than half the ring, with the ring identification built in, leaves a gap between a quarter and three quarters of the Casimir energy for k = L/2 + 1 (smaller for larger k), but costs d^{L/2} and is still worse than the fixed-support s = 2 program at the same L.

### Periodization argument (heuristic)

Let Pi_inf be the infinite-line program with the same ingredients (string variables on Z, the same PSD blocks at every separation, and Gamma(q), G(q) >= 0 for all real q), and let x be a feasible point of it whose connected two-block correlations C(r) decay. Define ring data by keeping every string supported in a window of length o(L) and replacing each two-block connected correlator by its periodization C_L(r) = sum_n C(r + nL). By Poisson summation the discrete Fourier transform of C_L at q = 2 pi m/L equals the continuous transform of C at the same q, so Gamma(q) and G(q) are unchanged at every ring momentum, and the linear equalities are preserved. The PSD blocks change by the tail tau(L) = sum_{|r| >= L/2} |C(r)|, and if Pi_inf has a strictly feasible point the repair costs O(tau(L)) in energy. Hence ell(L)/L <= ell_inf^{line} + O(tau(L)), with ell_inf^{line} <= e_0, and strictly below e_0 when the supports are finite (as the data show). The certified A_cert(L) then grows like (6/pi)(e_0 - ell_inf^{line}) L^2 - O(L^2 tau(L)). The data are consistent with this: at s = 1 the plateau is reached by L = 20 (fast decay of the optimal fake correlations), while at s = 2 the approach to the plateau is power-like, ell(L)/L = ell_inf - (pi A_relax/6)/L^2, as if the optimal fake state were itself critical with a slightly smaller Casimir coefficient. A proof would need a strictly feasible point of Pi_inf and control of the tail of the optimal solution, neither of which is established here.

### What this suggests for C and f

The bootstraps do not certify an upper bound on A at any fixed operator support, and the ingredients that would (supports growing with L) are exponentially expensive and give a certificate whose constant is set by how fast the offset e_0 - ell_inf falls with the support, not by C M f(d). The finite-L numbers are nonetheless tight in a useful sense: with h_q in the program, A_cert(L) at L = 12 is within 3 percent of the exact finite-size quantity (6L^2/pi)(e_0 - E_0(L)/L) for Heisenberg and XXZ, so a lattice proof that tracks the ring's energy-density mode at q = 2 pi/L would be working with the right object. All the bootstraps' own 1/L coefficients stay below pi M log_2 d, with the largest ratio (0.93) at Heisenberg, which is where the conjectured bound is saturated by the exact A.

## Reproducing

All scripts refuse to run outside a SLURM allocation. Dependencies beyond numpy and scipy: cvxpy (1.9.2 used) with the Clarabel and SCS solvers.

- `models.py`: the bonds (the C117 bonds and random seeds), known A values, and `gauge_M` for M.
- `ring_boot.py`: builds and solves one program. Example: `python3 ring_boot.py --model Heisenberg --L 16 --k 7 --s 2 --mloc 3 --out out/x.json`. Options: `--opset {full,onesite_h,h}`, `--pair-s`, `--no-pair-psd`, `--no-gram`, `--no-G` (with both, a purely local program), `--momenta all|m1,m2,...`, `--solver`. Runtimes on one Typhon core pair: seconds for s = 1 at L <= 48 with Clarabel; 5 to 40 minutes for s = 2, k = 6 or 7 at L <= 24; one to two hours for d = 3, k = 5 with SCS.
- `ed.py`: `--mode energies` (Lanczos ring energies) and `--mode feasible` (the exact-ground-state feasibility check of any program variant).
- `sweep.py` and `run_typhon.slurm`: run a task list in parallel in one node allocation, `sbatch --export=ALL,NPAR=24 --mem=180G run_typhon.slurm tasks.txt`. Memory is the binding resource for s = 2 (several GB per task).
- `tasks_main.txt`, `tasks_v3.txt`, `tasks_v4.txt`, `tasks_rd2.txt`: the task lists used, with outputs in `out/main`, `v3/out` and `v4/out` of the work directory. The first campaign ran an earlier revision of `ring_boot.py` with the same formulation but without the string pruning and the regularized Clarabel step; rerunning sample cases with the current file reproduces its s = 1 values to 1e-7. Its s = 2 runs used SCS and are superseded by `tasks_v4.txt`. Part of the first campaign was killed by the node memory limit (several d = 3 and k = 6 runs at large L are missing), and in `tasks_v4.txt` the s = 2, k = 6 full-operator runs at L >= 18 and the k = 7 runs at L = 22 were cancelled after more than two hours each; the tables above use only completed runs.
- `collect.sh` assembles the final result set, and `analyze.py` produces `summary.json` (every completed run, the ED data, M, the feasibility checks and the fits) and `analysis_final.txt` (the same as text tables), which are included here.

## Relation to the rest of the repository

This answers the question left open by C117 (REPORT.md section 9.79): fixed-support marginal bootstraps cannot see the Casimir term, and the natural ring-aware extension with low-momentum modes was the untested case. The result is that the energy-density mode h_q at q = 2 pi/L is the ingredient that lets a bootstrap track the finite-size energy, but no fixed-support bootstrap certifies A, so the bootstrap route reduces to the same missing crossover-free estimate as C114 and C116 (an offset that must vanish faster than 1/L^2). It does not change "Where things stand" in the top-level README.
