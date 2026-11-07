# C135: the LTI rate, a telescoping-certificate route to (K) at d = 2

Your folder is `../work_C135/`. Paths below are repository-relative, so prefix them with `../../` from your folder.

Read first: `specs/_Cd_target.md` (lessons 1 to 5), `specs/_bound_background.md` with all updates, `contrib/cd_free_route/RESULT.md` and `PROOF.md` (sections 1 and 3), `codex/work_C132/RESULT.md` and its `PROOF.md` (sections 1, 4, 7), `contrib/bootstrap_ring/README.md` (the fixed-support offset and its periodization heuristic), and `contrib/marginal_toolbox_map/MAP.md` sections 0, 1 and 7 together with `lti_probe.py` and `lti_probe_output.txt` in that folder. Literature to open with the browsing tool: Kull, Schuch, Dive, Navascués arXiv:2212.03014 (Section II F, Fig. 2, Section V and Fig. 5, the measured n^{-2} decay of the LTI error), Fawzi, Fawzi, Scalet arXiv:2305.06855 Appendix B (the proved O(1/l) rate, by tiling), and Blakaj, Wolf arXiv:2308.14585 (Omega_inf is not semialgebraic).

Rules of `_bound_background.md` apply. Tag every claim proved, numerical, heuristic or assumed. No probabilities, odds or percent-confidence in any file. Python (numpy, scipy, and cvxpy with Clarabel or SCS if already present; do not install anything) only inside your folder. No network commands except the browsing tool for literature. Plain prose without em-dashes. Keep PROGRESS.md and write RESULT.md and PROOF.md. Your sandbox cannot reach the cluster; write any computation that does not fit as a SLURM-ready script and say so.

## Setting

Local dimension d = 2. A bond h is a Hermitian operator on C^2 (x) C^2, M(h) its gauge-optimized half-spread (minimum over real a and Hermitian o of the half-spread of h + a 1 + o (x) 1 - 1 (x) o), E_L its L-ring ground energy and e_0 its energy density on Z. Omega_L and Omega_inf are the sets of nearest-neighbour marginals of translation-invariant states on the L ring and on Z. The target is

    (K)  every element of Omega_L is within trace distance K/L^2 of Omega_inf, for all L >= 2,

equivalently E_L >= e_0 L - K M(h)/L for every bond (cd_free_route PROOF section 1; C132 showed the dual-M norm equals the trace norm on the relevant slice). Only a 1/L^2 rate is of interest; O(1/L) is trivial.

Define LTI_n as the set of states rho_n on n sites with Tr_1 rho_n = Tr_n rho_n, and

    E_LTI(n) = min { Tr[(h (x) 1) rho_n] : rho_n in LTI_n }.

Facts recorded in MAP.md section 7, to be re-verified before use.

(F1) [Proved, elementary] For 2 <= n <= L the two-site marginals of LTI_n contain Omega_L, and they contain Omega_inf for every n. Hence E_LTI(n) <= min(e_0, min over L >= n of E_L/L).

(F2) [Proved, elementary] By SDP duality, E_LTI(n) = max { t : h (x) 1_{n-2} + X (x) 1 - 1 (x) X >= t 1 } over Hermitian X on n - 1 sites. Summing translates of this inequality around a ring of length L >= n telescopes X away and gives E_L >= t L.

(F3) [Proved, elementary] The uniform LTI rate

    (L-LTI)  e_0(h) - E_LTI(n) <= K M(h)/n^2   for every bond h and every n >= 2

implies (K) with the same K.

(F4) [Proved, elementary, check] The periodization construction of cd_free_route PROOF section 3, with any CPTP R satisfying R(rho_{l-1}) = rho_l (C132), uses only Tr_1 rho_l = Tr_l rho_l, so it applies to any LTI state and gives e_0 <= Tr[h sigma_2*(rho_n, R)].

(F5) [Numerical, probe, n <= 6] E_LTI(n) equals min over L >= n of E_L/L to solver accuracy for the Heisenberg, XX and critical TFI bonds. For a chiral bond (the Heisenberg bond twisted by pi/2 per site, minus 1/4) it equals a twisted-ring value instead, and for random bonds it is at most 3e-2 below the best untwisted ring at n >= 3, with the gap shrinking in n. Kull et al. measure e_0 - E_LTI(n) proportional to n^{-2} in six models (gapped ones faster), with different even and odd slopes.

The working hypothesis behind this run is "ring saturation": the LTI relaxation at level n is optimized by a (possibly twisted) ring state of length at least n, or by Omega_inf itself. If that is right, (L-LTI) and the twisted-ring form of (K) are the same statement, and the dual certificate X of (F2) is a concrete object to compute and to construct.

## Tasks

1. Verify (F1) to (F4), with complete proofs in PROOF.md. For (F4), check the support condition of the channel when rho_{l-1} is not full rank. State precisely what "twisted ring" means for a bond with an on-site symmetry g (that is, [h, g (x) g] = 0) and prove or refute that every n-site window marginal of the twisted ring ground state (or ground-space mixture) of length L >= n lies in LTI_n.

2. Ring saturation. Prove or refute, at d = 2: for every real bond with h = F h F (F the swap), E_LTI(n) = min(e_0, min over L >= n of E_L/L). Start from small n (n = 3, 4) where the SDP and its dual can be written in closed form for the XXZ family with fields. Then state the right general conjecture for bonds without reflection symmetry (allowing twisted rings), and test it numerically. A counterexample is as valuable as a proof, because it measures how much stronger (L-LTI) is than (K).

3. Dual certificates. Extract the optimal X of (F2) numerically for XX, Heisenberg and critical TFI at n = 4 to 8 and describe its structure. Is it, up to small corrections, a weighted sum of translated bond terms (an Anderson bound with tapered weights, which is sine-square-like and therefore exposed to lesson 2), or does it carry genuinely multi-site terms? For XX, try to construct X analytically from the free-fermion solution so that h + X (x) 1 - 1 (x) X >= (e_0 - c/n^2) 1 with explicit c. Such an X would be the first telescoping certificate with a Casimir-rate deficit, and its structure is what a general construction would have to reproduce.

4. Energy-weighted periodization. With rho_n the LTI optimizer and X its dual, complementary slackness gives (h (x) 1 + X (x) 1 - 1 (x) X - E_LTI(n)) rho_n = 0. Decide whether this identity, together with (F4) for the Petz map or the optimal channel, bounds Tr[h sigma_2*] - E_LTI(n) by K M/n^2 at d = 2. If an argument closes for a subclass (quasi-free, frustration-free, U(1)-symmetric), say exactly which. If it fails, give the precise obstruction, keeping in mind the C132 audit of the cd_free_route trace-norm obstruction.

5. Falsification computation (numerical, must be run). This is the deliverable even if Tasks 2 to 4 make no progress.

    (a) Directional gaps. By convex duality (cd_free_route PROOF Proposition 1), the worst trace distance of an element of Omega_L from Omega_inf equals sup over bonds of (e_0(h) - E_L(h)/L)/M(h). The exposed points of Omega_L are the ring ground-state marginals rho_L(h), so test these. Upper bounds on e_0 come from Petz periodization of a ring ground state (cd_free_route `run_checks.py`, rigorous, e_0 <= Tr[h sigma_2*]) at the largest feasible ring, or from a hand-written uniform MPS. Lower bounds on e_0 come from E_LTI(n). For each test bond and each L from 3 to 20 report the bracket

        L^2 (E_LTI(n) - E_L/L)/M  <=  Q_L(h)  <=  L^2 (e_0^up - E_L/L)/M,

    where the left side needs n > L and is therefore available only up to the largest solvable n (expect n = 8, perhaps 9 with a real formulation for real bonds and U(1) block reduction where it applies). A growing lower bound would refute (K) at that bond. A bounded upper bound is consistent with (K).

    (b) Trace distances of exposed points. For L = 3 to 8, also compute the actual trace distance in two ways. A lower bound is d_low(L, n) = min ||rho_L - Tr_{3..n} sigma_n||_1 over sigma_n in LTI_n with n > L. An upper bound is d_up(L) = min ||rho_L - tau||_1 over tau in the convex hull of a dictionary of points of Omega_inf: Petz-periodized marginals sigma_2* of rings of all lengths up to 20 for this bond and for the other test bonds, plus uniform-MPS marginals if you have them (an SDP in the hull weights). Report L^2 d_low and L^2 d_up, and check whether L^2 d_up stays bounded.

    (c) LTI rate. For n = 2 to 8 report n^2 (e_0^up - E_LTI(n))/M (an upper bound on the LTI constant) and n^2 (E_LTI(n + 2) - E_LTI(n))/M (a certified lower bound on it). A growing lower bound would refute (L-LTI) without touching (K).

    Test bonds: Heisenberg; XXZ at Delta = 0.5, 2, and Delta -> -1+ (-0.9, -0.99, -0.9999, the d = 2 long-crossover face of lesson 3); XX; TFIM at g = 1 and g = 1 + c/n for c in {-2, 0.5, 5} (the near-critical massive mode of lesson 2); XXZ with a DM_z term (twisted rings); the chiral bond of MAP section 7; and at least 30 random bonds radially projected to the M = 1 gauge-quotient sphere, both real and complex, including non-swap-symmetric ones. Then run a multistart local maximization over the M = 1 sphere (eight angles, as in `specs/C131_d2_cover.md` and `codex/work_C131/optimize_q.py`) of the lower bound in (c) at the largest feasible n, and of the bracket in (a) at L = 6, 8, 10. Report the maximizers and whether they drift toward a stratum as n or L grows.

    Numerical hygiene. Record solver status and residuals; Clarabel failed on the complex n = 7 LTI SDP in the probe, so use a real embedding, symmetry reduction (real bonds allow real rho_n; U(1) bonds allow magnetization blocks; reflection-symmetric bonds allow F-symmetric rho_n, but only after proving that this symmetrization preserves LTI), or SCS with tightened tolerances. Treat any quantity that relies on a non-converged solve as exploratory and say so.

6. State exactly what is proved, for which class and at d = 2, and the single remaining lemma if (L-LTI) is not closed.

## Falsification criteria, summarized

- A d = 2 bond with L^2 (E_LTI(n) - E_L/L)/M growing in L (n > L) refutes (K). Name the lesson it hits.
- A d = 2 bond with n^2 (E_LTI(n + 2) - E_LTI(n))/M growing in n refutes (L-LTI). It does not refute (K), and the run should then switch to measuring the gap between E_LTI(n) and the best (twisted) ring.
- A counterexample to ring saturation for a reflection-symmetric real bond means the LTI route is strictly stronger than (K) on that bond, and its size should be reported.

## What success means

In decreasing order of value. (i) A proof of (L-LTI), hence of (K), for all d = 2 bonds, or for a stratum of C131 not already covered by C129, C130 or C132. (ii) A proof of ring saturation for reflection-symmetric real d = 2 bonds, which turns Kull et al.'s n^{-2} data and any future LTI certificates into statements about (K) itself. (iii) An explicit telescoping certificate with deficit c M/n^2 for XX or for all quasi-free d = 2 chains, which would be the first rigorous Casimir-rate SOS certificate and a template for interacting chains. (iv) A clean numerical table from Task 5, with brackets on Q_L and on the LTI constant for every test bond, and either a bounded worst case or an explicit growing family. A run that delivers only (iv) is still useful, because the LTI lower bounds are certified (up to solver residuals), unlike the ring extrapolations used so far.
