# C134: The half-filled XXZ ring, a rigorous 1/L bound on the doubling defect

Your folder is `../work_C134/`. Paths below are repository-relative, so prefix them with `../../` from your folder. Read `specs/_Cd_target.md`, `specs/C129_d2_u1_bethe.md` and all of `codex/work_C129/` (RESULT.md, PROOF.md, CHECKS.md and the solver bethe_scan.py), `contrib/rigorous_toolbox_map/TOOLBOX.md` sections 4, 7 and 8, and `codex/work_C125/RESULT.md` (positive transfer, D(L) nonincreasing). Rules of `_bound_background.md` apply (tag claims proved, numerical, heuristic or assumed; no probabilities in files; Python only in your folder; browsing tool only for literature; plain prose without em-dashes; PROGRESS.md, RESULT.md, PROOF.md). Your sandbox cannot reach Typhon; write any large computation as a SLURM-ready script and say so.

## Setting

C129 reduced C(2) on the U(1) spin-1/2 class (XXZ plus DM_z plus field) to twist-optimized, field-free sector defects of XXZ, H = J sum (Sx Sx + Sy Sy + Delta Sz Sz), M = J/2 for |Delta| <= 1. It proved the free point, the region J_z <= -J, dilute sectors (n^3 <= L) and a band near the ferromagnetic point. The single missing piece, identified there, is half filling. Kozlowski arXiv:1508.05741 gives the 1/L term below half filling (Corollary 4.3), but at rho = 1/2 only convergence without a rate (Theorem 3.5), and the paper itself flags missing uniform energy control near half filling after Lemma 4.1.

Target. With F_L the twist-optimized ground energy of the half-filled sector (n = L/2 on even L, and the two neighbouring sectors on odd L, with the minimizing flux of C129),

    L [F_{2L} - 2 F_L] <= C J   for all L >= 2 and all Delta in (-1, 1],

with explicit C. A version uniform on [-1 + delta, 1 - delta] with explicit dependence on delta is already worth having, since C129 controls the band near -1 and Delta = 1 can be treated separately.

## Routes

Note first that only one inequality is needed. An upper bound on F_{2L} can be variational. The hard direction is a lower bound on F_L of Casimir form, F_L >= e_0 L - K J/L, which by the telescoping identity of `_Cd_target.md` is equivalent to the target up to constants.

(a) Exact finite-size integral equations. The Destri–de Vega and Klümper–Batchelor–Pearce nonlinear integral equations represent the half-filled ground state energy at finite L exactly, and the classical derivations of the 1/L term (de Vega–Woynarovich 1985, Hamer–Quispel–Batchelor 1987, Klümper–Batchelor–Pearce 1991) go through them. Decide whether any of these derivations has been made rigorous (search the literature, including Kozlowski's later work and work on the six-vertex model), and if not, which step fails and whether at half filling (no Fermi boundary, roots on the whole line) the analysis is simpler than in Kozlowski's setting. Try to prove the needed one-sided estimate.
(b) Finish Kozlowski at half filling. Read sections 3 and 4 of 1508.05741 and the remark after Lemma 4.1. Identify what replaces the Fermi-boundary analysis when the root density has full support, and whether a rate for Theorem 3.5 follows with constants uniform on compact Delta intervals.
(c) Bethe-free routes. The six-vertex transfer matrix commutes with H and is positive for the relevant Delta range (Yang–Yang, C125). Duminil-Copin, Kozlowski, Krachun, Manolescu and Tikhonovskaia arXiv:2012.11675 give sector free energies to O(1/N) per site. Decide whether their method, or a monotonicity of D(L) = F_{2L} - 2F_L in L (C125) combined with a bound at L of order the crossover, gives the 1/L rate.
(d) Delta = 1. Only an upper bound on the defect is needed, and the data decrease (4, 2.79, 2.62, 2.56 at L = 2, 4, 6, 8). Decide whether Delta = 1 follows from continuity of the bound from Delta < 1 together with monotonicity in Delta of a suitable quantity, or needs its own argument.

## Tasks

1. Literature: state precisely what is rigorous about the finite-size ground energy of the half-filled XXZ ring for -1 < Delta <= 1, with citations checked against the papers.
2. Prove the target on the largest Delta range you can, with explicit constants, and state exactly what remains.
3. Numerics: with the C129 Bethe solver, compute L [F_{2L} - 2F_L]/M at half filling for L up to several thousand on a Delta grid including -0.999, -0.99, -0.9, 0, 0.5, 0.9, 0.99, 1, and check any constant you prove. Cross-check with exact diagonalization at L <= 12.
