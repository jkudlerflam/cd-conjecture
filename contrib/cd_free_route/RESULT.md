# cd_free_route result

2026-11-06, Opus 5.5 agent. Derivations are in [PROOF.md](PROOF.md), and every numerical claim is reproduced by [run_checks.py](run_checks.py) (`python3 run_checks.py` in about 15 seconds, `--long` in about one minute; logs of the slower L = 8 runs are in `logs/`). No finite C(d) is proved. The route is reduced to one lemma, stated below, together with a proved obstruction to the obvious way of proving it.

## Route and why it avoids the lessons

The route is variant (a), E_L >= e_0 L - K(d) M/L, through periodization of ring states. It is not a permutation twist (C126), it uses no crossed-channel spectrum (C127, C128, cd_crossed_fock), and it is not the sine-square deformation.

- Lesson 1. It uses the full microscopic ring ground state, not crossed spectral data.
- Lesson 2. No deformed Hamiltonian appears. The trial state is built from the uniform ring state, and massive modes come out exponentially suppressed (TFIM, below).
- Lesson 3. No counting is involved. A decoupled slow copy enters weighted by its own bond scale (proved), and the long-crossover C106 chain stays bounded uniformly in its small parameter (numerical).
- Lesson 4. The bound is a variational upper bound on e_0, with no comparison of Gibbs operators and no seam pinching.
- Lesson 5. No convexity is assumed.

## Claims

[Proved] Kinematic reformulation (PROOF section 1). Variant (a) with constant K for all bonds at fixed d is equivalent to a statement about marginals. Every nearest-neighbour marginal of a translation-invariant L-ring state must lie within dual-M distance K/L^2 of the nearest-neighbour marginal of a translation-invariant state on Z. Classical marginals satisfy this with K = 0, by Markov extension.

[Proved] Quasi-free chains (section 2). Momentum box smoothing gives E_L >= e_0 L - (pi^2/(6L))|e_hop| for every quadratic chain, so A <= pi|e_hop|. For XX this reads E_L >= -csc(pi/L). It holds with equality on every even L (ED, L = 2 to 12) and gives the exact A = 1.

[Proved] Universal periodization channel (section 3). Iterating the Petz map of the ring marginal rho_l gives a translation-invariant state on Z. All of its nearest-neighbour marginals equal the (A, B_1) marginal sigma_2* of (id_A (x) R_{B->BC})(rho_l). Hence e_0 <= Tr[h sigma_2*] for every bond, and variant (a) follows from the key inequality delta_L = Tr[h sigma_2*] - E_L/L <= K M/L^2.

[Proved] Reduction (section 4). delta_L <= M ||(id (x) R)(rho_AB) - rho_ABC||_1, where the gauge cancels exactly.

[Proved obstruction] Section 5. Data processing plus Alicki–Fannes–Winter gives I(A:C|B) <= 2 eps log d + O(eps log 1/eps) for the trace-norm recovery error eps. In CFT that conditional mutual information is (c/3) pi^2/L^2, which carries c rather than A, and c is unbounded at fixed d (conditional constructions in STATUS.md). The trace-norm step therefore cannot give K(d), and any proof must use h itself. [Proved] The construction is velocity-weighted on decoupled stacks: delta is additive, and delta(eps h) = eps delta(h).

[Numerical] Values of L^2 delta_L/M are bounded and decrease with L on every test.

- Heisenberg: 4.05 to 3.33 for L = 6 to 20.
- XXZ at Delta = 0.5: 3.32 to 2.77 for L = 6 to 14.
- XXZ at Delta = 2: 3.29 to 2.20 for L = 6 to 14.
- SU(3) swap: 9.00, 7.49, 5.39, 4.93 at L = 3, 6, 9, 12.
- FM Heisenberg with the whole ground multiplet: 1.33 to 1.09.
- TFIM at g = 1.6: 0.31 to 0.02.
- C106 chain and a coupled fast/slow ladder: at most 3.3, uniformly in the small parameter down to 0.03 and 0.01.

[Numerical, heuristic] On critical chains L^2 delta_L tends to pi A/3. For XX, 2.109 at L = 20 against 2.094. For Heisenberg, 3.333 against 3.290. So the Petz state overshoots e_0 by exactly the Casimir deficit, and the route loses at most a factor 2.

[Assumed sufficient, open] Lemma P. For the ring ground state, Tr[h~((id (x) R)(rho_AB) - rho_ABC)] <= K(d) M/L^2 at l = L/2. This implies A <= (6/pi) K(d) M, and at d = 2 it would give the spin-1/2 Casimir statement. In CFT it is asymptotically as strong as the bound itself, up to the factor 2.

## Dead routes

- [Proved] Trace-norm or CMI closure of Lemma P. The obstruction is the one above.
- [Heuristic] Flux averaging for interacting U(1) chains, E_{2L} - 2E_L <= E_L(pi/2) + E_L(-pi/2) - 2E_L. This is exact for free fermions. In CFT it needs Luttinger K >= 1, so it fails at Heisenberg (K = 1/2), and neutral sectors are invisible to it.
- [Heuristic] Copy-Fourier product trial for interacting fermions. It collapses to the Hartree–Fock energy at cost O(M) per site, which is C126 territory.

No counterexample to a finite C(d) was found. No cluster jobs were needed, and all runs took under a minute locally.
