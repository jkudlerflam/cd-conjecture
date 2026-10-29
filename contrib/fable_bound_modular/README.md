# The modular route: A is the ground-state energy defect of the lattice entanglement Hamiltonian, and relative entropy alone cannot bound it

Author: Claude (Anthropic, Fable 5.1), 2026-10-29. Closed book apart from the files named in Section 8. Numerics are numpy and scipy only, in this folder, all run locally (nothing needed the cluster). Nothing was delegated to Codex: every analytic step here is short, and the one hard question (Section 6) is not self-contained enough to pose. Every claim is tagged proved, numerical, heuristic, CFT (standard conformal field theory, the infrared assumption of the background note) or assumed. Natural logarithms unless a base is written. Lattice spacing and hbar are one.

## Claim

Setting as in `specs/_bound_background.md`: a one-site translation-invariant nearest-neighbour chain with bond h on (C^d)^2, infrared a finite set of unitary CFT sectors (c_a, v_a), A = sum_a c_a v_a, B = sum_a c_a / v_a, C = sum_a c_a, gauge-optimized half-spread M, ground-state energy per bond e_0. For an interval I of ell sites with bonds x = 0, ..., ell-2 let

    w_x = (x+1)(ell-1-x)/ell,     K_w = sum_x w_x h_x,     W = sum_x w_x,

the lattice version of the interval entanglement Hamiltonian, and let tau_beta = exp(-beta K_w)/Z_w(beta) be its Gibbs state. Let rho_I be the ground state of the infinite chain reduced to I.

1. [Proved] Exact lattice identity. Tr[rho_I K_w] = e_0 W, so

       S(rho_I || tau_beta) = beta [ e_0 W - F_w(beta) ] - S(rho_I) >= 0,      F_w = -beta^{-1} log Z_w .

   The bond enters only through the linear pairing e_0 W - F_w(beta) = sum_x w_x Tr[ h (rho_2^{(0)} - rho_x^{(beta)}) ], the pairing the background note asks for. As beta -> infinity the bracket tends to the defect

       Delta(ell) = e_0 W - E_0(K_w) = max over states psi of  sum_x w_x ( e_0 - <h_x>_psi ) .

2. [CFT] Each term has a log ell coefficient fixed by the sector data:

       S(rho_I) = (C/3) log ell,
       log Z_w(beta) = -beta e_0 W + sum_a (c_a/6) (n_a + 1/n_a) log ell,     n_a = beta v_a / (2 pi),
       Delta(ell) = (A/(12 pi)) log ell,
       S(rho_I || tau_beta) = [ beta A/(12 pi) + pi B/(3 beta) - C/3 ] log ell = sum_a (c_a/6) (sqrt(n_a) - 1/sqrt(n_a))^2 log ell,

   all up to O(1). So A/(12 pi) is the log ell coefficient of the ground-state energy defect of the lattice entanglement Hamiltonian, equivalently of the large-beta slope of the relative entropy. The velocities enter because K_w = e_0 W + sum_a v_a H_w^a with H_w^a the velocity-free modular generator of sector a, whose ground-state energy -(c_a/(12 pi)) log ell is universal. Pointwise, the ground state of K_w (the lattice Boulware state) has local energy e_0 - (A/(24 pi)) ell^2/(x^2 (ell-x)^2) at bond x, that is e_0 - A/(24 pi x^2) at distance x from an edge.

3. [Numerical] Items 1 and 2 are checked exactly on the XX chain to ell = 2048 (Gaussian states), including a two-velocity chain (two decoupled XX chains with v = 1 and 1/2, d = 4), where the relative-entropy coefficient beta A/(12 pi) + pi B/(3 beta) - C/3 is reproduced to four digits, and by exact diagonalization on Heisenberg, XXZ (Delta = 0.5) and SU(3) swap chains. On every single-velocity chain S(rho_I || tau_beta) at beta = 2 pi / v is 0.013 and flat in ell, as the formula predicts.

4. [Proved] Consequence for the conjecture: given item 2 the statement A <= pi M log_2 d is equivalent to

       lim_ell Delta(ell)/log ell <= (M log_2 d)/12,     equivalently     lim_x x^2 ( e_0 - <h_x>_Boulware ) <= (M log_2 d)/24 .

   The Heisenberg chain has Delta(ell) = (1/24) log ell and M log_2 d /12 = 1/24, so the reformulated bound is saturated there too.

5. [Proved] Every positivity inequality among the natural lattice states (rho_I, its normalized powers rho_I^n / Tr rho_I^n, tau_beta for any beta, Gibbs states of open or periodic uniform chains) evaluates, in the CFT, to a Cauchy-Schwarz type statement: min over beta of S(rho_I || tau_beta) = (1/3)(sqrt(AB) - C) log ell >= 0 is (sum c_a v_a)(sum c_a/v_a) >= (sum c_a)^2, and the reversed and powered versions are of the form sum_a c_a (x_a - y_a)^2 >= 0. None of them compares A with M or d. The reason is structural: entanglement and Renyi entropies of ground states are velocity blind, thermal entropies carry 1/v_a, and the only velocity-positive quantity in the whole family is the K_w energy, which relative entropy can only bound from below.

6. [Proved] Layer-cake inequality. For any unimodal nonnegative weight, K_w = sum_k Delta_k H_{I_k} with nested intervals, so E_0(K_w) >= sum_k Delta_k E_0^open(n_k) and

       Delta(ell) <= sum_k Delta_k [ n_k e_0 - E_0^open(n_k) ] = (max w)(-f_b) + (pi A/24) sum_k Delta_k / n_k + o(1),

   with f_b the open-chain boundary energy in the chosen gauge. For the quadratic weight sum_k Delta_k / n_k -> 1/2, so Delta(ell) <= (ell/4)(-f_b) + pi A/48; comparing with (A/(12 pi)) log ell forces f_b < 0 strictly for every critical chain with A > 0, in every gauge (a small theorem at the physics-standard level, consistent with f_b = 1/2 - 2/pi for XX and the ED values -0.25 for Heisenberg, -0.20 for XXZ at Delta = 0.5, -0.48 for SU(3)). For a tent weight of height R and half-width R the two sides carry the same log R, with coefficients A/(12 pi) (true) against pi A/48 (layer cake), a loss of exactly pi^2/4, but the term (-f_b) R is O(M R) and dominates. Energy superadditivity therefore stops at the same place as every energy route in C98 and `contrib/fable_bound_proof`: a boundary energy of order M per unit of weight.

7. [Numerical, with a recalled exact fact] For the XX chain the top Schmidt vector of rho_I is exactly the ground state of K_w (the weights w_x are, up to normalization, the hopping amplitudes of the Eisler-Peschel matrix that commutes with the restricted correlation matrix), so its K_w energy equals E_0(K_w) to all digits. For Heisenberg, XXZ and SU(3) the top Schmidt vector is within one percent of E_0(K_w) in the defect at ell <= 9. The lattice entanglement Hamiltonian is a good variational object even with interactions, but this does not help bound its ground energy.

No bound of the form A <= C M log d with C independent of a crossover scale is proved here. What is new is the identification of A with a ground-state energy defect of a specific lattice operator (items 2 and 4), the exact relative-entropy identity that makes this defect a relative-entropy slope (item 1), and the demonstration that relative-entropy positivity by itself is velocity-neutral (item 5).

## Status

Item 1 and the inequalities of items 5 and 6 are proved at the level of elementary operator inequalities (positivity of relative entropy, superadditivity of ground-state energies, translation invariance). Item 2 is standard CFT (Renyi entropies of an interval, the Bisognano-Wichmann and Cardy-Tonni form of the modular Hamiltonian, the trace anomaly), not re-derived from the lattice. Item 4 is proved given item 2. Items 3 and 7 are numerics in `xx_modular.py`, `xx_weights.py`, `ed_modular.py` with outputs `*.out`. The conjecture itself is untouched except for the reformulation.

## 0. Conventions and the infrared assumption

[Assumed, same as the background note] E_0(L) = e_0 L - pi A/(6L) + o(1/L) on rings, the ground state is a unique translation-invariant state on the line with local expectations approaching their CFT values, and the lattice energy density decomposes in the infrared as

    h_x - e_0 = sum_a v_a eps_a(x) + (irrelevant operators and total derivatives),

where eps_a is the dimensionless (unit-velocity) energy density of sector a, eps_a = (T_a + Tbar_a)/(2 pi) in the normalization <T T> = c/(2 z^4). Relevant operators cannot appear with nonzero coefficient, since sum_x h_x is critical; derivatives sum to zero in sum_x h_x and contribute at order w' in K_w, which is subleading.

[Definition] rho_I = Tr_{I^c} |0><0| on ell consecutive sites. K_rho = -log rho_I. K_w = sum_x w_x h_x with the quadratic weight above (it vanishes at the virtual bonds x = -1 and x = ell-1). The Boulware state is the ground state of K_w; for x << ell it is the ground state of the half-chain Rindler Hamiltonian sum_x x h_x.

[CFT] Interval modular Hamiltonian: K_rho = 2 pi sum_a int dx w(x) eps_a(x) + c_0 with w(x) = x(ell-x)/ell and c_0 = S(rho_I) fixed by Tr rho_I = 1. In terms of the physical sector energies eps_a^phys = v_a eps_a, K_rho = 2 pi sum_a (1/v_a) int w eps_a^phys. The lattice operator K_w - e_0 W = sum_a v_a int w eps_a has the velocities in the numerator. For a single sector K_rho = (2 pi / v)(K_w - e_0 W) + c_0 up to lattice corrections; for several sectors no multiple of K_w equals K_rho.

## 1. The exact identity

[Proved, Proposition 1] For every x, Tr[rho_I h_x] = <h_x>_0 = e_0 by translation invariance of the line ground state, so Tr[rho_I K_w] = e_0 W. Relative entropy positivity, S(rho||tau) = -S(rho) - Tr[rho log tau] >= 0 with log tau_beta = -beta K_w - log Z_w(beta), gives

    S(rho_I) <= beta e_0 W + log Z_w(beta)      for every beta > 0,

and the identity in Claim 1. On a ring of length L the same holds with e_0 replaced by E_0(L)/L (`ed_modular.py` prints the residual "rho-defect" as zero to machine precision).

[Proved, Proposition 2] F(beta) = beta e_0 W + log Z_w(beta) = log Tr exp(-beta sum_x w_x (h_x - e_0)) is convex in beta, with F(0) = ell log d, F'(beta) = e_0 W - <K_w>_beta increasing from e_0 W - Tr K_w / d^ell to the defect Delta(ell) = e_0 W - E_0(K_w) as beta -> infinity. In particular Delta(ell) = lim F(beta)/beta and S(rho_I || tau_beta) = F(beta) - S(rho_I).

[Proved, Proposition 3 (variational form)] Delta(ell) = max_psi sum_x w_x (e_0 - <h_x>_psi), because E_0(K_w) = min_psi <K_w>_psi. Both the maximizer (the Boulware state) and the comparison (the line ground state, with e_0 - <h_x> = 0 at every bond) are explicit; the defect is the weighted local energy gain the Boulware state achieves by departing from the vacuum near the edges where w -> 0.

## 2. The CFT evaluation

[CFT] Per sector, write H_w^a = int w eps_a (unit velocity). From the Renyi entropies of an interval, Tr rho_I^n = const (ell/epsilon)^{-(c_a/6)(n - 1/n)} sector by sector, and K_rho = 2 pi sum_a H_w^a + c_0 with Tr exp(-K_rho) = 1 and c_0 = S(rho_I) = sum_a (c_a/3) log(ell/epsilon_a) + const, one gets

    log Tr exp(-2 pi n H_w^a) = n c_0^a - (c_a/6)(n - 1/n) log(ell/epsilon_a) = (c_a/6)(n + 1/n) log(ell/epsilon_a) + n s_a .

The n -> infinity limit gives the ground-state energy E_0(H_w^a) = -(c_a/(12 pi)) log(ell/epsilon_a) + const and the lowest eigenvalue of K_rho equal to (C/6) log ell = S/2, the Cardy-Tonni value. Since tau_beta = (x)_a exp(-beta v_a H_w^a)/Z_a is a product over sectors with n_a = beta v_a/(2 pi),

    log Z_w(beta) = -beta e_0 W + sum_a (c_a/6)(n_a + 1/n_a) log ell + O(1),
    E_0(K_w) = e_0 W - sum_a v_a (c_a/(12 pi)) log ell + O(1) = e_0 W - (A/(12 pi)) log ell + O(1).

[CFT, trace anomaly] For a general smooth weight the same quantity is the Casimir energy of a CFT with position-dependent velocity v_a w(x), that is in the metric -w(x)^2 dt^2 + dx^2, and the trace anomaly gives

    e_0 W - E_0(K_w) = (A/(24 pi)) int dx (w')^2 / w + O(1),

which for the quadratic weight is (A/(24 pi))(2 log ell - 4), for the tent min(x, 2R - x) is (A/(12 pi)) log R, and for the quartic x^2 (ell-x)^2/ell^3 is (A/(24 pi))(4/3) with no logarithm. The local form is e_0 - <h_x> = (A/(24 pi))(w'^2 - 2 w w'')/w^2, which is (A/(24 pi)) ell^2 /(x^2 (ell-x)^2) for the quadratic weight and A/(24 pi x^2) for the Rindler weight (the Minkowski vacuum is thermal at local temperature T(x) = 1/(2 pi x) relative to the Boulware state, and each sector's thermal energy density pi c_a T^2/(6 v_a) in units where eps_a^phys = v_a eps_a is c_a v_a /(24 pi x^2) once the physical energy density is restored).

[CFT] Entropies for comparison. S(rho_I) = (C/3) log ell, velocity blind. S(tau_beta) = sum_a (c_a/3)(1/n_a) log ell = (2 pi B/(3 beta)) log ell, inverse velocities (the local temperature 1/(beta w(x)) gives entropy density pi c_a/(3 v_a beta w(x)), integrated against dx). S(rho_I^n / Tr rho_I^n) = (C/(3n)) log ell, velocity blind. Tr[(rho_I^n / Tr rho_I^n) K_w] = e_0 W - (A/(12 pi))(1 - 1/n^2) log ell, velocity positive: the normalized powers of the true reduced state are the only states in the family besides the Boulware state whose K_w energy is below e_0 W by an A-weighted logarithm.

## 3. Numerical checks

[Numerical, exact Gaussian] `xx_modular.py` (output `xx_modular.out`), XX chain h = S^xS^x + S^yS^y, e_0 = -1/pi, v = 1, c = 1, M = 1/2, ell = 64 to 2048.

| quantity | measured log ell slope (ell = 1024 to 2048) | CFT |
|---|---|---|
| Delta(ell) = e_0 W - E_0(K_w) | 0.02587 (rising: 0.02526, 0.02547, 0.02564, 0.02576, 0.02587) | A/(12 pi) = 0.02653 |
| e_0 W - <chi| K_w |chi>, chi the top Schmidt vector | identical to the line above to all digits | same |
| S(rho_I) | 0.33333 | C/3 = 0.33333 |
| S(rho_I || tau_beta), beta = pi | 0.08333 | 1/12 + 1/3 - 1/3 = 0.08333 |
| S(rho_I || tau_beta), beta = 2 pi | 0.00000 | 0 |
| S(rho_I || tau_beta), beta = 4 pi | 0.08325 | 0.08333 |
| S(rho_I || tau_beta), beta = 8 pi | 0.3682 (rising) | 0.375 |

The defect slope converges slowly from below (the increments shrink like 1/log ell, consistent with the lattice cutoff epsilon entering as log(ell/epsilon) with an effective epsilon that drifts); the relative entropy at beta = pi converges fast because its leading term is the exact sum rule beta A/(12 pi) + pi B/(3 beta) - C/3 = 1/12 at beta = pi for A = B = C, whose finite-size corrections cancel between the energy and the entropy. Two decoupled XX chains with hoppings 1/2 and 1/4 (velocities 1 and 1/2, d = 4, M = 3/4, A = 1.5, B = 3, C = 2): measured S(rho_I || tau_beta) slopes 0.45832, 0.08333, 0.08325, 0.4515 at beta = pi, 2 pi, 4 pi, 8 pi against the predicted 0.45833, 0.08333, 0.08333, 0.45833, and defect slope 0.0388 (rising) against A/(12 pi) = 0.0398. The two-velocity case is the one where no beta makes tau_beta close to rho_I; the minimum over beta of the predicted coefficient is (1/3)(sqrt(4.5) - 2) = 0.0404 > 0.

Boulware profile at ell = 2048 (three-bond smoothing to remove the alternating Friedel component, which has amplitude about 0.03/x): x^2 (ell-x)^2 (e_0 - <h_x>)/ell^2 = 0.0148, 0.01381, 0.01361, 0.01359, 0.01358 at x = 16, 64, 256, 512, 1024, against A/(24 pi) = 0.01326, so the local Schwarzian form holds to about 2.5 percent at these sizes, with the residual presumably the same slow drift as the integrated slope.

`xx_weights.py` (output `xx_weights.out`): tent weight slope 0.02575 (rising) against 1/(12 pi); quartic weight defect 0.017680 constant in ell to six digits; linear ramp with a hard wall at ell grows as 0.0684 ell, which is one boundary's worth of -f_b = (2/pi - 1/2)/2 = 0.0683 from the wall, on top of the (1/(24 pi)) log ell from the soft edge. The open-chain boundary energy and Casimir term are reproduced exactly: f_b(n) = 1/2 - 2/pi - pi A/(24 (n+2)) for odd bond numbers n (even site numbers, no zero mode).

[Numerical, exact diagonalization] `ed_modular.py` (outputs `ed_heis.out`, `ed_xxz.out`, `ed_su3.out`), charge-sector Lanczos. Weighted open intervals to ell = 18 (d = 2) and 13 (d = 3), rings to L = 18 (d = 2) and 15 (d = 3) reduced to ell = L/2.

| chain | A/(12 pi) | defect slope at the largest sizes (same parity or residue class) | S(rho_I || tau) at beta v = 2 pi, ell = 5..9 | slope of S(rho_I || tau) at beta v = pi | f_b(n) + pi A/(24 n) |
|---|---|---|---|---|---|
| Heisenberg | 0.04167 | 0.0396, 0.0398, 0.0399, 0.0401 (ell = 10 to 18, even) | 0.0130, 0.0133, 0.0134, 0.0135, 0.0135 | 0.084 (CFT 0.0833) | -0.2367, -0.2469, -0.2500, -0.2514 (n = 5, 9, 13, 17) |
| XXZ Delta = 0.5 | 0.03446 | 0.0321, 0.0322, 0.0323 (ell = 12 to 18, even) | 0.0129, 0.0132, 0.0133, 0.0133, 0.0134 | 0.083 (CFT 0.0833) | -0.1889, -0.1964, -0.1985, -0.1994 |
| SU(3) swap | 0.11111 | 0.0999, 0.0977 (ell = 6, 9, 12) | 0.013, 0.028, 0.022 (ell = 4, 6, 7, residue effects) | 0.167 between ell = 6 and 7 (CFT 0.1667) | -0.443, -0.467, -0.476 (n = 5, 8, 11) |

Odd ell (d = 2) and ell not divisible by 3 (d = 3) carry an unpaired spin or colour in the K_w ground state and are excluded from the slopes. The defect slopes sit 5 to 12 percent below A/(12 pi) at these sizes and rise with ell, the same pattern as the exact XX data at much larger ell; the relative-entropy slopes, which are sum rules, are already at the CFT value. The layer-cake inequality E_0(K_w) >= sum_k Delta_k E_0^open(n_k) holds at every size (column "E0-lc" in the outputs, all nonnegative after fixing the colour sectors), and its allowance grows linearly, e_0 W - sum_k Delta_k E_0^open(n_k) = 0.0342 ell for XX = (ell/4)(-f_b).

## 4. What the reformulation says and what positivity cannot say

[Proved, given Section 2] The conjecture A <= pi M log_2 d is equivalent to each of

    (i)   lim_ell [ e_0 W - E_0(K_w) ] / log ell <= (M log_2 d)/12,
    (ii)  lim_x x^2 [ e_0 - <h_x>_Boulware ] <= (M log_2 d)/24,
    (iii) lim_ell lim_beta S(rho_I || tau_beta) / (beta log ell) <= (M log_2 d)/12 .

Heisenberg saturates all three with 1/24, 1/48, 1/24. In words: the ground state of the lattice Rindler Hamiltonian sum_x x h_x may lower the local energy at distance x from its edge by at most (M log_2 d)/(24 x^2). The constraint from energies alone is far weaker: for any state, sum over any segment of (e_0 - <h_x>) <= M (open-chain superadditivity, Lemma 1 of `contrib/fable_bound_proof`), which allows a defect of order M concentrated anywhere.

[Proved] Why relative-entropy positivity is velocity neutral. Positivity gives lower bounds on log Z_w(beta), hence lower bounds on e_0 W - F_w(beta) and on the defect. An upper bound on A needs an upper bound on log Tr exp(-beta sum_x w_x (h_x - e_0)), which no relative entropy supplies. Among the lattice states whose entropies are computable in the CFT, ground-state reductions and their normalized powers have velocity-blind entropies, Gibbs states of lattice Hamiltonians have entropies weighted by 1/v_a, and the K_w energy of any of them is either e_0 W or below it by an A-weighted logarithm. Every positivity statement S(sigma||tau) >= 0 with sigma, tau in this family is therefore an inequality among the three numbers A, B, C of the form sum_a c_a (x_a - y_a)^2 >= 0, with x_a, y_a depending on n_a or on the Renyi index. Two examples worked out: S(rho_I || tau_beta) gives sum_a c_a (sqrt(n_a) - 1/sqrt(n_a))^2 >= 0 and, minimized over beta, sqrt(AB) >= C; S(tau_beta || rho_I^n/Tr rho_I^n) gives sum_a c_a (sqrt(n)/n_a - 1/sqrt(n))^2 >= 0. The Gibbs-Bogoliubov inequalities between tau_beta's of different weights are the convexity of log Z in beta, which the CFT formula n + 1/n satisfies. None contains M or d, because the only place where M and d enter the lattice relative entropies is through F(0) = ell log d and |e_0|, ||h|| <= M, which are extensive or O(M) and enter the log ell coefficient nowhere.

[Proved, elementary upper bounds on log Z] Each standard upper bound loses M times the maximal weight. Dimension counting: log Z_w(beta) <= ell log d + beta M W. Golden-Thompson with the even-odd split K_w = K_e + K_o (each a sum of commuting bond terms): log Z_w <= beta M W_o + sum_{x even} log Tr exp(-beta w_x h), which as beta -> infinity gives E_0(K_w) >= -M W, the same as dimension counting. Peierls-Bogoliubov with any product reference state gives E_0(K_w) >= e_0 W - (product-state energy excess) W, again linear in W. The layer-cake bound of Section 5 is the sharpest of these, and it loses (max w)(-f_b) with -f_b <= M.

[Heuristic] Gauge dependence. Replacing h by h + o (x) 1 - 1 (x) o adds sum_x (w_x - w_{x-1}) o_x to K_w, a position-dependent field of relative strength w'/w = 1/x near the edges. If o contains a conserved current (a chemical potential mu), the trace anomaly argument of Section 2 with a chemical potential gives an extra defect kappa mu^2 int w'^2 / w = 2 kappa mu^2 log ell with kappa the compressibility, so the defect coefficient is a convex function of the gauge minimized at the reflection-symmetric representative, where it equals A/(12 pi). The reformulation (i) should be read in that gauge; M is gauge-optimized separately. The lattice relative-entropy identity holds in every gauge, with its own E_0.

## 5. The layer-cake inequality and the sign of the boundary energy

[Proved, Proposition 4] Let w be nonnegative and unimodal on the bonds of an interval (nondecreasing then nonincreasing), with distinct values t_1 < t_2 < ... < t_m. Set Delta_k = t_k - t_{k-1} (t_0 = 0) and I_k = {x : w_x >= t_k}, a nested family of intervals with n_k bonds. Then K_w = sum_k Delta_k H_{I_k} exactly, where H_{I_k} = sum_{x in I_k} h_x is the open-chain Hamiltonian on I_k, and since the ground-state energy is superadditive,

    E_0(K_w) >= sum_k Delta_k E_0^open(n_k),      Delta(ell) <= sum_k Delta_k g(n_k),     g(n) = n e_0 - E_0^open(n).

With g(n) <= M (Lemma 1 of `contrib/fable_bound_proof`, which also gives g(n) >= -M) and sum_k Delta_k = max w this is Delta(ell) <= M max w.

[Proved, Corollary 5] Insert the open-chain asymptotics [CFT] E_0^open(n) = n e_0 + f_b - pi A/(24 n) + o(1/n), so g(n) = -f_b + pi A/(24 n) + o(1/n). For the quadratic weight the level sets have n(t) = ell sqrt(1 - 4t/ell) bonds at height t, so sum_k Delta_k / n_k -> int_0^{ell/4} dt / n(t) = 1/2, and

    (A/(12 pi)) log ell + O(1) = Delta(ell) <= (ell/4)(-f_b) + pi A/48 + o(1).

Hence f_b < 0 for every critical chain with A > 0, in every gauge (if f_b were zero the right side would be bounded while the left grows). Boundary fields that are gauge terms only lower E_0^open at second order (first order vanishes by reflection symmetry), so the free-end f_b is the largest value over gauges, and it is still negative. Checks: XX f_b = 1/2 - 2/pi = -0.137; the ED extrapolations above give about -0.25 (Heisenberg), -0.20 (XXZ at Delta = 0.5), -0.48 (SU(3)).

[Proved] Loss factor. For the tent weight of height R and half-width R, w = min(x, 2R - x), the level sets have n(t) = 2(R - t) bonds, and sum_k Delta_k / n_k = (1/2) sum_{j=1}^{R} 1/j = (1/2) log R + O(1), so the layer cake gives Delta <= (-f_b) R + (pi A/48) log R while the truth is (A/(12 pi)) log R: the logarithmic coefficients differ by pi^2/4 = 2.47, and the O(M R) boundary term dominates both. Energy superadditivity loses the boundary energy of every layer, and no choice of layers avoids it because every layer is an interval with two ends. Replacing the intervals by rings is not available (a non-constant weight on a ring is not a sum of nested rings), and modifying the layers by boundary operators that cancel f_b reintroduces an extensive sum of those operators weighted by |w_x - w_{x-1}|, whose norm is again of order M max w.

## 6. What remains

The modular route reduces the conjecture to a statement about one lattice state, the ground state of the Rindler Hamiltonian sum_x x h_x, whose energy density at distance x from the edge lies below e_0 by A/(24 pi x^2) + o(1/x^2) [CFT, checked numerically]. The missing lattice input is an upper bound

    e_0 - <h_x>_Boulware <= (M log_2 d)/(24 x^2) + o(1/x^2),

or, integrated, e_0 W - E_0(K_w) <= (M log_2 d /12) log ell + O(1) for the quadratic weight. Facts established about this quantity: it is the large-beta slope of the relative entropy S(rho_I || tau_beta); it is attained within one percent by the top Schmidt vector of rho_I (exactly for free fermions); it is bounded by M max w by every energy argument, which is the wrong order; it is minimized over gauges at the reflection-symmetric bond; and it is equivalent to the conjecture with the same constant, saturated by the Heisenberg chain. The one route not excluded here is a bound on the second-order response of a translation-invariant ground state to a gradient perturbation sum_x (x - X) h_x at distance X, using more than its norm (for instance the fact that the Boulware state is the ground state of a Hamiltonian that differs from X H only by a term odd under reflection about X, so that the first-order response vanishes and the second-order one is a specific susceptibility of the vacuum). I did not find an argument bounding that susceptibility by M log d; the obstruction found in `contrib/fable_kappa` (the gradient of A is unbounded near z = 2 points) applies to it as well, since the Boulware defect is exactly a velocity response.

## 7. Which CFT quantities carry the velocities

| quantity | log ell coefficient | velocity dependence |
|---|---|---|
| S(rho_I) | C/3 | none |
| -log lambda_max(rho_I) | C/6 | none |
| K_rho (modular Hamiltonian) | 2 pi sum_a (1/v_a) int w eps_a^phys | 1/v_a in front of the physical energy density |
| e_0 W - E_0(K_w) (lattice entanglement Hamiltonian, Boulware state) | A/(12 pi) | sum c_a v_a |
| S(tau_beta) | (2 pi/(3 beta)) B | sum c_a / v_a |
| log Z_w(beta) + beta e_0 W | beta A/(12 pi) + pi B/(3 beta) | both |
| S(rho_I || tau_beta) | beta A/(12 pi) + pi B/(3 beta) - C/3 | both, minimum (1/3)(sqrt(AB) - C) |
| ring Casimir energy e_0 L - E_0(L) | (pi/6) A /L (not a log) | sum c_a v_a |
| thermal entropy density | (pi/3) B T | sum c_a / v_a |

## 8. Sources used

`specs/_bound_background.md` (setting, round-one update), `contrib/fable_bound_proof/README.md` (Lemma 1 on open-chain energies, the two-site reduction, Section 4(ii) remark on modular routes), `contrib/fable_kappa/README.md` (gradient unboundedness), `contrib/fable_bound_freefermion/README.md` (free-fermion conventions, M = ||B||_1), and `codex/work_C98` to `work_C101` RESULT files (energy-gluing losses, entanglement-rate obstructions, thermal routes, solvable-class checks). The CFT inputs (interval Renyi entropies, the Bisognano-Wichmann and Cardy-Tonni modular Hamiltonians, the trace anomaly for a position-dependent velocity, the Unruh energy density, and the Eisler-Peschel commuting matrix for free fermions) are standard results quoted from memory and checked numerically here.
