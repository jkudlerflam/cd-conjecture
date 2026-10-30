# The reconnection lemma: an independent attempt

2026-10-30. Natural logarithms throughout. Every claim carries one of the tags [proved], [numerical], [heuristic], [assumed]. [assumed] marks either an imported published theorem or the infrared CFT input; which one is said each time.

## Claim

The reconnection lemma of `specs/_reconnection_lemma.md` (C121's open lemma P) was not proved for interacting chains. What this folder establishes is the following.

1. [proved] An exact decomposition of the reconnection functional,

       G = R_L - 2 eps_L + eps_D - H_2(pi_L),

   where eps_L, eps_D are the pinching losses of the single and doubled rings and H_2(pi_L) is the Renyi-2 entropy of the normalized record distribution. Every term other than R_L lies in [0, c N log d] with c = 1, 2, 1 respectively.

2. [proved] Lemma P is equivalent, with explicit constants and for any fixed a > 0, to the all-temperature ratio bound R_L(beta) <= K(d) + kappa(d) M beta/L. It implies, for every ring length L, the ground-state statement

       L (E_2L - 2E_L) <= [b_1(d) + 3 log d] M / a.

   More sharply, at fixed L and fixed slab length tau, the cost per pinching interval satisfies

       tau Delta_L - 3 log d <= liminf G/N <= limsup G/N <= tau Delta_L + 2 log d,   Delta_L = E_2L - 2E_L,

   as N tends to infinity with beta = N tau. The cross terms that survive between consecutive pinchings therefore cost tau Delta_L per interval, up to O(log d), and nothing else. At tau = aL/M this is a (L Delta_L/M), the dimensionless Casimir defect.

3. [proved] Continuous pinching makes the gluing exact for every finite chain, stoquastic or not: as N tends to infinity at fixed beta, G tends to a limit G_inf <= 0. The price is a pinching loss growing like beta (E_Zeno - E_L), which is of order beta M, not beta M/L.

4. [proved] For stoquastic chains (in any finite Trotterization) the exact reconnection identity is R_L = -log E_{x,y}[exp(-DF(x,y))] up to a nonnegative correction, where x,y are independent seam column histories drawn from the L-ring measure and DF(x,y) = log[B_xx B_yy/(B_xy B_yx)] is a four-point reconnection free energy. Coarse records leave exactly the cross terms X = sum over distinct histories with equal records of B_xx B_yy, to be compared with Y = sum over distinct pairs of B_xy B_yx.

5. [proved, assuming the Bravyi-Vyalyi structure theorem] Translation-invariant commuting nearest-neighbour chains satisfy R_L <= log d for all L and beta, so lemma P holds for them with b_0 = log d, b_1 = 2 log d at every a.

6. [assumed infrared CFT and clustering; then proved] If the seam projector correlators factorize at Euclidean separations of order L (true in a unitary CFT with power-law corrections), then G = R_L + N log[sum p^4/(sum p^2)^2] + o(1) at fixed a, N, where p is the one-site marginal. So G <= R_L + o(1), and lemma P holds with b_1 = (pi a/4)(A/M) + o(1). Conversely item 2 forces b_1 >= (pi a/4)(A/M) - 3 log d. Under the infrared CFT, a uniform b_1(d) exists if and only if A/M is bounded at fixed d. The lemma is the target, not a step towards it.

7. [proved implication under the infrared CFT; heuristic about the route] The suggested level-spacing argument (interval tau ~ L/M, spacing ~ v/L) controls each sector separately, giving b_1 <= pi e a sum_a c_a via the proved Lieb-Robinson ceiling v_a <= 4eM (C109). It bounds the total central charge, which is not bounded in terms of d. It cannot close the gap.

The strongest new statement proved is item 2 together with item 1. The remaining gap is stated in the last section.

## Status

Lemma P for general interacting chains: open. Items 1 to 5 are proved, item 6 is proved conditional on the infrared CFT plus clustering, item 7 is an analysis of a proposed route. All candidate inequalities were checked numerically (XX by free fermions to L = 8192, Heisenberg, XXZ, XXZ in a field and SU(3) by full-spectrum ED, with larger Typhon runs submitted but still queued). What would upgrade the status: a proof of L(E_2L - 2E_L) <= C(d) M together with control of the thermal remainder, which is the original target in its zero-temperature form.

## How it was established and checked

- Author: Claude (Opus 5.5), working independently of codex/work_C124, for Jonah Kudler-Flam. No Codex or other external model was used.
- Inputs read: specs/_bound_background.md, specs/_reconnection_lemma.md, codex/work_C118, C121, C123 (RESULT and PROOF). The XX free-fermion history formulas are copied unchanged from codex/work_C123/checks.py into `xx_c123.py`.
- Checked by: the numerical checks below. No independent human or agent review.

## Argument

Notation as in C121 and C123. H_L is the L-ring, Z_L = Tr exp(-beta H_L), R_L = 2 log Z_L - log Z_2L, N >= 1 pinchings at spacing tau = beta/N, P_s the seam projectors in a fixed basis (d of them), q_s >= 0 the recorded history weights, Q_L = sum_s q_s, S = sum_s q_s^2, Q_D the doubled-ring channel trace with sites 0 and L pinched, G = log(S/Q_D), and N = max(1, ceil(beta M/(aL))).

### 1. The exact decomposition [proved]

Define pi_s = q_s/Q_L, H_2(pi) = -log sum_s pi_s^2, eps_L = log Z_L - log Q_L and eps_D = log Z_2L - log Q_D. Then S = Q_L^2 exp(-H_2(pi)), so

    G = 2 log Q_L - H_2(pi) - log Q_D
      = 2 log Z_L - 2 eps_L - H_2(pi) - log Z_2L + eps_D
      = R_L - 2 eps_L + eps_D - H_2(pi).                               (1.1)

C121 proved 0 <= eps_L <= N log d and 0 <= eps_D <= 2N log d. Since pi is supported on d^N records, 0 <= H_2(pi) <= N log d. Hence

    R_L - 3N log d <= G <= R_L + 2N log d,                             (1.2)

which is C121 (33a), now with each loss identified. The slack in C121's record Cauchy-Schwarz step is N log d - H_2(pi). The comparison of S with Q_D contains no information beyond R_L and these three bounded entropic terms.

### 2. Equivalence with the ratio bound and the per-interval cost [proved]

(a) Lemma P implies the ratio bound. Suppose G <= b_0 + b_1 N for all L, beta with N = max(1, ceil(beta M/(aL))). Then by (1.2) and N <= 1 + beta M/(aL),

    R_L <= b_0 + (b_1 + 3 log d) N <= b_0 + b_1 + 3 log d + [(b_1 + 3 log d)/a] M beta/L.    (2.1)

(b) The ratio bound implies lemma P. Suppose R_L <= K + kappa M beta/L. Since N >= beta M/(aL) in both branches of the max, beta M/L <= aN, and (1.2) gives

    G <= R_L + 2N log d <= K + (a kappa + 2 log d) N.                  (2.2)

So lemma P holds with b_0 = K, b_1 = a kappa + 2 log d. The two statements are equivalent up to these explicit constants, at every a.

(c) The necessary ground-state statement. Write R_L = beta Delta_L + 2 log W_L - log W_2L with W_n = Tr exp[-beta(H_n - E_n)] (C118 (1)). At fixed L, W_n tends to the ground degeneracy g_n in [1, d^n] as beta grows. Divide (2.1) by beta and let beta tend to infinity:

    L Delta_L <= (b_1 + 3 log d) M/a.                                   (2.3)

(d) Per-interval cost. Fix L and tau > 0, let beta = N tau and N tend to infinity. Then R_L/N = tau Delta_L + (2 log W_L - log W_2L)/N tends to tau Delta_L, since the W terms converge to finite positive limits. Dividing (1.2) by N,

    tau Delta_L - 3 log d <= liminf G/N <= limsup G/N <= tau Delta_L + 2 log d.        (2.4)

With the prescription N = ceil(beta M/(aL)) and beta to infinity at fixed L, tau tends to aL/M and (2.4) reads a L Delta_L/M + O(log d). C123 showed lim G/N equals the log ratio kappa_L(tau) of the top eigenvalues of its two temporal transfer maps, so (2.4) also pins that spectral rate: kappa_L(tau) = tau Delta_L + O(log d) for every finite chain. This is the precise answer to "which cross terms survive between consecutive pinchings": per interval they cost the reconnection energy tau Delta_L, and nothing that depends on the seam basis beyond O(log d).

(e) Consequence. A bound b_1(d) on the cost per interval is the same thing as a bound L Delta_L <= C(d) M on the finite-size ring defect, up to thermal remainders. Under the vacuum Casimir expansion Delta_L = pi A/(4L) + o(1/L) [assumed infrared CFT, C116], (2.3) gives A <= (4/(pi a))(b_1 + 3 log d) M, which is the conditional theorem of C121. No version of lemma P can be weaker than the target.

### 3. Why level spacing does not enter [proved under the infrared CFT; heuristic about the route]

The proposed per-interval estimate multiplies tau ~ aL/M by a level spacing ~ v/L, giving a v/M = O(1) by C109's v_a <= 4eM. By (2.4) the quantity that must be O(1) is tau Delta_L, not tau times a gap. Under the infrared CFT with decoupled sectors, L Delta_L tends to (pi/4) sum_a c_a v_a. Each sector contributes (pi a/4) c_a v_a/M <= pi e a c_a to the per-interval cost. Summing, the route proves b_1 <= pi e a sum_a c_a + O(log d), a bound by the total central charge. The background records constructions with c = n + 1 at fixed d and arbitrarily slow extra sectors, so sum_a c_a is not bounded by any f(d). The level spacing controls each sector's velocity, and the open problem is the velocity-weighted sector count itself. C123 Section 7 gave the related observation that a gap alone leaves the low-energy multiplicity uncontrolled.

### 4. Exact gluing in the continuous-pinching limit [proved]

Fix L and beta and let N tend to infinity. The numerator has C123's transfer form S_N = <u_S, T_S^(N-1) u_S> on two copies H_L (x) H_L, where the compression is J_S(X) = sum_s Pi_s X Pi_s with Pi_s = P_s (x) P_s. These Pi_s are orthogonal projections, not summing to the identity. For X in the range of J_S, X = sum_s Pi_s X Pi_s, so Pi_s {K, X} Pi_s = {Pi_s K Pi_s, X_s} with X_s = Pi_s X Pi_s. The finite-dimensional product-formula argument of C121 Section 6.2 then applies verbatim with D replaced by J_S, and with the initial operator J_S(I) = sum_s Pi_s. With K_S = H_L (x) 1 + 1 (x) H_L,

    lim_N S_N = sum_s Tr_{Pi_s} exp[-beta Pi_s K_S Pi_s] = sum_s z_ss^2,

where z_st = Tr exp[-beta H_open(s,t)] and H_open(s,t) is the open chain on the L-1 non-seam sites with the boundary fields <s|h|s> and <t|h|t> left by pinning the seam neighbours. Indeed P_s h_{L-1,0} P_s = (1 (x) <s|) h (1 (x) |s>) (x) P_s, an operator on site L-1 alone, and similarly for the bond (0,1); all other bonds commute with P_s. The same argument on the doubled ring, pinched at 0 and L with labels (s,t), splits it into two open chains with boundary labels (s,t) and (t,s):

    lim_N Q_D = sum_{s,t} z_st z_ts.

Every z_st > 0. Hence

    G_inf = log[sum_s z_ss^2] - log[sum_{s,t} z_st z_ts] <= 0,            (4.1)

with the diagonal terms reproducing the numerator exactly. This is the classical gluing proof, and it needs neither stoquasticity nor a gap: a continuously pinched seam is classical. The cost sits in eps_L, which tends to log Z_L - log sum_s z_ss. At fixed L and beta to infinity this is beta(E_Zeno - E_L) + O(1), where E_Zeno is the smallest ground energy of the H_open(s,s). C121 computed E_Zeno - E_L = 1/2 + tan(pi/(2L))/2 tending to M for XX. Exact gluing costs time-extensive energy of order M per unit imaginary time, while the target allows order M/L.

Item 4 is a sharp formulation of the obstruction. Sparse records leave the seam quantum for a time tau, and the reconnection then costs tau Delta_L (Section 2). Dense records make it classical at a cost of order M per unit time. Neither regime has a loss of order M/L per unit time unless L Delta_L is already known to be O(M).

### 5. Stoquastic chains: the exact cross terms [proved]

Take a stoquastic bond (nonpositive off-diagonal entries in a product basis) and any checkerboard Trotterization with K slices. A column history x lies in [d]^K, and the spatial transfer matrix B = (AC)^(L/2) between the seam column and itself is entrywise nonnegative. Then Z_L = Tr B and Z_2L = Tr B^2 for the Trotterized model.

(a) Reconnection identity. Let mu_x = B_xx/Z_L and DF(x,y) = log B_xx + log B_yy - log B_xy - log B_yx for pairs with B_xy B_yx > 0. Then

    exp(-R_L) = Z_2L/Z_L^2 = E_{mu (x) mu}[exp(-DF)] + Z_L^(-2) sum_{B_xx B_yy = 0} B_xy B_yx,

so R_L <= -log E_{mu (x) mu}[exp(-DF)] <= E_{mu (x) mu}[DF] by Jensen, with E[DF] = +infinity whenever some mu-typical pair has B_xy B_yx = 0. DF is a mixed second difference of the strip free energy F(x,y) = -log B_xy. It cancels every boundary term that depends on x or y alone and retains only the interaction between the two seam histories across width L. In the crossed channel it is a boundary-condition-changing free energy, heuristically of order v beta/L in a CFT [heuristic].

(b) Coarse records. Group histories by their values at N of the K slices, c(x) in [d]^N, and define the coherent record weights a_s = sum_{c(x)=s} B_xx (C121 eq. 20). Then

    sum_s a_s^2 = D + X,   X = sum_{x != y, c(x) = c(y)} B_xx B_yy,
    Z_2L        = D + Y,   Y = sum_{x != y} B_xy B_yx,       D = sum_x B_xx^2.        (5.1)

The diagonal D is common. The coherent reconnection inequality of C121 (32) holds with zero cost if and only if X <= Y, and with cost exp(c) if and only if D + X <= exp(c)(D + Y). These are all the surviving cross terms. When the records resolve every Trotter slice (N = K), X = 0 and the classical proof applies, but the record count d^K diverges with the Trotter number. This is item 4 in Trotter language.

(c) What positivity does and does not give. Nonnegativity gives D <= Z_2L, hence the classical bound for fine records. It gives no comparison between B_xx B_yy and B_xy B_yx for x != y. C121's example B = I already shows X can exceed Y by any factor. Bounding X/Y per interval is again equivalent, via (1.1) to (2.4), to bounding the reconnection energy tau Delta_L.

### 6. Commuting chains [proved, assuming Bravyi-Vyalyi]

[assumed published] Bravyi and Vyalyi (quant-ph/0308021) prove that commuting 2-local terms acting on a qudit induce a decomposition of its space as a direct sum over alpha of tensor products, with each term acting blockwise on its own tensor factor. For a translation-invariant chain whose bonds all commute, site j decomposes as the direct sum over alpha of H^L_alpha (x) H^R_alpha, and h_{j,j+1} is block diagonal in (alpha_j, alpha_{j+1}), acting on H^R_{alpha_j} (x) H^L_{alpha_{j+1}}.

[proved] Then H_L is block diagonal in the sector string (alpha_1, ..., alpha_L), and in each block it is a sum of commuting operators on disjoint factors. Hence Z_L = Tr T^L with the m x m entrywise positive matrix T_{alpha beta} = Tr exp(-beta h^{alpha beta}), m <= d. The classical inequality of C118 Section 3 gives (Tr T^L)^2 <= m Tr T^(2L), so R_L <= log m <= log d for all L and beta. By (2.2), lemma P holds with b_0 = log d and b_1 = 2 log d at every a. These chains have A = 0, so this class does not test the velocity-weighted bound.

### 7. The infrared regime [assumed infrared CFT and clustering; then proved]

Hypothesis (H). Fix a and N. As L tends to infinity along vacuum-compatible lengths at beta = N tau and tau = aL/M, the thermal correlators of any product of seam projectors at contour times pairwise separated by at least tau/2 factorize, <prod_k P_{s_k}(t_k)> = prod_k p_{s_k} + o(1), with p the one-site marginal, and the same holds for the two seams of the 2L ring at spatial separation L, with a marginal p' - p = o(1). In a unitary CFT the corrections are powers L^(-x) with x the smallest scaling dimension in the seam projectors. C123 proved (H) for XX with errors O(L^(-2)).

Derivation. Each record label appears twice on the reflected contour, so under (H) q_s/Z_L = prod_j p_{s_j}^2 + o(1), Q_L/Z_L = (sum p^2)^N + o(1) and Q_D/Z_2L = (sum p^2)^(2N) + o(1). Therefore eps_L = -N log sum p^2 + o(1), eps_D = 2 eps_L + o(1), H_2(pi) = -N log[sum p^4/(sum p^2)^2] + o(1), and by (1.1)

    G = R_L + N log[sum_s p_s^4/(sum_s p_s^2)^2] + o(1).                (7.1)

The logarithm is between -log d and 0 by Cauchy-Schwarz. For a maximally mixed marginal it is -log d, which is C123's XX formula G = R_L - N log 2 + O(L^(-2)).

Consequence. With R_L tending to a sum of torus ratios with large-argument slope (pi/4) A beta/L, (7.1) gives G <= (pi a/4)(A/M) N + O(1). Item 2 gives the converse. Under (H) and the infrared CFT, the best b_1 is (pi a/4)(A/M) + O(log d). A universal b_1(d) is the same as the target bound on A/M. Formula (7.1) is a fixed-N limit; uniformity as N grows with L is not established by this argument, as in C123 for XX.

Caveat [numerical]. The torus statement R_L <= (pi/4) A beta/L holds only asymptotically. For XX at L = 16, L Delta_L = 0.78983 exceeds pi/4 = 0.78540, and R_L exceeds pi beta/(4L) by 0.177 at beta = 40L. The excess decays like L^(-2).

## Numerical checks

Run from this folder with numpy and scipy. Local runtimes are on a laptop with four threads.

| Script | What it checks | Runtime | Output |
| :--- | :--- | ---: | :--- |
| `xx_run.py` | (1.1), (1.2), (7.1) and G <= R_L for XX at L = 8 to 8192, a = 0.25, 1, 2, N = 1, 2, 4; torus inequality scan | 5 s | `xx_results.json`, `xx_run.txt` |
| `ed_run.py local` | the same for XX, XXZ (Delta = -0.5, 0.5), Heisenberg, XXZ (Delta = 0.5) in field 0.4, SU(3); L = 4, 6 (SU(3) L = 3) | 4 s | `ed_local.json`, `ed_local.txt` |
| `ed_run.py typhon_qubit`, `typhon_su3` | qubit L = 8 (doubled 16), SU(3) L = 6 (doubled 12), SLURM scripts `run_typhon_*.slurm` | Typhon, queued | `ed_typhon_qubit.json`, `ed_typhon_su3.json` (pending) |
| `zeno_run.py` | convergence of G_N to G_inf and (4.1) | 3 s | `zeno_results.json`, `zeno_run.txt` |
| `stoquastic_run.py` | (5.1), the reconnection identity, the Trotter transfer against dense Trotter products | 30 s | `stoquastic_results.json`, `stoquastic_run.txt` |
| `commuting_run.py` | Section 6 on a random Bravyi-Vyalyi chain with d = 3, m = 2 | 2 min | `commuting_results.json`, `commuting_run.txt` |

`ed_core.py` holds the block ED (every eigenvalue retained, symmetry multiplicities used only for record-summed quantities). `xx_c123.py` is the copied C123 free-fermion code.

### Results

[numerical] Identity (1.1) holds to 4e-15 in all 99 local ED cases and to 2e-15 in all 54 XX cases. The pinching and entropy bounds of (1.2) hold everywhere. The explicit history sum agrees with independent channel propagation to 2e-14.

[numerical] Clustering formula (7.1) on XX: G - R_L + N log 2 is -9.9e-9 at L = 8192, a = 1, N = 4, and -6.5e-4 at L = 32. Every component converges separately, eps_L/N and H_2/N to log 2 and eps_D/N to 2 log 2, so all three losses saturate their maximal values in the infrared.

[numerical] Interacting ED, G - R_L at a = 1, N = 4 versus the clustering value from (7.1):

| Model | L | G - R_L | (7.1) value | L Delta_L/M |
| :--- | ---: | ---: | ---: | ---: |
| XX | 4, 6 | -2.820, -2.792 | -2.773 | 1.722, 1.636 |
| XXZ Delta = -0.5 | 4, 6 | -2.915, -2.833 | -2.773 | 1.066, 1.042 |
| XXZ Delta = 0.5 | 4, 6 | -2.860, -2.813 | -2.773 | 2.282, 2.146 |
| Heisenberg | 4, 6 | -3.086, -2.960 | -2.773 | 2.791, 2.618 |
| SU(3) | 3 | -5.000 | -4.394 | 4.183 |

[pending] The larger sizes (qubit L = 8 with doubled ring 16, largest block 12870; SU(3) L = 6 with doubled ring 12, largest block 34650) were submitted to Typhon as SLURM jobs 17530900 and 17530901 (`run_typhon_qubit.slurm`, `run_typhon_su3.slurm`). At the time of writing they were queued behind other jobs, with an estimated start of 2026-10-31. They write `ed_typhon_qubit.json` and `ed_typhon_su3.json` in the remote run folder; copy them here and add the a = 1, N = 4 rows to this table. No conclusion above depends on them.

The deviation from (7.1) is always negative and shrinks with L, consistent with power-law clustering corrections. The XXZ chain in a field is not vacuum-compatible at these sizes (its ground magnetization differs between L and 2L), and is kept only as a check of (1.1) with a non-uniform marginal; at a = 0.25 its marginal is p = (0.418, 0.582) and (1.1) still holds to machine precision.

[numerical] Candidate inequality G <= R_L. It holds in all 99 local ED cases (largest G - R_L is -0.654), in all 54 XX cases, and was not yet tested at the Typhon sizes. I have no proof. By (1.1) it is equivalent to Q_D/Z_2L >= sum_s (q_s/Z_L)^2, a comparison between the doubled ring's pinching survival and the collision probability of two independent L-ring records. A proof would sharpen (1.2) to R_L - 3N log d <= G <= R_L, which does not touch the gap.

[numerical] Continuous pinching (Section 4), fixed beta, N = 1 to 10 (SU(3) to 7):

| Model, L, beta | G_inf | E_Zeno - E_L | G_N at largest N | eps_L at largest N |
| :--- | ---: | ---: | ---: | ---: |
| XX, 4, 4 | -0.6931 | 0.7071 | -0.7747 | 1.42 |
| XX, 4, 8 | -0.6931 | 0.7071 | -0.9688 | 3.25 |
| Heisenberg, 4, 4 | -0.3205 | 0.6340 | -0.4225 | 1.50 |
| Heisenberg, 4, 8 | -0.0433 | 0.6340 | -0.2279 | 3.15 |
| SU(3), 3, 1.5 | -1.1315 | 2.0 | -1.3129 | 1.35 |
| SU(3), 3, 3 | -0.7301 | 2.0 | -0.9742 | 3.58 |

G_N decreases from positive values at N = 1 (0.76, 1.89, 2.29 at the larger beta, reproducing C121) through a minimum and approaches G_inf <= 0 from below, while eps_L grows towards beta(E_Zeno - E_L). The XX value of E_Zeno - E_L matches 1/2 + tan(pi/8)/2.

[numerical] Stoquastic Trotter network (Section 5), K = 12 slices, stoquastic XX and Heisenberg images, L = 4, 6, 4096 column histories. The Trotter transfer reproduces the dense Trotter traces to 1e-10, and the reconnection identity reproduces R_L to 1e-9. The Jensen bound E[DF] is infinite in every case, because the conserved magnetization forces B_xy = 0 for incompatible column pairs, so the Jensen step is useless without first resolving conserved charges. The cross-term ratio X/Y_c, with Y_c the part of Y with equal records, exceeds one in every case (1.004 to 8.5), so the surviving cross terms are not dominated record by record. The full ratio X/Y is below one except Heisenberg at beta = 6, L = 4, N = 1, where X/Y = 1.055 and log(sum_s a_s^2/Z_2L) = +0.045 > 0. At N = K, X = 0 exactly.

[numerical] Commuting chain: the bonds sharing a site commute to 1e-12, the transfer formula reproduces log Z_L to 6e-14, and R_L <= log 2 for L = 2, 3, 4 at beta = 0.5 to 32 (largest value 0.176).

## Remaining gap

The precise open statement is the all-temperature ratio bound

    R_L(beta) <= K(d) + kappa(d) M beta/L   for all L, beta,

or, through Section 2, any version of lemma P. Its zero-temperature content is the ring defect bound L(E_2L - 2E_L) <= C(d) M, which under the infrared CFT is A <= (4/pi) C(d) M. Nothing in the coarse-seam formulation reduces this: by (1.1) the reconnection functional differs from R_L only by three terms each bounded by 2N log d, and by (2.4) the cost per pinching interval is the reconnection energy tau Delta_L up to O(log d). Stoquasticity removes signs but leaves the cross terms X of (5.1), whose size relative to Y encodes the same energy. Dense pinching (Section 4) makes the seam classical and the gluing free, at a cost of order M per unit imaginary time instead of M/L.

A proof therefore has to bound the finite-size ring defect E_2L - 2E_L by O(M/L) using microscopic information, in a way that sees the velocity-weighted sum over infrared sectors rather than the sector count. In the stoquastic language of Section 5, it would be a lower bound on the reconnection Boltzmann factor B_xy B_yx/(B_xx B_yy) of typical seam histories, of the form exp[-K - kappa M beta/L], which is a statement about how much the two boundaries of a strip of width L interact. Three directions seem worth trying: (i) a variational upper bound on E_2L built from the L-ring ground state by spreading the cut over the whole ring, since stoquastic ground states are positive and superpositions of shifted cut products cannot cancel; (ii) a crossed-channel entropy bound, since R_L is the Renyi-2 entropy of the crossed transfer state whenever that state is positive, and its low-temperature coefficient is (pi/3) A; (iii) resolving the conserved charges first in Section 5, so that the Jensen bound becomes finite and DF can be bounded charge sector by charge sector.
