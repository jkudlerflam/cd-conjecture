# C137 proofs (Claude Opus agent run)

Conventions are those of C135 PROOF section 1. D_n(X) = X (x) 1 - 1 (x) X for X on n - 1 sites, M(h) is the gauge-optimized half-spread, e_0(h) the energy density on Z, E_LTI(n) the LTI value. Pauli matrices X, Y, Z; S = sigma/2. Sites of a window are 1, ..., n and bond j joins j and j + 1.

## 1. The window class and its place in the LTI dual

[Proved, definition and Lemma 1.1] Call X on n - 1 sites a window counterterm if

    X = sum_{j=1}^{n-2} x_j h_{j,j+1} + sum_{k=1}^{n-1} y_k,          (1.1)

with real x_j and Hermitian one-site operators y_k acting on site k. Then

    h_{12} (x) 1 + D_n(X) = sum_{j=1}^{n-1} c_j h_{j,j+1} + sum_{k=1}^{n} u_k =: W(c, u),     (1.2)

with c_1 = 1 + x_1, c_j = x_j - x_{j-1} (2 <= j <= n-2), c_{n-1} = -x_{n-2}, and u_1 = y_1, u_k = y_k - y_{k-1}, u_n = -y_{n-1}. The map (x, y) -> (c, u) is a bijection onto

    A_n = { (c, u) : c in R^{n-1}, sum_j c_j = 1, u_k Hermitian on C^d, sum_k u_k = 0 },

with inverse x_j = sum_{i<=j} c_i - 1 and y_k = sum_{i<=k} u_i. Proof: X (x) 1 contains h on bonds 1..n-2 and y_k on sites 1..n-1; 1 (x) X contains h on bonds 2..n-1 and y_k on sites 2..n. Collecting coefficients gives (1.2). The sums telescope to 1 and 0, and the partial sums invert the map. []

The admissible class is therefore: arbitrary real weights c_j of either sign summing to one, one-site operators u_k summing to zero (as operators on C^d), and both may depend on h and n. The scalar part of u_k is irrelevant (it sums to zero), so u_k may be taken traceless.

[Proved, Lemma 1.2, gauge representatives] Writing each bond in its own gauge representative h~^(j) = h + a_j + o_j (x) 1 - 1 (x) o_j and weighting it by c_j produces W(c, u) - (sum_j c_j a_j) 1 with u_k = c_k o_k - c_{k-1} o_{k-1} (c_0 = c_n = 0), and sum_k u_k = 0. Conversely, if every c_j is nonzero, every sum-zero u arises in this way (solve o_1 = u_1/c_1, o_k = (u_k + c_{k-1} o_{k-1})/c_k; the last equation u_n = -c_{n-1} o_{n-1} then holds because sum u = 0). So "each bond in its own gauge representative" and "one-site terms summing to zero" are the same class whenever no weight vanishes, and the latter is the cleaner definition. Define

    W_n(h) = sup_{(c,u) in A_n} lambda_min(W(c, u)).                (1.3)

[Proved, Proposition 1.3, restriction of the LTI dual] For every bond, every n >= 2 and every L >= n,

    W_n(h) <= E_LTI(n) <= min(e_0(h), E_L(h)/L).

If (W) holds, i.e. W_n(h) >= e_0(h) - K M(h)/n^2 for all h, n, then (T) holds with the same K, hence E_L >= e_0 L - K M/L and A <= (6/pi) K M. Proof: by Lemma 1.1 every window is h (x) 1 + D_n(X) for an admissible X, so it is feasible for the LTI dual of C135 (F2) with t = lambda_min. The rest is C135 (F1) to (F3). []

[Proved, Proposition 1.4, primal of the window problem] For every h and n,

    W_n(h) = min { e : rho a state on (C^d)^{(x) n}, Tr(rho h_{j,j+1}) = e for all j, rho_k = rho_1 for all k },     (1.4)

where rho_k is the one-site marginal at site k. Proof: f(c, u; rho) = Tr[rho W(c, u)] is affine in each argument and the state space is compact and convex, so Sion's minimax theorem gives sup_{A_n} min_rho f = min_rho sup_{A_n} f. For fixed rho, sup over c with sum c = 1 of sum_j c_j Tr(rho h_j) is +infinity unless all bond energies coincide, in which case it equals their common value e; sup over sum-zero u of sum_k Tr(u_k rho_k) is +infinity unless all one-site marginals coincide, in which case it is 0. The minimum on the right of (1.4) is attained by compactness and the feasible set is nonempty (any LTI state). []

So the window problem is the relaxation of LTI_n in which, instead of equality of the two (n-1)-site end marginals, one asks only that the bond energies of the given h and the one-site marginals be translation invariant. The constraint set depends on h through its bond energies only. For a universal (h-independent) window c, the same argument gives the kinematic statement: (W) with that fixed c and free u is equivalent to the inclusion of { sum_j c_j rho_{j,j+1} : rho on n sites with equal one-site marginals } in Omega_inf + (K/n^2) B_* (support functions, cd_free_route PROOF Proposition 1). An h-dependent optimal window is not of this kinematic form.

[Proved, Lemma 1.5, elementary properties]
(a) Gauge covariance: W_n(h + a + o (x) 1 - 1 (x) o) = W_n(h) + a, since the gauge term only shifts u by the sum-zero family u_k = (c_k - c_{k-1}) o. With u optimized, the window value is independent of the gauge representative.
(b) Monotonicity: W_{n+1}(h) >= W_n(h). Pad (c, u) by c_n = 0 and u_{n+1} = 0; the extra site decouples and lambda_min is unchanged. More generally, for a probability vector p of length m, the convolution of a window with p (a convex combination of translates) is a window of length n + m - 1 whose lambda_min is at least the old one, by concavity of lambda_min.
(c) Trivial rate: e_0 - W_n <= 2M/(n - 1). Take c_j = 1/(n-1), u from a centered optimal gauge representative h~ (spectrum in [a - M, a + M]). If psi is the open-chain ground state, the product of copies of psi on consecutive blocks of n sites, averaged over translations, is a TI state on Z with energy density at most (E_open + a + M)/n. Hence n e_0 <= E_open + a + M and E_open/(n-1) >= e_0 - (M + a - e_0)/(n-1) >= e_0 - 2M/(n-1), using e_0 >= a - M.
(d) Twist blindness: if g is unitary with [h, g (x) g] = 0 and h' = (1 (x) g) h (1 (x) g*), then W_n(h') = W_n(h), e_0(h') = e_0(h) and M(h') = M(h). The open-chain unitary V = (x)_k g^k (site k carries g^k) gives V h_{j,j+1} V* = (g^j (x) g^{j+1}) h (g^j (x) g^{j+1})* = (1 (x) g) h (1 (x) g*) = h' on bond j, using the symmetry, and V u_k V* = g^k u_k g^{-k}. First replace u by its average over the closure of the group generated by the global conjugation u_k -> g u_k g* (all k at once). This global conjugation is implemented by the symmetry (x)_k g of W(c, u), so it preserves the spectrum and sum u = 0, and concavity of lambda_min shows that the average is at least as good. The averaged u commutes with g site by site, so V W_h(c, u) V* = W_{h'}(c, u) with the same admissible u. Hence W_n(h') >= W_n(h), and the converse follows because h is the twist of h' by g* (h' has the same symmetry). On Z the same product automorphism maps any TI state to a state whose translation average has the same energy density, which gives e_0(h') <= e_0(h) and the converse by symmetry. For M, average an optimal gauge o over conjugation by powers of g (the symmetry (g (x) g) of h makes every conjugate optimal), so o commutes with g and the same o is a gauge for h'. The chiral bond of C135 is the pi/2 twist of the Heisenberg bond minus 1/4, so its window constant equals the Heisenberg one exactly. Open windows do not see twists, unlike rings and LTI (C135 section 3).

[Proved, Lemma 1.6, the C135 certificate is a window] C135 section 6 sets X_n = sum_j x_j h_{j,j+1} with x_j = sum_{k<=j} a_k - 1 and no one-site terms. This is (1.1) with y = 0, and Lemma 1.1 returns c_j = a_j = v_j v_{j+1}/s and u = 0. So the C135 XX certificate is the window c = a, u = 0 in A_n.

[Proved, Remark 1.7, relation to the sine-square route set aside in _Cd_target.md] For n = L, the sine-square window c_j = (2/L) sin^2(pi j/L), j = 1..L-1, is the normalized SSD Hamiltonian of the L ring, and sum over the L translates of it is the ring Hamiltonian H_L (the identity sum_x sin^2(pi(j - x)/L) = L/2). So the window lemma restricted to the SSD window is the set-aside route, and (W) is its version with h-dependent weights and one-site terms.

## 2. The trivial classes, K = 0

[Proved, Proposition 2.1] W_2(h) = max over gauge representatives of lambda_min(h~), and by Lemma 1.5(b) W_n >= W_2. Hence if some gauge representative satisfies h~ >= e_0(h) 1 ("gauge frustration-free"), then W_n(h) = e_0(h) for every n and (W) holds with K = 0, with the uniform window and the corresponding gauge. More generally any window with all c_j >= 0 then has lambda_min >= e_0.

[Proved, Corollary 2.2, classical bonds] If h is diagonal in a product basis, e_0(h) is the minimum mean weight of a cycle in the complete digraph on {1..d} with arc weights h(s, s'), and LP duality for the minimum mean cycle problem (Karp's potential) gives a function o with h(s, s') + o(s) - o(s') >= e_0 for all s, s'. Thus h~ = h + o (x) 1 - 1 (x) o >= e_0 and K = 0. The same holds for ferromagnetic XXZ with Delta <= -1 (the fully polarized product state is a bond ground state, e_0 = lambda_min(h)), for every bond with a product TI state annihilating all bond terms after a shift, and, by Lemma 1.5(d), for all twists of these.

## 3. Optimal windows: structure observed numerically, and an exact consequence

[Proved, Proposition 3.1, ring lower bound on the window deficit] For even n with L = n, Proposition 1.3 gives e_0 - W_n(h) >= e_0 - E_n(h)/n. So the window constant is never smaller than the ring Casimir constant at L = n. If the sine-square window is exact, in the sense that lambda_min(SSD window of length n) = E_n/n, then W_n(h) = E_LTI(n) = E_n/n and the SSD window is optimal. This holds for free fermions by Katsura's theorem. [Numerical] In RESULT.md section 2 it holds to solver precision for XX at every even n <= 16 and for critical TFIM (with the canonical site weights) at every n <= 16. For Heisenberg and XXZ at Delta = 0.5, -0.9 it holds at n = 4, 6 and fails by a small amount from n = 8 on (n^2 (E_n/n - W_n)/M between 1e-4 and 2.5e-3, with bundle gaps below 1e-10, so the failure is real), so for these chains the window is very slightly weaker than the ring.

## 4. Theorem Q: an h-independent window for every quadratic qubit bond

[Proved, setting] Call h quadratic if, after a common on-site unitary and a gauge change, it has the form

    h = a 1 + J_xx XX + J_xy XY + J_yx YX + J_yy YY + (z/2)(Z (x) 1 + 1 (x) Z).      (4.1)

This is the class of C135 section 5 (XY anisotropy, transverse field, z-directed DM, pairing; no ZZ). Unequal Z coefficients z_1 ZI + z_2 IZ are brought to the symmetric form by the gauge o = -(z_1 - z_2) Z/4, so (4.1) is no restriction within the class. Write h_cross = J_xx XX + J_xy XY + J_yx YX + J_yy YY.

Let v_j = sqrt(2/(n+1)) sin(pi j/(n+1)) for 1 <= j <= n, s = cos(pi/(n+1)), and c_j = v_j v_{j+1}/s.

[Proved, Theorem Q] For every quadratic h in the form (4.1) and every n >= 2,

    lambda_min( sum_{j=1}^{n-1} c_j h_{j,j+1} ) >= e_0(h) - (1 - s) ||h_cross||
                                               >= e_0(h) - [1 - cos(pi/(n+1))] M(h)
                                               >= e_0(h) - pi^2 M(h) / (2 (n+1)^2).

The window is the C135 sine window, independent of h, with no one-site terms (u = 0) in the representative (4.1). So (W) holds on the whole quadratic class with K = pi^2/2, by a fixed envelope, uniformly in gaps, fields and crossover lengths, and in particular uniformly along TFIM at g = 1 + c/n for all c.

Proof. Step 1 (one-site weights). The bond-induced weight of the one-site Z term at site k is (c_{k-1} + c_k)/2 (with c_0 = c_n = 0) times z. Since v is the lowest Dirichlet eigenvector, v_{k-1} + v_{k+1} = 2 s v_k for 1 <= k <= n (with v_0 = v_{n+1} = 0), by sin(theta(k-1)) + sin(theta(k+1)) = 2 cos(theta) sin(theta k). Therefore (c_{k-1} + c_k)/2 = v_k (v_{k-1} + v_{k+1})/(2s) = v_k^2, and

    sum_j c_j h_{j,j+1} = W' + ((1 - s)/s) sum_j v_j v_{j+1} h_cross,(j,j+1),
    W' := a 1 + sum_j v_j v_{j+1} h_cross,(j,j+1) + z sum_k v_k^2 Z_k.                   (4.2)

Step 2 (Lemma Q1: W' >= e_0(h) 1). Use Majorana operators a_j = (prod_{l<j} Z_l) X_j, b_j = (prod_{l<j} Z_l) Y_j on the open chain, gamma_j = (a_j, b_j)^T. They are Hermitian, square to one and anticommute. Direct multiplication gives

    Z_j = -i a_j b_j,  X_j X_{j+1} = -i b_j a_{j+1},  Y_j Y_{j+1} = i a_j b_{j+1},
    X_j Y_{j+1} = -i b_j b_{j+1},  Y_j X_{j+1} = i a_j a_{j+1}.

(For instance b_j a_{j+1} = Y_j Z_j X_{j+1} = i X_j X_{j+1}, the strings cancelling except Z_j.) Hence h_cross,(j,j+1) = i gamma_j^T B gamma_{j+1} with the real 2 x 2 matrix B, B_ba = -J_xx, B_ab = J_yy, B_bb = -J_xy, B_aa = J_yx, and z Z_k = (i/2) z gamma_k^T C gamma_k with C = [[0, -1], [1, 0]]. Define, for k in [-pi, pi), the two-component operator Gamma(k) = sum_j v_j e^{ikj} gamma_j and the Hermitian 2 x 2 symbol

    R(k) = (i/4) Ahat(k),   Ahat(k) = 2 z C + 2 B e^{-ik} - 2 B^T e^{ik}.

Expanding the Fourier integral, int dk/(2 pi) Gamma(k)* R(k) Gamma(k) = sum_{j,l} v_j v_l gamma_j^T R_{j-l} gamma_l with R_0 = (i/2) z C, R_{-1} = (i/2) B, R_1 = -(i/2) B^T. The on-site terms give z sum_k v_k^2 (i/2) gamma_k^T C gamma_k. The two off-diagonal terms combine, by anticommutation of Majoranas on different sites, into sum_j v_j v_{j+1} i gamma_j^T B gamma_{j+1}. Hence

    W' - a 1 = int_{-pi}^{pi} dk/(2 pi) Gamma(k)* R(k) Gamma(k).

Let R(k) = sum_{s=+,-} lambda_s(k) e_s e_s*. Then Gamma* R Gamma = sum_s lambda_s D_s* D_s with D_s = e_s* Gamma(k), and

    {D_s, D_s*} = sum_{p,q} conj(e_{s,p}) e_{s,q} sum_{j,l} v_j v_l e^{ik(j-l)} {gamma_{j,p}, gamma_{l,q}} = 2 sum_p |e_{s,p}|^2 sum_j v_j^2 = 2.

So 0 <= D_s* D_s <= 2 and Gamma* R Gamma >= 2 sum_s min(lambda_s(k), 0). If mu_1(k), mu_2(k) are the eigenvalues of i Ahat(k), then lambda_s = mu_s/4 and 2 min(mu/4, 0) = (mu - |mu|)/4. Moreover Tr i Ahat(k) = 4 Tr(B) sin k integrates to zero. Therefore

    W' >= a - (1/4) int dk/(2 pi) (|mu_1(k)| + |mu_2(k)|) =: e_*(h).

Finally e_*(h) >= e_0(h), because e_* is the energy density of a TI state on Z: the quasi-free state whose Majorana symbol is the negative spectral projection of i Ahat(k). Its existence as a TI spin state on Z, with even moments given by the interval Jordan-Wigner map and odd moments zero, is C135 PROOF section 5 (existence detail). Its energy density is computed from the same symbol, giving e_*. (In fact e_* = e_0, the standard free-fermion formula, but only e_* >= e_0 is used.)

Step 3. Since v_j > 0, sum_j v_j v_{j+1} h_cross,(j,j+1) >= -(sum_j v_j v_{j+1}) ||h_cross|| = -s ||h_cross||. With (4.2) and Step 2, lambda_min >= e_0 - (1 - s)||h_cross||.

Step 4 (||h_cross|| <= M). Let h~ = h + b + o (x) 1 - 1 (x) o be a centered optimal representative, ||h~|| = M. Total parity Z (x) Z commutes with h; averaging h~ with its conjugate by Z (x) Z gives another representative, of norm at most M by convexity, whose gauge o is even (in span{1, Z}). For this representative, the part odd under conjugation by Z (x) 1 is [h~ - (Z (x) 1) h~ (Z (x) 1)]/2 = h_cross, since the scalar, the gauge and the Z terms are even and the four cross terms are odd. Its norm is at most (||h~|| + ||h~||)/2 = M.

Step 5. 1 - cos x <= x^2/2. []

[Proved, checks] (i) XX: h_cross = (XX + YY)/4, ||h_cross|| = 1/2 = M, e_0 = -1/pi. Theorem Q gives deficit at most (1 - s)/2, while the sharper C135 bound uses lambda_min(sum v_j v_{j+1} h) >= -1/pi to get (1/s - 1)/pi; both are of order pi^2 M/n^2 and they agree in structure, because Lemma Q1 for XX is C135's positivity argument. (ii) n = 2: v = (1/sqrt 2, 1/sqrt 2), s = 1/2, c_1 = 1, and the claim is lambda_min(h) >= e_0 - ||h_cross||/2, which is checked numerically in check_quadratic.py. (iii) Pure on-site bond: h_cross = 0 and the window has zero deficit. (iv) Numerical check (check_quadratic.py, 40 random quadratic bonds including critical and near-critical TFIM, XY in a field and DM, n from 2 to 400, Majorana method cross-checked against exact diagonalization for n <= 10): max(e_0 - W') = -3.4e-6, max violation of Theorem Q = -2.1e-5 (no violation), max(||h_cross|| - M) = 6.7e-16 (the inequality of Step 4 is attained, for example at XX), and max n^2 (e_0 - lambda_min)/M = 1.67, below pi^2/2.

[Proved, Remark 4.1, what Theorem Q adds to C135] C135 proved the LTI rate for the quadratic class from the primal side (covariance smoothing) and gave an explicit certificate only for XX. Theorem Q is the dual statement on the whole class: an explicit telescoping certificate X_n = sum_j (sum_{i<=j} c_i - 1) h_{j,j+1}, built from translated bonds only, with deficit (1 - s) M. The sine window's one-site weights are automatically v_k^2, which is why no one-site counterterm is needed in the symmetric representative.

## 5. Lesson 2 for windows

[Proved, Proposition 5.1] On the quadratic class the fixed sine envelope is uniform in the mass. Theorem Q bounds the deficit by (1 - s)||h_cross|| for every quadratic bond, including TFIM at g = 1 + c/n for every real c and n, so n^2 (e_0 - W_n)/M <= n^2 [1 - cos(pi/(n+1))] < pi^2/2 along every near-critical sequence. No h-dependence of the weights is needed, and the only h-dependence of the certificate is through the bond itself (and the choice of the symmetric representative (4.1)).

[Proved, mechanism in the free case] In (4.2) the window is W' plus a correction. W' is a superposition of quadratic forms in normalized windowed modes, int dk/(2 pi) Gamma(k)* R(k) Gamma(k), and each such form is bounded below by the negative part of the infinite-volume symbol R(k) at the same momentum, whatever the gap. A massive mode is therefore never penalized in W'. The whole deficit comes from the mismatch between the bond normalization sum_j v_j v_{j+1} = s and the site normalization sum_k v_k^2 = 1, which multiplies the hopping (kinetic) part by s. Its cost is at most (1 - s) times the kinetic scale ||h_cross|| <= M, and it is of order pi^2 M/n^2 whether or not the mode is massive. Lesson 2 says that a local energy functional responds to massive modes at the same order as to the conformal ones, which is what this mechanism shows. It does not say that the response can exceed the kinetic scale, and in the free class it cannot. Lesson 2 is an obstruction to isolating A from a window, not to bounding the window deficit.

[Heuristic, gradient expansion and consistency with C105] For a smooth window c_j = (1/n) f(j/n), a local-density expansion gives e_0 - lambda_min ~ kappa sum_j (c_{j+1} - c_j)^2/c_j, with kappa a static response coefficient. For the sine-square window sum_j (c_{j+1} - c_j)^2/c_j -> 4 pi^2/n^2. In a CFT the SSD deficit is the ring Casimir term pi A/(6 n^2) (Proposition 3.1 with SSD exactness), so kappa_CFT = A/(24 pi), which is the smooth value 1/(24 pi) found for XX in C105 (filtered profile). Lesson 2 is the statement that massive modes add to kappa terms not proportional to their c v. The numerics of RESULT.md section 3 show that for the free TFIM these terms are smaller than the critical one and vanish as the correlation length becomes small compared to n (n^2 delta/M decreases roughly like xi/n at fixed mass).

## 6. Proof attempts for interacting d = 2 bonds, and where each stops

[Literature, read abstract] Gosset and Mozgunov (arXiv:1512.00088) improve Knabe's local gap threshold for frustration-free translation-invariant chains from 1/(n - 1) to 6/(n(n+1)), by a weighted (quadratic-in-position) sum of squares of subchain Hamiltonians. [Proved, obstruction] Their argument needs (i) frustration-freeness, so that every subchain Hamiltonian has ground energy exactly 0 and its square is controlled by its gap, and (ii) a local gap as input; its output is a gap, not a ground energy. For ground energies of frustration-free bonds the window statement is trivial with K = 0 (Proposition 2.1), so the method has nothing to add there. For a frustrated bond, shifting h to h - e_0 destroys positivity of the terms, and the cross terms of the squares are no longer controlled. A positive local decomposition of H_open - e_0 (n - 1) with O(M/n) boundary error would be needed to restart the argument, and that is (T) itself. So the Gosset–Mozgunov mechanism does not transfer.

[Proved, obstruction] Many-body IMS localization. For a partition of unity sum_alpha f_alpha(x)^2 = 1 on the ring, H_L = sum_alpha sum_x f_alpha(x)^2 h_{x,x+1} holds exactly, so localizing the Hamiltonian density costs nothing and gives back E_L >= sum_alpha lambda_min(window_alpha), which is the window route itself (Remark 1.7). The IMS error term arises only when one localizes a state or a one-body kinetic operator, and here all the content sits in the lower bound for one window. IMS therefore reduces to (W) and supplies no estimate for it.

[Proved, obstruction] Deformation from the free certificate. Let h_t = h_0 + t V with h_0 quadratic. With a fixed window of nonnegative weights, |lambda_min(W_{h_t}) - lambda_min(W_{h_0})| <= |t| ||V||, and |e_0(h_t) - e_0(h_0)| <= |t| ||V||. So the deficit at level n grows by at most 2|t| ||V||, and (W) holds on a ball of radius of order M/(n^2 ||V||) around the quadratic class. This radius shrinks with n, so it gives no uniform statement. The first-order variation of the deficit is t [omega_0(V) - sum_j c_j psi_0(V_j)], the difference between the infinite-volume and the window ground-state expectations of V, which is O(1/n^2) only when V is irrelevant or marginal at h_0. For a relevant V (for example the field direction at critical TFIM, or any direction leaving the critical surface) the deficit is a scaling function n^{-2} F(n^{1/nu} t) and perturbation theory in t converges only for t of order n^{-1/nu}. The free case escapes this because Theorem Q controls F exactly. In the interacting case, uniform control along a path needs a nonperturbative bound on the scaling function, which is what a constructive-RG argument near free points (C130) would have to supply. No such argument is available here.

[Numerical, Proposition-level conclusion] h-dependent weights. The optimization of Task 2 finds that h-dependence of the bond weights is not what makes windows work at d = 2. For every bond with a symmetry fixing the one-site terms (XX, Heisenberg, XXZ at all Delta tested except the gapped Delta = 2, chiral), the optimal c is the sine-square window at even n and the C135 sine window at odd n, with relative sup deviation below 5e-3, for all n computed; for XYZ near criticality the better fixed window is within 3e-3 in K of the optimum. For TFIM, including g = 1 + c/n for c in {-2, 0.5, 5}, the optimal c is the sine-square window at all n, and the only h-dependence is in the one-site terms u. On the seven random bonds completed (RESULT.md section 2) the optimal shape can differ from both fixed windows, but the better fixed window with optimized u is within 2e-3 in K of the optimum. So the weights do not track the scale at which massive modes appear; a fixed envelope with h-dependent one-site terms suffices on every tested bond.

## 7. What is proved, and the remaining lemma

[Proved] (W) holds with K = pi^2/2 on the quadratic class, by an h-independent window (Theorem Q), with K = 0 on gauge frustration-free bonds (Proposition 2.1, including all classical bonds and ferromagnetic XXZ with Delta <= -1), and with the same constants on all twists of these (Lemma 1.5(d)). The window class is a restriction of the LTI dual (Proposition 1.3), so on these classes the LTI rate, (K) and A <= (6/pi) K M follow. On the quadratic class this reproduces C135's LTI theorem from the dual side, with an explicit certificate.

[Assumed, the single remaining lemma] Fixed-envelope window lemma at d = 2. Let a_j = v_j v_{j+1}/s be the C135 sine window. There is K such that for every bond h on C^2 (x) C^2 and every n >= 2 there are one-site terms u_k with sum_k u_k = 0 such that

    lambda_min( sum_{j=1}^{n-1} a_j h_{j,j+1} + sum_{k=1}^n u_k ) >= e_0(h) - K M(h)/n^2.

Equivalently (Proposition 1.4 with c fixed), for every state rho on n qubits whose one-site marginals are all equal, the sine-weighted average sum_j a_j rho_{j,j+1} of its bond marginals lies within K/n^2 of Omega_inf in the dual-M norm. This is a kinematic statement, with no Hamiltonian and no h-dependence. Theorem Q proves it, with K = pi^2/2, in every quadratic direction h (for states with equal one-site marginals the gauge, and so the choice of representative (4.1), drops out). The missing directions are the non-quadratic ones: Z Z after rotation to the Jordan-Wigner frame, DM components transverse to the hopping plane, and transverse fields. The same lemma with the sine-square window in place of a_j is numerically equivalent (RESULT.md), but Theorem Q is proved only for a_j. Numerically the lemma holds on every bond tested, with K_sine_u at most about 1.95 (Heisenberg, n = 4) and decreasing towards the ring value with n on the main bonds.
