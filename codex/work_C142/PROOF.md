# C142 proof: the two U(1) strips

2026-11-07. Run by a Claude Opus agent, not by Codex. J=1 in all analytic sections (restore J by homogeneity). Tags: [Proved] is a mathematical deduction, [Proved, cited] uses a result of an earlier run that has had a referee pass (named), [Numerical] is floating point without enclosures, [Heuristic] has an unproved step, [Assumed] is a hypothesis.

Part I closes the repulsive strip 63/64 < Delta <= 1 and redoes C139 Section 8.2 line by line. Part II proves a version of (D) and (D1) near Delta = -1 that is weaker by a factor log(1/eta), and isolates what remains. Part III is the assembly.

# Part I. The repulsive side, 0 <= Delta <= 1

Notation of C139 Section 0 and C140: Delta = cos gamma, 0 < gamma <= pi/2, b = pi/gamma, a = pi sin gamma, ab = pi^2 sin(gamma)/gamma <= pi^2, q = 1 - 2gamma/pi, q_J = (pi-2gamma)/(2(pi-gamma)) < 1/2, r0(x) = 1/(2 gamma cosh(bx)), K the positive kernel, J = (I+K*)^(-1)K* >= 0 the hole kernel (even, unimodal, ||J||_1 = q_J), s_Q = r0 + J*(1_out s_Q), out = {|x| > Q}, h_Q = int_out s_Q. Ground-state roots x_1 < ... < x_n, mu = (1/L) sum delta_(x_j), sigma = sigma_L = a1 - K*mu, counting function z_L with z_L' = sigma, z_L(x_j) = I_j/L, I_j = j - (n+1)/2, nu = n/L, 1 <= n <= L/2.

## I.1 Where C140 lost the factor

[Proved] Lemma I.1 (sharp tail integral). For every Q >= 0,

    int_Q^inf r0(x) dx <= (pi/(2b)) r0(Q),

with equality at Q = 0, and the ratio of the two sides tends to 2/pi as Q -> infinity.

Proof. With t = e^(-bQ) in (0,1], int_Q^inf sech(bx) dx = (2/b) atan(t), so int_Q^inf r0 = atan(t)/(gamma b), while r0(Q) = t/(gamma(1+t^2)). The ratio is f(t)/b with f(t) = atan(t)(1+t^2)/t. Then f'(t) = [t - (1-t^2) atan(t)]/t^2 >= 0 on (0,1], since atan(t) <= t and 1-t^2 <= 1. So f is increasing, f(1) = pi/2, f(0+) = 1. End.

C140 Lemma 3.3 used atan(t) <= t <= 2t/(1+t^2), i.e. the constant 2 in place of f(t) <= pi/2. This factor 4/pi is the whole 6 percent shortfall.

[Proved, under (P): sigma >= 0 on [x_1, x_n]] Lemma I.2 (location of the outermost root).

    L r0(x_n) >= b(1-q_J)/pi = b/(2(pi-gamma)) >= b/(2pi) = 0.15915 b.

Proof. C140 Lemma 3.3 verbatim up to its last line (refereed in C140 REFEREE Section 1): by (F3), Corollary 3.2 and C139 Lemma 3.3, 1/(2L) <= int_(x_n)^inf sigma <= int_Q^inf s_Q = h_Q/2 <= int_Q^inf r0/(1-q_J), with Q = x_n. Lemma I.1 gives 1/(2L) <= (pi/(2b)) r0(Q)/(1-q_J). Finally 1-q_J = pi/(2(pi-gamma)). End.

The needed value is lambda_* b = ((2-ln2)/pi^2) b = 0.13241 b. Since 2pi(2-ln2) = 8.211 < pi^2 = 9.870, the bound exceeds it at every gamma, by a factor at least pi^2/(2pi(2-ln2)) = 1.2020 (at gamma -> 0). So (L*) holds, strictly, along the branch.

[Proved, under (P)] Theorem I.3 (uniform margin). For x in [x_1, x_n],

    sigma(x) >= theta_1(gamma) r0(x),   theta_1(gamma) = 1 - 2(2-ln2)(pi-gamma)/pi^2 >= 1 - 2(2-ln2)/pi = 0.16803.

Proof. C140 Theorem 3.4 gives sigma(x) >= r0(x) - (2-ln2)/(pi gamma L). Since r0(x) >= r0(x_n) on [x_1, x_n] (symmetry and monotonicity of r0), Lemma I.2 gives (2-ln2)/(pi gamma L) <= (2-ln2) r0(x) 2(pi-gamma)/(pi gamma b) = r0(x) 2(2-ln2)(pi-gamma)/pi^2, using pi gamma b = pi^2. End.

[Numerical, check_repulsive.json, 160 cases, Delta in {0.5, 0.9, 0.99, 0.999}, L from 4 to 256, seven fillings per L] L r0(x_n)/b >= 0.442 against the proved 0.159 to 0.239, and min sigma/r0 on [x_1, x_n] >= 0.974 against theta_1 between 0.168 and 0.45.

## I.2 The continuation argument on all of 0 < gamma <= pi/2

[Proved] Theorem I.4. For every 0 < gamma <= pi/2 (0 <= Delta < 1), every L >= 2 and 1 <= n <= L/2, the branch of real, strictly ordered critical points of the Yang-Yang functional continued from the free solution at gamma = pi/2 exists on [gamma, pi/2], is real-analytic and reflection symmetric, its Hessian is positive definite, and

    (i)  sigma_L >= theta_1 r0 on [x_1, x_n],
    (ii) int_(x > x_n) (sigma_L)_- dx <= (pi-2gamma)/(2pi L), and the same on x < -x_n.

Its Bethe energy is the sector ground energy F_L(n) at flux pi per bond, at both parities of L.

Proof. n = 1 is C140 Theorem 4.1 (x_1 = 0, sigma > 0 on R). Let n >= 2 and fix gamma_0 in (0, pi/2]. At gamma = pi/2, K = 0, sigma = a1 > 0, and the Hessian Y'' = 2pi L diag(sigma(x_j)) + Theta (C139 Lemma 2.2, Theta positive semidefinite by Bochner) is positive definite. By the implicit function theorem there is a unique real-analytic branch x(gamma) of critical points near pi/2. Let gamma_* be the infimum of the gamma_1 in [gamma_0, pi/2] such that the branch extends real-analytically, through critical points with positive definite Hessian, to (gamma_1, pi/2]. Suppose gamma_* >= gamma_0.

Step 1 ((P) along the branch). On the connected set (gamma_*, pi/2] the branch is continuous, strictly ordered, and reflection symmetric (the reflected family solves the same equations, agrees with x(gamma) at pi/2, and the coincidence set is open by local uniqueness and closed by continuity). The set where (P) holds is closed in (gamma_*, pi/2], since sigma_L(x; gamma) is jointly continuous and [x_1, x_n] moves continuously. Where (P) holds, Theorem I.3 gives sigma >= theta_1 r0 > 0 on the compact [x_1, x_n], and by continuity (P) holds strictly nearby, so the set is open. It contains pi/2, hence it is all of (gamma_*, pi/2], and (i) holds there.

Step 2 (limit point). The roots stay in the a priori ball of C139 Section 2.1, whose radius is bounded for gamma in [gamma_0, pi/2]. Take gamma_k decreasing to gamma_* with x(gamma_k) -> x_*. Then x_* is a critical point at gamma_* (grad Y is jointly continuous), with distinct coordinates (C139 Section 2.1), ordered and symmetric. By (i) and Lemma I.2, sigma(x_j(gamma_k)) >= theta_1 r0(x_n(gamma_k)) >= theta_1 b/(2(pi-gamma_k)L), which passes to the limit, so the Hessian at x_* is positive definite.

Step 3 (no change of branch). The implicit function theorem at (x_*, gamma_*) gives a neighbourhood U of x_*, an eps > 0, and a real-analytic x~ on (gamma_* - eps, gamma_* + eps) such that x~(gamma) is the only critical point in U at each such gamma. The set I = {gamma in (gamma_*, gamma_* + eps) : x(gamma) = x~(gamma)} contains gamma_k for k large (x(gamma_k) is in U), is closed (continuity), and is open (if x(gamma) = x~(gamma) is in U, then x(gamma') is in U for gamma' near gamma, hence equals x~(gamma')). So I = (gamma_*, gamma_* + eps), and x~ continues the branch analytically below gamma_*, with positive definite Hessian near gamma_*. This contradicts the definition of gamma_* unless gamma_* < gamma_0. Hence the branch exists on [gamma_0, pi/2], with (i).

(ii) From C140 Corollary 3.2, (sigma)_- <= phit/L on out, phit = (I-J_out)^(-1) psi. psi and J_out commute with the reflection, so phit is even and int_(y > x_n) phit = ||phit||_1/2 <= ||psi||_1/(2(1-q_J)) = q_J/(2(1-q_J)) = (pi-2gamma)/(2pi). (C140 stated the two-sided total 1/L.)

Identification. The facts used, all sign-free (C139 Section 1 item 2, C134 Section 2, C140 REFEREE Section 2): the sector Hamiltonian at flux pi per bond (total L pi, at either parity) has off-diagonal entries -1/2 in the S^z basis, so it is stoquastic and irreducible and F_L(n)(gamma) is a simple, real-analytic eigenvalue; along the analytic branch the coordinate Bethe vector is an analytic eigenvector with eigenvalue E(gamma) = Delta L/4 + sum e(x_j) wherever it is nonzero; at gamma = pi/2 it is the free-fermion Slater determinant of the symmetric Fermi sea, nonzero and the sector ground state. So E = F_L(n) near pi/2, and E and F_L(n) are both real-analytic on [gamma_0, pi/2], hence equal there by the identity theorem. End.

[Numerical] The integrated negative mass beyond x_n is zero to grid accuracy in all 160 cases (sigma_L > 0 on R), consistent with C140's sigma_L >= 0.952 r0 on R.

## I.3 Theorem C of C139, redone line by line with the proved sawtooth

Fix 0 < gamma <= pi/2, L >= 2, 1 <= n <= L/2, and the roots of Theorem I.4. Let Q = Q_L be the first x >= x_n with z_L(x) = nu/2 (Q_L > x_n because z_L(x_n) = nu/2 - 1/(2L); Q_L is finite for nu < 1/2 because z_L(inf) = c > nu/2; at even half filling Q_L = infinity is allowed). r = r_Q, Z = Z_Q, eps = eps_Q, m = m(Q), h = h(Q) are the restricted continuum objects of C139 Section 3, R = mu - 1_[-Q,Q] sigma dx, Cum R(x) = R((-inf, x]) = C_mu(x) - z_L(x).

[Proved, cited] Continuum inputs, all uniform in Q and gamma (C139 Lemmas 3.3 to 3.5, refereed in C140 REFEREE Sections 3 and 4): r0 <= r <= c_r r0 with c_r = 1 + 4/pi; |r'| <= b r; 1/(1+q) <= Z <= 1; Var_[-Q,Q](Z) <= q_J; ||(I+K_Q)^(-1)||_(inf->inf) <= (1+q_J)/(1-q_J) <= 3; |eps'| <= ab r; v = eps'/r has |v| <= ab and Var_[-Q,Q](v) <= V_v = (2+40/pi) ab; 0 <= e_bulk'' <= a2 = ab(1+q)^2/2 <= 2pi^2; eps(+-Q) = 0. At Q = infinity the same hold with r = r0, Z = 1/(1+q), v = ab tanh(bx).

[Proved, cited] Exact identity (C139 (5.1), sign-free, rederived in C140 REFEREE): with nu' = mu - 1_Q r dx, (I + K_Q) nu' = R on [-Q,Q] as measures, and

    F_L(n) - L[E(Q) + h(nu - m)] = L int eps dR = -L int_(-Q)^Q eps'(x) Cum R(x) dx.        (I.1)

The second form is integration by parts using eps(+-Q) = 0 and Cum R(+-Q) = 0.

### Lemma I.5 (sawtooth) [Proved]

    sup_[-Q,Q] |Cum R| <= s/L,   s = 1 - gamma/pi in [1/2, 1).

Proof. On (x_j, x_(j+1)), C_mu = j/L - nu/2 and z_L increases (Theorem I.4(i)) from I_j/L to I_j/L + 1/L, so Cum R decreases from j/L - nu/2 - I_j/L = 1/(2L) to -1/(2L). On (x_n, Q], Cum R = nu/2 - z_L(x) >= 0 by the definition of Q, and z_L(x) >= z_L(x_n) - int_(x_n)^x (sigma)_- >= nu/2 - 1/(2L) - (pi-2gamma)/(2pi L) by Theorem I.4(ii), so Cum R <= (1/2 + (pi-2gamma)/(2pi))/L = s/L. The left side is the mirror image. At Q = infinity the same argument applies on (x_n, inf). End.

[Numerical] L sup|Cum R| is between 0.4997 and 0.5000 in all 160 cases (no negative mass occurs numerically), against s in [0.67, 0.9997].

### Lemma I.6 (mismatch and cumulative error) [Proved]

    |nu - m| <= q_J s/L,
    ||G||_inf <= g/L,   g = (2+q) [(1+q_J)/(1-q_J)] s (1 + q q_J/2) <= 45/4,

where G(x) = nu'((-inf, x]) = C_mu(x) + nu/2 - int_(-Q)^x r.

Proof. Mismatch: nu - m = int Z dR (apply the self-adjoint (I+K_Q) to Z = (I+K_Q)^(-1)1), and int Z dR = -int Z' Cum R dx, so |nu - m| <= Var(Z) sup|Cum R|.

Cumulative error. Let k(y) = int_0^y K (odd, |k| <= q/2). Integrate (I+K_Q)nu' = R over [-Q, x]. Integrating by parts in u, with G(-Q) = 0 and G(Q) = nu - m,

    int_(-Q)^x (K_Q nu')(y) dy = int [k(x-u) - k(-Q-u)] dnu'(u)
                               = (K_Q G)(x) - (K_Q G)(-Q) + (nu - m)[k(x-Q) + k(2Q)].

So G + K_Q G = S + c_0 - (nu - m) psi on [-Q,Q], with S = Cum R, c_0 = (K_Q G)(-Q), psi(x) = k(x-Q) + k(2Q) in [0, k(2Q)] subset [0, q/2]. Hence G = g_1 + c_0 Z with g_1 = (I+K_Q)^(-1)[S - (nu-m)psi], and G(-Q) = 0 gives c_0 = -g_1(-Q)/Z(-Q), |c_0| Z(x) <= ||g_1|| (1+q). So ||G|| <= (2+q)||g_1|| and ||g_1|| <= [(1+q_J)/(1-q_J)] (s + q_J s q/2)/L. The value 45/4 = 3 * 3 * 5/4 is the gamma -> 0 limit, and each factor is monotone in gamma (constants_rep.py evaluates g(gamma); sup = 11.24999). End.

This replaces C139 8.2's ||G|| <= 27/(4L) (sawtooth 1/2, with the cruder psi <= q) and C140's unchecked rescaling.

[Numerical] L||G|| is between 0.500 and 0.532, L|nu - m| <= 0.057 (proved q_J s up to 0.49).

### Lemma I.7 (density) [Proved]

On [-Q,Q], sigma - r = -K*nu' = -int G(y) K'(x-y) dy (G = 0 below -Q and G = nu - m above Q, and K vanishes at infinity), so

    |sigma - r| <= ||K'||_1 ||G|| <= D/L,      D := 2g/(pi gamma)  (>= 2K(0) g = 2 cot(gamma) g/pi),
    |sigma' - r'| <= ||K''||_1 ||G|| <= D1/L,   D1 := 4 M_1 g,

with ||K'||_1 = 2K(0) = 2cot(gamma)/pi and ||K''||_1 = 4 max|K'| <= 4M_1, M_1 = (2/pi)(0.325/sin^2 gamma + 0.25/sin gamma) (C139 Lemma 5.0). D is chosen as an upper bound that does not vanish at gamma = pi/2, so that the ratio D1/D stays bounded there.

[Numerical] L sup|sigma - r| <= 0.92 against D up to 152 at Delta = 0.999, and L sup|sigma' - r'| <= 35 against D1 up to 4576.

### Lemma I.8 (quadrature) [Proved]

Assume L >= 4 gamma D = 8g/pi (so L >= 29 suffices at every gamma). Let A >= 0 solve r0(A) = 2D/L and call [-A,A] the core.

Core. On [-A,A] intersected with [-Q,Q], r >= r0 >= 2D/L, so sigma >= r - D/L >= r/2 > 0. Let x_(j1), ..., x_(j2) be the roots there (consecutive, z_L(x_j) = I_j/L) and f = eps o z_L^(-1). The composite trapezoid rule with step 1/L and the Peano kernel (t-a)(b-t)/2 <= 1/(8L^2) give

    |sum_(j1)^(j2) eps(x_j) - [eps(x_j1) + eps(x_j2)]/2 - L int_(x_j1)^(x_j2) eps sigma dx| <= V_f/(8L),

where V_f bounds the variation of f' = eps'/sigma = v g, g = r/sigma in (0, 2]. By the product rule for variation, Var(vg) <= ||v|| Var(g) + ||g|| Var(v) <= ab Var(g) + 2V_v. Writing r'sigma - r sigma' = r'(sigma - r) - r(sigma' - r') and using sigma^2 >= r^2/4, |r'| <= br and Lemma I.7, |g'| <= 4(bD + D1)/(L r). Since int_(-A)^A dx/r0 = (4gamma/b) sinh(bA) <= 2/(b r0(A)) = L/(bD),

    Var(g) <= 4 + 4 D1/(bD),   V_f = ab Var(g) + 2V_v <= 8ab + 80ab/pi + 4ab D1/(bD).

Here 4ab D1/(bD) = 4ab * 4M_1 g/((pi/gamma)(2g/(pi gamma))) = 8 ab gamma^2 M_1 = 8pi^2 gamma sin(gamma) M_1, and gamma sin(gamma) M_1 <= (2/pi)(0.325 gamma/sin gamma + 0.25 gamma) <= 0.575 (gamma/sin gamma <= pi/2 and gamma <= pi/2). So V_f <= 8pi^2 + 80pi + 45.4 <= 375.7, and the core contributes at most 46.97/L.

Ends. Integrating by parts on (x_j2, Q] with Cum R(x_j2^+) = 1/(2L) and eps(Q) = 0,

    sum_(j > j2) eps(x_j) - L int_(x_j2)^Q eps sigma dx + eps(x_j2)/2 = -L int_(x_j2)^Q eps' Cum R dx,

and similarly on the left (Cum R(x_j1^-) = -1/(2L)). By Lemma I.5 and |eps'| <= ab r, each end is at most s ab int_(x_j2)^Q r. If A < Q, there is no root in (x_j2, A], so int_(x_j2)^A r = G(x_j2^+) - G(A) <= 2g/L, and int_A^Q r <= c_r int_A^inf r0 <= c_r (pi/(2b)) r0(A) = c_r pi D/(bL) by Lemma I.1. If A >= Q, then j2 = n and int_(x_n)^Q r = G(x_n^+) - G(Q) <= 2g/L. With D/b = 2g/pi^2, both ends together are at most

    2 s ab g (2c_r/pi + 2)/L.

Empty core. If no root lies in [-A,A] (then A < Q), (I.1) gives |L int eps dR| <= s ab int_(-Q)^Q r <= s ab [2g + 2c_r pi D/b]/L, using int_(-A)^A r = G(-A) - G(A), which is below the ends bound. So no lower threshold on L beyond L >= 8g/pi is needed (C139 used L >= 77 to guarantee a nonempty core).

Curvature. |e_bulk(nu) - e_bulk(m) - h(m)(nu-m)| <= a2 (nu-m)^2/2 <= a2 (q_J s)^2/(2L^2). At Q = infinity there is no mismatch.

### Theorem I.9 (repulsive sector theorem, uniform through Heisenberg) [Proved]

For every 0 <= Delta < 1, L >= 2 and 0 <= n <= L/2,

    |F_L(n) - L e_bulk(Delta, n/L)| <= K_U/L,
    K_large(gamma) = V_f/8 + 2 s ab g (2c_r/pi + 2) + a2 (q_J s)^2/2  (for L >= 29),
    K_large(gamma) <= 46.96 + 765.50 + 2.47 <= 815 (analytic), sup over a 200001-point grid = 811.30,
    K_U = max(29^2, 815) = 841.

Proof. n = 0 is exact. For L >= 29 add the core, the two ends (or the empty-core bound) and the curvature, using F_L(n) - L e_bulk(nu) = [F_L - L(E(Q) + h(nu-m))] - L[e_bulk(nu) - e_bulk(m) - h(m)(nu - m)] with E(Q) = e_bulk(m), h = e_bulk'(m). For L < 29 use |F_L - L e_bulk| <= L (bond spectral diameter 1, e_bulk in the bond interval, as in C139 Section 5.4). The analytic majorant bounds each term separately: V_f/8 <= 375.7/8 = 46.96 (Lemma I.8), ends <= 2 pi^2 (45/4)(2c_r/pi + 2) = 765.50 (s <= 1, ab <= pi^2, g <= 45/4), curvature <= 2pi^2 (1/2)^2/2 = pi^2/4. constants_rep.py evaluates the exact gamma-dependent expression on 200001 points of (0, pi/2]; its supremum, 811.30, is approached as gamma -> 0 (where V_f/8 -> 43.42). Pointwise: K_large = 53.6 at Delta = 0, 376 at gamma = 0.5, 629 at Delta = 63/64, 811.3 at Delta -> 1. End.

[Numerical] L|F_L(n) - L e_bulk| <= 0.91 in all 160 cases, a factor 890 below K_large.

### Corollary I.10 (repulsive U(1) class including Heisenberg) [Proved]

For every bond of the U(1) spin-1/2 class with 0 <= J_z/J <= 1, every field, DM angle and L >= 2,

    L(E_(2L) - 2E_L) <= C_rep M,   C_rep = 5 * 815 + pi^2/4 = 4077.5 <= 4078.

(With the exact gamma-dependent K_large evaluated on the grid, 5 * 811.30 + pi^2/4 = 4058.96; the stated 4078 uses only the analytic majorant.)

Proof. Matched doubling (2L, 2n): for L >= 29, L[F_(2L)(2n) - 2F_L(n)] <= L[K/(2L) + 2K/L] = (5/2)K, i.e. 5K in units of M_0 = J/2. For L < 29 seam gluing gives at most 4L <= 112 in units of M_0 (C136 Section 9). The C129/C136 reductions (DM angle to flux, field cancellation in matched sectors, flux cost pi^2 M/4, M >= M_0 for |Delta| <= 1) are sign-free and refereed (C136 REFEREE, Assembly). At Delta = 1, for fixed L, E_L, E_(2L) and M are continuous in the bond, so the bound extends from Delta < 1. End.

This supersedes C139 Theorem A's 345600 on [0, 63/64] and C140's sketched 19935 + pi^2/4. The constant is larger than C139's conditional 2215 + pi^2/4 because the proved sawtooth is s <= 1 instead of the assumed 1/2. In C139 8.2 the factor 27 in the end integral was 4||G||L; the direct bound above uses 2||G||L.

### I.11 Which variables are uniform at Delta -> 1

[Proved, bookkeeping] All statements are in the rapidity x of C139, where the kernel width is sin gamma and the root density scale is b = pi/gamma. Every constant above is expressed through gamma-free combinations: L r0(x_n)/b, D/b, gamma D, ab, gamma sin(gamma) M_1, q, q_J, s. The rescaled rapidity lambda = x/gamma (XXX rapidities at gamma -> 0) makes r0 and K converge, but no step needs that limit: Delta = 1 itself is reached by finite-matrix continuity at fixed L.

### I.12 Consistency checks

* Delta = 0 [Proved]: K = J = 0, q = q_J = 0, s = 1/2, g = 2 * 1 * (1/2) = 1, sigma = r0 = a1 and Cum R is the exact sawtooth 1/(2L). K_large(pi/2) = 53.6, consistent with (coarser than) the exact free constant 8 - 4 sqrt 2.
* gamma -> 0 [Proved]: theta_1 -> 1 - 2(2-ln2)/pi = 0.168 > 0, L r0(x_n)/b >= 1/(2pi); no constant diverges.
* L = 2 [Numerical, ed_check.json]: the end-to-end ED defect L(F_(2L)(2n) - 2F_L(n))/M_0 is 3.987, 3.9987, 4.000 at Delta = 0.99, 0.999, 1 (L <= 6), below 4078 and equal to the known supremum 4 at Heisenberg.
* ED [Numerical]: Bethe energies equal ED sector energies at flux pi per bond to 2.3e-14 for all 1 <= n <= L/2, L = 2 to 12, Delta = 0.99, 0.999.
* Casimir and Delta = 1/2: unchanged from C139 (no new energy identity here).

# Part II. The ferromagnetic side, -1 < Delta <= 0

Notation of C136 and C138: eta = pi - gamma in (0, pi/2], Delta = -cos eta, T = -theta'/(2pi) >= 0 even and decreasing on x > 0, q = ||T||_1 = 1 - 2eta/pi, B = 1/(1-q) = pi/(2eta), T(0) = cot(eta)/pi, sigma_L = a1 + T*mu > 0, z_L globally increasing, Q = Q_L with z_L(Q_L) = nu/2 (Q_L = infinity at even half filling), r = r_Q = (I - T_Q)^(-1) a1 on [-Q,Q], R_L = mu - 1_[-Q,Q] sigma_L dx.

## II.1 Two exact facts

[Proved] Lemma II.1 (oscillation bound). For every f of bounded variation on [-Q,Q] (on R if Q = infinity),

    |int f dR_L| <= (1/L) Var_[-Q,Q](f).

Proof. z_L is increasing (sigma_L > 0) and the cells [I_j/L - 1/(2L), I_j/L + 1/(2L)] tile [-nu/2, nu/2] = z_L([-Q,Q]). Their x-images X_j tile [-Q,Q], with int_(X_j) sigma_L = 1/L and x_j in X_j. So int f dR_L = sum_j int_(X_j) [f(x_j) - f(x)] sigma_L(x) dx, and the j-th term is at most (1/L) osc_(X_j) f. Sum over j. End. (This is the mechanism of C138 Theorem A, there applied to f = eps_Q.)

[Proved] Lemma II.2 (representation). Let calR_Q = sum_(k>=1) T_Q^k, a positive symmetric kernel on [-Q,Q]. Then for x in [-Q,Q]

    sigma_L(x) - r(x) = int calR_Q(x,y) dR_L(y),   sigma_L'(x) - r'(x) = int d_x calR_Q(x,y) dR_L(y).

Proof. On [-Q,Q], sigma_L - r = T*mu - T_Q r = T_Q(mu - r dx) = T_Q R_L + T_Q(sigma_L - r), since mu is supported in [-Q,Q] and mu = R_L + sigma_L dx there. So (I - T_Q)(sigma_L - r) = T_Q R_L, and (I - T_Q)^(-1) T_Q = calR_Q. The kernel is smooth for eta > 0, so the x-derivative passes under the integral. (C138 Section 6 states the first identity in operator form.) End.

So (D) is a statement about the midpoint-quadrature error of the single function y -> calR_Q(x,y), and Lemma II.1 bounds it by the variation of that function.

## II.2 The variation of the restricted dressing kernel, uniform in Q

Let calR(x) be the full-line kernel, Fourier transform That/(1-That) (That(w) = sinh((pi-2eta)w/2)/sinh(pi w/2), C138 (2.1)).

[Proved] Lemma II.3. calR is even and nonincreasing on x > 0, so Var_R(calR) = 2 calR(0), and

    calR(0) = (1/pi) int_0^inf That/(1-That) dw <= (log(1/eta) + pi/2)/(pi eta).

Proof. calR = sum_(k>=1) T^(*k). T is even and nonincreasing in |x|, and convolutions of even unimodal integrable functions are even unimodal (Wintner), so each term, and the sum, is even and nonincreasing on x > 0. Its transform is nonnegative and integrable (see below), so calR(0) = sup calR = (1/(2pi)) int That/(1-That). For the bound write 1 - That = 2cosh(gamma w/2) sinh(eta w/2)/sinh(pi w/2) (C138 (2.1)), so

    That/(1-That) = sinh(A)/(2 cosh(A + eta w/2) sinh(eta w/2)),   A = (pi-2eta)w/2 >= 0.

Since sinh(A) <= e^A/2 and cosh(A + c) >= e^(A+c)/2, sinh(A)/cosh(A + eta w/2) <= e^(-eta w/2), so That/(1-That) <= 1/(e^(eta w) - 1). Also That is decreasing in |w| and u/(1-u) is increasing, so That/(1-That) <= q/(1-q) = B - 1. Split at w = 1:

    pi calR(0) <= (B - 1) + int_1^inf dw/(e^(eta w) - 1) = pi/(2eta) - 1 + (1/eta) log(1/(1 - e^(-eta))).

Since (1 - e^(-eta))/eta >= e^(-eta/2) (equivalent to sinh(eta/2) >= eta/2), log(1/(1-e^(-eta))) <= log(1/eta) + eta/2, and pi calR(0) <= [pi/2 + log(1/eta) - eta/2]/eta. End.

[Numerical, check_attractive_<eps>.json] calR(0) by quadrature is 0.239, 7.150, 1762.6 at 1+Delta = 0.5, 1e-2, 1e-6, against the bound 0.463, 7.93, 1830.3. The inverse Fourier transform agrees with the quadrature value of calR(0) to 1e-8 and is nonincreasing on [0, 6] within its quadrature error.

[Proved] Lemma II.4 (restricted kernel). For every Q in (0, infinity] and x in [-Q,Q],

    Var_(y in [-Q,Q]) calR_Q(x, y) <= 2(1+q) calR(0).

Proof. Fix y in [-Q,Q]. The function u_y = calR_Q(., y) solves u = T(. - y) + T_Q u on [-Q,Q], and U_y = calR(. - y) solves U = T(. - y) + T*U on R, hence U = T(. - y) + T_Q U + T*(1_out U) on [-Q,Q]. So w = U - u solves w = T_Q w + T*(1_out U_y), that is

    calR(x - y) - calR_Q(x, y) = w(x, y) = int_out K_out(x, u) calR(u - y) du,   K_out(x,u) = [(I - T_Q)^(-1) T(. - u)](x) >= 0.

For fixed x, y -> w(x,y) is a nonnegative combination of translates of calR, so Var_y w(x, .) <= Var(calR) int_out K_out(x,u) du = 2 calR(0) int_out K_out(x,u) du. And int_out K_out(x,u) du = [(I - T_Q)^(-1)(T*1_out)](x), with T*1_out = q - T_Q 1 on [-Q,Q], so it equals (I - T_Q)^(-1)(q - 1 + (I - T_Q)1) = 1 - (1-q) Z_Q(x) <= q, since Z_Q >= 1 (C136 (2.2)). Finally Var_y calR(x - .) <= 2calR(0). End.

The bound is uniform in the Fermi boundary, including the dilute limit and Q = infinity. The analogous bound for d_x calR_Q uses the same splitting. With M_1^T = max|T'| <= (2cos eta/pi)(0.32476/sin^2 eta + 0.25/sin eta) (T has the form of C139's K with sin 2gamma replaced by sin 2eta, so C139 Lemma 5.0 applies, and T' is unimodal on (0,infinity), so ||T''||_1 = 4M_1^T):

[Proved] Lemma II.5. Var_y d_x calR_Q(x, .) <= 4M_1^T + 4(2+q) T(0) calR(0).

Proof. d_x calR_Q(x,y) = calR'(x-y) - d_x w(x,y). First, Var(calR') = ||calR''||_1, and calR'' = T'' + T'*calR' gives ||calR''||_1 <= ||T''||_1 + ||T'||_1 ||calR'||_1 = 4M_1^T + 2T(0) * 2calR(0). Second, d_x w(x,y) = int_out d_x K_out(x,u) calR(u-y) du, so Var_y d_x w <= 2calR(0) int_out |d_x K_out(x,u)| du. From K_out = T(. - u) + T_Q K_out, d_x K_out(x,u) = T'(x-u) + int_(-Q)^Q T'(x-v) K_out(v,u) dv, hence int_out |d_x K_out(x,u)| du <= [T(Q-x) + T(Q+x)] + ||T'||_1 sup_v int_out K_out(v,u) du <= 2T(0) + 2T(0) q. Add. End.

## II.3 (D) and (D1) up to a logarithm

[Proved] Theorem II.6. For every -1 < Delta <= 0, L >= 2, 0 <= n <= L/2 and x in [-Q_L, Q_L] (all of R at even half filling),

    eta L |sigma_L(x) - r_(Q_L)(x)|     <= C_D(eta)  := 2(1+q) eta calR(0) <= (4/pi)(log(1/eta) + pi/2),
    eta^2 L |sigma_L'(x) - r_(Q_L)'(x)| <= C_D1(eta) := eta^2 [4M_1^T + 4(2+q) T(0) calR(0)]
                                         <= (8/pi)(0.32476 (eta/sin eta)^2 + 0.25 eta^2/sin eta) + (12/pi^2)(log(1/eta) + pi/2).

Proof. Lemmas II.1, II.2, II.4, II.5, II.3, with T(0) eta = eta cot(eta)/pi <= 1/pi and 2 + q <= 3. End.

This holds at every point of the occupied interval, not only on the core, and at every filling. Compared with what was available: C136's density bound is of order 1/(eta^2 L) (C138 2.4 item 2; numerically eta times it is 10.8 to 1061 at 1 + Delta = 1e-2 to 1e-6), while Theorem II.6 gives eta L |sigma - r| <= C_D(eta) = 3.87, 5.48, 7.00, 9.97 at 1 + Delta = 1e-2, 1e-3, 1e-4, 1e-6.

[Numerical, check_attractive_<eps>.json, 95 cases, 1+Delta from 0.5 to 1e-6, L from 16 to 1024, fillings 2/L to 1/2] eta L sup|sigma_L - r_Q| <= 0.317 and eta^2 L sup|sigma_L' - r_Q'| <= 0.216 (away from the Fermi point by eta/4), against C_D(eta) from 0.67 to 9.97 and C_D1(eta) from 1.46 to 10.35. Table in RESULT.md.

## II.4 What this gives, and what is left

[Proved, conditional on (W1), (W2), (W3) of C138 Section 4.2] C138 Theorem B uses (D) and (D1) only on its core K_c, with the threshold M = 4C_D, and its proof is unchanged if C_D, C_D1 are replaced by C_D(eta), C_D1(eta) (Theorem II.6 holds everywhere on [-Q,Q], a fortiori on K_c). So, granted (W1)-(W3) with the constants of C138 (C_W = 8, kappa = 7/10, C_e = 13/10),

    |F_L(n) - L e_bulk(Delta, n/L)| <= K_*(eta)/L,
    K_*(eta) <= (2pi/(9C_D))[C_D + C_D1(1+1/kappa)] + pi C_W/12 + 16pi C_D(1+2/kappa) + 16pi/3 + pi^2 C_e/8
             = 193.9 C_D(eta) + O(1),

with C_D = C_D(eta), C_D1 = C_D1(eta). This grows like (4/pi) 193.9 log(1/eta) = 247 log(1/eta). I have not refereed C138 Theorem B line by line; the conditional statement inherits its status.

[Assumed, the remaining estimate] The strip -1 < Delta < -1 + delta_0 is therefore not closed. What is left, in order of difficulty:

(D) without the logarithm. The loss is identified exactly. In Lemma II.1 each cell is charged the full oscillation of calR_Q(x, .), and calR(x - y) behaves like (1/(pi eta)) log(1/|x-y|) for eta < |x - y| < 1 (the symbol That/(1-That) is about 1/(eta |w|) for 1 < |w| < 1/eta). The near cells, |x - y| < eta, carry variation O(1/eta) with no logarithm; the log comes from the cells at distance between eta and 1, where the true per-cell error is second order. [Heuristic] A second-order (trapezoid) bound on those cells removes the log provided sigma_L >= r_Q/2 on the unit neighbourhood of x and the derivative term (1/L^2) int |calR'(x-y)| |sigma_L'(y)|/sigma_L(y)^2 dy over that neighbourhood is O(1/(eta L)). Inserting Theorem II.6, both hold once eta L r >= C log(1/eta) on that neighbourhood (the second term is then of order (log(1/eta))^2/(eta L lambda^2) with lambda = eta L r). On the shell C_0 <= eta L r <= C log(1/eta) the a priori control is not strong enough, and a bootstrap from Theorem II.6 does not close, because the second-order bound at x needs control on [x-1, x+1], which reaches into sparser regions. [Numerical] The true margin: eta L sup|sigma - r| <= 0.32 everywhere on [-Q_L, Q_L] in all sampled cases (C138: L to 8192, 1 + Delta to 1e-6; here 0.317 in 95 cases), against the proved C_D(eta) = 9.97 at 1 + Delta = 1e-6.

(W1)-(W3). These are statements about the restricted continuum equation alone (C138 Section 4.2) and are not addressed here beyond the following. [Proved] The domination calR_Q <= calR together with the representation r_0 - r_Q = int_out K_out(., u) r_0(u) du gives (W2) only at distance of order log(1/eta) from the Fermi point, since calR has an exponential tail with prefactor of order 1/eta; a uniform (W2) at distance 1 needs the Wiener-Hopf suppression of K_out near the edge, which this positivity argument does not see. [Numerical, C138 5.3] Var(v_Q)/eta <= 8.0, e_bulk''/eta^2 <= 4/pi, min r_Q/r_0 at distance >= 1 from the edge >= 0.762.

Why the window was not closed by Theorem A plus Theorem II.6. C138 Theorem A gives L|F - L e_bulk| <= (5pi/4) eta L, and the conditional bound above gives K_*(eta) of order log(1/eta). Their minimum is still unbounded as eta -> 0 (both are of order log(1/eta) at eta L of order log(1/eta)). So even granting (W1)-(W3), no uniform constant follows; the log must go.

# Part III. Assembly on J_z/J in [-1, 1]

[Proved] Covered regions, every field, DM angle and L >= 2, constant for L(E_(2L) - 2E_L) <= C M:

* 0 <= J_z/J <= 1: C_rep = 4078 (Corollary I.10; C139 Lemmas 3.1-3.5 and identity (5.1), C140 Lemmas 3.1, 3.2 and Theorem 3.4, all refereed in C140 REFEREE; Lemmas I.1-I.8 here, not yet refereed).
* -1 + delta <= J_z/J <= 0, any 0 < delta <= 1: C = 6*10^15 delta^(-25/3) (C136 Theorem, refereed in C136 REFEREE).
* J_z/J = -1: zero defect (C129 Section 8).
* Seam at J_z/J = 0: both theorems are stated on closed intervals containing Delta = 0 and use the same C129 reductions and the same M >= M_0, so on [-1+delta, 1] the constant is max(4078, 6*10^15 delta^(-25/3)) = 6*10^15 delta^(-25/3).
* Seam at J_z/J = 1 (Heisenberg): included in Corollary I.10 by finite-matrix continuity at fixed L.

[Assumed, open] -1 < J_z/J < -1 + delta for every fixed delta: not covered uniformly. Conditional on (D), (D1), (W1)-(W3) of C138 with the stated constants, C <= 628 there (C138 Theorem B), and the whole class would have the single constant C = 4078. Conditional on (W1)-(W3) only, the constant on that strip is 5K_*(eta) + pi^2/4, of order 5 * 247 log(1/eta) = 1235 log(1/eta) (Section II.4), which is length-uniform at each eta but diverges logarithmically as Delta -> -1.

So the full U(1) theorem with one explicit constant is proved on [-1 + delta, 1] for every delta > 0, with C = 6*10^15 delta^(-25/3) there and C = 4078 on [0,1], and on the single point -1; the ferromagnetic strip remains open, with the missing inequality isolated as the logarithm-free (D) on the shell, plus (W1)-(W3).
