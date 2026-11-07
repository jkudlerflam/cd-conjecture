# cd_crossed_fock: derivations

2026-11-06. Opus 5.5 agent. Conventions are those of `specs/_bound_background.md`, `specs/_Cd_target.md` and `codex/work_C125/PROOF.md` (cited as C125). Lattice spacing and hbar are one, logarithms are natural. A tag at the start of a paragraph applies to everything up to the next tag. "Numerical" always means a computation in this folder, reproduced (in reduced form) by `run_checks.py`.

## 0. Summary of the logic

The route of `specs/C127_crossed_majorization.md` asks for a bosonic (Fock) majorant of the crossed partition function Z~(L) = Tr (T_beta/lambda_0)^L of a positive spatial transfer operator. Sections 1 and 2 set up the objects and record what a majorant buys. Section 3 makes the single-insertion sector precise and proves a symmetry selection rule for it. Section 4 proves the majorization for XX with explicit constants (the free-fermion check of the route). Section 5 states the conditional theorem that the route needs (a trace-level domination of Z~ by the Fock space of single-insertion levels) and shows it implies the target with A <= (d^2 - 1)(12.95 + 1.91/gamma_0) M, gamma_0 an infrared input. Section 6 gives the obstruction: on the Heisenberg chain the level-by-level (min-max) form of that domination is false, because of a crossed composite that the SU(2) selection rule makes invisible to every single-site insertion, and two-insertion states are far from additive. Section 7 collects the numerical methods and checks.

## 1. The crossed formulation and what a majorant gives

### 1.1 Identities

[Proved] Let T = T_beta >= 0 be a positive transfer operator with eigenvalues lambda_0 >= lambda_1 >= ... and Tr T^L = Z(L, beta) for the even L considered. Put eps_j = -log lambda_j and Z~(L) = sum_j exp(-L (eps_j - eps_0)) = Z(L,beta)/lambda_0^L. With rho_L = T^L / Tr T^L,

    log Z~(L) = -log lambda_max(rho_L) = S_min(rho_L),
    R_L = log[Z(L)^2/Z(2L)] = S_2(rho_L),
    log Z~(L) <= R_L <= 2 log Z~(L).                                    (1.1)

The last line uses 1 <= Z~(2L) <= Z~(L) (the vacuum term is one and every other term decreases). Also log Z~(L) is convex and nonincreasing in L, and

    log Z~(L) = sum_{k>=0} 2^{-k-1} R_{2^k L} + lim_n 2^{-n} log Z~(2^n L).      (1.2)

[Proved] Since lambda_0^L >= <psi|T^L|psi> for every unit vector psi, every column state gives a rigorous upper bound log Z~(L) <= log Z(L) - log <psi|T^L|psi>; the crossed vacuum saturates it. Physical boundary states psi have <psi|T^L|psi> exponentially small in beta relative to lambda_0^L (a boundary free energy), so test vectors useful for the target must be nonlocal in imaginary time. Equivalently lambda_0 = lim_n Z(nL)^{1/(nL)} and the target is the statement that gluing n rings of length L into one costs at most K + kappa M beta/L per ring in log Z. This is the reconnection formulation of section 9.86 of REPORT.md and is recorded only to fix language.

### 1.2 The majorant and its normalization

[Proved] Write F(x) = -sum_{k>=1} log(1 - exp(-2 pi k x)). The majorant of the spec is log Maj(L) = K + S F(L/(V' beta)). By the modular property of the Dedekind eta function,

    F(x) = pi/(12x) + (1/2) log x - pi x/12 + F(1/x),   F(x) <= pi/(12 x),          (1.3)

so log Maj <= K + S pi V' beta/(12 L). One c = 1 boson of velocity v has two oscillators per |k| (both signs of crossed momentum) and its oscillator part is exactly 2F(L/(v beta)), with log of it equal to pi v beta/(6L) + O(log) as L/beta -> 0, which is the CFT value (pi/6) c v beta/L. So S counts oscillators per |k| with both signs and S = 2 is one boson (normalization check of the spec).

[Proved] If log Z~(L) <= K + kappa M beta/L for all even L and beta, then letting beta -> infinity at fixed L (Z~ is then dominated by the ground states, log Z~(L)/beta -> L e_0 - E_L, where e_0 is the bulk ground-energy density) gives the Casimir lower bound

    L (L e_0 - E_L) <= kappa M   for every even L,                                (1.4)

and on the CFT subsequence L(L e_0 - E_L) -> pi A/6, so A <= (6/pi) kappa M. A majorant with (S, V') gives kappa M = S pi V'/12 and A <= S V'/2.

[Proved] Necessary conditions on any majorant. (a) By (1.4) at beta -> infinity, S pi V'/12 >= D_max := sup_L D(L), D(L) = L(L e_0 - E_L). (b) As L -> infinity at fixed beta, Z~(L) -> g_0, the multiplicity of the top crossed eigenvalue, so K >= log g_0. (c) If S V'/2 equals the infrared value A exactly, the zero-mode factor of the CFT (of order (beta/L)^{1/2} per compact boson) and the (S/2) log x term of (1.3) force K to grow like log(beta/L); a finite K needs S V'/2 > A strictly. These three facts fix what "smallest S, V', K" can mean: V' is bounded below by the lattice deficit D at the smallest ring, not by the CFT.

## 2. The checkerboard transfer

### 2.1 Construction and positivity

[Proved, C125 section 2.1, reimplemented] With 2m imaginary-time slices of width beta/(2m), Trotter step delta = beta/m and gate g = exp(-delta h') for the bond h' in its reflection-positive frame, the column space is the space of periodic histories (C^d)^{(x)2m} and

    K_e = G^{(x) m} on the slice pairs (0,1),(2,3),...,  K_o = P^T K_e P on (1,2),...,(2m-1,0),
    G[(s_j, s_{j+1}), (t_j, t_{j+1})] = g[(s_{j+1}, t_{j+1}), (s_j, t_j)],
    Z_m(L, beta) = Tr (K_e K_o)^{L/2},  U_sym = K_e^{1/2} K_o K_e^{1/2} >= 0,           (2.1)

with P the cyclic shift by one slice. G >= 0 whenever -h' = sum_alpha J_alpha A_alpha (x) conj(A_alpha) + B (x) 1 + 1 (x) conj(B) + c with J_alpha >= 0: the exponential of a sum of such terms is a norm limit of products of gates with positive realignment, and the realignment of a product of such gates is again positive because (A (x) conj A)(A' (x) conj A') = (AA') (x) conj(AA'). (C125 states the version without the B terms. Here the term delta (B (x) 1 + 1 (x) conj B) of -delta h' exponentiates to exp(delta B) (x) conj exp(delta B), whose realignment is a rank-one positive matrix.) The crossed eigenvalues are u_j = exp(-2 eps_j).

[Numerical] `cflib.py` builds K_e^{1/2}, K_o as sparse Kronecker products and checks (2.1) against the direct Trotter trace Tr[(exp(-delta H_e) exp(-delta H_o))^m] on rings L = 2, 4, 6 for every model below, to 1e-14 (`run_checks.py`, check T1). G is PSD in every case.

### 2.2 The d = 2 class

[Proved] For d = 2 the sufficient condition above covers, in the frame where it holds, all bonds

    h' = -Jx SxSx + Jy SySy - Jz SzSz - (B (x) 1 + 1 (x) B)/2,   B = gx Sx + gz Sz,   Jx, Jy, Jz >= 0,      (2.2)

because conj(Sx) = Sx, conj(Sz) = Sz, conj(Sy) = -Sy and B is real, so -h' = sum_a J_a S^a (x) conj(S^a) + (B (x) 1 + 1 (x) conj B)/2. In the original frame (a pi rotation about y on alternate sites) these are the antiferromagnetic XYZ chains h = Jx SxSx + Jy SySy + Jz SzSz with fields that are staggered in Sx and Sz; the transverse-field Ising chain H = -sum sx sx - g sum sz is already of the form (2.2) with no rotation. A uniform field along any axis in the original XXZ frame is not of the form (2.2). The covered critical families are the XXZ line 0 <= Delta <= 1 (c = 1, A = pi sin(gamma)/(2 gamma) with Delta = cos gamma, at J = 1) and the Ising lines of the XY chain in a field (c = 1/2). Models used here (`models.py`): XX, XXZ at Delta = 0.5, Heisenberg, XXZ at Delta = 2 (gapped Neel), TFIM at g = 1 (c = 1/2, v = 2, M = sqrt 2), XY on its Ising line with gamma = 3 (v = 0.75), XYZ (1, 0.5, 0.3) (gapped), XXZ Delta = 0.5 with staggered field 0.3 (gapped). M is the gauge-optimized half spread, computed by `gauge_M` (minimization over one-site gauge terms o (x) 1 - 1 (x) o).

### 2.3 Symmetries of the crossed transfer

[Proved] (i) Crossed momentum. The rotation R2 of the column by two slices (one Trotter step, imaginary time beta/m) maps even pairs to even pairs and odd pairs to odd pairs, so it commutes with K_e^{1/2}, K_o and U_sym. Its eigenvalues are exp(2 pi i n/m), crossed momentum p = 2 pi n/beta modulo 2 pi m/beta. (ii) Charges. If h' conserves Sz_x - Sz_{x+1} (positive-frame XXZ, which conserves the staggered magnetization), every gate obeys s_j - s_{j+1} = t_j - t_{j+1} on its pair (Sz units), so Q(s) = sum_{j even}(s_j - s_{j+1}) is conserved by K_e and by K_o (on odd pairs Q_o = -Q_e because the total telescopes to zero). Q is the net Sz world-line winding through the bond. If h' changes spins only in pairs (XYZ, TFIM, no gx), the parity of the number of flips at even slices is conserved instead. (iii) On-site symmetries. Suppose h' is real and (u (x) conj u) h' (u (x) conj u)^dagger = h' for all u in a group G. Let V(u) = (x)_j w_j on the column with w_j = conj(u) on even slices and w_j = u on odd slices. Writing out G[(s_j, s_{j+1}),(t_j, t_{j+1})] = g[(s_{j+1}, t_{j+1}),(s_j, t_j)], the conjugate (conj u (x) u) G (conj u (x) u)^dagger on an even pair is the realignment of (u (x) conj u) g (u (x) conj u)^dagger = g, so V commutes with K_e; on odd pairs the roles of u and conj u are exchanged and the complex conjugate of the hypothesis (g is real) gives the same for K_o. Hence V commutes with K_e^{1/2}, K_o and U_sym. For the Heisenberg chain in the positive frame, -h' = sum_a S^a (x) conj(S^a) is a multiple of the projector on sum_a |a a> plus a constant, invariant under u (x) conj u for every u in SU(2). [Numerical] For a random u in SU(2), ||[V, K_e]|| and ||[V, K_o]|| vanish to machine precision at m = 3, while the uniform choice w_j = u does not commute (check I3). [Numerical] The decomposition into (charge, momentum) sectors reproduces all traces (check T1), so the charges are conserved as claimed.

## 3. The single-insertion sector

### 3.1 Insertion states

[Proved] Insert a one-site operator O at the column that sits between K_e and K_o, at slice j; on the column it acts on the slice-j coordinate (as O or O^T according to the parity of j, which is immaterial for spans). Write Ohat_j for this column operator and psi_0 for the top eigenvector of U_sym. For two insertions on the same column parity, at columns 1 and 2r + 1 of an infinite ring,

    lim_L Z(Ohat, Ohat')/Z = <psi_0| K_e^{1/2} Ohat K_o K_e^{1/2} U_sym^{r-1} K_e^{1/2} Ohat' K_o K_e^{1/2} |psi_0> / u_0^{r+1}.      (3.1)

To see this, write the limit as <L_0|A|R_0> with R_0 = K_e^{1/2} psi_0 and L_0 = K_e^{-1/2} psi_0 the right and left top eigenvectors of K_e K_o, use (K_e K_o)^{r-1} K_e^{1/2} = K_e^{1/2} U_sym^{r-1}, and note that the inverse square root cancels against the leftmost K_e. So the right insertion state is

    |Y_{O,j}> = K_e^{1/2} Ohat_j K_o K_e^{1/2} |psi_0>,   |Y_1> = u_0^{1/2} |psi_0>,        (3.2)

and the crossed spectral support of the insertion is that of the expansion of (3.2) in the eigenbasis of U_sym. In fact the identity holds exactly at finite L: with X_O = K_e^{1/2} Ohat K_o K_e^{1/2},

    Tr[K_e Ohat K_o (K_e K_o)^{r-1} K_e Ohat' K_o (K_e K_o)^s] = Tr[X_O U_sym^{r-1} X_{O'} U_sym^s],

because K_o (K_e K_o)^{k} K_e^{1/2} = K_o K_e^{1/2} U_sym^{k}. Both sides involve only K_e^{1/2}, never its inverse, so (3.1) is well defined even though G has eigenvalues of order delta.

### 3.2 Span, multiplicity and the selection rule

[Proved] Let W_n be the span of the momentum-n components of |Y_{O,j}> over traceless O and over the two slice parities j in {0, 1} (the rotation R2 relates all other slices). Then dim W_n <= 2(d^2 - 1) at finite Trotter number. In continuous time only one insertion time per frequency exists and the bound is d^2 - 1; at finite m the second parity is a second, inequivalent insertion point inside one Trotter step (in effect it adds an operator localized on the neighbouring bond). Charge selection: Sz insertions keep Q, S+ and S- insertions change it by one unit, so for the d = 2 U(1) models W_n lives in Q in {0, +1, -1} with at most 2 dimensions in each.

[Proved] Selection rule. Suppose the bond has an on-site symmetry group G as in 2.3(iii) and the top eigenvector psi_0 is unique. Then psi_0 is invariant (a one-dimensional representation of a connected compact simple group is trivial; for U(1) it is a phase, which does not change the argument), V Ohat_{O,j} V^dagger = Ohat_{w_j O w_j^dagger,j}, and the base vector K_o K_e^{1/2} psi_0 is invariant, so O -> |Y_{O,j}> intertwines the action O -> w_j O w_j^dagger on the traceless matrices sl(d) (the adjoint representation, or its complex conjugate, which is equivalent for SU(2)) with V. Hence every W_n lies in the isotypic components of V that occur in the adjoint representation of G on sl(d). For the Heisenberg chain (G = SU(2), d = 2) single-site insertions reach only crossed states of total spin 1. Every crossed state of spin 0 other than the vacuum, and every state of spin 2 or more, is orthogonal to all W_n: it is a composite, reachable only by two or more insertions. The same argument with G = U(1) gives the charge selection above.

[Numerical] psi_0 is unique (gap in the (Q, n) = (0, 0) sector) in every computed case, and in the Heisenberg (0, 0) sector the spin-0 levels have overlap with W below 1e-30 (`run_checks.py`, check I1).

## 4. The XX majorization with explicit constants

[Proved] Theorem 4.1. Let h = SxSx + SySy (M = 1/2, v = 1, A = 1). Fix X > 1/2, put V1 = X/asinh(X), take r > 1 + 6/(pi^2 X V1) and V' = r V1. Then for every even L >= 2 and every beta > 0

    log Z~_XX(L) <= K(X, r) + 2 F(L/(V' beta)),
    K(X, r) = sup_{x > 0} Phi(x),
    Phi(x) = 2F(x) - 2F(x/r) + log(theta(x)^2 + theta_1(x)^2/2) + 4 (2X)^{-2} + 1/(pi X V1 x),          (4.1)

with theta(x) = sum_Q exp(-pi x Q^2), theta_1(x) = sum_Q exp(-pi x (Q^2 + Q)). For X = 2 and r = 1.5 this gives V' = 2.078 (V'/M = 4.16) and K = 2.04; for X = 2, r = 2 it gives V'/M = 5.54, K = 1.44. Consequently log Z~_XX(L) <= K + (pi/6) V' beta/L and the majorization route reproduces A_XX <= V' with S = 2, the normalization of one boson.

Proof. Start from the exact crossed product of C125 (eq. (6) there, proved in its section 3, and re-verified here against exact diagonalization, check X1),

    Z~(L) = (1/2)[A_+ + A_-] + 2 exp(-L delta_beta) B,
    A_+- = prod_{n>=0} (1 +- exp(-L kappa_n))^4,  B = prod_{n>=1} (1 + exp(-L mu_n))^4,
    kappa_n = asinh(nu_n), nu_n = (2n+1) pi/beta;  mu_n = asinh(2 n pi/beta);  delta_beta > 0.

Since A_- <= 1 <= A_+ and exp(-L delta_beta) <= 1, Z~ <= A_+ + 2B.

Step 1 (low frequencies). asinh is concave with asinh(0) = 0, so asinh(nu)/nu decreases and asinh(nu) >= nu/V1 for 0 <= nu <= X. Hence every factor with nu_n <= X obeys 1 + exp(-L kappa_n) <= 1 + q^{n+1/2}, q = exp(-2 pi x), x = L/(V1 beta). Extending the product to all n only adds factors >= 1. By the Jacobi triple product, prod_{n>=1}(1 - q^n)(1 + q^{n-1/2})^2 = sum_Q q^{Q^2/2} = theta(x), so

    prod_{n>=0} (1 + q^{n+1/2})^4 = theta(x)^2 exp(2F(x)).

For the periodic frequencies, prod_{n>=1}(1 - q^n)(1 + q^n)(1 + q^{n-1}) = sum_Q q^{Q(Q+1)/2} = theta_1(x) and prod(1 + q^{n-1}) = 2 prod_{n>=1}(1 + q^n), so prod_{n>=1}(1 + q^n)^4 = theta_1(x)^2 exp(2F(x))/4.

Step 2 (high frequencies). asinh(nu) = log(nu + sqrt(nu^2 + 1)) >= log(2 nu), so for the factors with nu > X, log(1 + exp(-L kappa)) <= exp(-L kappa) <= (2 nu)^{-L}. The frequencies have spacing 2 pi/beta and (2nu)^{-L} decreases, so

    4 sum_{nu_n > X} (2 nu_n)^{-L} <= 4 (2X)^{-L} + 4 (beta/2pi) int_X^infinity (2 nu)^{-L} d nu
                                  = 4 (2X)^{-L} + (beta/pi) (2X)^{1-L}/(L - 1).

For 2X > 1 the function L (2X)^{-L}/(L - 1) decreases on L >= 2 (the ratio of consecutive values is (L^2 - 1)/(2X L^2) < 1), so with L >= 2 the right side is at most h := 4(2X)^{-2} + (beta/L)/(pi X) = 4(2X)^{-2} + 1/(pi X V1 x). The same bound holds for the periodic frequencies.

Step 3. Combining, Z~ <= exp(2F(x) + h) (theta(x)^2 + theta_1(x)^2/2), so log Z~ - 2F(x/r) <= Phi(x), and x/r = L/(V' beta). Phi is finite on (0, infinity). As x -> 0, (1.3) gives 2F(x) - 2F(x/r) = -(pi/6)(r - 1)/x + log r + O(x), while log(theta^2 + theta_1^2/2) = log(3/(2x)) + o(1) by Poisson summation and the last term is 1/(pi X V1 x); so Phi -> -infinity exactly when (pi/6)(r - 1) > 1/(pi X V1), which is the hypothesis on r. As x -> infinity, Phi -> log 3 + 4(2X)^{-2}. QED.

[Numerical] K(X, r) is the maximum of an explicit smooth function of one variable. `xx_theorem.py` evaluates it on 2 x 10^5 logarithmically spaced points in [1e-4, 1e3] (derivative of Phi in log x below 0.08 near the maximum, grid step 8e-5, so the grid error is below 1e-5), and checks (4.1) against the exact log Z~ on 70 values of beta in [0.05, 4096] and 19 even L in [2, 4096]: the smallest margin is 1.39 for (X, r) = (2, 1.5) and 1.26 for (2, 2), both at L = 2. The optimum over all shapes found by direct scan of the exact product is V'/M = 4.08 at K = log 2 and S = 2 (`xx_exact_majorant.py`), against the necessary value 12 D(2)/(2 pi M) = 2.78 from section 1.2(a).

[Proved, consistency] XX limit. As beta -> infinity at fixed L/beta small, 2F(L/(V' beta)) ~ pi V' beta/(6L) with V' > 1 = v, matching log Z~_XX ~ (pi/6) beta/L (one boson, c = 1, v = 1); the theorem cannot reach V' = 1 because the zero modes theta^2 cost log(beta/L) (section 1.2(c)).

## 5. What the route needs: a conditional theorem

[Proved] Theorem 5.1. Work in continuous imaginary time (m = infinity) for a positive-transfer chain with V = 4eM. Suppose there is a family of one-particle levels eps_{k,a}, k in Z (Matsubara frequency omega_k = 2 pi k/beta), a = 1, ..., S_1, such that

    (H1) eps_{k,a} >= alpha(|omega_k|) for k != 0, alpha(omega) = omega/V for omega <= V, log(omega/(4M)) for omega > V;
    (H2) Z~(L) <= exp(K_0) prod_{k,a} (1 - exp(-L eps_{k,a}))^{-1} for every even L >= 2;
    (H3) eps_{0,a} >= gamma_0/(beta M) for some gamma_0 > 0.

Then for every even L >= 2 and every beta > 0

    log Z~(L) <= K_0 + 2.32 S_1 e^{-2} + S_1 [2 pi e/3 + 9.28/(pi e) + 1/gamma_0] M beta/L,              (5.1)

and hence, by section 1.2, A <= (6/pi) S_1 [2 pi e/3 + 9.28/(pi e) + 1/gamma_0] M = S_1 [4e + 55.7/(pi^2 e) + 6/(pi gamma_0)] M, numerically S_1 [12.95 + 1.91/gamma_0] M.

Proof. Split the logarithm of the product in (H2) into k = 0, 0 < |omega_k| <= V and |omega_k| > V. For 0 < |omega_k| <= V, (H1) and monotonicity give -log(1 - exp(-L eps)) <= -log(1 - exp(-L |omega_k|/V)); summing over both signs of k and extending to all k >= 1 gives at most 2 S_1 F(L/(V beta)) <= 2 S_1 pi V beta/(12 L) = S_1 (2 pi e/3) M beta/L by (1.3). For |omega_k| > V, (H1) gives exp(-L eps) <= (4M/|omega_k|)^L <= e^{-L} <= e^{-2}, so -log(1 - y) <= y/(1 - e^{-2}) <= 1.16 y, and sum_{omega_k > V} (4M/omega_k)^L <= e^{-L} + (beta/2 pi) int_V^infinity (4M/omega)^L d omega = e^{-L} + (2 M beta/pi) e^{1-L}/(L - 1) <= e^{-2} + (4/(pi e)) M beta/L, using that L e^{1-L}/(L-1) decreases on L >= 2. With both signs this contributes at most 2.32 S_1 [e^{-2} + (4/(pi e)) M beta/L]. For k = 0, -log(1 - e^{-t}) <= log(1 + 1/t) <= 1/t with t = L gamma_0/(beta M). Adding the three pieces gives (5.1). QED.

[Proved conditionally, C125 section 4] (H1) holds for the single-insertion levels that are detected by a reflection-positive insertion pair, conditional (as in C125) on the positive spectral representation of the reflected insertion correlator, by the Laplace-transform Lieb-Robinson bound |G_AB(i omega, x)| <= 2||A|| ||B|| V exp(-omega x/V)/[omega(V - omega)] (0 < omega < V) and its logarithmic high-frequency version. In continuous time S_1 = d^2 - 1 by section 3.2. So the route gives A <= (d^2 - 1)[12.95 + 1.91/gamma_0] M, a polynomial C(d), provided (H2) and (H3) hold.

[Heuristic] (H3) is an infrared input. For U(1) sectors it can be replaced by the zero-mode factor of the spec, sum_Q exp(-Q^2 L/(2 rho_s beta)) with sum rho_s <= 4M (Kohn budget), which is harmless; a neutral sector with an arbitrarily light zero-frequency primary is not controlled by anything proved here.

[Proved] (H2) is the whole difficulty. It implies the target with explicit constants (Theorem 5.1), so it is at least as strong as the target. Its natural strengthenings are false or unproved, as section 6 shows.

## 6. The obstruction

### 6.1 Level-by-level domination fails on the Heisenberg chain

Write e_0 = 0 <= e_1 <= e_2 <= ... for the crossed excitation energies eps_j - eps_0 (with multiplicity) and f_0 = 0 <= f_1 <= ... for the energies of the bosonic Fock space built on the single-insertion levels (all multisets, energies added). (H2) with K_0 = 0 would follow from the level-by-level statement e_k >= f_k for all k (equivalently N(E) <= N_F(E) for every E), and that statement is what any proof by an operator inequality H~ - eps_0 >= Phi^dagger H_Fock Phi (Phi an isometry) would give, by min-max.

[Proved, given two numerical facts stated in the proof] Proposition 6.1. For the antiferromagnetic Heisenberg chain, at any beta and Trotter number at which the lowest spin-0 excitation eps_s of the crossed transfer lies below twice the lowest single-insertion level eps_t, the level-by-level domination is false, globally and in the (Q, n) = (0, 0) sector.

Proof. By the selection rule of section 3.2 every single-insertion level has spin 1. The lowest one is the spin-1 multiplet at zero crossed momentum (the staggered magnetization), with energy eps_t and multiplicity 3. [Numerical fact 1: the next single-insertion levels, at crossed momentum +-1, lie above eps_s in all computed cases.] So the Fock space has, below min(2 eps_t, next level), only the vacuum and those three one-particle states. [Numerical fact 2: the state at eps_s is a spin singlet; it has no partner level in the Q = +-1, +-2 sectors at the same crossed momentum, while the spin-2 composite at 0.625 (m = 6, beta = 8) does.] The true crossed spectrum therefore has, in addition, the spin-0 state at eps_s < 2 eps_t, orthogonal to every single-insertion state, and e_4 = eps_s < f_4 = 2 eps_t. QED.

[Numerical] The hypothesis holds from beta = 4 on, at every Trotter number computed. Values (m = 7; the m-dependence at fixed beta is of order (beta/m)^2):

| beta | eps_t | eps_s | eps_s/eps_t | lambda_crit = eps_s/(2 eps_t) |
| --: | --: | --: | --: | --: |
| 2 | 0.788 | above 2 eps_t | | 1 (no violation) |
| 4 | 0.401 | 0.758 | 1.89 | 0.945 |
| 6 | 0.259 | 0.467 | 1.80 | 0.902 |
| 8 | 0.183 | 0.321 | 1.75 | 0.877 |

At beta = 8 the ratio eps_s/eps_t is 1.84, 1.80, 1.77, 1.75, 1.74 for m = 4, 5, 6, 7, 8 (m = 8: eps_t = 0.190, eps_s = 0.332, lambda_crit = 0.871); at beta = 6 and m = 8 it is 1.80 (lambda_crit = 0.899).

[Assumed published, Affleck, Gepner, Schulz and Ziman 1989] In the infrared, eps_t and eps_s are the zero-momentum crossed levels of the spin-1 and spin-0 components of the SU(2)_1 primary g of dimension 1/2, with logarithmic corrections of opposite sign, x_t = 1/2 - g/4 and x_s = 1/2 + 3g/4 with g of order 1/log(beta). Hence eps_s/eps_t -> 1 and lambda_crit -> 1/2 as beta -> infinity: the violation persists and grows in the continuum limit. The data above decrease monotonically toward this limit.

[Proved] Consequently no isometry Phi can satisfy H~ - eps_0 >= Phi^dagger H_F Phi with H_F the free Fock Hamiltonian of the exact single-insertion levels. Any proof of (H2) must either work at the level of traces, paying for composites with the large overcounting of the Fock space at higher energy, or degrade the one-particle energies.

### 6.2 The degraded statement, and how far it holds

[Numerical] Define lambda_crit = min_k e_k/f_k over the enumerated levels (energies up to a cutoff of about three times the lowest one-particle level). Level-by-level domination holds with one-particle energies multiplied by any lambda <= lambda_crit, which only increases V' in the majorant by 1/lambda. On all eight d = 2 models (m <= 7, beta <= 8, and Heisenberg at m = 8), the global lambda_crit is at least 0.871 (the Heisenberg singlet at beta = 8, m = 8 is the worst case; every other model has global lambda_crit = 1 except staggered-field XXZ, 0.917 at beta = 8). The sector-resolved version (each (Q, n) separately) is at least 0.84 wherever the Fock enumeration is complete; in three high-temperature cases (beta = 2: XX at m = 6, 7 and XY at m = 6) a sector has a true level below the enumeration cutoff (E about 3.4 to 4) and no Fock level there, so lambda_sector is undetermined. The trace-level domination (H2) holds with K_0 = 0 in every computed case, with a margin log Z_F - log Z~ that grows with beta at fixed L (for Heisenberg at L = 2 it is 2.5 at beta = 4 and 11 at beta = 8, m = 5).

[Heuristic, supported by 6.1] In the infrared the Heisenberg composite forces lambda <= 1/2, so the weakest level-by-level statement that survives on actual lattice chains is

    (H2_lambda)  e_k >= lambda(d) f_k for all k, with a fixed lambda(d) in (0, 1/2],

which implies (H2) with the one-particle energies scaled by lambda, hence the target with V replaced by V/lambda. (H2_lambda) is the minimal extra input identified here. It is a crossed binding-energy bound: a crossed state that no single insertion reaches must cost at least lambda times the energy of the cheapest multi-particle Fock state in its quantum numbers.

### 6.3 Two-insertion states are not additive

[Numerical] `two_insertion.py` builds |ab> = sum Y_a(j_a + 2 r_1) Y_b(j_b + 2 r_2) |base> with Fourier phases, removes the vacuum and the single-insertion span of the target sector, and records its crossed spectral weights. At m = 7, beta = 6 (lowest one-particle levels 0.26 to 0.6):

- like charges or like signs of frequency: the weight sits well above eps_a + eps_b (Heisenberg Q = 2 from two Q = 1 triplet levels at 0.259: dominant level 0.927, 1.8 times the sum; XX Q = 2: 0.969 against 0.525). This is the crossed image of the quadratic growth of charged zero-mode energies.
- opposite charges or opposite frequencies: most of the weight sits far below the sum (Heisenberg (Q, n) = (1, 1) with (-1, -1): 67 percent of the weight below 0.9 (eps_a + eps_b), dominant level the spin-0 composite at 0.467 against the sum 1.188; XX two n = 1 bosons into n = 2: 99 percent below, at 1.629 against 1.938, from the concavity of asinh).

So crossed energies of multi-insertion states are neither additive nor sub- or super-additive in a uniform way. The Fock count survives because each true low state is reached by many insertion labels, not because energies add.

### 6.4 Why the composites are the general obstruction

[Heuristic] In a free n-boson crossed theory realized with fewer than n independent single-site operators, the oscillators outside the span of the operators' current components are reached only as composites, for example as the currents in the operator product of a vertex operator with its conjugate, at energy |omega|/v. These composites have exactly the cone energies allowed by (H1), and their number grows with n at fixed d. The single-insertion count (H1 with S_1 = d^2 - 1) therefore does not by itself bound the number of crossed species: the C125 lesson reappears inside the lattice framework, as the composite sector. Whether lattice locality at fixed d bounds composites (for instance through (H2_lambda)) is exactly as open as the target. The Heisenberg singlet shows that composites are present and light already in the simplest interacting chain, and that any such bound must lose at least a factor 2 in the zero-mode energies.

[Proved, lattice fact] One lattice-specific mechanism that a future proof could use: two insertions at the same site and coinciding times collapse to a single insertion of the product operator, which is again in M_d. Composites that are not single insertions therefore come from insertions at finite time separation. Turning this into a quantitative binding bound would require a control of the crossed transfer at imaginary-time separations of order 1/M that is not available here.

## 7. Numerical methods and consistency checks

[Numerical] Spectra. Full crossed spectra of U_sym, resolved by charge and crossed momentum, for m = 2..8 (XX), 2..7 (other U(1) models and the Z_2 models; SLURM scripts for m = 8, 9 are in the folder, not run), at beta in {0.5, 1, 2, 3, 4, 6, 8, 12, 16} (`run_spectra.py`, `run_spectra_model.py`). Single- and two-insertion analyses with eigenvectors for m = 4..7 at beta in {2, 4, 6, 8} (`analysis_insertion.py`, `two_insertion.py`).

[Numerical] Continuous time. eps_0(beta) by Lanczos at m up to 11 with tensor-structured matvecs (`eps0_lanczos.py`) and Richardson extrapolation in delta^2 with three correction terms (`infinite_m.py`). Check against XX's exact p_c: errors 6e-14, 2e-12, 1e-10, 1.5e-9, 6e-8, 8e-7 at beta = 1, 2, 3, 4, 6, 8. Combined with exact ring spectra for L <= 14 this gives log Z~ at m = infinity, which converges from above in m like (beta/m)^2 (Heisenberg, beta = 4, L = 4: 0.829, 0.756, 0.714, 0.688 for m = 4..7 against 0.614 at m = infinity; the differences times m^2 are 3.4 to 3.6).

[Numerical] Zero temperature (`zeroT.py`). Exact e_0 (Bethe ansatz for XXZ, free fermions for the XY family; for the gapped XYZ and staggered-field models e_0 = E_18/18, exponentially accurate) and Lanczos ground energies of even rings L <= 18. D(L) = L(L e_0 - E_L) is largest at L = 2 and decreases monotonically in L for all eight models, to pi A/6 for the critical ones (XX: D(18) = 1.0509 M against pi A/6 = 0.5236 = 1.0472 M).

[Numerical] Consistency checks: the finite-m XX crossed levels converge to C125's exact crossed spectrum like m^{-2} (at beta = 4 the lowest (Q, n) = (0, 1) level is 1.4059, 1.4155, 1.4218 for m = 6, 7, 8 against 2 asinh(pi/4) = 1.4425, with m^2 times the difference equal to 1.32 for all three; the lowest (Q, n) = (1, 0) level converges to delta_beta = 0.4164 in the same way, which identifies the odd-charge sectors with C125's periodic sector); the XX exact product against exact diagonalization (L <= 10, beta in {0.5, 2, 8}: error below 1e-6 at L = 2, where the product converges slowly, and below 3e-11 for L >= 4); the XX m = infinity pipeline (Richardson plus asymptotics) reproduces the exact-product value V'_min/M = 4.08 at S = 2, K = log 2; the one-boson normalization (L/beta) log Z~ -> pi A/6 (XX: 0.526 at L = 16, beta = 8192, against 0.5236); crossed traces equal direct Trotter traces to 1e-14.
