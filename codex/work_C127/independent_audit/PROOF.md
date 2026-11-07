# C127 proof and limits

2026-11-06. Natural logarithms and the bond convention of the shared target are used. A status tag applies to the following displayed calculation. All lattice lengths below are even. General positive-spectrum statements use L >= 4, where C125 constructs the common exact spectrum. The explicit XX formulas also cover L = 2.

## 1. Scope of the result

[Proved] A dimension bound on the range of the single-insertion map is not a bound on the number of crossed eigenstates detected by that map. The missing operation is taking the spectral cyclic span. The distinction occurs in the actual XX chain, not only in an abstract spectral construction.

[Proved] The exact XX excitation trace obeys, for all beta > 0 and even L >= 2,

    log Ztilde_XX(L) <= log 3 + (2 beta/pi) L/(L^2-1)
                     <= log 3 + 8 beta/(3 pi L).                 (1)

It also obeys the proposed scalar boson-product bound with

    S = 6, V' = 2e = 4eM, K = Kstar = log 3 + 8/(3e).            (2)

This is a statement about the full spin trace, including its parity and zero-momentum sectors.

[Assumed or unresolved] No universal product bound for the interacting positive-transfer class, and no lattice counterexample to that scalar bound with arbitrary dimension-dependent constants, is established here. In particular, the failure of a momentum-preserving comparison below must not be reported as a failure of the scalar trace bound. Section 5 proves that existence of the latter with unrestricted constants is equivalent to the requested K + kappa M beta/L estimate.

## 2. The single-insertion map

### 2.1 Rank, spectral support, and reconstruction

[Assumed] For this subsection let the crossed representation include a ground vector Omega, a nonnegative excitation generator Hc = Htilde - eps_0, and a unitary time-circle rotation U(tau) commuting with Hc, with U(beta) = 1 and U(tau) Omega = Omega. Suppose it also reconstructs single-site insertions as linear maps B(O) on a domain containing Omega. These hypotheses include more than the compactness construction of eigenvalues in C125. That construction alone does not reconstruct insertions or continuous rotations. If there are several ground vectors, the argument applies to each chosen ground vector separately; summing over them can multiply the rank by their number.

[Proved under this realization] Write omega_k = 2 pi k/beta and let P_k be the Fourier projection of U. For k != 0 define

    Phi_k(O) = beta^(1/2) P_k B(O) Omega,
    S_k = ran Phi_k.

At k = 0 first subtract the vacuum expectation from B(O). Since B(1) Omega = Omega, the identity is in the kernel of these maps. Thus

    dim S_k <= d^2-1.                                         (3)

One may choose a traceless basis of M_d to parametrize the quotient by the identity. The count is at most d^2-1 at each signed frequency, hence at most 2(d^2-1) for the pair k,-k. This counts source vectors, not eigenstates or oscillator species.

[Proved] In general Hc S_k is not contained in S_k. The invariant space detected by the sources is instead

    C_k = closure span { E_Hc(B) v : v in S_k, B Borel }.

For an isolated energy E the detected subspace E_Hc({E}) S_k has dimension at most d^2-1. There is no bound of that size on the sum over different E. With discrete spectrum, one vector with nonzero projections onto n distinct energies generates an n-dimensional spectral cyclic space. With continuous spectrum it can generate an infinite-dimensional space. Compressing Hc to S_k gives at most d^2-1 Rayleigh-Ritz values, but does not give its spectrum on C_k or bound the heat trace there.

### 2.2 What the C125 estimate implies

[Proved, using the inherited C109 estimate] At nonzero bosonic Matsubara frequency the KMS commutator identity is

    G_AB(i omega,x) = i integral_0^infinity
                         exp(-omega t) <[A_x(t),B_0]> dt.

With V = 4eM, the inherited bound is

    ||[A_x(t),B_0]|| <= 2||A||||B|| min(1,exp(-x+Vt)).

For 0 < omega < V, split the integral at x/V. Its two parts give

    |G_AB(i omega,x)|
      <= 2||A||||B|| { [exp(-omega x/V)-exp(-x)]/(V-omega)
                       + exp(-omega x/V)/omega }
      <= 2||A||||B|| V/[omega(V-omega)] exp(-omega x/V).          (4)

No beta-dependent constant has been suppressed. For omega = V the bound is 2||A||||B||(x+1)exp(-x)/V. The same conclusions for negative frequencies use adjoints and |omega|.

[Proved conditional spectral implication] Suppose the reflected insertion is compatible with the physical correlation function, so that, in the normalization above,

    <Phi_k(O), exp(-x Hc) Phi_k(O)> = G_reflected(O,O)(i omega_k,x)

and reflection preserves the operator norm. The left side is a positive Laplace transform. Equation (4), with A the reflected adjoint of O, gives the constant

    C_O,k = 2||O||^2 V/[|omega_k|(V-|omega_k|)].

If its spectral measure had mass m > 0 in [0,a] with a < |omega_k|/V, then

    m exp(-a x) <= C_O,k exp(-|omega_k| x/V)

on arbitrarily large even x, a contradiction. Exhausting the interval below |omega_k|/V proves

    E_Hc([0,|omega_k|/V)) Phi_k(O) = 0,
    Hc restricted to C_k >= |omega_k|/V,
    0 < |omega_k| < V.                                        (5)

The conclusion is also valid for quadratic forms on S_k. It does not make S_k invariant. A polynomial prefactor at omega = V still excludes energies below one.

[Assumed or unresolved] Positivity of the spectrum of T is not by itself the reflected-insertion identity used in (5). Ordinary A,B Matsubara correlators can have signed spectral residues. One must construct the reflected source pair, or prove a corresponding support theorem. C125 states its spectral implication conditionally for this reason. Neither the finite source rank nor the numerical similarity transform below repairs that reconstruction gap.

[Proved, inherited optimization] Retaining the exponential Lieb-Robinson parameter mu gives the asymptotic detected-sector lower rate

    alpha(omega) = omega/(4eM),           0 < omega <= 4eM,
                   log(omega/(4M)),      omega >= 4eM.            (6)

The matching case can have a polynomial prefactor in x. Equation (6) is not a linear cone at high frequency. It also gives no lower bound at omega = 0.

### 2.3 An actual lattice witness for the rank distinction

[Proved from the C125 exact XX representation] Take h = Sx Sx + Sy Sy, M = 1/2. There are four fermionic crossed levels at each positive odd Matsubara frequency nu = (2j+1) pi/beta, two of each sign of momentum, with energy asinh(nu). The first Fock sector has even total occupation. For fixed bosonic k = 1, pairs with momenta

    +(2j+1) pi/beta and -(2j-1) pi/beta, j >= 1,

have the same total p = 2 pi/beta and distinct energies

    E_j = asinh((2j+1) pi/beta) + asinh((2j-1) pi/beta).          (7)

Their number below any fixed positive energy cutoff grows without bound as beta grows. This is not an unbounded-species claim. Even a fixed-species Fock space has many multiparticle states at a given total momentum.

[Proved local detection] The XX density Matsubara bubble is a sum over fermionic nu of G_x(i nu) G_-x(i nu-i omega), with

    |G_x(i nu)| = exp(-|x| asinh|nu|)/sqrt(1+nu^2).

At even x its signed coefficient for the exponent

    asinh|nu| + asinh|nu-omega|

is, up to the common convention for the sign of the connected correlator,

    sign[nu(nu-omega)] /
      [beta sqrt((1+nu^2)(1+(nu-omega)^2))].

At omega = 2 pi/beta the terms nu = (2j+1) pi/beta, j >= 1, and their reflected partners have the exponents (7), nonzero coefficients, and the same sign. Different j give different exponents, so these terms do not cancel. A single local Sz frequency component therefore detects infinitely many distinct decay rates. A finite-dimensional invariant single-insertion space could only produce finitely many distinct rates. This proves the failure of that proposed invariance on the lattice without assuming positivity of the density bubble's individual residues.

[Numerical] At E <= 0.5 the displayed j >= 1 family has 4, 20, and 80 distinct energies for beta = 100, 500, and 2000, respectively. The one-source frequency map still has rank one. At m = 4, beta = 2, a single Sz source in the actual checkerboard transfer has 3, 8, 8, and 10 distinct resolved energies in a nonzero rotation sector for XX, Delta = 0.5, Heisenberg, and Delta = 2. These are finite-Trotter findings, not extrapolations to the quantum limit.

## 3. Multi-insertion states and the missing energy comparison

[Proved] The rank of a source space does not imply an injective or surjective Fock construction for its spectral cyclic completion. For several insertions, time ordering, relations in the local algebra, and the transfer evolution all enter. Even if products of local insertions are dense, density is only completeness. It gives no additive lower bound on the energy of a product state and no ordering of the eigenvalues by occupation numbers.

[Assumed sufficient comparison] One sufficient extra input would be a comparison of ordered eigenvalues, including multiplicities, between Hc and a reference whose states are labelled by

    (neutral sector, charges Q_a, occupations n_(k,s)),
    E_ref = sum_a Q_a^2/(2 rho_a beta)
            + sum_(k>=1,s<=S) n_(k,s) alpha(2 pi k/beta),

with a finite neutral multiplicity g <= exp(K0). Here S counts both momentum signs. A form-theoretic sufficient condition is an isometric embedding of the physical crossed space into the reference space such that Hc >= W* H_ref W on the relevant form domain. The min-max principle then gives the ordered-eigenvalue comparison and hence heat-trace domination. For infinite spaces compact resolvents and the indicated summability must also hold. This comparison is an additional assumption, not a consequence of (3)-(6).

[Proved location of the interaction problem] Source two-point functions determine spectral measures of vectors, not the energy of their products. Connected higher-point functions and interaction binding energies are the missing information. Attractive binding can lower product-state energies; null relations can remove states; additional charge or neutral sectors can fail to be generated by one insertion. A proof has to control these effects in an aggregate trace or eigenvalue comparison. The source rank alone controls none of them.

[Proved limitation of Gibbs variational comparison] The exact variational formula is

    log Tr exp(-L Hc) = sup_sigma { S(sigma) - L Tr(sigma Hc) }.

Substitution of a quasi-free trial state gives a lower bound on this logarithm. It does not give the upper bound wanted here. For H_lambda = H_ref + lambda W,

    d/dlambda log Tr exp(-L H_lambda) = -L <W>_lambda.

An upper bound therefore requires a uniform estimate of this interacting expectation, or a stronger spectral comparison. It is not supplied by matching the two-point functions of a trial state.

[Heuristic, not a proof] A polymer expansion could supply an upper bound if its connected activities had summable weights with the necessary suppression per temporal excitation, uniformly in beta/L. The usual small beta M expansion only controls a high-temperature region. In the required region beta/L can be arbitrarily large, and no uniform connected-activity estimate was derived. Assigning an independent decay penalty to each insertion would assume the missing additive-energy comparison.

[Proved limitation of the band argument] Even granting (6) on all detected cyclic spaces, energies below E < 1 can only involve single-insertion frequencies |omega| <= V E. This is a bound on which frequency labels occur, not on the rank of the spectral projection within one label. Equation (7) supplies that distinction already on XX. For general multiple-insertion histories there are several time arguments, and total frequency can hide arbitrarily large opposite frequencies. Counting band-limited single functions does not count these histories or control their energy without an additional many-insertion estimate.

[Proved] A bound N(E) <= exp(K+c beta E) gives, by integration by parts,

    Tr exp(-L Hc) = L integral_0^infinity exp(-LE) N(E) dE
                  <= exp(K) L/(L-c beta),

only when L > c beta. It does not control the opposite aspect-ratio regime. The bosonic occupation structure supplies a different counting law. If log Z_ref(s) <= K+B beta/s, then

    N_ref(E) <= exp(sE) Z_ref(s)
              <= exp(K+2 sqrt(B beta E))

after minimizing over s > 0. When a comparison is only available for s >= L0, this optimization is restricted accordingly. Square-root counting is the required additional structure, not a consequence of a linear count.

## 4. Zero modes and the logarithmic high-frequency tail

### 4.0 A lattice lower bound on K

[Proved] Consider the diagonal ferromagnetic d-state bond h(s,t) = -J delta_(s,t), J > 0. Its optimized half-spread is M = J/2. The original representative is optimal because swapping the two sites reverses any telescoping one-site gauge, and averaging the two representatives cannot increase the spectral spread. Its spatial transfer matrix has eigenvalues exp(beta J)+d-1 and exp(beta J)-1, the latter with multiplicity d-1, all positive. Therefore

    Ztilde(L,beta) = 1+(d-1) r_beta^L,
    r_beta = [exp(beta J)-1]/[exp(beta J)+d-1].

Choose even L asymptotic to exp(beta J/2). Then beta/L tends to zero, L(1-r_beta) tends to zero, and Ztilde tends to d. For every fixed S and V' proportional to M the oscillator product tends to one. Any uniform scalar majorant exp(K) times that product must therefore have

    K(d) >= log d.                                             (8a)

This is a lattice failure of K = 0 and a necessary lower bound on the additive constant. It is not a failure of the proposed majorant when K is allowed. There are no oscillator excitations to blame in this example; the surviving weight is the classical zero-frequency sector.

### 4.1 Charge sectors

[Assumed] Suppose the separate charge factor is bounded by independent integer Gaussian sums

    Z0 <= g product_a theta(rho_a beta/L),
    theta(a) = sum_(Q in Z) exp(-Q^2/(2a)),
    sum_a rho_a <= 4M.

This requires a finite-size lower bound on the charge energies as well as the normalization and conservation assumptions of the inherited Kohn budget. The budget alone does not prove a Gaussian spectrum. Neutral primary sectors are not controlled by that budget.

[Proved under these assumptions] For a > 0, Q^2 >= Q at integer Q >= 1 gives

    theta(a) <= 1 + 2/[exp(1/(2a))-1] <= 1+4a,
    log theta(a) <= 4a.

Consequently, independently of the number of charges,

    log Z0 <= log g + 4(beta/L) sum_a rho_a
            <= log g + 16 M beta/L.                            (8)

This estimate remains valid for limits of finitely many factors with summable stiffness. It avoids assuming a dimension-only count of U(1) sectors. At large a an individual Gaussian sum grows as sqrt(a); it cannot in general be bounded by a constant alone. Light states with energy tending to zero also prevent a K = 0 oscillator comparison at large L unless they have their own factor. A finite multiplicity g can be absorbed into K.

### 4.2 A conditional reference using the correct ultraviolet rate

[Proved conditional on the Fock comparison in Section 3] The logarithmic high-frequency rate does not by itself prevent a beta/L heat-trace bound. Let

    B_alpha(L,beta) = product_(k>=1)
                          (1-exp(-L alpha(2 pi k/beta)))^(-S).

The function F(w) = -log(1-exp(-L alpha(w))) is decreasing. Right-endpoint integration gives

    log B_alpha <= S beta/(2 pi) integral_0^infinity F(w) dw.

In the low-frequency part substitute w = V u. In the high-frequency part substitute w = 4M exp(u). Since V = 4eM, the matching point is u = 1. Thus

    integral F <= V integral_0^1 -log(1-exp(-Lu)) du
                  + 4M integral_1^infinity exp(u)
                                       [-log(1-exp(-Lu))] du
               <= V pi^2/(6L)
                  + 4M exp(-(L-1))/[(L-1)(1-exp(-L))].

For L >= 2 the second term is at most 8M/[e(1-exp(-2))L]. Therefore

    log B_alpha <= S c_alpha M beta/L,
    c_alpha = e pi/3 + 4/[pi e(1-exp(-2))].                      (9)

The infinite product converges for L > 1, as expected from its k^(-L) tail. Together with (8), the assumed comparison would imply

    log Ztilde <= log g + [16+S c_alpha] M beta/L.               (10)

This is a conditional theorem with explicit constants. What is unproved is the lattice Fock comparison, not control of the ultraviolet tail after that comparison.

## 5. The scalar product inequality is equivalent to the target form

[Proved] Define

    B_S(t) = product_(k>=1) (1-exp(-k/t))^(-S),
    b = -log(1-exp(-1)).

For f(u) = -log(1-exp(-u)), monotonicity gives

    log B_S(t) <= S t integral_0^infinity f(u) du
                = S pi^2 t/6.                                (11)

There is no positive uncontrolled remainder in this upper bound. Conversely, for t >= 1, retain k <= floor(t). Since floor(t) >= t/2,

    log B_S(t) >= S floor(t) b >= S b t/2.                      (12)

The integral in (11) is sum_(r>=1) 1/r^2, obtained by expanding f and integrating positive terms.

[Proved equivalence with changes of constants] Set x = M beta/L and V' = c M, with c > 0 fixed. A product estimate gives

    log Ztilde <= K + log B_S(c x/(2 pi))
                <= K + (S pi c/12) x.                         (13)

Conversely, suppose log Ztilde <= K0 + kappa x uniformly. Choose an integer

    S >= 4 pi kappa/(b c),
    K = K0 + 2 pi kappa/c.

For t = c x/(2 pi) >= 1, (12) dominates kappa x. For t < 1, the added constant in K dominates kappa x. This proves the product estimate for all x. The same conclusion holds with a separate zero factor satisfying (8), by modifying kappa. Thus, when K and S are freely dimension-dependent, deciding the scalar boson-product question would decide the original target. A failed derivation from insertion rank is not a refutation of either statement.

[Proved normalization] A linear reference with S = 2 and V' = v has leading logarithm pi v beta/(6L), the oscillator contribution of one c = 1 boson. For general S its coefficient is S pi V'/12, hence A <= S V'/2 if zero sectors contribute only subleading terms in the infrared coefficient. If a coarse bound such as (8) is used instead, its linear coefficient must be included in the resulting A bound. It cannot silently be dropped.

[Proved relation to the Casimir and doubling targets] In a common positive crossed representation the physical bulk pressure is -eps_0(beta). Hence

    lim_(beta->infinity) log Ztilde(L,beta)/beta = L e_0-E_L.

The estimate K+kappa M beta/L gives E_L >= L e_0-kappa M/L. Positivity also gives E_L <= L e_0, since Ztilde >= 1. Writing delta_L = L e_0-E_L therefore gives

    0 <= delta_L <= kappa M/L,
    E_(2L)-2E_L = 2 delta_L-delta_(2L) <= 2 kappa M/L.

On an even vacuum-compatible CFT sequence, A <= 6 kappa M/pi. The choice of sequence matters; this argument uses the even rings on which the positive transfer reproduces the physical trace, and a CFT vacuum expansion on that sequence when A is invoked. Small lengths outside a chosen L0 are handled by the inherited product-gluing bound.

## 6. An all-temperature XX proof

[Proved, using C125's exact product] Put a_n = asinh((2n+1) pi/beta), n >= 0, and b_n = asinh(2n pi/beta), n >= 1. In the notation of C125,

    Ztilde = (P_o^+ + P_o^-)/2 + 2 exp(-L delta_beta) P_e^+,
    P_o^+ = product_(n>=0) (1+exp(-L a_n))^4,
    P_o^- = product_(n>=0) (1-exp(-L a_n))^4,
    P_e^+ = product_(n>=1) (1+exp(-L b_n))^4,
    delta_beta > 0.

Term by term P_o^- <= P_o^+ and P_e^+ <= P_o^+. Thus

    Ztilde <= 3 P_o^+.                                        (14)

The factor three accounts for the first sector and the two shifted-sector copies. It bounds their partition sums; it does not assert that all three are degenerate vacua.

[Proved] Define f_L(z) = log(1+exp(-L asinh z)), z >= 0. Let g(u) = log(1+exp(-Lu)). Then g' <= 0, g'' >= 0, while asinh has nonpositive second derivative for z >= 0. The chain rule gives f_L'' >= 0. Midpoint convexity on consecutive intervals of width 2 pi/beta therefore gives

    sum_(n>=0) f_L((2n+1) pi/beta)
      <= beta/(2 pi) integral_0^infinity f_L(z) dz.

With z = sinh u and log(1+y) <= y,

    integral f_L(z) dz
      = integral_0^infinity cosh(u) log(1+exp(-Lu)) du
      <= integral_0^infinity cosh(u) exp(-Lu) du
      = L/(L^2-1), L > 1.

Combining this with (14) proves (1). In particular no global linear crossed cone has been assumed. The logarithmic large-z dispersion is retained through the integral.

[Proved] To get (2), let y = beta/L and t = e y/pi, so the oscillator exponent is 2 pi k L/(2e beta). If t >= 1, (12) with S = 6 gives

    log B_6(t) >= 3 b t >= (3/e)t = 3y/pi >= 8y/(3 pi),

where b >= 1/e follows from -log(1-z) >= z. If t < 1, then 8y/(3 pi) < 8/(3e), absorbed by Kstar. Equations (1) and (14) prove the claimed full trace bound with S = 2(d^2-1) at d = 2. This does not establish the same S or K for interacting d = 2 chains.

[Proved consistency checks] As beta tends to zero at fixed even L, the exact XX excitation trace tends to one, consistent with (1). For fixed beta and L tending to infinity, it also tends to one; Kstar is slack. At beta/L tending to infinity in the infrared, A = 1 gives log Ztilde with leading coefficient pi beta/(6L). The coefficient 8/(3 pi) in (1) exceeds pi/6, as required. Equation (1) gives A <= 16/pi^2, consistent with A = 1 but not sharp. The parity factor in (14) and the positive delta_beta were retained, so no physical spin sector was discarded.

## 7. Which majorizations are actually refuted

[Proved lattice failure of a stronger statement] A momentum-preserving comparison with oscillators of energy at least |p|/V' fails for XX for every finite V'. C125 supplies locally detected states

    p_k = 2 pi k/beta,
    E_k = asinh(pi/beta)+asinh((2k-1) pi/beta).

At fixed beta, E_k/p_k tends to zero. A Fock reference with additive energies |p|/V' has total energy at least |total p|/V' by the triangle inequality, even after adding zero-momentum charge states with nonnegative energy. For sufficiently large k, E_k < |p_k|/V'. No fixed finite multiplicity factor exp(K) repairs a heat-trace comparison separately in that momentum sector as L tends to infinity: the smaller actual decay exponent wins. The reference sector has finite trace at some fixed L and a lowest energy at least |p_k|/V', so its large-L decay cannot evade this argument through its multiplicities.

[Proved scope distinction] This counterexample refutes a comparison preserving momentum with a global linear reference energy. It does not refute the sum over momenta. In fact Section 6 proves a scalar linear-boson product bound for the very same XX lattice chain. It also does not refute a proposed comparison restricted to low-energy states and augmented by a controlled ultraviolet factor.

[Assumed sufficient repair, not a minimal theorem] Using the nonlinear ultraviolet rate (6), controlling zero sectors as in (8), and proving a many-insertion ordered-energy comparison as in Section 3 would suffice by (9)-(10). For XX, the exact fermionic factorization supplies a concrete replacement, leading to (1). For general positive-transfer chains no lattice criterion guaranteeing that comparison, and no uniquely minimal additional hypothesis, is established. Calling the comparison a consequence of fixed local dimension would repeat the rank error of Section 2. No actual lattice counterexample to the unrestricted scalar trace target has been constructed.

## 8. Numerical interpretation

[Numerical] checks.py constructs the positive checkerboard operator directly. In the alternating y-rotated spin basis, g = exp(-beta h'/m) has positive realignment. Let Ke and Ko be the two column matrices of C125 and set

    U = sqrt(Ke) Ko sqrt(Ke), T = sqrt(U).

Translation by two time-history entries has order m and commutes with Ke, Ko, and U. Orbit Fourier bases give blocks labelled by k modulo m, with representative p = 2 pi k/beta. An eigenvalue lambda gives excitation energy -log(lambda/lambda_max)/2. A one-entry time shift exchanges the two checkerboard layers and is not asserted to commute with this U. Finite m therefore resolves discrete, aliased Matsubara momenta, not an independent continuum frequency for every integer k.

[Numerical] A diagonal Sz insertion D_k on one column sublattice is transported through C = Ke Ko = sqrt(Ke) U invsqrt(Ke) as B_k = invsqrt(Ke) D_k sqrt(Ke). Applying B_k to the top eigenvector and projecting onto U's eigenspaces measures the spectral cyclic support of that one source. B_k need not be Hermitian. Its overlap test does not assume that the physical density correlator has positive spectral residues. Only well-conditioned beta >= 2, m >= 3 cases are used for this diagnostic; overlap weights below 1e-8 and eigenvalues below 1e-13 times the top one are omitted.

[Numerical] The physical L = 4 Trotter trace is independently contracted in the original spin Hilbert space. Its relative agreement with the spatial construction is at most 4.03e-12. Trotter error is separately measured against exact L = 4 diagonalization. The largest column dimension in the local run is 256. A low-frequency cone test on the full finite-step spectrum has a negative margin -0.07247 at Heisenberg, beta = 8, m = 2. The continuum estimate (4) is not asserted for a finite checkerboard step, or for arbitrary undetected eigenstates. At beta = 8 even m = 4 has substantial Trotter errors, recorded explicitly in CHECKS.md. Full momentum-sector arrays and all trace comparisons are saved, with a compact summary in CHECKS.md.

[Assumed or unresolved] Small Trotter number does not certify a continuum-temperature inequality. In particular beta = 8, m = 1 is coarse. Larger Trotter jobs are prepared for SLURM, not run locally or submitted from this sandbox.

[Assumed work record] This proof and its scripts are isolated in independent_audit because the parent checks.py and output schema changed concurrently during verification. Parent-folder edits from the other run are preserved. All statistics cited here are reproduced by this subfolder's scripts, with no imports from the parent run.
