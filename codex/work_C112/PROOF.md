# C112 derivations

## 1. Which deficit is being measured

[Definition] On a periodic lattice of L sites let H = sum_x h_x, H_w = sum_x w_x h_x, W = sum_x w_x, and E_L = E_0(H). The finite-ring adaptation deficit used in the computations is

    D_L(w) = (E_L/L) W - E_0(H_w).

[Proved] A translation-averaged ground-state density matrix of H has bond expectation E_L/L. Its expectation of H_w is (E_L/L)W. The variational principle gives D_L(w) >= 0, including when the ground state is degenerate. For constant weights D_L is zero.

[Definition] The infinite-volume-reference deficit on the same finite ring is

    Def_infty,L(w) = e_infty W - E_0(H_w)
                  = D_L(w) + (e_infty-E_L/L)W.

[Proved under the stated Casimir expansion] The last term is pi A W/(6L^2) + o(W/L^2). Consequently Def_infty,L does not vanish for constant weights. It cannot satisfy a pure gradient bound on finite rings without a finite-size term. The specification's exact translation-invariance statement uses E_L/L on a finite ring, not e_infty. All reported ring ratios use D_L. The exact gapped counterexample below has E_L/L = e_infty, so this distinction does not affect that counterexample.

## 2. Stress-tensor normalization and all infrared sectors

[Assumed] In this section the modulated operator is the stress tensor itself, with no extra derivative operator or microscopic contribution at the same order. The infrared is a product of unitary vacuum CFT sectors with charges c_a and positive velocities v_a. Set A = sum_a c_a v_a. This assumption about the modulated density is stronger than identifying only the uniform Hamiltonian with its CFT limit.

[Proved within CFT] For a smooth positive periodic w define

    ell_w = integral_0^L dx / w(x),
    b = L/ell_w,
    f'(x) = b/w(x),
    G(w) = integral_0^L dx w'(x)^2/w(x).

The map f is an orientation-preserving diffeomorphism of a circle of length L. Its Schwarzian is

    {f,x} = f'''/f' - (3/2)(f''/f')^2
          = -w''/w + (1/2)(w'/w)^2.

Thus integral w{f,x} dx = G(w)/2, since the integral of w'' vanishes. For one physical chiral stress tensor the transformation law gives

    U(f) [v integral w T_ch dx] U(f)^(-1)
      = vb integral T_ch dy - cv/(24 pi) integral w{f,x} dx.

There are two chiralities. Adding their anomalies gives -cv G/(24 pi). The uniform vacuum energy is -pi cv/(6L), so the rescaled vacuum energy is -pi cv/(6 ell_w). Summing sectors gives

    E_0(H_w) - e_infty W = -pi A/(6 ell_w) - A G(w)/(24 pi),       (2.1)
    D_L(w) = A G(w)/(24 pi)
             + (pi A/6)[1/ell_w - W/L^2].                         (2.2)

These are continuum identities in the vacuum representation. Lattice matching and its remainders remain separate assumptions. If left and right velocities differ, each chirality contributes c_ch v_ch/(48 pi) to the coefficient of G. With the convention A = (1/2)sum_ch c_ch v_ch, the answer is unchanged. No on-site symmetry or conserved charge is used.

[Proved within CFT] The stress-tensor coefficient is therefore

    kappa = 1/(24 pi) = 0.0132629119243246... .                    (2.3)

[Source and convention check] The chiral transformation law and conjugation used above are equations (5.19a) and (5.20) of [Per Moosavi, Inhomogeneous Conformal Field Theory Out of Equilibrium](https://link.springer.com/article/10.1007/s00023-021-01118-0), Lemma 5.3. The integration, ground-energy specialization, and subtraction in (2.1)-(2.2) are written out here to fix the normalization. They also agree with C105's integrated coefficient A/(12 pi) when G has leading term 2 log L. This does not justify applying (2.1) to every microscopic energy-density gauge.

### 2.1 Cosine, including the finite-ring cancellation

[Proved within CFT] Let q = 2 pi n/L, n a positive integer, and w = 1 + epsilon cos(qx), |epsilon| < 1. Direct integration gives

    W = L,
    ell_w = L/sqrt(1-epsilon^2),
    G(w) = L q^2 [1-sqrt(1-epsilon^2)].

For the third identity put u = qx. The average of epsilon^2 sin(u)^2/[1+epsilon cos(u)] is 1-sqrt(1-epsilon^2), obtained by polynomial division and the average of 1/[1+epsilon cos(u)]. Substitution gives

    D_L(w) = (pi A/(6L))(n^2-1)[1-sqrt(1-epsilon^2)]
           = (pi A epsilon^2/(12L))(n^2-1) + O(epsilon^4/L).       (2.4)

Consequently the ratio to the continuum gradient functional is

    D_L/G = [A/(24 pi)](1-1/n^2).                               (2.5)

The fundamental n = 1 cosine has zero deficit in the CFT for every |epsilon| < 1. Increasing L at fixed n does not remove this factor. To extract kappa without a Casimir correction one needs 1 << n << L, or a localized bump with 1 << sigma << L. Alternatively one compares directly with (2.2).

[Proved within CFT, independent second-order check] In one sector,

    H_0 = (2 pi v/L)(L_0 + Lbar_0 - c/12),
    V = (pi v/L)(L_n + L_(-n) + Lbar_n + Lbar_(-n)),
    H_w = H_0 + epsilon V.

The Virasoro algebra and vacuum conditions give

    ||L_(-n)|0>||^2 = <0|[L_n,L_(-n)]|0>
                    = (c/12)n(n^2-1),
    excitation energy = 2 pi v n/L.

Second-order perturbation theory gives a positive deficit equal to

    epsilon^2 (pi v/L)^2
      * [2(c/12)n(n^2-1)] / (2 pi v n/L)
      = pi cv epsilon^2(n^2-1)/(12L).

Summing sectors reproduces (2.4). The velocity is in the numerator because the two matrix elements contribute v^2 and the excitation denominator contributes v. The n = 1 norm vanishes by the global conformal vacuum condition.

### 2.2 Gaussian bump

[Proved within CFT] For w = 1 + epsilon f and a localized f on a large ring, expand

    ell_w = L - epsilon integral f + epsilon^2 integral f^2 + O(epsilon^3),
    1/ell_w - W/L^2
      = -epsilon^2 integral (f-fbar)^2/L^2 + O(epsilon^3),
    G(w) = epsilon^2 integral f'^2 + O(epsilon^3).

Thus

    D_L = A epsilon^2/(24 pi)
            [integral f'^2 - (2 pi/L)^2 integral(f-fbar)^2]
          + O(epsilon^3).                                      (2.6)

For f(x) = exp[-x^2/(2 sigma^2)] on the line,

    integral f'^2 dx
      = sigma^(-4) integral x^2 exp[-x^2/sigma^2] dx
      = sqrt(pi)/(2 sigma),
    Def = A epsilon^2/(48 sqrt(pi) sigma) + O(epsilon^3/sigma).    (2.7)

The ring term is relatively of order sigma^2/L^2. The numerical bump is periodized by summing Gaussian images and has sigma = L/10. Its finite-ring term is retained, not discarded. The cosine and bump have the same kappa.

[Proved, consistency checks] Equation (2.2) vanishes for constant weights. Its second-order expression is nonnegative by the periodic Poincare inequality, with equality for the constant and first harmonics. It is additive over decoupled sectors, homogeneous in velocities and weights, and scales as inverse length for a stretched profile. All four checks pass.

## 3. Exact XX response

[Definition] Use h_x = S_x^x S_(x+1)^x + S_x^y S_(x+1)^y, equivalently hopping amplitude 1/2. Its dispersion is epsilon_k = cos k, e_infty = -1/pi, M = 1/2, and A = 1. Center the Fourier phase on the bond. The energy-density Fourier component has one-particle matrix element

    <k+q|h_q|k> = cos(k+q/2).

[Proved] The second-order cosine deficit is the particle-hole sum

    D/epsilon^2 = (1/4) sum_(k occupied, s=+1,-1, k+sq empty)
                   cos(k+sq/2)^2 / [cos(k+sq)-cos k].            (3.1)

The two orientations have equal sums. Set s_q = q/2, 0 < q < pi. At a Fermi edge, write t = k+q/2-3pi/2, with -s_q < t < s_q. The squared matrix element is sin(t)^2 and the excitation energy is 2 sin(s_q) cos(t). Replacing sum_k/L by integral dk/(2pi) gives

    D/(L epsilon^2)
      = [log(sec(s_q)+tan(s_q))-sin(s_q)]/[4 pi sin(s_q)].        (3.2)

Indeed integral sin(t)^2/cos(t) dt = log(sec t+tan t)-sin t. Expanding the numerator as s_q^3/3 + O(s_q^5) yields

    D/(L epsilon^2) = q^2/(48 pi) + O(q^4).

Since G/L = epsilon^2 q^2/2 + O(epsilon^4 q^2), this fixes kappa to 1/(24 pi) independently of the CFT calculation. This is a thermodynamic response with many wavelengths. Equation (2.4), not its n -> infinity limit, applies to a fixed harmonic on a finite ring.

[Numerical] The finite-ring check diagonalizes the one-particle weighted hopping matrix, with antiperiodic fermions for L divisible by four and the half-filled spin-chain vacuum. The uniform energy used in the subtraction is the exact -csc(pi/L). At epsilon = 0.02 the ratios D_lattice/D_CFT for n = 2 are 1.00807963, 1.00201099, 1.00050219 at L = 32,64,128. For the Gaussian they are 1.01258964, 1.00311337, 1.00077630. The n = 1 deficit is zero up to subtraction error. These are exact free-fermion reductions evaluated in floating point, not a many-body truncation.

## 4. An exact gapped counterexample to the stated lattice lemma

[Proved] Fix any d >= 2. On one site define

    a = diag(-delta, delta, ..., delta),  delta > 0,
    X = |0><1| + |1><0|,
    h = (a tensor I + I tensor a)/2
          + t(X tensor I - I tensor X).                        (4.1)

The discrete divergence cancels on a uniform ring, so H = sum_i a_i. Its unique ground state is |0> on every site, its gap is 2 delta, its bulk energy is -delta, and it has no conformal infrared sector. Thus A = 0 for every delta > 0.

[Proved] The gauge-optimized scale is M = delta. Removing the divergence leaves a bond with spectral half-spread delta. Conversely, the expectations in |00> and |11> differ by 2 delta in every gauge, because a discrete divergence has zero expectation in |u,u>. The spectral spread is at least the difference of these expectations. This supplies the matching lower bound. This argument applies without change for d = 2,3,4.

[Proved] For positive periodic weights put

    bar_w_i = (w_i+w_(i-1))/2,
    g_i = w_i-w_(i-1).

The weighted chain is a sum of independent one-site matrices,

    H_w = sum_i [bar_w_i a_i + t g_i X_i].

The two active levels have eigenvalues +/-sqrt(delta^2 bar_w_i^2+t^2 g_i^2). Every spectator has energy delta bar_w_i, above the negative eigenvalue. Therefore the exact deficit is

    D_L(w) = sum_i [sqrt(delta^2 bar_w_i^2+t^2 g_i^2)-delta bar_w_i]
           = sum_i t^2 g_i^2 /
              [sqrt(delta^2 bar_w_i^2+t^2 g_i^2)+delta bar_w_i].  (4.2)

This is an analytic many-body spectrum, valid for arbitrary length and all positive weights. The second expression avoids numerical cancellation.

[Proved] For a slowly stretched positive profile, with delta and t fixed,

    D_L = (t^2/(2delta)) sum_i g_i^2/bar_w_i
             + O(sum_i t^4 g_i^4/(delta^3 bar_w_i^3)),
    D_L / [M sum_i (w_(i+1)-w_i)^2/w_i]
          -> t^2/(2delta^2).                                  (4.3)

The shifted denominators differ by a relative amount tending to zero. For every proposed finite constant at fixed d, first choose t/delta large enough and then stretch the profile far enough. The fourth-order remainder is negligible for that bond, while the ratio exceeds the proposed constant. Thus a remainder that is o(G) for each fixed bond cannot rescue the statement.

[Proved] There is no finite supremum over the bonds in the literal lemma at any d >= 2. This is already true by taking t -> infinity at fixed delta, with H, M, and A unchanged. Taking t fixed and delta -> 0 gives another sequence, whose raw half-spread sqrt(delta^2+4t^2) stays bounded. Even replacing M by this raw scale does not cure the latter sequence, since D_L/(raw_scale G) ~ t/(4delta). The M-minimizing representative has t = 0 and deficit zero. The counterexample therefore refutes the all-representatives statement, not a proposed theorem with a specified minimizing representative.

[Numerical] For t = 1, delta = 0.001 and L = 4096, the scalar sum (4.2) gives ratios 499999.999912 for a cosine and 499999.999980 for a Gaussian. These numbers hold for each of d = 2,3,4. The same formula was checked against independent full-space diagonalization at (d,L) = (2,4),(3,3),(4,3), with energy errors at most 1.2e-16. The 4096-site calculation is only a scalar sum of known one-site energies, not a large diagonalization.

## 5. The XXZ ferromagnetic endpoint

### 5.1 The standard representative

[Assumed, exact XXZ data recorded in the campaign] For -1 < Delta <= 1,

    gamma = arccos Delta,
    A = v = pi sin(gamma)/(2gamma),  M = 1/2.

Thus A tends to zero as sqrt(1+Delta)/sqrt(2) when Delta decreases to -1, even though its bond derivative diverges. The conventional stress-tensor prediction for the gradient coefficient tends to zero.

[Proved] At Delta = -1 on an even ring, a staggered rotation turns the conventional bond into -S_i dot S_(i+1). Each bond is bounded below by -1/4, and every fully symmetric spin state attains that value on every bond. Consequently every positive w has E_0(H_w) = -W/4 and D_L(w) = 0. For Delta < -1, each standard XXZ bond is bounded below by Delta/4 and the fully polarized product states attain it, so D_L(w) is again zero. This checks the endpoint and the gapped polarized side without relying on ED.

### 5.2 A discrete divergence produces a divergent long-wavelength response

[Proved] Add t(S_i^z-S_(i+1)^z) to the XXZ bond at Delta = -1. The uniform Hamiltonian is unchanged, and M remains 1/2 by definition of the gauge quotient. After the staggered rotation the weighted perturbation is a ferromagnetic exchange modulation plus the field

    t sum_i (w_i-w_(i-1)) S_i^z.

Take even L, q = 2pi/L and w_i = 1+epsilon cos[q(i+1/2)]. Let |S,m> denote the fully symmetric multiplet, S = L/2. The exchange modulation annihilates this multiplet after its constant expectation is subtracted. The field has zero matrix element inside the multiplet because sum_i(w_i-w_(i-1)) = 0 and P S_i^z P = S_total^z/L.

[Proved] Write g_i = cos[q(i+1/2)]-cos[q(i-1/2)], so sum g_i^2 = L(1-cos q). Within |S,m>,

    <(S_i^z)^2> = 1/4,
    <S_i^z S_j^z> = (m^2-L/4)/[L(L-1)],  i != j.

The norm of the field-created state is therefore

    ||t sum_i g_i S_i^z |S,m>||^2
      = t^2 [(L^2/4-m^2)/(L(L-1))] sum_i g_i^2.

This state has excitation energy 1-cos q. To check the denominator, use

    H_F - E_F = (1/2)sum_i(1-P_(i,i+1)).

Acting on S_q^z|S,m>, each permutation exchanges the two Fourier coefficients while acting trivially on |S,m>. Summing the differences yields (1-cos q)S_q^z|S,m>. This proves the denominator without a dispersion assumption.

[Proved] Degenerate second-order perturbation theory is diagonal in m and lowers the energy most at m = 0. It gives the actual ground-energy expansion at fixed L,

    D_L = epsilon^2 t^2 L^2/[4(L-1)] + O_L(epsilon^4),
    lim_(epsilon -> 0) D_L/(M G_discrete)
      = t^2 L/[4(L-1)sin(pi/L)^2]
      ~ t^2 L^2/(4pi^2).                                     (5.1)

The energy is even in epsilon because translation by L/2 reverses it. A sequence epsilon_L small enough for the fixed-L expansion proves divergence along smooth positive profiles. The underlying endpoint has quadratic magnons and no nonzero linear-velocity CFT coefficient. Its uniform ground energy is exactly -L/4, and the CFT coefficient approaching from Delta > -1 tends to zero. There is no divergence of A.

[Numerical] At t = 1 and epsilon = 0.001, the measured ratios at L = 6,8,10,12 are 1.19999907, 1.95097646, 2.90892016, 4.07131480. The second-order values are 1.20000000, 1.95097918, 2.90892665, 4.07132815. Halving epsilon from 0.002 to 0.001 reduces the departure by a factor close to four. The calculation uses the balanced sector, which the preceding effective-multiplet calculation proves is the minimizer for sufficiently small epsilon at each fixed even L.

### 5.3 A fixed-amplitude counterexample without an exchange of limits

[Proved] The same fixed bond also fails the bound when epsilon is held fixed in (0,1) and L tends to infinity. This can be seen from a product-state variational bound, used here to lower-bound the deficit and refute its proposed upper bound. In the rotated ferromagnetic variables choose spin-coherent directions

    n_i = (sqrt(1-b^2 sin(qi)^2), 0, b sin(qi)),  0 < b < 1.

The field contribution to the trial energy is

    t sum_i (w_i-w_(i-1)) n_i^z/2
      = -t b epsilon L sin(q/2)/2.

The exchange energy above -W/4 is sum_i w_i |n_(i+1)-n_i|^2/8. The smooth interpolation obeys |n'(x)|^2 <= b^2 q^2/(1-b^2). Cauchy-Schwarz on each unit interval bounds the exchange cost by (1+epsilon)L b^2 q^2/[8(1-b^2)]. Therefore

    D_L >= t b epsilon L sin(q/2)/2
             - (1+epsilon)L b^2 q^2/[8(1-b^2)].                 (5.2)

At q = 2pi/L the first term tends to pi t b epsilon/2 and the second tends to zero. In contrast,

    G_discrete <= 2L epsilon^2 sin(q/2)^2/(1-epsilon) = O(1/L).

Thus D_L/(M G_discrete) diverges at least linearly in L for this fixed positive-amplitude smooth cosine. The raw half-spread of the fixed bond is finite as well, so replacing M by that raw scale would not fix this example. An o(G) remainder for the fixed bond cannot hide the failure. The use of product states is only a variational lower bound, not a claim that the weighted ground state is a product state.

[Heuristic, not needed for either counterexample] In the critical XXZ phase the same divergence couples to the Luttinger density. Completing the square in the density energy pi v rho^2/(2K) gives an extra coefficient t^2 K/(2pi v) multiplying G. With the usual XXZ K = pi/[2(pi-gamma)], this diverges as t^2/[4(1+Delta)] while A tends to zero. This is the continuum explanation of the soft density response. The exact lattice endpoint calculation (5.1) and gapped example (4.2) do not require that dictionary or an exchange of its infrared limits.

## 6. How the random-bond scale is checked

[Proved] The quotient norm has the finite-dimensional dual lower bound

    M >= Tr(hD)

for every Hermitian D satisfying Tr D = 0, Tr_1 D = Tr_2 D, and ||D||_1 <= 1. Indeed Tr(hD) is unchanged by scalar and divergence additions, and is at most the centered operator norm times ||D||_1. Every explicit gauge representative supplies the upper bound given by its half-spread.

[Numerical method] The code minimizes a smooth approximation to the half-spread over all traceless Hermitian one-site gauge generators, using increasing spectral inverse scales. From its two extremal spectral density matrices it constructs D = (rho_plus-rho_minus)/2. Put r = Tr_2 D - Tr_1 D and replace

    D -> D - (r tensor I-I tensor r)/(2d).

The difference of partial traces of the subtracted operator is r, so the result satisfies the dual gauge constraint. Dividing by its trace norm gives a feasible dual matrix. The reported M_lower and M_upper bracket the optimum up to floating-point error. They are not rigorous interval-arithmetic enclosures. The reported ratio interval uses both endpoints. No numerical optimization is used for the M values in the analytic counterexamples or the standard spin and swap chains.

## 7. Scope of the conclusion

[Proved] The literal smooth-modulation lemma is false because it combines a gauge-invariant right-hand scale with a gauge-dependent deficit and permits the representatives (4.1) and (5.1). Both failures have a mechanism and a known uniform spectrum. Neither is a counterexample to a bound on A.

[Assumed, unresolved] A repaired statement must specify the modulated density, explain its relation to the representative attaining M, and control additional response terms when matching it to the stress tensor. On rings it must also use the finite-size prescription in Section 1. The present work does not prove or disprove the repaired statement for a specified minimizing representative.

[Numerical limitation] Small tests of standard XXZ and SU(d) representatives and of nine generic random bonds show finite ratios, but the sizes and sample counts cannot test uniformity over critical, gapped, and nonconformal infrared limits. The requested moderate-size unrestricted survey is prepared as a SLURM job and was not executed because Typhon connections timed out. No absence of a counterexample is inferred from these small tests.
