# C116. Ring traces, spatial transfer rank, and modular data

2026-10-29. Natural logarithms throughout. Claims labeled proved below either follow from finite-dimensional identities or are explicitly conditional on the stated infrared input. The requested universal bound remains open.

## 1. Setting and infrared assumptions

[Assumed] Let H_L = sum_j h_(j,j+1) be the periodic, one-site translation-invariant nearest-neighbor chain of the background specification. The oriented bond h need not be invariant under exchanging its two sites. The local dimension is d. Its gauge-optimized scale is

    M = inf_(o,alpha) ||h + o tensor I - I tensor o + alpha I||.

Choose a representative of norm m <= M + epsilon. Coboundaries telescope on every ring. A scalar shift adds alpha L to H_L. Every bulk-canceling combination below is invariant under both changes, so epsilon can be sent to zero after proving a norm estimate.

[Assumed] Write e_inf for the infinite-chain ground energy density and

    g_L(beta) = log Tr exp(-beta H_L),
    p(beta) = lim_(N -> infinity) g_N(beta)/N.

The pressure limit exists for bounded finite-range interactions. We use a vacuum-compatible sequence of L closed under doubling. Full traces, including all charge sectors on each selected ring, are used. This prescription selects ring lengths, not a projection onto the ground state or one charge sector.

[Assumed] The full infrared trace has a torus scaling limit at fixed tau > 0,

    g_L(tau L) + tau e_inf L^2 -> F(tau),                         (1)
    F(tau) = (pi A/6) tau + log g_vac + o(tau) as tau -> infinity,
    A = sum_a c_a v_a.

This is stronger input than convergence of the ground energy alone. It requires control of the trace over excited states. No theorem establishing (1) for all microscopic chains is supplied here. Marginal corrections may slow its convergence.

[Proved within the assumed CFT] For an actual product of nonchiral sectors,

    exp F(tau) = product_a Z_a(v_a tau),
    Z_a(y) = Tr_(H_a) exp[-2 pi y (L_0 + Lbar_0 - c_a/12)].       (2)

The vacuum term gives pi c_a v_a tau/6. Positive excitation energies contribute a sublinear term at large tau, assuming a finite torus trace and finite vacuum multiplicity. A general infrared Hilbert space can impose selection rules between sectors. Then a sum of character products replaces the literal product in (2). The large-tau slope in (1) still follows if the conformal vacuum is present. A finite number of character labels is not a finite number of states. Each character includes infinitely many descendants.

[Proved normalization caveat] No spatial-reflection assumption is needed. When a nonchiral sector has unequal left and right speeds, its contribution is

    A_a = c_a (v_(a,L) + v_(a,R))/2.

In the dual thermal coefficient the corresponding contribution is

    B_a = c_a (1/v_(a,L) + 1/v_(a,R))/2,

not c_a divided by the average speed. Equivalently, resolve chiral components and assign each its half-central-charge weight. The simpler notation B = sum_a c_a/v_a below applies when both speeds in each sector coincide, or when that chiral convention is used.

[Proved within modular CFT] If the sector partition functions have their usual modular transformations with a vacuum in the dual channel, then

    F(tau) = pi B/(6 tau) + O(1) + exponentially small terms

as tau tends to zero. For a product of modular-invariant theories with unique vacua the O(1) term vanishes. The associated infinite-volume low-temperature pressure is

    p(beta) = -beta e_inf + pi B/(6 beta) + o(beta^-1).            (3)

The lattice itself has not been assumed modular invariant at finite L or beta.

## 2. Which ratios retain the Casimir term

[Proved from (1)] Direct substitution gives

    g_L(tau L) - g_(2L)(2 tau L)
        = 3 tau e_inf L^2 + o(1).                               (4)

The two torus terms cancel. Thus Z(L,beta)/Z(2L,2 beta) both retains the bulk and loses the desired torus function. Bulk-normalizing each factor makes its scaling limit one, which still loses A.

[Proved from (1)] Coefficients 4 and -1 repair this cancellation:

    Q_L(beta) = 4 g_L(beta) - g_(2L)(2 beta),
    Q_L(tau L) -> 3 F(tau),
    A = (2/pi) lim_(tau -> infinity) lim_(L -> infinity)
                         Q_L(tau L)/tau.                         (5)

More generally, for fixed positive s_i,t_i and coefficients a_i,

    sum_i a_i g_(s_i L)(t_i beta)

cancels the bulk if sum_i a_i s_i t_i = 0. Its large-tau vacuum slope is

    (pi A/6) sum_i a_i t_i/s_i.

Integer lengths and vacuum compatibility must be enforced.

[Proved from (1)] The equal-temperature comparison is better suited to eliminating the thermal bulk:

    R_L(beta) = 2 g_L(beta) - g_(2L)(beta),
    R_L(tau L) -> r(tau) = 2 F(tau) - F(tau/2),
    r(tau) = (pi A/4) tau + log g_vac + o(tau),
    A = (4/pi) lim_(tau -> infinity) lim_(L -> infinity)
                         R_L(tau L)/tau.                         (6)

At small tau the term pi B/(6 tau) cancels from r. At every lattice L and beta, all terms proportional to L at that beta cancel from R, including nonconformal contributions to the pressure.

[Proved from (1) and (3)] Another useful bulk subtraction is

    Xi_L(beta) = g_L(beta) - L p(beta),
    Xi_L(tau L) -> F(tau) - pi B/(6 tau),
    A = (6/pi) lim_(tau -> infinity) lim_(L -> infinity)
                         Xi_L(tau L)/tau.                        (7)

The pressure subtraction is a limit at fixed beta before beta = tau L is inserted. Interchanging these operations without (3) is unjustified.

[Proved, finite-dimensional] At fixed L, write

    W_L(beta) = Tr exp[-beta(H_L - E_L)].

Then

    R_L(beta) = beta [E_(2L) - 2 E_L]
                  + 2 log W_L(beta) - log W_(2L)(beta).          (8)

All energies and degeneracies enter W. Sending beta to infinity at fixed L gives

    R_L(beta)/beta -> E_(2L) - 2 E_L.

Under the vacuum Casimir expansion, sending L to infinity afterwards recovers A through (6). This establishes the sequential zero-temperature extraction without needing full trace convergence at fixed tau. Equation (1) is needed for the alternative fixed-aspect-ratio order requested here.

## 3. A same-shape upper bound has an additional slow-sector obstruction

[Proved within the known XX scaling limit] Consider two decoupled XX chains per site, with the combined bond

    h_epsilon = [h_XX^(1) + epsilon h_XX^(2)]/(1+epsilon),
    0 < epsilon <= 1.

Here h_XX = Sx tensor Sx + Sy tensor Sy. The settled free-chain normalization gives

    d = 4,       M = 1/2,
    v_1 = 1/(1+epsilon),       v_2 = epsilon/(1+epsilon),
    A = v_1 + v_2 = 1,
    B = 1/v_1 + 1/v_2 = (1+epsilon)(1+1/epsilon).                 (9)

The optimized M also follows by using the product extremal bond eigenstates, whose one-site marginals are maximally mixed, so coboundaries cannot reduce their energy separation.

For each fixed tau, take the infrared limit at fixed epsilon first. The slow sector then has

    F_XX(v_2 tau) = pi/(6 v_2 tau) + o(1)

as epsilon tends to zero, by modular vacuum dominance. Consequently the scaling limit of Q_L(tau L) diverges at fixed d, M, A, and tau.

[Proved consequence] No upper estimate of the form

    Q_L(beta) <= K(d) + k M f(d) beta/L

can hold universally at all scales, or even for its scaling limit, with finite K(d), k and f(d). There is also an elementary obstruction at beta = 0, where Q_L(0) = 2 L log d. The slow-sector example shows why merely restricting to the infrared scaling regime does not repair this particular form. It is not a counterexample to any bound on A.

[Proved distinction] For R, the slow-sector contribution is

    2 F_XX(v_2 tau) - F_XX(v_2 tau/2) -> 0

in this example. Its leading reciprocal-velocity terms cancel. This is a reason to use R or Xi rather than Q for a new upper bound. It does not establish such a bound.

## 4. Crossover-free inequalities that do hold

[Proved, trace perturbation bound] For Hermitian K,V on the same space, min-max eigenvalue monotonicity gives

    exp(-beta ||V||) Tr exp(-beta K)
      <= Tr exp[-beta(K+V)]
      <= exp(beta ||V||) Tr exp(-beta K).                         (10)

One need not assert operator monotonicity of the exponential. Equivalently, interpolate the free energy and bound the expectation of V by its norm.

[Proved] Compare two L-rings to one 2L-ring by removing the two closing bonds and adding the two crossed bonds. In a gauge representative of norm m the difference has norm at most 4m. Equation (10) and gauge optimization give

    |R_L(beta)| <= 4 beta M.                                     (11)

This is the full-spectrum counterpart of the ground-energy gluing estimate already recorded in C98. It is valid at every beta >= 0 and ordinary ring length L >= 3, without a CFT or symmetry assumption.

[Proved] Compare n disjoint L-rings with one nL-ring. There are n removed and n inserted bonds, so

    |g_(nL)(beta) - n g_L(beta)| <= 2 n beta M.

Divide by n and let n tend to infinity at fixed L,beta. This yields

    |Xi_L(beta)| <= 2 beta M.                                    (12)

Both (11) and (12) are crossover-free. They are too weak for the target. For example, (11) at beta = tau L gives |R|/tau <= 4ML, whose right side diverges in the scaling limit. Retaining the entire trace has not improved the seam norm.

[Proved] Splitting an even ring into even and odd bond layers H_even and H_odd, Golden-Thompson followed by Hilbert-Schmidt Cauchy-Schwarz gives

    Z_L(beta)
      <= Tr[e^(-beta H_even) e^(-beta H_odd)]
      <= {Tr e^(-2 beta H_even) Tr e^(-2 beta H_odd)}^(1/2)
      = [tr_(d^2) e^(-2 beta h)]^(L/2).                          (13)

This bound holds for an oriented, nonsymmetric h. It is exact at beta = 0. If h_min is the lowest eigenvalue of the selected bond, then

    g_L(beta) + beta e_inf L
      <= beta L (e_inf - h_min) + O(L)                          (14)

at low temperature. The coefficient e_inf - h_min is nonnegative. If it were zero, H_L >= L h_min would forbid a positive vacuum Casimir deficit. Thus this two-layer comparison has unavoidable extensive slack on a chain with A > 0. Upper bounds on two different traces also cannot be subtracted to obtain an upper bound on their difference.

[Proved direction check] Peierls-Bogoliubov applied to a trial state rho reads

    log Z_L(beta) >= S(rho) - beta Tr(rho H_L).

Restricting an infinite ground state to an interval and closing its seam therefore supplies a lower bound on log Z, not the needed upper bound. Reverse comparisons controlled solely by the seam norm reproduce (11) or (12). This statement excludes these elementary uses of the inequality, not every possible choice of an auxiliary Hamiltonian.

## 5. Why d^2 is not the spatial transfer dimension

[Proved, finite matrix lemma] Suppose a spatial transfer representation has the form

    Z_L(beta) = Tr T_beta^L

with finite dimension D and T_beta positive semidefinite. If its eigenvalues are lambda_j >= 0, Cauchy-Schwarz gives

    [sum_j lambda_j^L]^2 <= D sum_j lambda_j^(2L),
    0 <= R_L(beta) <= log D.                                    (15)

The lower inequality follows by expanding the square. The same conclusion holds for an even power of a Hermitian transfer matrix, with the length convention adjusted if the transfer advances two sites.

[Proved obstruction] Under (1) with A > 0, equation (6) makes r(tau) unbounded as tau tends to infinity. Therefore a temperature-independent D, in particular D <= d^2, is incompatible with (15). Keeping all excited states does not remove this obstruction. Their transfer eigenvalues are what force the effective rank to grow.

[Proved stronger obstruction, positivity unnecessary] Suppose instead that T_beta is an arbitrary finite matrix, its traces reproduce the positive partition functions on all sufficiently large rings, and their thermodynamic limit exists. Let rho_beta be its spectral radius. Then

    limsup_(N -> infinity) |Tr T_beta^N|^(1/N) = rho_beta.

To see this, expand the trace as a sum of eigenvalue powers, including algebraic multiplicity. The upper bound is immediate. Along a subsequence, the phases of the finitely many maximal-modulus eigenvalues simultaneously approach one; their sum then has magnitude bounded below by a positive multiple of rho_beta^N. Lower-modulus terms vanish relative to it. The existence of the pressure limit consequently gives p(beta) = log rho_beta. Hence

    Z_L(beta) <= D rho_beta^L,
    Xi_L(beta) <= log D.                                         (16)

The same argument works on a fixed blocked-length subsequence. This lemma concerns an exact, length-independent transfer matrix, not an L-dependent fit or a bound using an incorrect bulk pressure.

[Proved explicit XX lower bound] C98 records, for the chosen XX normalization,

    E_4 = -sqrt(2),       e_inf = -1/pi,
    delta_4 = 4 e_inf - E_4 = sqrt(2) - 4/pi > 0.

The elementary entropy ceiling implies

    p(beta) <= -beta e_inf + log 2,
    g_4(beta) >= -beta E_4.

Substitution into (16) gives the necessary condition

    log D_beta >= beta (sqrt(2)-4/pi) - 4 log 2.                 (17)

Thus even allowing a different finite transfer matrix at each temperature, its dimension must grow at least exponentially in beta if it is to reproduce every ring trace. At beta = 32 the right side is 1.738579842, already larger than log 4 = 1.386294361. This is an analytic contradiction to D <= d^2 for XX, independent of a CFT limit and independent of spatial transfer positivity. Equation (17) is only a necessary condition; no exact finite-dimensional representation is asserted to exist.

[Proved tensor-network accounting] A two-site imaginary-time gate has operator Schmidt rank at most d^2 across the spatial cut. With n Trotter slices, the cut crosses a number of independent gate legs proportional to n. The uncompressed history space therefore has dimension d^O(n), for example at most d^(2n) in a usual even/odd construction. The local d^2 bound is not a bound on the whole history space. The Trotter limit sends n to infinity.

[Assumed published construction] Pozsgay's QTM discussion, Section 3 and equations (3.10)-(3.11), constructs the transfer matrix on the auxiliary Trotter spaces, not the original one-site space. See [Generalized Gibbs Ensemble for Heisenberg Spin Chains](https://arxiv.org/html/1304.5374v3). Faulmann, Gohmann and Kozlowski analyze the conformal low-temperature part of this spectrum in the infinite-Trotter limit for XXZ, with explicit qualifications on the part of the spectrum controlled. See [Low-temperature spectrum of the quantum transfer matrix of the XXZ chain in the massless regime](https://arxiv.org/abs/2305.06679). Neither result supplies a universal temperature-independent dimension bound.

[Proved scope] A claim of dimension d^2 per physical site could instead mean a total dimension growing like d^(2L). In that interpretation (15) gives R <= 2 L log d and again loses the required scaling limit. A useful new result would have to control the effective rank at separation L, not just the elementary tensor size.

## 6. Reflection positivity and the sign of R

[Proved matrix counterexample] Operator positivity of a gate does not survive a spatial rearrangement of its indices. For h = Z tensor Z, the gate exp(-beta h) is positive. In its spatial realignment, the subspace indexed by diagonal one-site matrix units has matrix

    [[exp(-beta), exp(beta)],
     [exp(beta), exp(-beta)]].

Its eigenvalues are 2 cosh beta and -2 sinh beta. The other two eigenvalues of this realignment are zero. Thus a rotation of the gate cannot itself justify the positive transfer matrix assumed in (15). This example does not exclude a model-specific change of reflection or a blocked representation.

[Assumed published theorem] Jaffe and Janssens, Theorem 5.2, characterize reflection positivity of reflection-invariant spin Hamiltonians by a positivity condition on the matrix of couplings across the reflection plane. Arbitrary Hermitian nearest-neighbor bonds do not automatically satisfy these hypotheses. See [Characterization of Reflection Positivity: Majoranas and Spins](https://arxiv.org/html/1506.04197v2#S5.SS4). Consequently chessboard estimates cannot be invoked for the task's full class merely from locality.

[Proved explicit even-ring sign counterexample] The twisted XX bond is

    h_theta = [exp(i theta) S^+ tensor S^-
                  + exp(-i theta) S^- tensor S^+]/2,
    theta = pi/4.

It has d = 2, M = 1/2 and infinite-chain A = 1. A site-dependent rotation removes theta on the line, but the ring retains flux L theta. The Jordan-Wigner single-particle dispersion is cos(k+theta), with antiperiodic momenta at even fermion number and periodic momenta at odd number.

For L = 4, the even-number grid shifted by theta has energies 0,-1,0,1. The lowest even-number many-body energy is -1. The odd-number grid has energies +1/sqrt(2),-1/sqrt(2),-1/sqrt(2),+1/sqrt(2); its allowed minimum is -1/sqrt(2). Thus E_4 = -1. For L = 8 the total flux is 2 pi, so E_8 = -csc(pi/8). It follows that

    lim_(beta -> infinity) R_4(beta)/beta
      = 2 - csc(pi/8) = -0.61312592975... < 0.                   (18)

Hence R >= 0 is not an all-ring theorem even for an even ring in a critical nearest-neighbor chain. Vacuum-compatible lengths here include multiples of eight. This does not affect extraction along that sequence, and it does not contradict the requested upper bound.

[Proved scope] Reflection inequalities can still be useful in their proper class. Positivity alone gives no upper effective-rank budget for an infinite-dimensional transfer matrix. C98 already shows why qualitative reflection energy inequalities plus an O(M) seam ceiling do not control the Casimir coefficient. That argument is not repeated here.

## 7. Massive spectators do not contribute to the extracted slope

[Proved, exact factorization] For decoupled critical and massive systems, log Z, R, Q and Xi add. It suffices to examine the massive factor.

[Assumed massive-sector input] Suppose that, along compatible rings, the massive factor has a bounded vacuum multiplicity g_m, a uniform gap Delta > 0 above that vacuum multiplet, and

    E_m(L) = e_m L + O(exp(-L/xi)).

Exponential splittings inside a finite vacuum multiplet can be included without changing the result. These hypotheses are more than a bare assertion that the infinite system is gapped. In particular, a frustrated incompatible ring is not covered.

[Proved under this input] If the spectator local dimension is d_m, its excited trace obeys

    0 <= Tr_exc exp[-beta(H_m(L)-E_m(L))]
         <= d_m^L exp(-beta Delta).                              (19)

Set beta = tau L. For tau Delta > log d_m, the right side tends to zero. The ground correction beta O(exp(-L/xi)) also tends to zero. For the 2L ring at the same beta one needs tau Delta > 2 log d_m. Therefore, for every sufficiently large fixed tau,

    R_m(L,tau L) -> log g_m,
    lim_(tau -> infinity) lim_(L -> infinity) R_m(L,tau L)/tau = 0. (20)

This proof uses no density-of-states approximation. It is sufficient for the ordered slope limit. It does not claim uniform convergence as Delta tends to zero. A near-critical massive sector may affect preasymptotic torus traces, but not (20) at a fixed Hamiltonian. For coupled sectors the infrared decoupling and massive finite-size behavior remain part of the assumed scaling description.

[Proved exactly for a test spectator] Take the ferromagnetic classical Ising chain h = -J Z tensor Z, J > 0. Its periodic trace is

    Z_m(L,beta) = (2 cosh beta J)^L + (2 sinh beta J)^L.

With x = tanh(beta J)^L in [0,1],

    R_m = log[(1+x)^2/(1+x^2)],
    0 <= R_m <= log 2.                                          (21)

The bounds follow from (1+x)^2 >= 1+x^2 and (1+x)^2 <= 2(1+x^2). For beta = tau L and fixed J,tau > 0, x tends to one, so R_m tends to log 2. Its contribution to the slope is zero. Tensoring this factor onto XX adds no A. This is a full-spectrum check, with the domain-wall states included, rather than a ground-energy subtraction.

[Proved caution] Suppression of massive modes is an infrared property of the limit, not an all-scale statement that massive modes never enter R. That distinction is required by the round-six update. A local modulation response need not have the suppression in (19)-(20).

## 8. The missing inequality and the bound it would imply

[Proved conditional implication] Suppose one proves for every chain in the task's class and all L,beta

    R_L(beta) <= K(d) + k M f(d) beta/L,                          (22)

with finite, explicit K(d), k and f(d), independent of crossover lengths. Taking the ordered limits in (6) gives

    pi A/4 <= k M f(d),
    A <= (4k/pi) M f(d).                                        (23)

Thus C = 4k/pi is explicit if (22) is established. More generally a remainder whose scaling-limit value is o(tau) is harmless. The corresponding Xi estimate with slope k would give C = 6k/pi. These are sufficient conditions, not newly proved lattice bounds.

[Proved restriction] An all-scale version of the sharp conjecture with f(d) = log_2 d and k = pi^2/4 fails, even allowing the constant K(d). At fixed L, beta -> infinity removes K(d) and would imply

    E_(2L)-2E_L <= pi^2 M log_2 d/(4L).

C98's Heisenberg L = 4 countercheck gives

    E_8 - 2E_4 = 0.348906591063...,
    pi^2 M/(4L) = pi^2/32 = 0.308425137534....

A larger coefficient, a length-dependent slope correction vanishing as L grows, or a different comparison is needed. This does not disprove the asymptotic conjecture or the relaxed target with any explicit C and f.

[Open] No estimate of the form (22), or another crossover-free estimate sufficient for (23), was proved here. Equations (11)-(14) are the unconditional bounds obtained, and none implies finite C M f(d). The finite transfer-rank assumption is false; symmetry-specific reflection theorems do not cover the class; and modular covariance alone supplies no microscopic effective-rank budget. This is a limitation of the methods investigated, not a no-go theorem for every partition-function proof.
