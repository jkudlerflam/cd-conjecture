# C106 proof and limitations

## 1. Setting and status

[Assumed] The microscopic hypotheses, the definition of the gauge-invariant bond scale M, and the infrared CFT identification are those of ../../specs/_bound_background.md, including both updates. Lattice spacing, hbar, and Boltzmann's constant are one. Logarithms without a subscript are natural. Write

    A = sum_a c_a v_a,
    N(u) = sum_{a: v_a >= u} c_a.

There are finitely many sectors in any individual chain, c_a > 0, and v_a > 0. For unequal left and right velocities, the entropy estimates below still hold with the average velocity, because (1/v_R + 1/v_L)/2 >= 1/[(v_R + v_L)/2]. The explicit example has equal left and right speeds.

[Assumed inherited dynamical identification] Each propagating sector has a local real-time realization to which the Lieb-Robinson estimate applies. C69 then supplies the conservative cap

    v_a <= V = C_LR M,       C_LR = 6 e.

The results below do not improve that constant. M = 0 gives trivial bulk dynamics and no positive-velocity sector, so take M > 0.

[Proved main negative result] There is a one-site translation-invariant nearest-neighbor chain with d = 4 for which v_a/M can be held at a fixed positive value while the temperature interval of fixed-fraction infrared entropy matching shrinks to zero. This disproves the proposed uniform crossover lemma, including a lower bound proportional to v_a^2/M. It does not disprove a bound on A or N(u).

[Open] No universal interacting bound on N(u) is established here. Sections 2, 3, and 7 state sufficient additional assumptions and prove their consequences. Section 4 gives the microscopic counterexample to the natural thermal assumption. All numerical work is a small check of explicit formulas.

## 2. Counting identities and the logarithm

[Proved] The layer-cake identity is

    A = integral_0^V N(u) du.                         (2.1)

For each sector, integral_0^V 1_{u <= v_a} du = v_a. Multiply by c_a and sum. Values at the finitely many endpoints do not affect the integral.

[Proved] If, for 0 < delta <= 1,

    N(u) <= g(d) (M/u)^(1-delta),

then

    A <= [C_LR^delta/delta] M g(d).                  (2.2)

This is the direct integral of u^(delta-1). If N(u) <= g(d) for every u > 0, there is no logarithmic loss: A <= C_LR M g(d). Such a uniform estimate also implies sum_a c_a <= g(d), by taking u below the smallest velocity. It would conflict with the conditional fixed-d constructions in the supplied background.

[Proved] The borderline tail

    N(u) <= g(d) M/u

only implies, for a finite set of sectors with v_min > 0,

    A <= M g(d) [1 + log(V/v_min)].                 (2.3)

On (0,v_min), use N(u) = N(v_min) <= M g(d)/v_min; integrate M g(d)/u above v_min. There is no uniform bound when v_min is allowed to approach zero. In contrast, a logarithmic tail N(u) <= g(d)[1 + log(V/u)] is integrable and gives A <= 2 V g(d). These are distinct meanings of a logarithmic counting estimate.

## 3. What the entropy route would need

### 3.1 A sufficient window assumption

[Assumed window lemma E_p] Fix constants alpha > 0, kappa > 0, and 1 <= p < 2. For every velocity window

    I(w) = {a: w <= v_a < 2w},
    T_w = kappa M (w/M)^p,

assume the thermodynamic entropy density satisfies

    s(T_w) >= (alpha pi T_w/3) sum_{a in I(w)} c_a/v_a.   (3.1)

The upper endpoint of the last window can be included. It is enough to impose (3.1) on the windows used in the proof. This is a simultaneous lower bound on the selected sectors in the actual microscopic Gibbs state. Individual informal crossover assignments do not by themselves imply it. Other sectors can interact with the selected ones outside their common infrared regime, and there is no exact sector decomposition of microscopic entropy in general.

[Proved] Write C(w) = sum_{a in I(w)} c_a. Since v_a < 2w,

    sum_{a in I(w)} c_a/v_a >= C(w)/(2w).

For an L-site state, S <= L log d. Divide by L and take the thermodynamic limit. Thus s(T_w) <= log d, and (3.1) gives

    C(w) <= K (M/w)^(p-1),
    K = 6 log d/(alpha pi kappa).                   (3.2)

Using a velocity window, rather than the entire set v_a >= w, gains one power of w/M relative to the estimate that replaces every velocity by V.

[Proved conditional tail] For 1 < p < 2, cover [u,V] by [2^j u,2^(j+1)u), with the final window truncated. Sum (3.2):

    N(u) <= K (M/u)^(p-1) sum_{j>=0} 2^[-j(p-1)]
         = K (M/u)^(p-1)/(1-2^(1-p)).               (3.3)

This has the requested exponent 1-delta = p-1, or delta = 2-p > 0. For p = 1, the same argument gives

    N(u) <= K [1 + log_2(V/u)],                     (3.4)

with an integer ceiling available at the endpoints. The inequality displayed follows because the number of windows is floor(log_2(V/u))+1. A logarithmic divergence is compatible with unbounded total central charge.

[Proved conditional A bound] There is a convenient direct constant that also includes p = 1. Let

    w_j = V/2^(j+1),
    I_j = {a: w_j < v_a <= 2w_j},       j = 0,1,...,
    C_j = sum_{a in I_j} c_a.

The windows partition all positive velocities. By (3.2),

    A <= 2 sum_j w_j C_j
      <= 2 K M^(p-1) sum_j w_j^(2-p)
       = [2^(p-1) K C_LR^(2-p)/(1-2^(p-2))] M.    (3.5)

At p = 1 this is A <= 12 C_LR M log d/(alpha pi kappa). The sum converges because p < 2. This establishes the implication of E_p, not E_p itself.

### 3.2 The quadratic borderline can work with entropy increments

[Assumed increment lemma Q] On the same velocity windows define

    T_j = kappa w_j^2/M.

Suppose that the actual entropy increments obey

    s(T_j)-s(T_j/2)
      >= (alpha pi T_j/6) sum_{a in I_j} c_a/v_a.   (3.6)

For example, it suffices that the specific heat per site obey

    C_V(T)/T >= (alpha pi/3) sum_{a in I_j} c_a/v_a

throughout T in [T_j/2,T_j]. This is a stronger and different assertion from the one-temperature lower estimate in (3.1).

[Proved conditional theorem] Since v_a <= 2w_j, (3.6) implies

    s(T_j)-s(T_j/2) >= alpha pi kappa w_j C_j/(12M).

The temperature intervals are disjoint, since T_(j+1) = T_j/4. Canonical entropy is nondecreasing in temperature: in finite volume its derivative is Var(H)/(L T^3) >= 0. The same monotonicity holds in the thermodynamic limit. Hence the sum of these disjoint entropy increments is at most log d, even if there is residual zero-temperature entropy. Therefore

    sum_j w_j C_j <= 12 M log d/(alpha pi kappa),
    A <= 24 M log d/(alpha pi kappa).               (3.7)

One can formulate the proof entirely with finite increments, without differentiability of the thermodynamic entropy. The quadratic crossover exponent is sufficient under Q, although it is only borderline under E_2.

[Proved limitation] Subtracting two lower bounds on s(T) does not give a lower bound on their difference. Thus Q does not follow from E_2 by subtraction. Nor does applying entropy capacity independently to each velocity window allow one to add the resulting capacities. In (3.7), disjoint increments spend the single total entropy budget only once.

## 4. An explicit fast-sector crossover counterexample

### 4.1 Microscopic Hamiltonian

[Definition] At each site put two complex fermion modes, with local Fock dimension d = 4. Let f_x be their two-component column and define

    H_epsilon = sum_x [f_x^dagger B f_(x+1)
                        + f_(x+1)^dagger B^dagger f_x
                        + f_x^dagger X f_x],
    B = -i t sigma_z,
    X = epsilon sigma_x - a epsilon 1,
    t > 0,       a > 1,       epsilon > 0.          (4.1)

All terms are even in fermions. Ordering the two modes within each site and applying the Jordan-Wigner map gives a Hermitian nearest-neighbor bond on C^4 tensor C^4. The site parity string in a hopping term lies within those two sites, so the range remains one bond and translation by one four-state site is a symmetry. The thermodynamic calculations do not depend on the boundary parity twist on a finite ring.

[Proved] With f_x proportional to exp(ikx), the Bloch matrix and bands are

    h_epsilon(k) = 2t sin(k) sigma_z + epsilon sigma_x - a epsilon 1,
    E_+(k) = -a epsilon + sqrt(4t^2 sin^2(k)+epsilon^2),
    E_-(k) = -a epsilon - sqrt(4t^2 sin^2(k)+epsilon^2).  (4.2)

Assume epsilon sqrt(a^2-1) < 2t. The lower band is strictly negative. The upper band has four simple zeros per Brillouin zone, determined by

    |sin k_F| = epsilon sqrt(a^2-1)/(2t).

The magnitudes of the four derivatives are equal:

    v_epsilon = 2t sqrt(1-a^(-2))
                   sqrt[1-epsilon^2(a^2-1)/(4t^2)].  (4.3)

Each pair of opposite-moving complex Fermi modes gives a c = 1 CFT. Thus this chain has two c = 1 sectors, both of speed v_epsilon, and

    A_epsilon = 2 v_epsilon,
    B_IR := sum_a c_a/v_a = 2/v_epsilon.             (4.4)

This is the elementary free-fermion infrared limit at fixed epsilon > 0. The limit epsilon -> 0 is taken only after defining the data of each chain. At epsilon = 0 the two original orbitals each have speed 2t. The noncommutation of the infrared and epsilon limits is the mechanism.

### 4.2 The gauge-optimized bond norm

[Proved] Split the on-site term equally between the adjacent bonds. On one bond the four-mode one-particle matrix is

    K_b = [ X/2     B   ]
          [ B^dag   X/2 ]
        = -a epsilon/2 + Q,
    Q = t tau_y tensor sigma_z + (epsilon/2) 1 tensor sigma_x.

The two displayed nonconstant matrices anticommute, and each squares to the identity times its coefficient squared. Consequently

    Q^2 = E_b^2 1,
    E_b = sqrt(t^2+epsilon^2/4).

The one-particle eigenvalues are -a epsilon/2 + E_b twice and -a epsilon/2 - E_b twice. For a epsilon < 2 E_b, there are two positive and two negative eigenvalues. The maximum many-body energy fills the two positive modes; the minimum fills the two negative modes. Their energies are

    lambda_max = -a epsilon + 2 E_b,
    lambda_min = -a epsilon - 2 E_b.

The bare bond half-spread is 2 E_b.

[Proved gauge lower bound] The one-particle projectors of the two extremal Slater states are

    P_+ = (1+Q/E_b)/2,
    P_- = (1-Q/E_b)/2.

Both diagonal site blocks of each projector are identical:

    (P_+)_11 = (P_+)_22 = 1/2 + epsilon sigma_x/(4E_b),
    (P_-)_11 = (P_-)_22 = 1/2 - epsilon sigma_x/(4E_b).

A number-conserving Gaussian reduced state is fixed by its one-particle covariance matrix, so each extremal state has equal left and right one-site marginals. Its expectation of o tensor 1 - 1 tensor o is zero for every Hermitian one-site o. An added scalar cancels between the two extremal expectations. Their energy difference, 4 E_b, is therefore unchanged by every allowed bond gauge. The spectral spread of any gauge representative is at least this difference. The bare representative attains it, proving

    M_epsilon = 2 E_b = sqrt(4t^2+epsilon^2).       (4.5)

The Gaussian marginal statement is also valid in the site tensor-product representation. These extremal states have definite fermion number; odd one-site moments vanish, and the even local moments agree with the fermionic reduced states.

### 4.3 Holding v/M fixed, rather than just bounded below

[Proved] Fix any r in (0,1). We can choose a = a(epsilon) so that

    v_epsilon/M_epsilon = r                       (4.6)

for every sufficiently small epsilon > 0. To see this without a limiting estimate, put x = a^2-1 and D = 4t^2+epsilon^2. Equations (4.3) and (4.5) give

    v_epsilon^2 = (4t^2 x - epsilon^2 x^2)/(1+x).

Equation (4.6) is therefore the quadratic equation

    epsilon^2 x^2 + (r^2 D-4t^2)x + r^2 D = 0.      (4.7)

At epsilon = 0 it has the simple positive root

    x_0 = r^2/(1-r^2),       a_0 = 1/sqrt(1-r^2).

The derivative with respect to x at that root is -4t^2(1-r^2), which is nonzero. Thus there is a positive smooth root x(epsilon) near x_0. Both inequalities used in Sections 4.1 and 4.2 hold for small epsilon. This proves exact fixed r, not merely convergence to r.

[Proved normalization] Multiplying each H_epsilon by 1/M_epsilon fixes M = 1 and v = r. Temperatures scale by the same factor. Hence an arbitrarily small crossover in this family cannot be attributed to an overall shrinking Hamiltonian scale.

### 4.4 Entropy in the two temperature regimes

[Definition] For a single fermionic energy level with E/T = y, its entropy is

    eta(y) = log(1+exp(-|y|)) + |y|/(1+exp(|y|)).

It is even, nonnegative, and decays exponentially up to a linear factor. The thermodynamic entropy of (4.1) is exactly

    s_epsilon(T) = integral_(-pi)^pi dk/(2pi)
          [eta(E_+(k)/T) + eta(E_-(k)/T)].           (4.8)

[Proved normalization integral] For y > 0, integrate the convergent exponential series for the two terms. Each integral is sum_(n>=1) (-1)^(n+1)/n^2 = pi^2/12. Thus

    integral_0^infinity eta(y) dy = pi^2/6,
    integral_(-infinity)^infinity eta(y) dy = pi^2/3.  (4.9)

[Proved fixed-epsilon infrared limit] At any fixed epsilon > 0, localize the integral to the four simple zeros of E_+. In each neighborhood use y = E_+(k)/T, so dk = T dy/(v_epsilon+o(1)). The complement and E_- have a nonzero gap from zero energy and give exponentially small contributions. Equation (4.9) gives

    lim_(T->0) s_epsilon(T)/T
       = 4 [1/(2pi v_epsilon)] [pi^2/3]
       = 2pi/(3 v_epsilon).                       (4.10)

This establishes the CFT coefficient directly from the lattice spectrum.

[Proved intermediate-temperature limit] Now let epsilon/T -> 0 and T/t -> 0, keeping a(epsilon) bounded. The only contributing neighborhoods are those of sin k = 0, namely k = 0 and pi on the Brillouin circle. In either neighborhood put q = 2t(k-k_0)/T. Then

    E_+(k)/T -> |q|,
    E_-(k)/T -> -|q|,
    dk = [T/(2t)] dq.

There are two neighborhoods and two bands. Equation (4.9) gives

    s_epsilon(T)/T ->
       2 [1/(4pi t)] integral_R 2 eta(q) dq
       = pi/(3t).                                (4.11)

Here is a direct control of this rescaling limit. Outside fixed small neighborhoods of the two zeros, |sin k| is bounded below and both level entropies are exponentially small in t/T. Within the neighborhoods, |sin(k_0+z)| >= b|z| for some b > 0. Since epsilon/T -> 0 and a stays bounded, for sufficiently large |q| both |E_+|/T and |E_-|/T are bounded below by b|q| minus a fixed constant. The exponentially decreasing kernel then provides an integrable bound independent of the small parameters. On each finite q interval the rescaled energies converge uniformly. Dominated convergence proves (4.11).

[Proved ratio] Combining (4.6), M_epsilon -> 2t, and (4.11),

    s_epsilon(T) / [2pi T/(3 v_epsilon)] -> r

whenever epsilon/T -> 0 and T/t -> 0.              (4.12)

The infrared entropy prediction overestimates the intermediate-temperature entropy by an arbitrarily large fixed factor if r is chosen small. The speed v/M = r remains fixed and strictly positive while epsilon is taken to zero.

### 4.5 The crossover lemma is false

[Definition] For fixed 0 < alpha < 1, define the one-sided thermal matching endpoint

    T_alpha(H) = sup {tau > 0:
        s_H(T) >= alpha (pi T/3) sum_a c_a/v_a
        for every 0 < T <= tau}.                  (4.13)

The set is nonempty for each member of this example by (4.10). Requiring a two-sided relative-error CFT approximation is stronger, so any failure of this one-sided definition also obstructs that requirement.

[Proved counterexample] Choose 0 < r < alpha and use the exactly fixed-r family. Set T_epsilon = sqrt(epsilon t). Equation (4.12) shows that, for all sufficiently small epsilon,

    s_epsilon(T_epsilon)
        < alpha 2pi T_epsilon/(3 v_epsilon).

Therefore

    T_alpha(H_epsilon) <= sqrt(epsilon t),
    T_alpha(H_epsilon)/M_epsilon -> 0,             (4.14)

while d = 4 and v_epsilon/M_epsilon = r are unchanged. There is no strictly positive universal function F_alpha(d,r) such that

    T_alpha(H) >= M F_alpha(d,v/M)

for all chains in this setting. In particular, for any fixed p and any positive kappa depending only on d and alpha, the bound T_alpha >= kappa M(v/M)^p fails.

[Proved stronger scaling statement] One can replace the square-root upper estimate in (4.14) by O(epsilon). Fix K > 0 and set T = K epsilon. Near each of the two valleys use q = 2t(k-k_0)/epsilon. The same dominated-convergence argument yields

    lim_(epsilon->0) s_epsilon(K epsilon)/(K epsilon)
       = 1/(2pi t K) integral_R dq
          {eta[(sqrt(1+q^2)-a_0)/K]
           + eta[(sqrt(1+q^2)+a_0)/K]}.

Next send K -> infinity, rescale q = K y, and apply (4.9). The last expression tends to pi/(3t). Thus for r < alpha there is a finite K, depending on r and alpha but not epsilon, for which the entropy inequality already fails at T = K epsilon for every sufficiently small epsilon. This proves T_alpha <= K epsilon. No numerical choice of K is needed for the proof.

[Proved scope] This disproves a uniform low-temperature window. It does not exclude an isolated matching temperature elsewhere, an entropy estimate with a sufficiently small velocity-dependent prefactor, or another observable with no thermal window. In this example all infrared sectors have the same speed, so failure of simultaneous matching is not caused by combining widely different speeds. If a lemma only lower-bounds the total entropy by one of the two equal sectors, choose r < alpha/2 to obtain the same obstruction.

### 4.6 Why locality and a spectral derivative bound do not repair it

[Proved] The Bloch derivatives obey

    ||h'_epsilon(k)|| <= 2t,
    ||h''_epsilon(k)|| <= 2t.

Nevertheless, for E_+(k) = -a epsilon + sqrt(epsilon^2+4t^2 sin^2 k),

    E_+''(0) = 4t^2/epsilon.                       (4.15)

At a Fermi point the divergent term is

    E_+''(k_F) = 4t^2/(a^3 epsilon) + O(epsilon),

at each of the equivalent Fermi points near the two valleys. To verify the first-valley formula, write R = sqrt(epsilon^2+4t^2 sin^2 k). Then

    R'' = 4t^2 cos(2k)/R
           - 16t^4 sin^2(k) cos^2(k)/R^3,

and substitute R(k_F) = a epsilon. The leading terms combine to 4t^2/(a^3 epsilon).

A bound on derivatives of the matrix-valued Bloch Hamiltonian does not bound derivatives of its individual eigenvalues near an avoided crossing. Eigenvector mixing introduces the small band separation in a denominator. Consequently, a scalar-band argument suggesting a linear window of width v^2/M does not generalize even to two free bands. A Lieb-Robinson upper cone also supplies no lower bound on this separation. The example has uniformly bounded microscopic couplings and obeys a uniform Lieb-Robinson estimate throughout.

[Proved consistency] Along the fixed-r family,

    N(u) = 2 for 0 < u <= v_epsilon, and 0 above v_epsilon,
    A/M = 2r < 2.

Thus the example satisfies the already proved free bound A <= 2M. It is a counterexample to the crossover route, not to fast-mode counting or the conjecture. It has no growing number of modes.

### 4.7 The universal versions of E_p and Q also fail

[Proved] For any prescribed positive alpha and kappa, the two-sector family disproves universal E_p for each p in [1,2), and universal Q. This is stronger than merely observing that a low-temperature-window proof of those lemmas fails.

For E_p, first choose a small fixed r, take w = theta v with a fixed theta in (1/2,1), and set T_w = kappa M(theta r)^p. Both sectors lie in I(w). Send epsilon to zero with v/M = r. At this positive temperature the entropy converges to that of two uncoupled bands with speed 2t near their zeros. Then let r tend to zero. Since T_w/t tends to zero, the same rescaling used in (4.11) shows

    s(T_w) / [2pi T_w/(3v)] = r [1+o(1)]

in this order of limits. Choose r small enough that this ratio is below alpha, and then choose a sufficiently small positive epsilon. This gives a member of the microscopic family violating (3.1). Both choices are finite; the limits establish their existence.

For Q, choose r arbitrarily small along values for which the common speed lies strictly inside one of the dyadic windows based at V = C_LR M. Then w_j/M is comparable to r, and T_j/M = kappa (w_j/M)^2 tends to zero as r tends to zero. First send epsilon to zero at fixed r. The entropy increment is asymptotic, as r subsequently tends to zero, to

    [pi/(3t)] (T_j/2).

Its CFT prediction in (3.6), before multiplication by alpha, is

    [2pi/(3v)] (T_j/2).

Their ratio is again r [1+o(1)]. Thus Q fails for sufficiently small fixed r and then sufficiently small positive epsilon. The conditional theorems in Section 3 remain valid for any restricted class in which their premises can be established, but neither premise is a theorem for the full class in the question.

## 5. Why slow convergence alone is not the argument

[Proved distinction] A long RG flow, a small gap in an extra spectator system, or slow logarithmic convergence does not alone disprove a one-sided entropy bound for surviving sectors. The total entropy might remain larger than the infrared expression. Section 4 instead computes the entropy and proves that it is smaller by the factor r in an explicit range of temperatures. This sign is what the counting argument needs.

[Proved distinction] The failure is not a large interaction correction and does not use a marginally irrelevant coupling. Each member is a solvable quadratic chain. The low-energy carrier pockets have width O(epsilon), depth O(epsilon), and a fixed edge slope. The upper and lower bands away from these pockets instead carry the larger speed 2t. The small avoided-crossing gap permits these features with fixed local dimension and norm.

## 6. What finite local operator dimension does and does not count

[Proved linear algebra] The complex vector space of operators on a fixed block of k sites has dimension d^(2k). If a specified linear scaling map from this space onto a space W of independent infrared fields is surjective, then

    dim W <= d^(2k).                               (6.1)

The same conclusion holds if one exhibits a collection of lattice operators whose images are independent. It does not follow from the weaker statement that each field has nonzero overlap with some local observable.

[Proved rank obstruction] For example, one lattice observable may have the scaling expansion

    O -> b_1 T_1 + ... + b_n T_n,

with every b_a nonzero and mutually orthogonal stress tensors T_a. This is one vector in an n-dimensional field space, not a surjection onto that space. Its correlation function can contain n distinct contributions. Resolving these contributions spectrally does not manufacture n independent vectors in the original one-dimensional source space.

[Proved locality qualification] Translating and smearing a local observable over longer distances enlarges its support. Time filtering also enlarges its effective spatial support. To distinguish frequencies separated by Delta omega, a Fourier filter typically uses time of order 1/Delta omega; the Lieb-Robinson estimate confines its tails outside a radius of order M/Delta omega. It is an upper locality bound for a chosen filter, not a guarantee of a uniformly short filter. For cones at frequencies v_a q, Delta omega vanishes as q -> 0 even if the velocities are distinct. A field visible after a scaling limit therefore need not be selectable within a fixed microscopic block. This observation explains the missing step; it is not a theorem excluding every other local construction.

[Proved central-charge obstruction] Counting one stress tensor per sector counts the number of sectors, not N(u). At the chiral vacuum level,

    ||L_(-2)|0>||^2 = <0|[L_2,L_(-2)]|0> = c/2.

The corresponding vector space has dimension one, regardless of the value of c. A dimension bound contains no upper bound on this norm. At any fixed Virasoro descendant level, the number of ordered monomials L_(-n_1)...L_(-n_j)|0>, with n_i >= 2, is bounded by a partition count independent of c. Null vectors can only reduce it. Thus Virasoro symmetry alone supplies no lower bound proportional to c on the number of such fields at fixed level.

[Assumed standard CFT example] A unitary sector need not have a spin-one chiral current at all. The Ising CFT, c = 1/2, has no weight-one current. The claim that each sector requires currents cannot be imposed on the general class in the question. A theorem relating c to a larger set of low-dimension fields would be an additional CFT input, beyond the existence of a stress tensor.

[Proved limitation] These statements do not construct a large-c lattice theory with few low-dimension fields. They show why the stated operator dimension argument has two missing hypotheses: enough independent continuum fields per unit c, and a uniformly bounded microscopic block that spans them.

## 7. A sufficient operator-content lemma, with its required support scale

[Assumed operator lemma O] For each threshold u in (0,V], suppose there is a space W_u of independent infrared fields of the sectors with v_a >= u and constants eta > 0 and k(d,u/M) such that

    dim W_u >= eta N(u),

and the scaling map from operators on a fixed block of k(d,u/M) sites is onto W_u. All quantum numbers needed to identify the fields are included in this map. This assumption entails independent access, not just nonzero overlaps.

[Proved conditional implication] Equation (6.1) immediately gives

    N(u) <= eta^(-1) d^[2k(d,u/M)].                (7.1)

To obtain an integrable tail directly from this estimate, it suffices that

    k(d,u/M) <= K(d) + beta(d) log(V/u),
    gamma(d) := 2 beta(d) log d < 1.               (7.2)

Indeed,

    N(u) <= eta^(-1) d^[2K(d)] (V/u)^gamma(d),
    A <= V d^[2K(d)]/[eta(1-gamma(d))].             (7.3)

Since V = 6 e M, this has the requested form with an explicit constant. A ceiling in k can be absorbed into K(d). If eta depends on d, that dependence is absorbed into f(d).

[Proved implication for weaker support control] A bound k proportional to M/u only gives a bound exponential in M/u through (7.1), which is not integrable at zero. A constant k independent of u would imply a bound on the entire central charge and would encounter the conditional fixed-d constructions in the background. Logarithmic growth with the coefficient in (7.2) is a concrete sufficient middle ground, but it is not supplied by a Lieb-Robinson upper cone.

[Open] Neither the lower bound dim W_u >= eta N(u) nor the support bound (7.2) has been proved for the general interacting class. If one restricts to a complete Gaussian infrared with a faithful microscopic realization of its individual currents, the first condition is natural, but the second and the independence of the images still require proof. The existing free-chain counting and stiffness results in the supplied reports already handle the microscopic quasi-free class by stronger algebraic methods. They are not rederived here.

## 8. Checks, references, and work record

[Numerical] checks.py uses NumPy and SciPy for a 16 by 16 bond diagonalization and one-dimensional entropy quadratures. It verifies the bond half-spread, equality of the extremal one-site marginals, exact fixed v/M = 1/4, the small-temperature CFT coefficient, and failure of the alpha = 1/2 entropy lower estimate at T = 10 epsilon. CHECKS.md records the output. These are small checks of an explicit analytic counterexample, not heavy many-body numerics.

[Proved consistency checks] Units are consistent in T_w = kappa M(w/M)^p and T_j = kappa w_j^2/M. Positive rescaling of H rescales M, v, and T together. The example has c_total = 2 and A/M = 2r, so it respects the inherited free bound. Its entropy tends to 2 log 2 = log 4 at high temperature, by (4.8). The fixed-epsilon infrared limit and the intermediate-temperature limit have different slopes, as required for the counterexample. The geometric sums in Section 3 converge only in their stated exponent ranges.

[Proved source record] The project background, REPORT.md sections 9.68 to 9.75, C69 RESULT.md, C100 RESULT.md, and contrib/fable_bound_freefermion/README.md supply the inherited results and conventions. The counterexample and conditional counting statements above are derived here. A primary-literature lookup located the lattice-operator matching method of Zou, Milsted, and Vidal, [Conformal fields and operator product expansion in critical quantum spin chains](https://arxiv.org/abs/1901.06439). Its abstract describes a variational construction demonstrated on Ising, rather than supplying the uniform operator-support lemma needed here. Their [Conformal data and renormalization group flow in critical quantum spin chains](https://arxiv.org/abs/1710.05397) describes numerical lattice-to-CFT identification and a tricritical-Ising-to-Ising flow. Neither paper is used as a premise for the proof or for a claim that no stronger theorem exists.

[Proved work record] The first script launch used an unavailable python executable; subsequent runs use python3. The first numerical assertion exposed an incorrect expression for v^2 in the root-finding check. It was corrected to (4t^2 x - epsilon^2 x^2)/(1+x), which also appears in (4.7), before the final checks passed. No failed check is used as evidence. No files outside work_C106 were changed. No network command, cluster job, commit, or push was run. The project STATE.md is outside the authorized folder, so the resumable handoff is in local STATE.md and PROGRESS.md.
