# C109 proofs and missing inputs

## 1. Scope and conventions

[Assumed] The Hamiltonian is H_L=sum_i h_(i,i+1) on a periodic chain with local dimension d, one-site translation invariance and nearest-neighbor interactions. The bond scale is the gauge-optimized half spectral spread M from _bound_background.md. The infrared contains finitely many unitary CFT factors, with c_a the nonchiral Virasoro central charge and v_a>0. Put

    A_neutral = sum_(neutral a) c_a v_a,
    B_p = sum_a c_a v_a^p.

Lengths are measured in physical lattice sites and hbar=1. Equal left and right velocities are used for the displayed correlators. With unequal velocities, replace each c v^p by c(v_R^p+v_L^p)/2. The target uses p=1.

[Assumed, stronger matching where stated] Using a Lieb-Robinson bound to constrain a CFT velocity requires a real-time scaling limit of local observables, with nonzero local stress response. A finite-size Casimir formula alone is not this assumption. The conditional velocity statements below explicitly use this dynamical matching. This is the usual physical meaning of a local CFT infrared, but it must be part of a rigorous theorem.

[Proved] For every epsilon>0 choose a centered gauge representative g with ||g||<=M+epsilon. Coboundaries telescope on rings and in the infinite-volume local dynamics. All lattice estimates below hold with M+epsilon; taking epsilon down to zero gives the displayed constants. No assumption that the optimizing gauge is attained is needed.

[Inherited] REPORT 9.68 through 9.77 already establish the failures of a bond-gradient norm bound, cutoff-independent thermal matching, elementary ring gluing and a general charge-stiffness counting theorem. C104 bounds full WZW factors associated to microscopic semisimple charges by (4 pi/3)d(d^2-1)M, conditional on its stiffness dictionary. An independent coset stress tensor is not included by that result. C100 already distinguishes the spatial and temporal current correlators used in Section 5 below.

## 2. What a discrete symmetry seam actually bounds

### 2.1 The lattice estimate

[Assumed] Let u(g) be an on-site unitary symmetry with an invariant bond representative, [h,u(g) tensor u(g)]=0. For a finite or compact group, average a near-optimal gauge representative under u tensor u. Convexity of the operator norm preserves ||g||<=M+epsilon, and the averaged coboundary comes from an invariant one-site operator. Thus the symmetric gauge class has the same infimum M.

[Proved] Insert the symmetry seam on the bond (L,1), replacing g by

    g^(g) = (1 tensor u(g)) g (1 tensor u(g)^dagger).

Only one term changes. The variational principle and unitary invariance of the norm give

    |delta E_g(L)| := |E_g(L)-E_1(L)|
       <= ||g^(g)-g|| <= 2M.                            (2.1)

An invariant coboundary still telescopes in the twisted ring, so (2.1) is compatible with the definition of M. A general defect without an on-site symmetry construction need not admit this estimate or a tension-free seam.

### 2.2 Matching the seam to a CFT

[Assumed] The seam has the topological symmetry-defect scaling limit, with no residual defect tension, and the chosen periodic lengths carry the untwisted conformal vacuum. Suppose that its lowest state has weights (h_a(g), hbar_a(g)) in each factor. Set x_a(g)=h_a(g)+hbar_a(g). For a product action with no correlated projection excluding the tensor product of twisted vacua,

    delta E_g(L) = (2 pi/L) sum_a v_a x_a(g) + o(L^-1).
                                                               (2.2)

If the symmetry permutes factors or a projection changes the ground state, use the actual twisted vacuum energy instead. Do not automatically add independent dimensions in that case. For unequal chiral velocities the summand is v_R h+v_L hbar.

[Proved from (2.1)-(2.2)] The direct comparison is

    sum_a v_a x_a(g) <= M L/pi + o(1).                  (2.3)

The right side diverges in the limit needed to identify the CFT coefficient. Positivity of x_a(g) does not remove L. The missing lattice estimate is a bound on the rescaled energy shift, for example

    limsup_(L->infinity) L delta E_g(L) <= K(d) M,      (S_g)

on compatible vacuum lengths. It is sufficient to prove this limit bound; a bound for every small L is unnecessary. No proof of (S_g) is supplied by the seam norm.

[Proved, limitation of the scalar estimates] For any b>0, the sequence delta_b(L)=min(2M,b/L) satisfies the same norm bound and has lim L delta_b(L)=b. These sequences are not asserted to be lattice spectra. They show that (2.1) and existence of the limit alone do not constrain its value.

[Proved, invisibility example] Tensor a gapped spin with an independent critical Ising spin at each site, so d=4. Let the microscopic Z_2 act only on the gapped spin and take its Hamiltonian to be a sum of Z fields. The twisted bond then equals the original bond, delta E_g(L)=0 for every L, while the Ising factor contributes cv>0. Thus availability of a microscopic discrete symmetry alone does not ensure detection of all neutral stress tensors. An emergent continuum symmetry also need not have a one-bond lattice seam.

### 2.3 Parafermion dimension accounting, conditional on the seam dictionary

[Assumed, standard parafermion data] For the usual Z_N parafermion theory,

    c_N = 2(N-1)/(N+2),
    x_k = k(N-k)/[N(N+2)],   0<=k<N.                  (2.4)

The disorder operator has chiral weight x_k/2. The formulas follow from the order/disorder spectrum in [Fateev and Pugai, Sections 2 and 3](https://arxiv.org/html/0909.3347). For the application here, additionally assume that the lattice g^k seam has this disorder operator as its twisted vacuum. Naming the infrared parafermion theory does not by itself verify this seam identification for an arbitrary lattice realization.

[Proved arithmetic] The generator twist gives c_N/x_1=2N. A stronger twist removes that growth. For k=floor(N/2),

    c_N/x_k = 2N(N-1)/[k(N-k)] <= 8.                 (2.5)

For N=2m this is 8-4/m, and for N=2m+1 it is 8-4/(m+1). Also,

    sum_(k=0)^(N-1) k(N-k) = N(N^2-1)/6,
    x_average = (N^2-1)/[6N(N+2)],
    c_N/x_average = 12N/(N+1) < 12.                 (2.6)

[Proved, conditional implication] Consider a finite group whose action on each parafermion factor surjects onto that factor's Z_N, with product twisted-vacuum matching as in (2.2). Uniform group averaging is then uniform on each cyclic image. Hence, with X_g=sum_a v_a x_a(g),

    A_neutral <= 12 |G|^-1 sum_g X_g.

If the still-unproved lattice input is

    limsup L |G|^-1 sum_g delta E_g(L) <= K(d) M,

then

    A_neutral <= (6/pi) K(d) M.                       (2.7)

This separates two questions. For these parafermion products, twist dimensions can cover the central charges with a constant independent of the parafermion orders and the number of factors. The missing step is the 1/L lattice seam budget, not a lower bound on the generator's dimension. Group averages here are finite arithmetic means, not statistical assertions.

[Proved checks] N=2 gives c=1/2, x_1=1/8 and c/x_1=4. N=3 gives c=4/5, x_1=x_2=2/15 and c/x_1=6. These are the Ising and three-state Potts cases. General minimal models and general cosets are not covered merely by (2.4).

### 2.4 Why distributing a discrete twist is not a stiffness proof

[Proved structural distinction] Choose a Hermitian logarithm q of u(g). A discrete symmetry gives invariance of h under u(g) tensor u(g), but does not give

    [h,q tensor 1+1 tensor q]=0.

Consequently a site-dependent interpolation exp(i theta_j q) generally rotates the bulk interaction even when neighboring theta_j differ by only O(1/L). The absolute rotation angle matters. The continuous-charge gauge transformation that leaves only angle differences on the bonds is unavailable.

[Proved example] For the transverse Ising model take u=Z and q=Z/2. The XX bond is invariant under Z tensor Z, but [XX,Z tensor 1+1 tensor Z] is nonzero. A uniform intermediate angle changes XX into a linear combination involving YY and XY+YX. It is not a symmetry. Replacing a single Z_2 seam by a smooth U(1) flux therefore changes the Hamiltonian in the bulk.

[Heuristic limitation] This does not rule out a new variational construction using discrete defects. It rules out borrowing the continuous Kohn-stiffness argument without supplying the missing continuous symmetry or another cancellation mechanism.

## 3. A universal discrete twist is precisely ring gluing

[Proved lattice identity] Take n independent identical copies of the chain, with local Hilbert space dimension d^n, and let the on-site symmetry cyclically permute the copies. The untwisted Hamiltonian has ground energy n E_0(L). Its twisted boundary bonds join the end of copy s to the beginning of copy s+1. Relabeling the nL tensor factors gives one ring of length nL, with the same bond and physical orientation. Therefore

    delta E_cycle(L) = E_0(nL)-n E_0(L).               (3.1)

This is an equality of finite-dimensional Hamiltonians up to a tensor-factor permutation. It uses neither conformal invariance nor an orbifold projection. There is no assertion that the same relabeling is an on-site unitary in the original replica coordinate.

[Proved norm estimate] There are n changed seam bonds, each bounded by 2M. Thus

    |delta E_cycle(L)| <= 2n M.                       (3.2)

[Proved from the assumed Casimir asymptotics] For fixed integer n and compatible vacuum lengths,

    E_0(nL)-nE_0(L)
       = [pi A/(6L)](n-1/n)+o(L^-1).                (3.3)

Comparing with the CFT seam formula gives the full twist dimension in each sector,

    x_cycle,a = (c_a/12)(n-1/n).

This also follows from the Schwarzian under the n-sheeted uniformization. The replica construction in [Calabrese and Cardy, Section III](https://arxiv.org/html/hep-th/0405152) uses a single-sheet weight c(1-1/n^2)/24; multiplying by n yields the replica chiral weight c(n-1/n)/24 used here.

[Proved consequence] Replica twists see every neutral sector without any microscopic internal symmetry in the original chain. They do not give a new lattice bound. In particular,

    lim L [E_0(2L)-2E_0(L)] = pi A/4.                (3.4)

The desired A<=C M f(d) is equivalent, under the Casimir assumption, to the upper bound C pi M f(d)/4 on this limit. This is C98's doubling problem. For larger fixed n, the factor n in the norm estimate cancels the leading n in the twist dimension and leaves the same factor L. Taking n large before verifying a uniform Casimir remainder is not a legitimate improvement.

## 4. The crossover-free bounds that do follow from locality

### 4.1 An explicit velocity ceiling

[Proved lattice estimate] Write m=M+epsilon for the norm of a centered representative. For A and B supported on single sites at distance r>0, use the standard commutator recursion

    C_B(X,t) <= C_B(X,0)
       +2 sum_(Z crossing X) ||g_Z|| integral_0^|t| C_B(Z,s) ds,
    C_B(X,t)=sup_(A in algebra(X)) ||[A(t),B]||/||A||.

The recursion is [Nachtergaele and Sims, equations (2.30)-(2.32)](https://arxiv.org/html/1004.2086). An initial singleton has two crossing bonds. A bond also has two crossing nearest-neighbor bonds, because its own internal interaction is removed in the local interaction picture. Each step has norm at most m. A path reaching the other site needs at least r bonds. Consequently

    ||[A(t),B]||
       <= 2||A||||B|| sum_(n>=r) (4m|t|)^n/n!
       <= 2||A||||B|| exp[-mu r+4m exp(mu)|t|].       (4.1)

The second inequality follows by multiplying the n-th term by exp(mu(n-r))>=1 and extending the sum to n>=0. Minimize exp(mu)/mu over mu>0, obtaining mu=1. Thus an admissible velocity is

    v_LR = 4 e M.                                    (4.2)

The same exponential cone holds for fixed finite supports with support-dependent prefactors and bounded distance shifts. No local-dimension factor enters this estimate. The constant is deliberately not claimed optimal.

[Proved conditional implication] Under the dynamical matching in Section 1, each v_a<=v_LR. To see the mechanism, suppose v_a>v_LR and smear a nonzero sector stress commutator in a small scaled neighborhood of x=v_a t outside the LR cone. The lattice commutator is exponentially small in the scaling length by (4.1), including after polynomial field renormalization. The CFT stress commutator has a nonzero distribution supported on its own light ray, with coefficient proportional to c_a>0. A test function with the appropriate nonzero derivative detects it. These two scaling limits contradict each other. A static Casimir assumption alone does not justify this step.

### 4.2 Single factors and a known number of factors

[Assumed standard CFT classification; proved bounds] A unitary Virasoro minimal model has c=1-6/[m(m+1)]<1, m>=3. Therefore an individual such factor, in any chain satisfying dynamical matching, obeys

    c v < 4 e M.                                     (4.3)

In particular an Ising factor obeys cv<=2 e M, a tricritical Ising factor cv<=(14e/5)M, and a three-state Potts factor cv<=(16e/5)M. A standard Z_N parafermion factor has c_N<2, so

    c_N v < 8 e M.                                   (4.4)

These estimates use neither a microscopic symmetry nor a lower bound on v. The standard central-charge formulas are reviewed in [Ginsparg, Applied Conformal Field Theory](https://arxiv.org/abs/hep-th/9108028); the parafermion data also appear in the source for (2.4).

[Proved scope] For r_min minimal-model factors and r_pf other standard parafermion factors, count disjoint factors once and obtain

    A_neutral <= 4 e M (r_min+2r_pf).                 (4.5)

This is crossover-free, but it is not a function of d until the number or total central charge of the counted factors is bounded in terms of d. Such a multiplicity theorem is not assumed or proved. In particular, identifying each infrared factor with a separate tensor factor of one microscopic site is invalid. The background explicitly warns that central-charge multiplicity can be hidden in slow modes.

[Inherited proved subclass] For a full quasi-free chain, including neutral Majorana factors, the campaign already proves A<=2M, stronger than (4.3)-(4.5). The present velocity argument is useful for interacting realizations of a specified bounded-c factor. It does not extend the free theorem to interacting products with uncontrolled multiplicity.

## 5. Vacuum current correlations and the two missing bounds

### 5.1 Normalization and velocity powers

[Assumed CFT stress normalization; proved calculation] Let

    <T_a(z)T_a(0)> = c_a/(2 z^4),
    h_a = v_a(T_a+Tbar_a)/(2 pi),
    j_a = v_a^2(T_a-Tbar_a)/(2 pi).

The expression for j follows by differentiating right and left moving h_a and imposing partial_t h_a+partial_x j_a=0. Independent chiralities have zero connected cross-correlator in the plane vacuum. Therefore at nonzero equal-time spatial separation,

    <h(x)h(0)>_stress = B_2/(4 pi^2 x^4),
    <j(x)j(0)>_stress = B_4/(4 pi^2 x^4).             (5.1)

At one spatial point and Euclidean time separation tau>0,

    <j(0,tau)j(0,0)>_stress = (sum_a c_a)/(4 pi^2 tau^4).
                                                               (5.2)

These are stress contributions. A raw lattice density can contain derivative improvements or additional oscillating operators. Identifying its full leading tail with (5.1) is a further matching statement, not automatic from the spectrum alone. Equations (5.1)-(5.2) agree with C100, Section 3. In particular, the spatial energy-current coefficient is B_4, not B_2. Under H->lambda H, h scales as lambda, j as lambda^2, v as lambda, and these powers are required by dimensional consistency.

### 5.2 Local norm does not bound the tail coefficient

[Proved lattice identities] With h_i a centered representative satisfying ||h_i||<=M,

    j_i = i[h_(i-1),h_i],
    d h_i/dt = j_i-j_(i+1),
    ||j_i|| <= 2M^2.

For any state, Cauchy-Schwarz for the connected covariance and Var(O)<=||O||^2 give

    |Cov(h_x,h_0)| <= M^2,
    |Cov(j_x,j_0)| <= 4M^4.                           (5.3)

[Proved limitation] Even grant a clean lower matching factor alpha>0 at every r>=R. Equations (5.1)-(5.3) then give at a chosen such r

    B_2 <= 4 pi^2 M^2 r^4/alpha,
    B_4 <= 16 pi^2 M^4 r^4/alpha.                     (5.4)

The best permitted r is a controlled matching distance. Sending r to infinity does not constrain B_p. No uniform R(d) follows from the infrared assumption.

[Proved covariance countermodel] This issue persists if one additionally imposes a nonnegative static structure factor. For any R>0 and K>0 define

    C_R(x) = K/[1+(x/R)^2]^2.

Then C_R(0)=K, |C_R(x)|<=K and

    Fourier[C_R](q) = (pi K R/2)(1+R|q|)exp(-R|q|) >=0,
    lim_(x->infinity) x^4 C_R(x)=K R^4.

Restricting this positive-definite kernel to integer sites preserves positive definiteness. Its tail coefficient is unbounded at fixed local variance. This is a counterexample to the proposed inference from bounded covariance and spectral positivity, not a construction of a critical spin-chain ground state or a counterexample to the target bound.

### 5.3 There is no velocity floor, even at d=2 with one Ising sector

[Proved exact lattice example] For 0<gamma<=1, consider

    H_gamma = -sum_i [(1+gamma) X_i X_(i+1)/2
                     +(1-gamma) Y_i Y_(i+1)/2+Z_i],

with the symmetric bond

    h_gamma = -(1+gamma) X tensor X/2
              -(1-gamma) Y tensor Y/2
              -(Z tensor 1+1 tensor Z)/2.            (5.5)

In the even two-spin parity block the matrix is [[-1,-gamma],[-gamma,1]], while in the odd block it is [[0,-1],[-1,0]]. Thus its eigenvalues are

    {-sqrt(1+gamma^2),-1,1,sqrt(1+gamma^2)}.

The physical swap P satisfies P h_gamma P=h_gamma and sends every coboundary o tensor 1-1 tensor o to its negative. Half spectral spread is convex and invariant under conjugation and scalar shifts. Averaging a gauge representative with its swapped representative returns h_gamma up to a scalar, so no coboundary lowers its spread. Hence

    M_gamma = sqrt(1+gamma^2).                        (5.6)

[Proved free spectrum] The Jordan-Wigner transformation turns (5.5) into hopping of amplitude one, pairing of amplitude gamma, and chemical potential of magnitude two. In a Nambu basis the Bloch matrix is, up to a constant unitary and sign convention,

    K(k)=2[(1-cos k) sigma_z+gamma sin k sigma_y].

Squaring this 2 by 2 matrix gives

    epsilon(k)=2 sqrt[(1-cos k)^2+gamma^2 sin^2 k].    (5.7)

There is a single particle-hole invariant gapless point k=0. For every fixed gamma>0 it is one massless real fermion, hence the Ising CFT with c=1/2. Expanding (5.7),

    epsilon(k)^2 = 4 gamma^2 k^2
                    +(1-4gamma^2/3)k^4+O(k^6),
    v_gamma = 2 gamma,
    A_gamma = gamma,
    v_gamma/M_gamma = 2gamma/sqrt(1+gamma^2) -> 0.    (5.8)

The XY dispersion and its normalization are also given in [Najafi, Rajabpour and Viti, equation (6)](https://link.aps.org/accepted/10.1103/PhysRevB.97.205103). The argument here supplies the bond optimization and the limiting ratio independently.

[Proved interpretation] For gamma>0 the pairing breaks the on-site fermion-number U(1), leaving its Z_2 parity. The Ising sector is neutral in the sense relevant to the charge-stiffness route. At gamma=0 the dispersion becomes quadratic, but that endpoint is not needed as a CFT example. Every member used in (5.8) has gamma>0, is critical, and has d=2. Thus a positive lower bound v>=m(d)M is false already for a single neutral sector with bounded central charge. Meanwhile A_gamma/M_gamma tends to zero, so the example does not violate the proposed upper bound on A.

[Proved crossover check] The linear regime requires |k| much smaller than gamma as gamma tends to zero, by comparing the k^2 and k^4 terms in epsilon^2. Thus the same example exhibits a diverging length for Ising matching. This is consistent with, and separate from, the stronger finite-speed crossover counterexample already in C106.

### 5.4 Even a moment bound needs more than an upper speed limit

[Proved algebra] If a velocity floor v_min existed, then

    A <= B_2/v_min,
    A <= B_4/v_min^3.

Section 5.3 rules out a universal positive v_min/M. A ceiling v<=v_LR gives inequalities in the wrong direction for this purpose. A bound on total central charge would allow A^2<=(sum c_a)B_2, but the required central-charge bound is not available.

[Proved algebraic countermodel] Take n abstract independent Ising factors with c_a=1/2 and v_a=n^(-1/2). Then

    B_2=1/2,   B_4=1/(2n),   max v_a<=1,
    A=sqrt(n)/2 -> infinity.

These arrays prove that neither a second nor a fourth velocity moment, together with a velocity ceiling, controls the first moment without a multiplicity or slow-sector constraint. They are not claimed to have fixed-d lattice realizations with bounded M. The distinction matters.

[Proved exact alternative accounting] Define the central-charge weighted number of sectors above a speed threshold,

    N(u)=sum_a c_a 1_(v_a>=u).

For a finite set of sectors, integrate each indicator to obtain

    A = integral_0^(v_LR) N(u) du.                   (5.9)

If B_p<=K_p M^p, then N(u)<=K_p(M/u)^p. For p=2 or 4 this upper envelope is nonintegrable at zero. A new integrable bound for N(u), or a direct integrated bound for it, would replace the false velocity-floor assumption. Equation (5.9) states the missing aggregate slow-sector control without proposing a bound on sum c_a itself.

## 6. Can neutral sectors be made charged at bounded cost?

### 6.1 Keeping a minimal model unchanged cannot create a U(1)

[Proved CFT obstruction] A nonzero unitary holomorphic U(1) current has positive level k and Sugawara stress tensor T_J=:JJ:/(2k) with central charge one. In a unitary theory the commuting remainder has nonnegative central charge. Therefore the total chiral c is at least one. A unitary minimal model with c<1 admits no such current. A continuous internal symmetry in a unitary two-dimensional CFT would require chiral conserved currents. Thus a deformation that retains the same minimal-model CFT cannot make it carry a nontrivial continuous charge. An extension or change of the infrared theory is needed.

[Proved spectator limitation] Adding a decoupled charged chain supplies a charge and a stiffness, but its generator acts as the identity on the original neutral factor. The original factor is still invisible. A low-energy extension has to mix the operator algebras in a specified way.

### 6.2 A parafermion extension exists in the continuum, with extra requirements

[Assumed standard parafermion construction; proved dimension check] The usual chiral extension of Z_k parafermions with a compact boson gives SU(2)_k currents

    J^+ = sqrt(k) psi_+ exp(+i sqrt(2/k) phi),
    J^- = sqrt(k) psi_- exp(-i sqrt(2/k) phi),
    J^3 = i sqrt(k/2) partial phi.

With <phi(z)phi(0)>=-log z, the vertex has weight 1/k and psi has weight 1-1/k. Hence J^+ and J^- have weight one. The central charges satisfy

    2(k-1)/(k+2)+1 = 3k/(k+2).

The general parafermion/boson realization is given in [Ginsparg, equation (9.34)](https://arxiv.org/pdf/hep-th/9108028), with roots normalized to length squared two. The displayed SU(2) weights and central-charge identity follow directly in that normalization.

[Proved limitation] Having weight one is not enough for a microscopic on-site charge. The combined fields must belong to the local extended operator algebra, and their chiral propagation requires the parafermion and boson velocities to agree. A tensor product at arbitrary radius, modular invariant or relative speed does not automatically satisfy these requirements. Even a valid continuum WZW extension does not specify a finite-dimensional nearest-neighbor lattice Hamiltonian realizing it.

### 6.3 A sufficient lifting theorem, with all required controls stated

[Proved conditional reduction] Suppose one could associate to every original chain a new one-site nearest-neighbor chain H' such that, in the same physical length and time units,

    d' <= g(d),
    M' <= K(d) M,
    A'_charged,WZW >= eta(d) A_neutral,
    eta(d)>0.

Assume the counted full WZW sectors of H' have microscopic on-site semisimple currents and satisfy the dictionary required by C104. Its established bound then gives

    A_neutral <= [4 pi K(d)/(3 eta(d))]
                   g(d)[g(d)^2-1] M.                (6.1)

For a purely U(1) lift, C104 alone is insufficient: one also needs its missing quantitative stiffness-weight comparison. Replacing neutrality by a nonzero but arbitrarily small charge coupling does not solve the task.

[Proved reason a bounded perturbation is insufficient] If two bonds differ in norm by at most delta, then on a ring |E'_0(L)-E_0(L)|<=L delta and their bulk energy densities differ by at most delta. Thus their finite-size Casimir estimators can differ by as much as

    (6L/pi) |(e'_0-e_0)L-(E'_0(L)-E_0(L))|
       <= 12 delta L^2/pi.

This bound says nothing about uniform control of A' relative to A. The endpoint example in Section 5.3 also has A tending to zero under a small parameter change. In a deformation argument, one must actually prove the third line of the lifting hypothesis, not infer it from a bounded change in M or continuity of finite-ring spectra. C102 has already ruled out general convexity and subadditivity shortcuts for A.

[Heuristic, unresolved] There is no general lattice construction here satisfying these requirements for arbitrary neutral minimal-model or coset factors. Coset notation identifies a continuum commutant stress tensor, not a bounded-cost inverse operation on arbitrary microscopic chains. No impossibility theorem for every possible lifting construction is claimed.

## 7. What has and has not been settled

[Proved, conditional on dynamical CFT matching] Locality gives a finite, crossover-free bound for a specified single bounded-central-charge neutral factor, including (4.3)-(4.4), and for a specified finite number through (4.5). The inherited free theorem remains stronger on the quasi-free subclass.

[Proved route obstructions] The one-bond twist norm has the wrong size dependence. A universal replica twist is exactly ring gluing. Local current norms do not bound infrared tail coefficients, the spatial current measures B_4, and a universal velocity floor is false even at d=2 on the Ising line. Bounded-cost deformation of the bond does not by itself compare infrared Casimir coefficients.

[Assumed, sufficient missing inputs] A complete general proof could follow from any one of the following additional results, with its stated matching requirements: (i) the averaged 1/L seam budget plus twist coverage in Section 2, or the universal replica estimate (3.4); (ii) a scale-free stress spectral estimate together with integrable slow-sector counting as in (5.9); (iii) a lattice lifting theorem with all three controls in Section 6.3. None is established here for all chains.

[Proved, status of the task] No bound A_neutral<=C M f(d) for an arbitrary collection of neutral or coset sectors is proved in this run, and no counterexample to that bound is constructed. The failures above are failures of specified deductions, not a claim that the desired theorem is false.

[Numerical, verification scope] checks.py tests the finite-dimensional seam identities, the optimized XY bond spectrum and its dispersion, the parafermion arithmetic, and a positive covariance kernel. The largest matrix diagonalized is 64 by 64. These checks support normalization and algebra only, not a lattice-to-CFT theorem. No heavy numerics or cluster jobs are needed.
