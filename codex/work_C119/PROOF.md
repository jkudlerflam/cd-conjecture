# C119. Full ring ratios and what they count

2026-10-29. Natural logarithms are used. The microscopic assumptions and M are those of ../../specs/_bound_background.md. The universal interacting bound remains open. Results conditional on an infrared description are labeled accordingly.

## 1. Exact fermion traces

[Definition] Write

    g_N(beta) = log Tr exp(-beta H_N),
    R_L(beta) = 2 g_L(beta) - g_(2L)(beta),       tau = beta/L.

[Proved] For a number-conserving quadratic chain with m orbitals, let epsilon_b(k) be its bands and set

    P_(N,alpha) = product_(j,b) [1 + exp(-beta epsilon_b(k_j))],
    S_(N,alpha) = product_(j,b) [1 - exp(-beta epsilon_b(k_j))],
    k_j = 2 pi (j+alpha)/N,       j = 0,...,N-1.

The grand fermion trace with fixed twist alpha is P. The trace with total fermion parity inserted is S. Projection gives (P+S)/2 for even number and (P-S)/2 for odd number. For the local spin chain obtained by the sitewise Jordan-Wigner map, including its closing bond,

    Z_N = [P_(N,1/2) + S_(N,1/2)
                     + P_(N,0) - S_(N,0)]/2.                 (1)

Indeed, even total number has antiperiodic fermion boundary conditions and odd total number has periodic boundary conditions. This includes every number sector. Scalar shifts multiply (1) by exp(-beta aN), which cancels from R. Substitution in the definition of R computes the ratio exactly at every size and temperature.

[Proved distinction] A tensor product of independent spin XX chains has an independent parity sum for each chain. Its Z is a product of expressions (1). A multi-orbital fermion chain mapped to one local spin chain has one total-parity sum. These are different torus theories even when their velocities and central charges agree. C106 belongs to the second class. Replacing its trace by two independent spin XX traces gives the wrong order-one term.

[Proved] Stable evaluation uses

    log P = -beta sum min(epsilon,0)
                         + sum log(1+exp(-beta |epsilon|)),
    S/P = (-1)^(number of negative epsilon)
                         product tanh(beta |epsilon|/2).

The zero-mode case follows by continuity. The code uses log1p and expm1 for projection. Direct spin diagonalization checks the boundary convention, not just the bands.

## 2. XX and independent stacks, including the first lattice correction

[Definition] Normalize h_XX = v(Sx Sx + Sy Sy), v>0. Then epsilon(k)=v cos k, d=2, M=v/2, e_inf=-v/pi, A=v, using the settled conventions of C116. Use L divisible by four for the following expansion. Formula (1) applies at other lengths too.

[Proved] Put y=beta v/L, q_r=exp(-2pi y r), and define

    T_+(y) = exp(pi y/6) product_(r=1/2,3/2,...) (1+q_r)^4,
    T_-(y) = exp(pi y/6) product_(r=1/2,3/2,...) (1-q_r)^4,
    T_R(y) = exp(-pi y/3) product_(r=1,2,...) (1+q_r)^4,
    Z_XX(y) = [T_+(y)+T_-(y)]/2 + 2 T_R(y),
    F_XX(y) = log Z_XX(y).

The AP filled sea has E_AP=-v csc(pi/L), even parity, and excitation energies v sin(2pi r/L), r half-integer, with multiplicity four. The P sea has E_P=-v cot(pi/L), two zero modes, and four excitations per positive integer r away from the band extrema. The two zero modes leave two states after the odd projection. Taking L to infinity in (1) gives these products.

[Proved] Define Theta(s)=sum_(n in Z) exp(-pi s n^2). Equivalently,

    Z_K(y) = Theta(y/(2K)) Theta(2Ky) / eta(iy)^2,              (2)
    Z_XX(y) = Z_(K=1)(y),
    eta(iy) = exp(-pi y/12) product_(n>=1)(1-exp(-2pi yn)).

Expanding the oscillator products and charge and winding sums gives (2), whose dimensions are Delta=n^2/(4K)+K w^2+N+Nbar. Poisson summation and the eta transformation give F_K(y)=F_K(1/y). Hence

    r_XX(y) = 2 F_XX(y)-F_XX(y/2).                             (3)

[Proved] The first positive XX dimension is 1/4, with multiplicity two. Expanding the theta sums and their logarithms gives

    r_XX(y) = pi y/4 - 2 exp(-pi y/4)
                          + 6 exp(-pi y/2)
                          + O(exp(-3pi y/4)),       y -> infinity,

    r_XX(y) = 4 exp(-pi/(2y)) - 6 exp(-pi/y)
                          + O(exp(-3pi/(2y))),      y -> 0.    (4)

For the second line use 2F(1/y)-F(2/y). The pi/(6y) thermal terms cancel identically. The small-y contribution is exponentially small, not proportional to 1/v.

[Proved] The first lattice correction is explicit. Define

    D_+(y) = 7pi^3 y/360
                  + (16pi^3 y/3) sum_(r half-integer>0) r^3 q_r/(1+q_r),
    D_-(y) = 7pi^3 y/360
                  - (16pi^3 y/3) sum_(r half-integer>0) r^3 q_r/(1-q_r),
    D_R(y) = -pi^3 y/45
                  + (16pi^3 y/3) sum_(r integer>0) r^3 q_r/(1+q_r),
    F_2(y) = {[T_+ D_+ + T_- D_-]/2 + 2 T_R D_R}/Z_XX(y).

Indeed beta v sin(2pi r/L)=2pi yr-(4pi^3 y r^3)/(3L^2)+O(L^-4). Expanding each logarithmic product gives the sums. The expansions of csc and cot give the vacuum terms. Exponential convergence at fixed y>0 justifies termwise expansion, with O_y(L^-4) remainder. Thus

    g_L(beta) = beta v L/pi + F_XX(y) + F_2(y)/L^2 + O_y(L^-4),
    R_L(beta) = r_XX(y)
              + [2F_2(y)-F_2(y/2)/4]/L^2 + O_y(L^-4).          (5)

The remainder constants need not stay bounded as y tends to zero. The finite products (1) have no such restriction.

[Proved check] The exact even-ring vacuum defect is

    E_(2L)-2E_L = v [2 csc(pi/L)-csc(pi/(2L))]
               = pi v/(4L) + 7pi^3 v/(192L^3) + O(L^-5).     (6)

Taking large y in (5) gives the same coefficient 7pi^3 y/192. Even XX cannot obey its sharp asymptotic slope at every finite L with only a beta-independent additive constant.

[Proved] For independent spin XX chains,

    d=2^m,       M=(sum_a v_a)/2,       A=sum_a v_a,
    R_L(beta)=sum_a R_L^XX(beta v_a),
    r(tau)=sum_a r_XX(v_a tau).

R_L^XX on the right uses unit velocity. The corrections (5) also add. Additivity of M follows from product bond extremizers with equal one-site marginals, as in C116.

[Proved scale qualification] At beta=tau L with tau fixed, a sector with v tau<<1 contributes the second line of (4). At beta=s L/v it instead contributes r_XX(s), independently of v. For example r_XX(1)=0.6326341534. It is not small just because the sector is slow. If multiplying L by the velocity is intended, the dimensionally specified choice beta=sLv/M^2 has y=s(v/M)^2 and is suppressed. The parameter y settles the scale question.

## 3. Several CFT sectors and the limitation of slow-mode suppression

[Assumed infrared input] Use C116's full-trace scaling assumption and vacuum-compatible rings. For an actual product,

    F(tau)=sum_a F_a(v_a tau),
    F_a(y)=pi c_a y/6 + log W_a(y),
    W_a(y)=sum_Delta g_(a,Delta) exp(-2pi y Delta).

For selection rules between sectors, a sum over allowed tuples replaces the product. No factorization is assumed in that case.

[Proved within this input] The full scaling answer is

    r(tau)=pi A tau/4 + 2log W(tau)-log W(tau/2),              (7)

where gaps in W are velocity-weighted dimensions. If the first positive weighted dimension is omega, with multiplicity divided by the vacuum multiplicity equal to b, then

    r(tau)=pi A tau/4+log g_0-b exp(-pi omega tau)
              + O(exp[-pi tau min(2omega,omega_2)]).

Coincident exponents are combined in the remainder. The leading excitation correction comes from the larger ring.

[Proved within modular product CFT] For a fixed sector with unique vacuum and first positive dimension Delta_*>0,

    2F_a(y)-F_a(y/2)
       =2b_* exp(-2pi Delta_*/y)
          +O(exp[-2pi min(2Delta_*,Delta_2)/y]),      y->0.     (8)

A finite vacuum multiplicity leaves log g_0. Unequal left and right speeds require separate chiral variables and half-central-charge weights, as in C116. Equation (8) concerns equal-speed modular sectors.

[Proved warning] Delta_* and its multiplicity belong to the theory. Knowing c and v alone does not make (8) uniform over varying CFTs. A vanishing velocity and a shrinking conformal gap can compensate. An unbounded number of small terms also requires a summability estimate. Both issues survive cancellation of the extensive thermal pressure.

[Proved continuum obstruction] At fixed y>0, the compact boson (2) has, as K->infinity,

    F_K(y)=(1/2)log(2K/y)-2log eta(iy)+o(1),
    r_K(y)=(1/2)log(K/y)-4log eta(iy)+2log eta(iy/2)+o(1).     (9)

This follows by Poisson summing the charge theta function while the winding sum tends to one. It diverges with c=1 and fixed v. Many light charge primaries, not many stress tensors, cause the divergence. This is not a fixed-d nearest-neighbor construction. A microscopic argument must constrain this spectral information, or only claim an ordered slope bound.

## 4. C106, parity selection and curvature corrections

[Assumed settled model data] C106 gives

    H(k)=2t sin k sigma_z + epsilon sigma_x-a epsilon I,
    epsilon_+/-=-a epsilon +/- sqrt(4t^2 sin^2 k+epsilon^2),
    k_0=arcsin[epsilon sqrt(a^2-1)/(2t)],
    v=2t sqrt(1-a^-2) sqrt[1-epsilon^2(a^2-1)/(4t^2)],
    M=sqrt(4t^2+epsilon^2),       d=4,       A=2v.

Here a>1 and epsilon sqrt(a^2-1)<2t. These also give a epsilon<M, the condition used for the bond norm. The four zeros are +/-k_0 and pi+/-k_0. Substituting both bands into (1) computes R exactly, with one parity constraint for both orbitals.

[Proved] Suppose k_0/(2pi) is rational. Take even L with Lk_0/(2pi) integer. The vacuum phases are half-integer on the AP grid and integer on the P grid. Linearizing the four crossings gives

    Z_C106^CFT(y)=[T_+(y)^2+T_-(y)^2]/2+8T_R(y)^2,
    F_C106(y)=log Z_C106^CFT(y),
    r_C106(y)=2F_C106(y)-F_C106(y/2),       y=v tau.            (10)

There are eight excitation factors per positive momentum and four periodic zero modes. Half of their sixteen occupations survive projection, giving eight. This is not Z_XX(y)^2. Its first positive dimension is 1/2 with multiplicity eight, so

    r_C106(y)=pi y/2-8 exp(-pi y/2)+O(exp(-pi y)).

The simultaneous spin-structure sum is modular invariant. Its small-y leading term is 16 exp(-pi/y). At epsilon=0 and L divisible by four the same linear theory applies with v=2t.

[Proved phase caveat] For incommensurate k_0, the phases {Lk_F/(2pi)} enter the trace and need not converge along all even lengths. Keep them in (1), or specify a convergent vacuum-compatible subsequence. Central charges and velocities do not determine the order-one answer without this information. Negative R values occur in the unaligned numerical checks. The bound in section 6 requires no commensurability.

[Proved massive-band contribution] At fixed epsilon>0, the lower C106 band is filled and has hole gap (a+1)epsilon. Its vacuum energy is an analytic periodic function with nearest complex singularities at distance asinh(epsilon/(2t)) from the real axis. Its bulk contribution cancels from R, its remaining vacuum contribution is exponentially small in L, and its hole trace is bounded by L exp[-beta(a+1)epsilon]. Its filled occupation is even on even rings, so it does not change the low-energy parity constraint. In the double scaling epsilon L=O(1), neither suppression is uniform and the lower band must remain in (1).

[Proved] The first curvature correction to (10) is also explicit. At k_0 set

    D(k)=epsilon^2+4t^2 sin^2 k,       E=a epsilon,
    kappa=D''/(2E)-(D')^2/(4E^3),
    lambda=D'''/(2E)-3D'D''/(4E^3)+3(D')^3/(8E^5).             (11)

The derivatives are at k_0. Orient q in the increasing-energy direction at each crossing. Then the local expansion is vq+kappa q^2/2+lambda q^3/6. Reflection gives the same kappa and oriented lambda at all four crossings. In each term of (10), use s=+1 or -1 and r positive half-integer for NS, and s=+1 and r positive integer for R. Put w_r=s q_r/(1+s q_r). Its log has L^-2 coefficient

    D_(alpha,s)=-(4pi^3 tau lambda/3) B_4(alpha)
       -(4tau lambda (2pi)^3/3) sum_r r^3 w_r
       +tau^2 kappa^2 (2pi)^4 sum_r r^4 w_r(1-w_r),            (12)

where B_4(1/2)=7/240 for NS and B_4(0)=-1/30 for R. The weighted mean of (12), with weights T_+^2/2, T_-^2/2 and 8T_R^2 divided by (10), is F_2,C106(tau). Hence, for fixed parameters and fixed tau>0 on the aligned grids,

    R_L(tau L)=r_C106(v tau)
       +[2F_2,C106(tau)-F_2,C106(tau/2)/4]/L^2+O(L^-4).       (13)

[Proved derivation] The particle and hole energies at positive q are vq +/- kappa q^2/2+lambda q^3/6. The terms of order 1/L in their two logarithms cancel. The cubic terms give -tau lambda(2pi r)^3 w_r/3, and the quadratic terms give tau^2 kappa^2(2pi r)^4 w_r(1-w_r)/4. Multiply by four crossings. Euler-Maclaurin gives vacuum energy term (pi^3/(3L^3)) sum_F lambda_F B_4(alpha), yielding the first term of (12) after multiplication by -beta. Two crossings with lambda=-v and kappa=0 reproduce section 2. The aligned particle-hole sums eliminate odd inverse powers of L. Unaligned phases can have order 1/L corrections and are not covered by (13).

[Proved nonuniformity] At fixed a>1 as epsilon->0, kappa=O(epsilon^-1) and lambda=O(epsilon^-2). Thus (13) is not uniform even at fixed v/M. This is the finite-size version of C106's thermal crossover obstruction. It does not make the exact ratio diverge. In the double scaling epsilon L/(2t)=u, near either zero of sin k,

    (L/(2t)) epsilon_+/-=-a u +/- sqrt[(2pi r)^2+u^2]+O(L^-2).

The exact products retain this crossover. The regimes u<<1 and u>>1 must not be interchanged in a linear CFT expansion. Section 6 bounds the ratio through the entire crossover.

## 5. Massive modes

[Proved, free lattice model] Take H_m(k)=2t sin k sigma_z+m sigma_x, the a=0 specialization of (1). For m>0 it has a unique filled-band vacuum on compatible rings and no massless sector. A fixed gapped free filled-band energy is analytic in a strip |Im k|<gamma. Its N-point quadrature error is O(N exp(-gamma' N)) for any gamma'<gamma. For this model one can take any gamma<asinh(m/(2t)).

[Proved] Excitation factors obey

    sum log(1+exp(-beta E_j))
                       <= (number of modes) exp(-beta m).

For even rings the AP vacuum has allowed parity. The P odd sector must contain an excitation because the filled lower band has even occupation. At beta=tau L with fixed positive m,tau,

    R_m=O(L^2 exp(-gamma' L))+O(L exp(-tau m L)) -> 0.         (14)

This deliberately nonoptimal bound holds for every fixed tau>0 in the free model. It is not uniform as m->0.

[Proved continuum check] For one complex Dirac mode with E(k)=sqrt(m^2+v^2 k^2) and momentum twist alpha, Poisson summation of the filled-sea energy gives

    delta E_N(alpha)=(2m/pi) sum_(n>=1)
                      cos(2pi n alpha) K_1(nmN/v)/n.          (15)

The Fourier transform at nonzero separation x of E(k) is -m K_1(m|x|/v)/(pi |x|), and the filled sea contributes with opposite sign. At m->0 and alpha=1/2, (15) gives -pi v/(6N), the c=1 Casimir value. At large mN/v it is of order sqrt(mv/N) exp(-mN/v). The vacuum part of R is beta[delta E_(2L)-2delta E_L]. Its thermal part is obtained from 2 sum log(1+exp(-beta E_j)), including any spin projection afterward. Both cones and their shared projection enter the lattice example.

[Proved interpretation] Modes with mL/v of order one contribute order one at beta of order L. This need not mean extra infrared sectors. At fixed Hamiltonian the contribution disappears in the infrared limit. Finite vacuum degeneracy instead leaves log g_m, as in C116's Ising spectator, and changes no large-tau slope.

## 6. A uniform lattice estimate for the requested free families

[Proved matrix quadrature lemma] Let H(k)=X+B exp(ik)+B^dagger exp(-ik) be a finite Hermitian Bloch matrix. Define

    phi_beta(x)=log[2 cosh(beta x/2)],
    G(k)=tr phi_beta(H(k)),
    I_H=integral_0^(2pi) ||H''(k)||_1 dk.

Then

    integral |G''(k)| dk <= beta I_H.                         (16)

To prove it without differentiating singular band eigenvectors, write

    G''=tr[phi_beta'(H) H'']+Q.

Q is the Hessian of the convex spectral trace. In an eigenbasis its diagonal coefficients are phi_beta''(epsilon_i), and its off-diagonal coefficients are divided differences of the increasing function phi_beta'. All are nonnegative, including at degeneracies by continuity. Since |phi_beta'|<=beta/2, the negative part of G'' is bounded by (beta/2)||H''||_1. Periodicity gives integral G''=0, so its absolute integral is twice its negative-part integral. This proves (16), avoiding C106's divergent individual-band curvature.

[Proved scalar quadrature step] For Ghat_n=(1/(2pi)) integral G(k) exp(-ink)dk, two integrations by parts give

    |Ghat_n| <= ||G''||_(L1)/(2pi n^2),       n != 0.

Only Fourier modes divisible by N survive a twisted N-point sum. Therefore, for every alpha,

    |sum_j G(2pi(j+alpha)/N)-N Ghat_0|
       <= [pi/(6N)] ||G''||_(L1)
       <= pi beta I_H/(6N).                                  (17)

The identity log(1+exp(-beta x))=phi_beta(x)-beta x/2 shows that the grand trace has the same error for N>=2. The trace of H has only modes 0,+1,-1, whose nonconstant modes sum to zero.

[Proved parity step] Suppose N is even and the band multiset at k equals that at -k. The AP grid consists of distinct k,-k pairs, so the number of negative AP energies is even. A zero gives S_AP=0. Hence

    (1/2)P_AP <= Z_N <= P_AP+P_P.

If p(beta) is the grand-trace pressure and a_N=pi beta I_H/(6N), (17) gives

    Np-a_N-log 2 <= log Z_N <= Np+a_N+log 2.

Subtract the lower bound at 2L from twice the upper bound at L:

    R_L(beta)<=3log 2+(5pi/12) I_H beta/L,       L even >=4.    (18)

No low-temperature or commensurability assumption occurs. A negative S_P is harmless. Evenness of the band multiset, rather than reflection invariance of the bond, is the condition used here. Equation (18) is not asserted for arbitrary paired or asymmetric free spin chains.

[Proved applications] For XX, I_H=4v=8M. For C106,

    H''(k)=-2t sin k sigma_z,       I_H=16t<=8M.

The a=0 massive model has the same I_H and M=sqrt(4t^2+m^2). They all obey

    R_L(beta)<=3log 2+(10pi/3)M beta/L.                        (19)

Independent spin XX stacks add their bounds and their M, giving

    R_L(beta)<=3log d+(10pi/3)M beta/L.                        (20)

The weaker intercept 3log d also covers C106 and its massive comparison. Thus k=10pi/3, f(d)=1 is a crossover-free estimate for all requested free families, uniform in velocities, avoided-crossing scale and masses. The coefficient is not sharp. It does not establish an interacting bound. For the untwisted grand fermion trace the parity intercept is absent.

## 7. Attempts to make R large at fixed d and M

[Proved obstruction in the requested families] Equation (20) forbids R_L(sL/M) from diverging at fixed s, even for parameters depending on L. Stacks spend an additive M budget. The C106 avoided crossing only changes how the finite answer crosses between spectral regimes. Small masses give a finite crossover contribution. None produces an unbounded ratio.

[Assumed published XXZ dictionary] For h=SxSx+SySy+Delta SzSz with Delta=cos(gamma), the compact-boson parameters in (2) are

    K=pi/[2(pi-gamma)],       v=pi sin(gamma)/(2gamma),

with M=1/2 for -1<Delta<=1. The dictionary is stated in [Collura, Essler and Groha, equations (16)-(17)](https://arxiv.org/html/1706.07939v1#S2). Applying it here assumes the full infrared torus description, not merely correlation-function scaling.

[Proved within this dictionary] The direct microscopic radius attempt is gamma=pi-eta, eta->0. Then K~pi/(2eta), v~eta/2 and Kv->pi/4. Apply modular covariance before taking this limit. At beta=tau L,

    r_XXZ(tau) -> 2log Theta(2/(pi tau))
                         -log Theta(4/(pi tau)).               (21)

The dual oscillator and winding corrections vanish. The charge sum retains the displayed theta function. It is finite at fixed tau and behaves as (1/2)log(pi tau) at large tau. The order of limits is L->infinity at fixed eta,tau, then eta->0, then, if desired, tau->infinity. This calculation does not replace the endpoint ferromagnet by a CFT.

[Proved consequence] Equation (21) is not evidence for a large A, because A=v->0. It disproves the unrestricted claim that every slow sector is negligible at fixed tau when its primary dimensions may shrink. Independent slow XX sectors do become negligible, because their dimensions stay fixed.

[Proved further distinction] At tau=1/v, the same XXZ sequence has r~(1/2)log K by (9), although A tau=1. Thus r<=pi A tau/4 plus a fixed constant uniform over Hamiltonians is false. The target with M instead of A is not contradicted, since M tau grows as 1/v. An ordered o(tau) remainder is the correct sharp statement about the slope.

[Open] No nearest-neighbor chain with fixed d,M and R_L(sL/M) unbounded was constructed. Equation (20) proves this impossible in the requested free families, but gives no impossibility theorem for interacting chains. Equation (9) identifies another potential mechanism, an unbounded light-primary tower at fixed aspect ratio, even with one fast CFT. Charged realizations face the campaign's stiffness budget; neutral sectors have no general substitute. The stated assumptions do not settle this construction problem.

## 8. Constants supported by the data and the missing seam lemma

[Proved distinction] There are three different estimates.

1. [Proved, requested free families] The all-even-length estimate (20) has K(d)=3log d, k=10pi/3, f(d)=1. It does not use the CFT assumption.
2. [Proved, free infrared slope] The inherited theorem A<=2M for all nearest-neighbor quasi-free chains gives lim_(tau->infinity) r(tau)/tau<=pi M/2 on compatible sequences. Thus k=pi/2, f(d)=1 is sharp for the free asymptotic slope. XX stacks saturate it. The stronger pointwise scaling inequality r_XX(y)<=pi y/4 was checked numerically here, not proved here.
3. [Heuristic, interacting envelope] The strongest specified campaign envelope is A<=pi M H_(d-1), where H_(d-1)=sum_(j=1)^(d-1)1/j. Its ratio version is

       limsup_(tau->infinity)
            [limsup_(L->infinity) R_L(tau L)]/tau
                              <=(pi^2/4)M H_(d-1).            (22)

This chooses k=pi^2/4, f(d)=H_(d-1). Heisenberg forces this constant at d=2 if f(2)=1. The growth of f at larger d is not established as optimal. Finite solvable data do not determine a unique sharp function of d. The weaker log_2 d envelope remains consistent.

[Proved restriction] Equation (22) must not be advertised as an all-scale inequality with this coefficient. C116's Heisenberg L=4 defect exceeds it; (6) shows the same issue for the sharp free coefficient. Massive crossover corrections and shrinking primary gaps further forbid uniform substitution of the Casimir term for the trace. A remainder b_h(tau)=o(tau), after the L limit, suffices and may depend on the fixed Hamiltonian. There is no assumption of a uniform onset scale.

[Proved exact seam representation] Put two L-rings and one reconnected 2L-ring on the same Hilbert space. Let H_0 be the decoupled pair, H_1 the reconnected ring, and

    V=H_1-H_0,       H_lambda=H_0+lambda V,
    V(s)=exp(sH_0)V exp(-sH_0),
    S_L(beta)=< T exp[-integral_0^beta V(s)ds] >_(H_0).

Dyson's identity inside the trace and free-energy differentiation give

    S_L(beta)=Z_(2L)(beta)/Z_L(beta)^2=exp[-R_L(beta)],
    R_L(beta)=beta integral_0^1 <V>_(H_lambda,beta) d lambda.   (23)

Brackets denote normalized thermal traces. No commutativity or reflection positivity is assumed. V has two removed and two inserted bonds. Its norm bound is only 4M and loses the separation L.

[Assumed minimal sufficient seam lemma] An improvement of the integrated response would suffice:

    limsup_(L->infinity) L integral_0^1
                       <V>_(H_lambda,tau L) d lambda
       <= k M f(d)+b_h(tau)/tau,
    b_h(tau)/tau -> 0 as tau -> infinity.                     (24)

Equivalently the normalized Euclidean reconnection insertion needs a lower bound exp[-k M f(d)tau-b_h(tau)-o_L(1)]. An all-scale strengthening is exp[-K(d)-k M f(d)beta/L]. Equation (24) allows finite-size and massive corrections to disappear in their proper order and bounds the integrated response, not each interpolating Hamiltonian separately. It avoids assuming a sign for individual seam correlations. With the C116 infrared trace assumption it implies A<=(4k/pi)Mf(d).

[Proved conditional transfer formulation] If an exact positive spatial transfer operator exists on a chosen blocked sequence and its relevant powers are trace class, let rho=T_beta^L/Tr T_beta^L. Then R_L=-log Tr rho^2. The needed effective Renyi-2 rank estimate is

    (Tr T_beta^L)^2/Tr T_beta^(2L)
                      <= exp[K(d)+k M f(d)beta/L],            (25)

or its asymptotic version (24). This concerns separation L in the entire imaginary-time history space. The rank of a single gate does not imply it. Arbitrary Hermitian bonds need not yield a positive spatial transfer operator, so (23)-(24) are the formulation for the full class. These statements isolate the missing inequality; locality alone has not proved it.

## 9. Checks and scope

[Numerical] checks.py compares the complete C106 spin trace with diagonalization at dimensions at most 256, with maximum log-trace error 4.98e-14. It checks (5), (13), slow-sector suppression, parity selection, massive suppression, radius growth and (21). A deterministic scan of 360 C106 length/parameter/aspect-ratio combinations obeys (19). These validate the explicit calculations, not a universal interacting estimate.

[Work record] All Python ran in work_C119. The small run took less than one second. An extended scan is prepared in run_typhon.slurm, and the Python --cluster option requires SLURM_JOB_ID. The sandbox cannot reach Typhon, no job was submitted, and no large diagonalization was substituted locally. No independent referee or sub-agent was used. Project-level files were not changed because the task restricts writes to this folder.
