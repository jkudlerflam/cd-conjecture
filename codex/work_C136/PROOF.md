# C136 proof and remaining estimates

2026-11-07. J=1 in the analytic sections. Restore J by homogeneity. Claims labeled proved are mathematical deductions, not numerical certifications. The final uniform theorem requested in C136 is not established here.

## 1. Conventions and inherited finite-volume identification

[Proved, inherited from C129 and the refereed C134 Section 2] Set

    J = sqrt(J_perp^2 + D_DM^2),       Delta = J_z/J,
    F_L(n) = min_Phi e_L(n,Phi).

For J>0 and -1<Delta<=0, F_L(n) is the sector ground energy of

    K_L = sum_x [-Sx_x Sx_(x+1) - Sy_x Sy_(x+1)
                                      + Delta Sz_x Sz_(x+1)].

The ring has L oriented summands, including two at L=2. The minimizing flux in the positive transverse convention is L pi modulo 2 pi. All integer lengths are allowed. Spin flip reduces to 0<=n<=L/2. We always double to (2L,2n), including when L is odd. This differs from C134's odd-length half-filled doubling to (2L,L).

For 1<=n<=L/2, let gamma=arccos Delta, alpha=pi-gamma, and

    p(x) = 2 atan[tanh(x) cot(gamma/2)],
    theta(x) = 2 atan[tanh(x) cot(gamma)],
    I_j = j-(n+1)/2.

The unique real solution of

    L p(x_j) - sum_k theta(x_j-x_k) = 2 pi I_j

has x_j=-x_(n+1-j), and

    F_L(n) = Delta L/4 + sum_j e(x_j),
    e(x) = -sin(gamma)^2/[2(sinh(x)^2+sin(gamma/2)^2)].       (1.1)

The identification already proved in C134 is for every n in this range, not just half filling. Its strict Yang-Yang Hessian, sign-adapted coercivity decomposition, free-point coordinate Bethe vector, analytic characteristic polynomial, and Perron-Frobenius simplicity apply unchanged. The referee's repairs are understood here. In particular, coercivity uses a ray decomposition adapted to the signs of the sorted coordinates. The set of parameters where the analytic eigenvalue is the ground eigenvalue is nonempty, closed, and open by isolation of the simple ground eigenvalue. At multiple contacts the missing hopping terms are sums of adjacent-pair contact terms; the scalar exchange factors satisfy each such condition.

Define

    a1(x) = p'(x)/(2 pi) = sin(gamma)/[pi(cosh(2x)-Delta)],
    K(x) = theta'(x)/(2 pi),       T(x)=-K(x)>=0,
    q = integral T = 2 gamma/pi-1,       B = 1/(1-q),
    a = pi sin(gamma),                   b = pi/gamma.

Then e=-a a1. T is even and decreasing on x>0, q<1, and a1 is positive, even, and decreasing on x>0. The same statements include T=0 at Delta=0.

## 2. Finite Fermi intervals and their density

[Proved] For Q>=0, let T_Q act by integration on [-Q,Q]. On both L1 and L-infinity its norm is at most q. Therefore

    (I-T_Q)^(-1) = sum_(k>=0) T_Q^k,
    ||(I-T_Q)^(-1)|| <= B.                                  (2.1)

This inverse is positivity preserving. Its norm is independent of Q, although it is not uniform as Delta tends to -1.

Let

    r_Q = (I-T_Q)^(-1) a1,       Z_Q = (I-T_Q)^(-1) 1,
    m(Q) = integral_(-Q)^Q r_Q(x) dx.

Extend the root density by zero outside its interval, writing r_Q^occ. Its counting density on the whole line is

    s_Q(x)=a1(x)+T*r_Q^occ(x).

It is strictly positive everywhere and equals r_Q inside the interval. Termwise comparison of the positive Neumann series with the full-line series gives

    a1 <= r_Q <= r_infinity = 1/[2 gamma cosh(bx)]            (2.2)

on [-Q,Q], and 1<=Z_Q<=B. Both r_Q and Z_Q are even and nonincreasing on the positive half interval. To see the latter assertion without assuming a rearrangement theorem, extend an even decreasing function by zero. Write it by its level sets as a positive superposition of centered interval indicators. Convolution of an even decreasing kernel with each such indicator is even decreasing, since its derivative at x>0 is T(x+c)-T(x-c)<=0. Every Neumann term has this property.

[Proved] Differentiation of the integral equation and self-adjointness give

    m'(Q)=2 r_Q(Q) Z_Q(Q)>0.                                (2.3)

In detail, the partial derivative of r_Q at fixed x solves

    (I-T_Q) partial_Q r_Q
        = r_Q(Q)[T(x-Q)+T(x+Q)].

On integrating and moving the inverse onto 1, the two boundary integrals each equal Z_Q(Q)-1. Adding the derivative of the integration limits gives (2.3). Thus m is strictly increasing, m(0)=0, and monotone convergence of the positive Neumann series gives m(infinity)=1/2. Every 0<nu<1/2 defines a unique finite Q with m(Q)=nu.

The following explicit upper bound on Q will be useful. Subtract the restricted and full-line equations and integrate their positive difference. The inside deficit is at most Bq times the full-line outside mass. Consequently

    1/2-m(Q) <= B integral_(|x|>Q) r_infinity(x) dx
              <= (2B/pi) exp(-bQ).                          (2.4)

The final inequality uses sech u<=2 exp(-u) and gamma b=pi. The direction of (2.4) matters: a lower bound on the missing density gives an upper bound on Q.

## 3. A cumulative-root bound uniform in filling

[Proved] Put nu=n/L, mu=L^(-1) sum_j delta_(x_j), and

    z_L(x)=p(x)/(2pi)-L^(-1)sum_j theta(x-x_j)/(2pi),
    sigma_L(x)=z_L'(x)=a1(x)+T*mu(x)>0.

Its range is (-c,c), with

    2c = alpha/pi + q nu,
    2c-nu = (alpha/pi)(1-2nu).                               (3.1)

Let C_mu(x)=mu((-infinity,x])-nu/2. If P_nu(t) denotes projection to [-nu/2,nu/2], midpoint counting implies

    ||C_mu-P_nu(z_L)||_infinity <= 1/(2L).                  (3.2)

The bound includes both one-sided limits at atoms. Outside the occupied counting interval the difference is zero. For the continuum density of mass nu, define C_nu analogously and z_nu(x)=integral_0^x s_Q(y)dy. Then

    C_nu=P_nu(z_nu).

The two measures have the same mass. Integration by parts, with vanishing difference of cumulative functions at both infinities, yields

    z_L-z_nu = T*(C_mu-C_nu).

Since projection to an interval is 1-Lipschitz, (3.2) gives

    ||C_mu-C_nu||_infinity <= B/(2L),
    ||z_L-z_nu||_infinity <= qB/(2L).                       (3.3)

This also holds at nu=1/2 with Q=infinity. No inverse endpoint density and no lower bound on nu occurs. It is stronger than applying C134's full-line discrepancy below half filling, which would compare measures of different masses.

For example, integration by parts once more gives

    ||sigma_L-s_Q||_infinity <= B ||T'||_1/(2L),
    ||sigma_L'-s_Q'||_infinity <= B ||T''||_1/(2L).          (3.4)

These are density-uniform but anisotropy-dependent. They do not, by themselves, prove an O(1/L) energy error. Pairing a smooth energy with (3.3) by a first-order bound only yields O(1) total energy error.

## 4. An exact finite-boundary quadrature identity

[Proved] Suppose nu<1/2. Define the finite boundary Q_L by

    z_L(Q_L)=nu/2.

It exists by (3.1). Every root lies strictly inside [-Q_L,Q_L]. Define the signed quadrature measure

    R_L = mu - 1_[-Q_L,Q_L] sigma_L(x) dx.

Its total mass is zero. Its cumulative function has sup norm at most 1/(2L), vanishes outside the interval, and is the midpoint sawtooth in the coordinate t=z_L(x).

In this section r,Z,m refer to the continuum objects at Q=Q_L, so m need not equal nu. The following identities are exact:

    (I-T_Q)(mu-r dx)=R_L,
    nu-m = integral Z R_L.                                 (4.1)

These follow by subtracting sigma_L=a1+T_Q mu from r=a1+T_Q r, restricted to the interval, then using self-adjointness. Because Z is even decreasing and 1<=Z<=B,

    |nu-m| <= Var(Z)/(2L) <= (B-1)/L <= B/L.                (4.2)

There is no integration-by-parts boundary term because the cumulative function of R_L vanishes at both endpoints.

Set

    h(Q) = -a r_Q(Q)/Z_Q(Q),
    epsilon_Q = (I-T_Q)^(-1)(e-h(Q)),
    E(Q) = Delta/4 + integral_(-Q)^Q e(x) r_Q(x) dx.

Then epsilon_Q(+-Q)=0 and

    F_L(n)-L[E(Q)+h(Q)(nu-m(Q))]
                          = L integral epsilon_Q R_L.       (4.3)

Thus the correct finite-boundary identity contains a chemical-potential subtraction and a density mismatch. Omitting these terms would leave a first-order endpoint error.

Under t=z_L(x), the right side of (4.3) is ordinary midpoint quadrature on [-nu/2,nu/2]. If

    V_L = Var_[-Q,Q](epsilon_Q'/sigma_L),

then the Peano-kernel estimate gives

    |L integral epsilon_Q R_L| <= V_L/(8L).                 (4.4)

There are no outside strips in (4.4). No clipping is needed at a finite Q. Clipping in a bound uniform as Q tends to infinity would introduce interior derivative jumps and a tail error, both of which would have to be estimated.

## 5. Bulk curvature and an explicit theorem below half filling

[Proved] Define the bulk function e_bulk(nu)=E(Q) at m(Q)=nu. This definition will also be identified as the finite-ring limit below. Differentiate E and use self-adjointness as in (2.3). It gives

    dE/dQ = 2 r_Q(Q) [(I-T_Q)^(-1)e](Q),
    e_bulk'(nu)=h(Q).

Because epsilon_Q vanishes at the endpoints, integration by parts in its x derivative gives

    (I-T_Q) epsilon_Q'=e'.                                 (5.1)

For x>0, e'(x)>0. On odd functions the operator T_Q reduces on [0,Q] to the nonnegative kernel T(x-y)-T(x+y), whose norm is at most q. Therefore epsilon_Q'(x)>0 for 0<x<=Q. Moreover |a1'|<=2a1 for Delta<=0, so positivity of the inverse gives

    |epsilon_Q'| <= 2a r_Q.                                (5.2)

Differentiating epsilon_Q at fixed x with respect to Q produces no endpoint source, because epsilon_Q(+-Q)=0. Hence partial_Q epsilon_Q=-h'(Q)Z_Q. Differentiating its boundary value gives

    h'(Q)=epsilon_Q'(Q)/Z_Q(Q),
    0 <= e_bulk''(nu)
       = epsilon_Q'(Q)/[2 r_Q(Q) Z_Q(Q)^2] <= a.             (5.3)

These statements extend by limits at nu=0 and nu=1/2. In particular h(0)=-(1+Delta), while h(Q) tends to zero as Q tends to infinity by (2.2). For any two densities u,v in [0,1/2], integrating (5.3) twice gives

    0 <= e_bulk(u)-e_bulk(v)-h(v)(u-v) <= a(u-v)^2/2.        (5.4)

Combining (4.2)-(4.4) with (5.4) proves the finite identity estimate

    |F_L(n)-L e_bulk(n/L)| <= [V_L/8+a B^2/2]/L.            (5.5)

The positive side of the error only needs V_L/8. This identifies a concrete variation estimate that would suffice to finish the density problem.

[Proved, new explicit restricted theorem] Fix 0<delta<=1 and 0<eta<1/2. For every -1+delta<=Delta<=0, all L>=2, and all 0<=n/L<=1/2-eta, a single explicit constant K_(delta,eta) satisfies

    |F_L(n)-J L e_bulk(Delta,n/L)| <= K_(delta,eta) J/L,
    L[F_(2L)(2n)-2F_L(n)] <= (5/2)K_(delta,eta) J.          (5.6)

Here e_bulk is dimensionless. One admissible constant is as follows:

    s0 = sqrt(delta(2-delta)),
    B0 = pi/[2 arccos(1-delta)],
    R0 = log[4 B0/(pi eta)],
    m0 = s0/[pi(cosh(2R0)+1)],
    S1 = 2/pi + 1/(pi s0^3),
    U1 = 2 B0,                 U2 = 8+4 B0/(pi s0),
    V0 = 2 R0 [U2/m0 + U1 S1/m0^2],
    K_(delta,eta) = max{(2B0/eta+1)^2, V0/8+pi B0^2/2}.     (5.7)

All constants are uniform down to zero density. No dilute asymptotic or fixed-particle-number assumption is used.

To prove (5.7), first suppose L>=2B0/eta. By (4.2), m(Q_L)<=1/2-eta/2. Equation (2.4), b>=1, and B<=B0 imply Q_L<=R0. On this interval

    sigma_L>=a1>=m0.                                      (5.8)

For the derivative bounds, write s=sin gamma and d(x)=sinh(x)^2+s^2. Then T(0)=|Delta|/(pi s), |d'|/d<=2/s^2, and

    ||T'||_infinity <= 2/(pi s0^3),
    ||T'||_1 = 2T(0) <= 2/(pi s0).

Also |a1'|<=2a1 and |a1''|<=8a1. The latter follows by differentiating log a1: its first derivative has magnitude at most 2, and its second at most 4 on Delta<=0. Since a1<=1/pi and mu has mass <=1/2,

    ||sigma_L'||_infinity <= S1,
    ||epsilon_Q'||_infinity <= U1,
    ||epsilon_Q''||_infinity
        <= ||e''||_infinity+||T'||_1 ||epsilon_Q'||_infinity
        <= U2.                                             (5.9)

The epsilon' bound uses (5.1), a<=pi, and the Neumann norm. Differentiating (5.1) with respect to x gives the second inequality in (5.9). Integrating the elementary derivative of epsilon'/sigma over an interval of length at most 2R0 now gives V_L<=V0. This proves the large-L branch by (5.5).

For the small-L branch, the bond spectral diameter is 1. The bulk e_bulk(nu) belongs to this same spectral interval: for rational nu<1/2 this follows from the just-proved large-L convergence along multiples of its denominator, and for irrational nu by continuity. Thus |F_L-L e_bulk|<=L. When L<2B0/eta, this is bounded by (2B0/eta+1)^2/L. The sector n=0 is exactly extensive. Finally apply the two-sided bound to the matched densities at L and 2L to get the factor 5/2 in (5.6).

[Proved, checks of the analytic formulas] At Delta=0, T=0, r_Q=a1, Z=1, m(Q_L)=nu, and (4.3) is the free midpoint identity. The curvature formula gives e_bulk''(nu)=pi sin(pi nu), agreeing with e_bulk=-sin(pi nu)/pi. The limit nu=0 gives the one-magnon chemical potential -(1+Delta), matching C129's saturation field. Every quantity in (5.7) is finite for delta,eta>0. Its divergence as eta tends to zero is not a uniform half-filled theorem.

## 6. Uniformity through the half-filled density endpoint

[Proved, new all-filling theorem] The dependence on eta in (5.7) can be removed. For fixed 0<delta<=1 define the following positive constants. Their size is not optimized.

    alpha0 = arccos(1-delta),       s0 = sqrt(delta(2-delta)),
    B0 = pi/(2alpha0),
    beta = 1/2,
    Mbeta = sqrt(2) sin(pi/4-alpha0/2),
    C0 = 4/[pi s0(1-Mbeta)],
    H1 = 2/s0^2,                   H2 = 8/s0^4+4/s0^2,
    c = (2/3) log max{1,(16/3)C0},
    t0 = s0/[16 pi cosh(c+1)^2],
    kappa = min{1/2,
                s0/[cosh(2(c+1))+1],
                (t0/2) exp(-2c)},
    P0 = (8/3)C0,      P1=P0 H1,      P2=P0 H2,
    W = [P2+4P0+2P1+H1(P1+2P0)]/kappa,
    Vv = pi[4+(4/3)W+(32/3)C0 H1/kappa],
    D = 1+3B0/(pi s0),              D1=(3/2)B0 H2,
    Vclip = 2Vv+8pi(H1+D1/D)+8pi,
    K_all(delta) = max{(4pi D/kappa+1)^2,
                      Vclip/8+pi[128D^2/kappa^2
                                  +72D B0/kappa+B0^2/2]}.

Then, for every Delta in [-1+delta,0], every L>=2, and every 0<=n<=L/2,

    |F_L(n)-J L e_bulk(Delta,n/L)| <= K_all(delta) J/L,
    L[F_(2L)(2n)-2F_L(n)] <= (5/2)K_all(delta) J.            (6.1)

The bulk energy is dimensionless. The same constant covers both parities and the whole dilute limit. Equations (5.6)-(5.7) remain useful smaller constants away from half filling. The proof of (6.1) follows.

### 6.1 An exponentially decaying resolvent majorant

[Proved] Let R_infinity=sum_(j>=1) T^(*j), the full-line positive resolvent kernel. With beta=1/2, its exponentially weighted convolution norm is controlled because

    integral exp(beta x)T(x)dx
       = sin[(2gamma-pi)beta/2]/sin(pi beta/2)
       = sqrt(2)sin(gamma/2-pi/4) <= Mbeta < 1.              (6.2)

This is the analytic continuation of the elementary kernel transform in C134 Section 3; it also follows by the same rectangular-contour integral. All integrals converge since beta<2. Moreover

    T(x)<=T(0)sech(x)^2<=4T(0)exp(-2|x|),
    T(0)=|Delta|/[pi sin(gamma)]<=1/(pi s0).

Apply the L-infinity/L1 convolution inequality to exp(beta x)T^(*j), then sum its geometric bound. Evenness gives

    R_infinity(x)<=C0 exp(-beta|x|).                        (6.3)

For any restricted resolvent R_Q, positivity of every Neumann term gives R_Q(x,y)<=R_infinity(x-y). Differentiating only the first kernel in each term, and using |T'|<=H1 T and |T''|<=H2 T, also gives

    |partial_x^j R_Q(x,y)|<=Hj C0 exp(-beta|x-y|),
                    j=1,2.                                (6.4)

The same bounds hold for the kernel carrying a source outside [-Q,Q] to a point inside it. There are no endpoint differentiations in (6.4).

### 6.2 Comparison to the full-line density

[Proved] Write r0=r_infinity and d_Q=r0-r_Q on [-Q,Q]. The restricted equation gives

    d_Q=(I-T_Q)^(-1) T_(outside Q) r0.

The positive source-to-interval kernel is bounded by R_infinity. Since r0(y)<=exp(-b|y|)/gamma, (6.3)-(6.4) yield, for |x|<=Q and j=0,1,2,

    |d_Q^(j)(x)|
      <= [2C0 Hj/(gamma(b+beta))]
                              exp(-(b+beta)Q) cosh(beta x),
    |d_Q^(j)(x)|/r0(x)
      <= Pj exp[-(3/2)(Q-|x|)],                            (6.5)

where H0=1. The second line uses r0(x)>=exp(-b|x|)/(2gamma), b>=1, and b+beta>=3/2. The definition of c ensures d_Q/r0<=1/2 on |x|<=Q-c whenever this core is nonempty.

We need a lower bound at the boundary as well. The resulting global comparison is

    kappa r0(x) <= r_Q(x) <= r0(x),       |x|<=Q.            (6.6)

Here are the details, including uniformity near the free point. If |Delta|<=1/16, then s=sin gamma>=9/10 and gamma-pi/2=arcsin|Delta|<=2|Delta|. The exact moment in (6.2) is at most sqrt(2)/16<1/10. Thus the pointwise version of C0 is at most

    4/[(16)(3)(9/10)(9/10)] < 3/16.

The first line of (6.5), before replacing b+beta by 3/2, then gives d_Q/r0<=4(3/16)/(3/2)=1/2 throughout the interval. At Delta=0 the assertion also follows directly from T=0.

For |Delta|>=1/16, the kernel satisfies

    T(y)>=s0/[16pi cosh(y)^2].

If Q>=c+1, integrate r_Q(Q)>=integral T(Q-y)r_Q(y)dy just over y in [Q-c-1,Q-c]. This interval lies in the core where r_Q>=r0/2, and r0(y)>=r0(Q). Hence r_Q(Q)>=t0 r0(Q)/2. On the remaining boundary layer, monotonicity of r_Q and

    r0(Q)/r0(x)>=exp[-b(Q-x)]>=exp(-2c)

give the third term in the definition of kappa. If Q<c+1, use r_Q>=a1, a1>=s0/[pi(cosh(2(c+1))+1)], and r0<=1/pi. This gives the second term in kappa. These cases prove (6.6), with constants independent of Q and density.

### 6.3 Variation of the dressed velocity

[Proved] Let v_Q=epsilon_Q'/r_Q. Equation (5.2) gives |v_Q|<=2a. We next prove Var(v_Q)<=Vv independently of Q.

Set ell=r0'/r0=-b tanh(bx). Its derivative has total integral at most 2b<=4, |ell|<=2, and |ell'|<=4. Also |r_Q'|<=H1 r_Q, by differentiating r_Q=a1+T_Q r_Q and using |a1'|<=2a1 and |T'|<=H1 T. The identity

    r_Q'/r_Q-ell = (-d_Q'+ell d_Q)/r_Q

can now be differentiated explicitly. Equations (6.5)-(6.6) give

    |(r_Q'/r_Q-ell)'|
       <= W exp[-(3/2)(Q-|x|)].

Indeed the numerator derivative contributes P2+4P0+2P1, and differentiating 1/r_Q contributes H1(P1+2P0), all divided by kappa. Integrating the exponential over [-Q,Q] gives at most 4/3. Therefore

    Var(r_Q'/r_Q)<=4+(4/3)W.                               (6.7)

To handle the chemical-potential term, differentiate Z_Q=1+T_Q Z_Q and integrate by parts in its integration variable. This gives

    Z_Q'(x)=Z_Q(Q)[R_Q(x,-Q)-R_Q(x,Q)].

Use (6.3)-(6.4) and |h(Q)|Z_Q(Q)=a r_Q(Q)<=a r0(Q). Combining r0(Q)/r0(x)<=2exp[-b(Q-|x|)] with (6.6) gives

    |h Z_Q'/r_Q| <= (4a C0/kappa) exp[-(3/2)(Q-|x|)],
    |(h Z_Q'/r_Q)'|
        <= (8a C0 H1/kappa) exp[-(3/2)(Q-|x|)].             (6.8)

In the derivative bound, one H1 comes from the resolvent derivative and one from |r_Q'/r_Q|. Finally

    v_Q=-a r_Q'/r_Q-h Z_Q'/r_Q.

Equations (6.7)-(6.8), their integrals, and a<=pi give Var(v_Q)<=Vv. This proves the bounded variation needed before inserting the finite counting density. The estimates persist at Q=infinity, where h=0 and r_Q=r0.

### 6.4 Finite-root errors, clipping, and tails

[Proved] Return to Q=Q_L and m=m(Q_L). The comparison in (3.3) is with the continuum density of mass nu. Restricted densities are pointwise ordered as their support interval grows, by their positive Neumann series. Consequently the sup difference of their uncentered cumulative functions is at most their mass difference. Equations (3.3) and (4.2) therefore imply, for the cumulative function G of mu-r_Q^occ dx,

    ||G||_infinity <= 3B0/(2L),
    |G(+infinity)|=|nu-m|<=B0/L.                            (6.9)

Let H be the cumulative function of sigma_L 1_[-Q,Q] dx-r_Q^occ dx. Since H=G-cumulative(R_L),

    ||H||_infinity<=2B0/L,
    H(+infinity)=nu-m.

Integrating T against G on the whole line has no boundary term even when G(+infinity) is nonzero, because T decays. On the occupied interval,

    ||sigma_L-r_Q||_infinity <= D/L,
    ||sigma_L'-r_Q'||_infinity <= D1/L.                     (6.10)

For the first bound use ||T'||_1<=2/(pi s0); for the second use ||T''||_1<=H2 q<=H2. The added 1 in D makes the free-point cutoff regular.

Suppose L>=4gamma D/kappa and choose A>=0 by

    r0(A)=2D/(kappa L),
    cosh(bA)=kappa L/(4gamma D).

Set A_eff=min(A,Q), and clip epsilon_Q to epsilon_Q(A_eff) outside [-A_eff,A_eff], but only within its occupied interval. Call the clipped function g. If Q<=A, no clipping is performed. On the core, (6.6) and (6.10) imply sigma_L>=r_Q/2. With d=sigma_L-r_Q,

    |(r_Q/sigma_L)'|
       <=4(H1 D+D1)/(L r_Q).

Moreover the explicit sech integral and (6.6) give

    integral_(-A_eff)^(A_eff) dx/(L r_Q(x))
       <=4gamma sinh(bA)/(L kappa b)<=1/(bD)<=1/D.

For f(t)=g(z_L^(-1)(t)), its derivative in the core is v_Q r_Q/sigma_L. Variation of a product, Var(v_Q)<=Vv, and |v_Q|<=2a give

    Var(f')<=2Vv+8pi(H1+D1/D)+8pi=Vclip.                   (6.11)

The last 8pi bounds the two jumps to zero at clipping points. If no clipping occurs these jumps are absent. The midpoint Peano estimate applies also when a cell crosses a kink, with its jump counted in the derivative variation measure.

If A<Q, the full-line tail mass is at most 4r0(A)/b. Using (6.9) at both boundaries, including one-sided limits at atoms, gives

    mu(|x|>A) <= integral_(|x|>A) r_Q^occ +4B0/L,
    integral_(A<|x|<=Q) sigma_L
               <= integral_(|x|>A) r_Q^occ +5B0/L.

The constants include the possible total mass difference nu-m. Thus the sum of the two tail masses is at most

    [16D/(kappa b)+9B0]/L.

Because epsilon_Q(Q)=0 and |epsilon_Q'|<=2a r_Q,

    |epsilon_Q(x)-epsilon_Q(A)|
       <=2a integral_A^Q r_Q <=4a r0(A)/b
       =8aD/(kappa bL),             A<=x<=Q.

The two tails therefore cost at most

    L |integral (epsilon_Q-g) R_L|
       <= a[128D^2/kappa^2+72D B0/kappa]/L.                (6.12)

If A>=Q this error is zero. Combining (4.3), midpoint quadrature with (6.11), (6.12), and the curvature bound (5.4) proves the large-L branch of K_all. The density correction is at most pi B0^2/(2L).

At exact half filling Q_L=infinity, R_L is the full-line midpoint discrepancy, m=nu=1/2, and h=0. The same clipped identities hold by integrability or by taking the preceding bounds to the full-line limit. The cumulative estimates and tail calculation remain valid. The density correction vanishes.

For the small-size branch use |F_L-L e_bulk|<=L as in Section 5, using C134 for the half-filled bulk identification. If the cutoff condition fails, then L<4pi D/kappa. This gives the first entry of K_all. Matched-density doubling proves the second bound in (6.1). This completes Task 1 for fixed delta>0, including all dilute sequences, without relying on the limited C129 dilute patch.

### 6.5 A simpler closed constant

[Proved] The following coarser version avoids the auxiliary constants in a theorem statement:

    K_all(delta) <= 10^15 / delta^(25/3).                   (6.13)

Here are explicit majorizations. Since alpha0>=sqrt(2delta), s0>=sqrt(delta), and

    1-Mbeta = 1-cos(alpha0/2)+sin(alpha0/2)
             >=sqrt(delta/2),

we have B0<=9/(8sqrt(delta)), C0<=2/delta, H1<=2/delta, H2<=12/delta^2, D<=3/delta, and D1<=21/delta^(5/2). Also exp(c)<=(11/delta)^(2/3). The third entry in the definition of kappa is smaller than its second entry by exp(-2c)/(16pi), and is smaller than 1/2. Hence kappa is that third entry. Using cosh(c+1)<=exp(c+1), e^2<8, pi<22/7, and 11^(8/3)<600 gives

    1/kappa <= 32pi e^2 11^(8/3) delta^(-19/6)
               < 500000 delta^(-19/6).

The Pj bounds are P0<=6/delta, P1<=12/delta^2, P2<=72/delta^3. It follows that

    W <= 84000000 delta^(-37/6),
    C0 H1/kappa <= 2000000 delta^(-31/6),
    Vv <= 430000000 delta^(-37/6),
    Vclip <= 900000000 delta^(-37/6).

In the last line use D>=1 and delta<=1. The dominant large-size term is bounded by

    128pi D^2/kappa^2
       <= (2.88 * 10^14) pi delta^(-25/3).

The other terms contribute at most

    [112500000+(121500000)pi+2] delta^(-25/3).

Their sum is below 10^15 delta^(-25/3). The small-size term is bounded by (6000000pi+1)^2 delta^(-25/3), which is also below that ceiling. This proves (6.13).

## 7. Where half filling entered C134

[Proved, dependency audit]

1. Its full-line inverse and explicit sech density in Section 3 describe mass 1/2. At a smaller density the occupied density solves a restricted equation. Full-line dressing alone uses the wrong bulk reference.
2. Its Section 4 residual mu-sigma has an O(1/L) cumulative norm only because the unused strips have width zero for even half filling, or at most 1/(2L) for the odd neighboring sector. At fixed nu<1/2, (3.1) gives strips of fixed positive width. Equation (3.3) replaces this step with projection onto the occupied interval.
3. Its derivative estimates compare to a globally smooth sech density. Here the occupied density has endpoint jumps. The smooth comparison is s_Q, not r_Q^occ. Restricting a convolution and integrating by parts without the endpoint terms would be wrong.
4. Its energy identity has no chemical-potential or moving-boundary terms. Equations (4.1)-(4.3) are the replacements, including the mismatch m(Q_L)-nu.
5. Its cutoff uses epsilon=-a r_infinity, the explicit tanh logarithmic derivative, and an explicit integral of 1/r_infinity. These make the clipped derivative variation independent of log L. A finite-interval dressed energy is not -a r_Q after imposing epsilon(Q)=0; it is -a r_Q-h Z_Q. The same cancellation cannot be asserted without proof.
6. Its tail masses are compared to the full-line density with mass 1/2. The occupied finite interval must instead be retained, or any clipped part must be estimated together with the Fermi boundaries.
7. Its odd lower bound discards the outside strips with a favorable sign, and its even upper bound uses their absence. A matched odd-density doubling requires both bounds at nu=(L-1)/(2L). C134's doubled sector n=L is different.
8. Its small-size spectral argument survives once the density-dependent bulk energy has been identified. Its doubling algebra survives only when the same bulk density is used at both lengths.

[Proved] Section 6 supplies the missing clipped variation and tail bounds uniformly in Q_L. This is the extra work beyond the compact finite-boundary estimate (5.7). C134 alone covers the half-filled boundary, not the neighboring sequences with L/2-n unbounded and (L/2-n)/L tending to zero.

[Proved, limitation of the originally proposed dilute patch] The C129 bound 2 pi^2 n^3/L is constant only when n=O(L^(1/3)). For example n=floor(L^(1/2)) has nu->0 but n^3/L diverges. A theorem for each fixed nu>=rho0>0 plus this dilute bound would leave that sequence uncovered. The restricted theorem (5.6) does cover it at fixed delta, because its constant is uniform down to nu=0.

## 8. The ferromagnetic crossover

[Proved, tracked failure of the C134 constants] Put epsilon=1+Delta and alpha=pi-gamma. Then epsilon=1-cos alpha, so alpha is asymptotic to sqrt(2 epsilon). For the actual C134 constants,

    B = pi/(2alpha),
    D = 1+cot(alpha)/alpha,
    D1 = B q [8/sin(alpha)^4+4/sin(alpha)^2],
    a = pi sin alpha,       b=pi/(pi-alpha).

As alpha tends to zero,

    D ~ alpha^(-2),       D1 ~ 4pi alpha^(-5),
    L_cut = 4gamma D ~ 4pi alpha^(-2),
    K_small=(4gamma D+1)^2 ~ 16pi^2 alpha^(-4),
    K_large=a[3b/2+D1/(2D)+32D^2/b+8DB]
                         ~ 32pi alpha^(-3).                 (7.1)

These are expansions of explicit formulas, not energy asymptotics. The compact majorant in C134 replaces several useful factors by constants and grows as epsilon^(-5/2); the pointwise K has dominant small-size term epsilon^(-2). Even dropping that term after L>=L_cut leaves K_large of order epsilon^(-3/2).

C129's band is epsilon L^2<=A, or L<=sqrt(A/epsilon). The C134 cutoff requires L of order epsilon^(-1), a parametrically larger scale. Even optimizing its small-size branch with the band at L of that latter order gives a ceiling of order epsilon^(-1), not a constant. The large-size energy bound is K_large/L; multiplying the defect by L cancels this 1/L. There is no remaining L factor that suppresses K_large as L increases. Thus retaining both delta and L in the existing estimates does not solve Task 2.

[Proved] The divergence of the full-line inverse norm is real in the norm being used: (I-T*)1=(1-q)1, so the L-infinity inverse has norm exactly B. Taking larger L cannot change that operator norm. Energy cancellation may avoid paying this norm, but an assertion that the same Neumann constant becomes uniform beyond the crossover is false.

[Assumed, exact missing crossover estimate] One sufficient new input is explicit A,K_* such that, for every -1<Delta<=0, every 0<=n<=L/2, and all L with (1+Delta)L^2>=A,

    |F_L(n)-J L e_bulk(Delta,n/L)| <= K_* J/L.              (7.2)

Then the matched-sector defect has ceiling 5K_* in units of M_0=J/2 beyond the crossover, and C129's band supplies A below it. This is stronger than necessary. The directly sufficient one-sided input is

    L[F_(2L)(2n)-2F_L(n)]/M_0 <= K_d
                   whenever (1+Delta)L^2>=A.                (7.3)

Neither (7.2) nor (7.3) is proved here. In particular the norm estimates of Section 3 and the a B^2 contribution in (5.5) still diverge. Better quadrature alone would not remove every loss in that estimate.

## 9. Assembly with the field, DM angle, and bond scale

[Proved, inherited reductions spelled out] The physical bond is

    h_bond = J_perp(XY) + D_DM(Sx Sy-Sy Sx) + J_z Sz Sz
                                  -(B_field/2)(Sz tensor 1+1 tensor Sz).

The rotation exp(-i theta sum_x x Sz_x), theta=atan2(D_DM,J_perp), moves its hopping phase to flux Phi=Ltheta. The doubled ring carries 2Phi. If n_* minimizes the length-L energy including the field, its doubled sector has exactly twice the field energy, so the field cancels from the comparison of numerators. No even-length assumption or differentiability in the field is used.

At a minimizing flux, every ground vector has zero current expectation. The exact uniform-gauge trial identity is

    H(a+u)-H(a)=(cos u-1)T_hop(a)+sin u I(a).

Using ||T_hop||<=JL/2 and choosing a flux difference of magnitude at most pi gives

    0<=e_L(n,Phi)-F_L(n)<=JL[1-cos(pi/L)]/2.

Only the doubled ring pays this cost. Therefore

    L(E_(2L)-2E_L)
       <= L[F_(2L)(2n_*)-2F_L(n_*)] + pi^2 M(B_field)/4.    (8.1)

The optimized physical bond scale is half the range of the four numbers

    J_z/4 +/- B_field/2,        -J_z/4 +/- J/2.

Indeed averaging a telescoping correction o tensor 1-1 tensor o over z rotations cannot enlarge spectral diameter. The remaining correction t(Sz tensor 1-1 tensor Sz) replaces +/-J/2 by +/-sqrt(J^2/4+t^2), which only enlarges the range. Hence t=0 minimizes it. In the present anisotropy interval M_0=J/2 and M(B_field)>=M_0. One must not set M(B_field)=J/2 at arbitrary field.

[Proved, restricted assembled theorem] If Delta in [-1+delta,0] and the selected minimizing sector, after spin flip, obeys n_*/L<=1/2-eta, then

    L(E_(2L)-2E_L) <= [5K_(delta,eta)+pi^2/4] M(B_field).   (8.2)

This follows from (5.6), M_0<=M, and (8.1). Independently, (1+Delta)L^2<=A gives the constant A+pi^2/4 for all fillings. The endpoint Delta=-1 has zero physical defect, and Delta=0 has the sharper C129 constant 8-4sqrt(2)+pi^2/4. Saturated sectors have zero physical defect. When J=0 the bond is classical and repetition of a minimizing cyclic spin word gives a nonpositive defect.

[Proved, assembled all-filling compact theorem] By (6.1), the sector restriction can be removed entirely at fixed delta>0. For every bond in the stated U(1) class with -1+delta<=Delta<=0, every field and DM angle, and every integer L>=2,

    L(E_(2L)-2E_L) <= [5K_all(delta)+pi^2/4] M(B_field).     (9.1)

This is a theorem on each compact attractive anisotropy interval with all fillings, all fields, and all lengths. It supplies an explicit, very large constant. It does not give a common constant as delta tends to zero.

In particular, (6.13) gives the single closed choice

    C(delta)=6 * 10^15 / delta^(25/3).                      (9.2)

No numerical estimate enters this constant. The smaller expression in (9.1) is available when useful.

[Proved, conditional full assembly] If (7.2) were established, the requested theorem would hold with

    C=max{A,5K_*}+pi^2/4.                                  (8.3)

Under (7.3), replace 5K_* by K_d. This is a conditional formula, not an explicit constant proved for the full class. All statements use every L>=2, the two-oriented-bond L=2 convention, and the physical M including the field. At L=2 the two DM terms cancel in the physical sum, consistent with the flux reduction. No vacuum-compatible subsequence is introduced.

[Proved, an explicit region expressed in bond parameters] The sector restriction in (8.2) can be eliminated in a high-field region without assuming an asymptotic magnetization formula. For each fixed delta>0, all bonds with

    -1+delta <= Delta <= 0,
    |B_field| >= (15/16) J(1+Delta)

satisfy the original all-length bound with

    C = max{28,5K_(delta,1/4)} + pi^2/4.                    (8.4)

This includes a nonempty interacting unsaturated interval below |B_field|=J(1+Delta), as well as the already solved saturated region. It is not the arbitrary-field theorem requested in the task.

Here is the sector-placement proof. Let t=-Delta in [0,1), let A=sum eta_x eta_(x+1), and let T_hop have J=1. C129's local positive-operator inequality at Delta=-1 says T_hop-A>=-n in sector n. Hence

    K_L(Delta)|_n = Delta L/4+t n+T_hop-t A
       = Delta L/4+t n+(1-t)T_hop+t(T_hop-A)
       >= Delta L/4+(1+Delta)T_hop.

The free-sector energy therefore gives

    F_L(n) >= J Delta L/4
              -J(1+Delta) sin(pi n/L)/sin(pi/L).            (8.5)

For B_field>=0, subtract the polarized-state energy. The physical sector energy is at least

    B_field n-J(1+Delta) sin(pi n/L)/sin(pi/L).

For L>=8 and n/L>=1/4, monotonicity of sin x/x and of x/sin x implies

    sin(pi n/L)/[n sin(pi/L)]
       <= [2sqrt(2)/pi] [pi/(8 sin(pi/8))]
       = cos(pi/8) < 15/16.

The last strict inequality follows from cos(pi/8)=sqrt(2+sqrt(2))/2, by squaring twice. Thus these sectors have energy strictly above the polarized trial energy whenever the stated high-field condition holds and Delta>-1. Spin flip treats negative field. Every minimizing sector then has n_*/L<1/4 at L>=8, so (8.2) applies with eta=1/4. For 2<=L<=7, take two copies of an L-ring sector ground state and reconnect the two seams. Each changed bond expectation costs at most its spectral diameter 2M_0, hence the sector defect is at most 4M_0. Multiplying by L gives at most 28M_0. Equation (8.1) finishes (8.4). The endpoints Delta=-1 and J=0 are already covered separately.

[Assumed, unresolved] Task 1 and the compact-interval assembly are proved. Task 2 and the requested uniform assembly through Delta=-1 remain unresolved. No finite C for the whole interval [-1,0] and all fields has been proved in this run. The optional positive-Delta task is not pursued because its prerequisite has not been met.
