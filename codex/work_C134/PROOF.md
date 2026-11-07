# C134 proof

2026-11-07. Status tags refer to mathematical claims, not to numerical certification. The proof below establishes a compact-parameter theorem on the attractive half of the massless interval. No independent referee has checked this new argument.

## 1. Statement and constants

[Proved] Use S=sigma/2, J>0, and C129's oriented ring convention at L=2. Let F_L be the twist minimum in sector n=floor(L/2). Spin flip identifies the other neighboring sector at odd L. Then, for every -1<Delta<=0, there is the following explicit constant, valid for every integer L>=2.

Set

    gamma = arccos Delta,             pi/2 <= gamma < pi,
    q = 2gamma/pi - 1,                B = 1/(1-q),
    V1 = 2|cot gamma|/pi,
    V2 = q[8/sin(gamma)^4 + 4/sin(gamma)^2],
    D = 1+B V1,                       D1 = B V2,
    b = pi/gamma,                     a = pi sin gamma,
    K = max{(4gamma D+1)^2,
             a[3b/2 + D1/(2D) + 32D^2/b + 8DB]}.

Here a is dimensionless. All energy calculations below set J=1 and restore it by homogeneity. The conclusions are

    F_L >= L e_0 - K J/L                         (all L),
    |F_L-L e_0| <= K J/L                         (even L),
    L(F_(2L)-2F_L) <= (5K/2) J                   (all L).       (1.1)

Thus the requested constant in units of J is C=5K/2. In units of M=J/2, the ceiling is 5K. The constants are deliberately coarse.

[Proved] A single closed expression works uniformly on [-1+delta,0], for any 0<delta<=1. Define

    alpha = arccos(1-delta),       s = sqrt(delta(2-delta)),
    Bbar = pi/(2alpha),            Dbar = 1+1/(alpha s),
    D1bar = Bbar(8/s^4+4/s^2),
    Kdelta = (4pi Dbar+1)^2
               + pi[3 + D1bar/2 + 32Dbar^2 + 8Dbar Bbar].

Then Cdelta=5Kdelta/2 is an explicit uniform constant in (1.1). In particular, the simpler choice

    Cdelta = 3000/delta^3                                    (1.2)

works. This covers an interval of interacting models at all sizes. It is not uniform as delta tends to zero and does not cover the positive-anisotropy interval. Section 11 separately covers Delta=1/2.

## 2. Flux and finite Bethe ground state, including odd lengths

[Proved, inherited from C129] The twist minimum equals the sector ground energy of the periodic stoquastic Hamiltonian

    K_L = sum_x [-Sx_x Sx_(x+1) - Sy_x Sy_(x+1)
                                  + Delta Sz_x Sz_(x+1)].     (2.1)

The physical antiferromagnetic flux is Lpi modulo 2pi. The configuration exchange graph is connected in every nontrivial sector. Perron-Frobenius gives a simple ground eigenvalue at all real Delta. There is no even-length restriction in this statement.

[Proved] The following finite Bethe construction on -1<Delta<=0 also has no parity restriction. For 1<=n<=L/2 define

    I_j = j-(n+1)/2,
    p(x) = 2 arctan[tanh(x) cot(gamma/2)],
    theta(x) = 2 arctan[tanh(x) cot(gamma)].

The equations are

    L p(x_j) - sum_k theta(x_j-x_k) = 2pi I_j.                (2.2)

Let P'=p and T'=theta, with P(0)=T(0)=0. Their Yang-Yang functional is

    Y(x) = L sum_j P(x_j) - sum_(j<k) T(x_j-x_k)
                                      - 2pi sum_j I_j x_j.

Since p'>0 and theta'<=0, its Hessian satisfies

    u^T Y'' u = L sum_j p'(x_j)u_j^2
                    - sum_(j<k) theta'(x_j-x_k)(u_j-u_k)^2 > 0

for every nonzero u. Existence requires coercivity, not just this Hessian test. The asymptotes are

    p(+infinity)=pi-gamma,
    theta(+infinity)=pi-2gamma,
    P(x)=(pi-gamma)|x|+O(1),
    T(x)=(pi-2gamma)|x|+O(1).

The O(1) terms are bounded for fixed gamma. Sorting coordinates increasingly cannot increase Y, by the rearrangement inequality for its last term. Every sorted vector is a nonnegative sum of rays that move the largest k coordinates together to the right, or the smallest k together to the left. On a right-moving ray, the recession slope is

    k[L(pi-gamma) -(pi-2gamma)(n-k) - pi(n-k)]
       = k(pi-gamma)[L-2n+2k] > 0.

The left-moving slope is the same. There are finitely many generating rays, so the recession function is positive on the unit sphere in the sorted cone. Hence Y is coercive and has one minimizer. Equations (2.2) hold there. Subtracting two equations shows that unequal ordered I's require unequal ordered roots. Uniqueness and reflection give x_j=-x_(n+1-j). The positive Hessian and the implicit function theorem give a real analytic root family for pi/2<=gamma<pi.

[Proved, coordinate Bethe verification] Put z_j=exp(i p(x_j)). The candidate energy is

    E_L = Delta L/4 + sum_j e(x_j),
    e(x) = -sin(gamma)^2/[2(sinh(x)^2+sin(gamma/2)^2)]
         = -Delta-cos p(x).                                  (2.3)

In coordinates y_1<...<y_n, use a sum of plane waves over permutations of the z_j. On adjacent particles, the local eigenvalue equation requires the exchange amplitude

    A_(...,k,j,...)/A_(...,j,k,...)
       = -[1+z_j z_k+2Delta z_k]/[1+z_j z_k+2Delta z_j].       (2.4)

This follows by comparing the missing two hopping terms at contact with the interaction Delta. More explicitly, the contact condition is

    [1+zw+2Delta w] A_12 + [1+zw+2Delta z] A_21 = 0.

The elementary rapidity identity

    [1+z_j z_k+2Delta z_j]/[1+z_j z_k+2Delta z_k]
          = exp(i theta(x_j-x_k))

makes the periodic condition z_j^L equal to (-1)^(n-1) times the product of these phases. This is (2.2), since exp(2pi i I_j)=(-1)^(n-1). Neither this boundary condition nor the local calculation assumes even L. Multiple contacts are handled by the same scalar exchange factors, which commute. At n=1 the plane wave verifies the equation directly, including the two-oriented-bond L=2 case.

At gamma=pi/2, the plane-wave sum is a nonzero Slater determinant occupying the n momenta centered at zero. It is the free sector ground state. All denominators in (2.4) are nonzero near that point, since at the free point every pair sum of occupied momenta lies strictly between -pi and pi. Thus the construction is a nonzero eigenvector on a neighborhood of the free point.

There is no need to assume that this particular normalization of the Bethe vector remains nonzero at every gamma. The finite-dimensional characteristic polynomial

    det[K_L(cos gamma)|_n - E_L(gamma) 1]

is real analytic on the connected interval. It vanishes on the neighborhood just constructed, and therefore everywhere by the identity theorem. The energy (2.3) is a continuous eigenvalue branch throughout the interval. It cannot leave the simple ground branch: a first departure would require a degeneracy of the ground eigenvalue, contradicting Perron-Frobenius. This proves E_L=F_L(n). It also supplies the odd-length identification needed below, without importing an odd-length version of Kozlowski's theorem.

## 3. Full-line dressing

[Proved] Define

    a1(x)=p'(x)/(2pi),       K(x)=theta'(x)/(2pi)
          = sin(2gamma)/[2pi(sinh(x)^2+sin(gamma)^2)].

Here K is a kernel, distinct from the scalar constant K in Section 1. Where confusion is possible, write the kernel as mathcal K. For pi/2<=gamma<pi,

    mathcal K<=0,             ||mathcal K||_1=q<1,
    integral a1 = 1-gamma/pi.

Consequently I+mathcal K* is invertible by its convergent Neumann series on bounded functions and on finite signed measures, with inverse norm at most B.

The full-line solution is

    rho(x)=1/[2gamma cosh(pi x/gamma)],
    (I+mathcal K*)rho=a1,             integral rho=1/2.        (3.1)

For an explicit check, with Fourier convention fhat(w)=integral exp(iwx)f(x)dx,

    a1hat(w)=sinh[(pi-gamma)w/2]/sinh(pi w/2),
    Khat(w)=sinh[(pi-2gamma)w/2]/sinh(pi w/2).

Their ratio a1hat/(1+Khat) is 1/[2cosh(gamma w/2)], whose inverse transform is (3.1). The values at w=0 follow by continuous limits. These transforms also follow by integrating the elementary hyperbolic kernels around a rectangle of height i pi.

Since e=-pi sin(gamma) a1, its dressed version is

    epsilon=(I+mathcal K*)^(-1)e=-a rho,
    a=pi sin gamma,           b=pi/gamma.                     (3.2)

Define

    e_infinity=Delta/4+integral e(x)rho(x)dx.

The estimates below prove convergence to e_infinity directly, so its identification with the half-filled bulk energy e_0 does not require an assumed finite-size expansion.

## 4. A quantitative cumulative root estimate

[Proved] From now on n=floor(L/2), nu=n/L. Define the finite root measure and centered counting function

    mu=(1/L)sum_j delta_(x_j),
    z(x)=p(x)/(2pi) -(1/L)sum_j theta(x-x_j)/(2pi),
    sigma=z'=a1-mathcal K*mu.

The negative sign of the kernel gives sigma>0 on the entire real line. Thus z is a global increasing coordinate, and the roots are its equally spaced midpoint nodes z(x_j)=I_j/L. Let z(+infinity)=c and z(-infinity)=-c. Direct evaluation gives

    2c-nu = (1-gamma/pi)(1-2nu)=:w.

For even L, w=0 and nu=1/2. For odd L, 0<=w<=1/(2L). The central interval [-nu/2,nu/2] consists of n cells of width 1/L, with the roots at their midpoints. The full range [-c,c] adds two strips of total width w.

Let G(x) be the cumulative signed mass of mu-rho dx from -infinity to x. Similarly let H be the cumulative mass of sigma dx-rho dx and R that of mu-sigma dx. Counting midpoint jumps gives

    ||R||_infinity <= 1/(2L)+w/2 <= 1/L.                    (4.1)

This bound includes values just before and just after every root and both infinite tails. Integrating sigma-rho=-mathcal K*(mu-rho) yields

    H=-mathcal K*G,             R=G-H=(I+mathcal K*)G.

Fubini is valid for these finite signed measures and the integrable kernel. The Neumann bound gives

    ||G||_infinity <= B/L,       ||H||_infinity <= qB/L.       (4.2)

In particular the centered counting function obeys

    |z(x)-integral_0^x rho(t)dt| <= 2qB/L.                  (4.3)

Equation (4.3) is a quantitative real-axis replacement for the rate-free convergence needed here. It is not an all-order expansion or a complex-strip estimate.

Differentiate the convolution once and twice, integrating the signed measure against its cumulative function. Then

    ||sigma-rho||_infinity <= B||mathcal K'||_1/L <= D/L,
    ||sigma'-rho'||_infinity <= B||mathcal K''||_1/L <= D1/L. (4.4)

For the first derivative, |mathcal K| decreases on x>0, so ||mathcal K'||_1=2|cot gamma|/pi=V1. For the second, let d(x)=sinh(x)^2+s^2, s=sin gamma. The inequalities

    d(x)>=s^2 cosh(x)^2,
    |d'(x)|<=2cosh(x)^2,       |d''(x)|<=4cosh(x)^2

give

    |mathcal K''(x)| <= |mathcal K(x)|(8/s^4+4/s^2).

Integrating proves the V2 bound in Section 1. The added 1 in D keeps the cutoff regular at the free point.

## 5. Exact dressed quadrature identity

[Proved] By the definition of sigma and (3.1),

    mu-sigma dx=(I+mathcal K*)(mu-rho dx).

The kernel is even. Pairing this equation with epsilon and using (3.2) gives the exact finite identity

    F_L-L e_infinity
      = L integral epsilon(x)[mu(dx)-sigma(x)dx].            (5.1)

Thus it is the quadrature error of the dressed energy in the finite counting coordinate. There is no asymptotic replacement of the roots in (5.1).

A direct midpoint estimate for epsilon composed with z^(-1) is not justified in the tails. At fixed finite L, sigma has a different tail exponent from rho. The following cutoff removes that issue at a cost of order 1/L in total energy.

## 6. Cutoff and bounded variation estimate

[Proved] Suppose L>=4gamma D and choose A>=0 by

    rho(A)=2D/L,
    cosh(bA)=L/(4gamma D).

On [-A,A], (4.4) gives sigma>=rho/2. Clip the dressed energy to its boundary value:

    g(x)=epsilon(x)             if |x|<=A,
         epsilon(A)             if |x|>A.

It is continuous, nonpositive, and constant in both tails. Set f(t)=g(z^(-1)(t)) on [-c,c], extended by its limiting values at the endpoints. In the central rapidity interval,

    f'(z(x)) = epsilon'(x)/sigma(x)
              = a b tanh(bx) rho(x)/sigma(x).

Write d=sigma-rho. Since |d|<=D/L and |d'|<=D1/L,

    |(rho/sigma)'|
       = |rho' d-rho d'|/sigma^2
       <= 4(bD+D1)/(L rho).                                 (6.1)

The derivative of tanh supplies the other term. Therefore

    integral_(-A)^A |d/dx [epsilon'/sigma]| dx
       <= integral_(-A)^A [2ab^2 sech(bx)^2
                           +4ab(bD+D1)/(L rho(x))] dx.

Here

    integral_(-A)^A dx/[L rho(x)]
      = 4gamma sinh(bA)/(bL) <= 1/(bD).

The first contribution is at most 4ab, the second at most 4ab+4aD1/D. The two jumps of f' to zero at the clipping points contribute at most 4ab. Monotone reparametrization preserves total variation. Consequently

    Var(f') <= 12ab+4aD1/D.                                 (6.2)

For A=0, f is constant and the same bound holds. This is the step that avoids a logarithmic loss from the growing rapidity interval.

## 7. Midpoint quadrature and tails

[Proved] If a continuous function has derivative of bounded variation, midpoint quadrature with cell width h has error bounded by h^2 Var(f')/8. On a single cell the Peano kernel for the functional h f(midpoint)-integral f has magnitude at most h^2/8. Integrating it against the signed derivative measure df' proves the claim; summing cells adds the variation measures without an extra factor equal to the number of cells.

Applying this with h=1/L on [-nu/2,nu/2] gives

    |sum_j g(x_j)-L integral_(-nu/2)^(nu/2) f(t)dt|
       <= a[3b/2+D1/(2D)]/L.                                (7.1)

For even L this is the entire counting interval. For odd L,

    sum_j g(x_j)-L integral sigma(x)g(x)dx
      = [the midpoint error] - L integral_(outer strips) f(t)dt.

Because f<=0, the last term is nonnegative. Hence (7.1) is a lower bound on the full clipped discrepancy for both parities. Its absolute-value form is needed only for even lengths.

It remains to bound the cost of clipping. The hyperbolic density obeys, for A>=0,

    rho(A+t)<=2rho(A)exp(-bt),      t>=0,
    rho({|x|>A})<=4rho(A)/b=8D/(bL).                         (7.2)

Both mu and sigma dx have total mass at most 1/2:

    integral mu=nu,
    integral sigma=(1-gamma/pi)+q nu <=1/2.

Using (4.2) at the two tail boundaries therefore gives

    mu({|x|>A}) <= rho({|x|>A})+2B/L,
    integral_(|x|>A) sigma <= rho({|x|>A})+2qB/L.

Endpoint atoms may be included by one-sided limits. Since q<=1, the sum of these two tail masses is at most

    [16D/b+4B]/L.                                           (7.3)

Moreover |epsilon-g|<=a rho(A)=2aD/L, and epsilon-g vanishes on the core. Thus

    L |integral (epsilon-g)(mu-sigma dx)|
       <= a[32D^2/b+8DB]/L.                                (7.4)

Combining (5.1), (7.1), and (7.4) proves the large-L part of (1.1), with

    Klarge=a[3b/2+D1/(2D)+32D^2/b+8DB].

The same argument proves convergence F_L/L -> e_infinity, also at odd L, although only the lower estimate is required there. For example the omitted outer strips have total width O(1/L), g is bounded, and the remaining terms vanish after dividing by L. We henceforth write e_0=e_infinity.

## 8. Small sizes, uniformity and doubling

[Proved] Every bond in (2.1) has spectral diameter 1 for |Delta|<=1. Every finite-ring energy density lies in its bond spectral interval, including sector minima. The limiting density e_0 lies in the same interval. Thus

    |F_L-L e_0|<=L.

If L<4gamma D this is at most (4gamma D+1)^2/L. Taking the maximum with Klarge proves the first two statements in (1.1) at every L>=2. Restoring J multiplies all energies by J.

The doubled length 2L is even, so its upper estimate and the length-L lower estimate yield

    F_(2L)-2F_L <= KJ/(2L)+2KJ/L = 5KJ/(2L).

The different sectors at odd L cause no density mismatch in this argument, since the bulk reference is the half-filled e_0 at both lengths.

For the compact bound, gamma<=pi, sin gamma>=s, B<=Bbar, D<=Dbar, D1<=D1bar, 1<=b<=2, and a<=pi. Also D>=1. These inequalities majorize each term in Section 1 and prove the Kdelta expression there. No parameter-dependent unspecified threshold remains.

For the simpler (1.2), use alpha>=sqrt(2delta), s>=sqrt(delta), Bbar<=9/[8sqrt(delta)], Dbar<=2/delta, and D1bar<=27/[2delta^(5/2)]. The first bound follows from 1-cos(alpha)<=alpha^2/2. Since delta<=1,

    Kdelta <= [(8pi+1)^2 + 623pi/4]/delta^3
            < 1200/delta^3.

The last inequality follows, for example, by pi<22/7, which bounds the bracket by 114949/98<1200. Multiplication by 5/2 proves (1.2).

## 9. Free point, endpoint, and restricted bands

[Proved, C129 plus the odd half-filled specialization] At Delta=0,

    F_L=-J csc(pi/L)                           for even L,
    F_L=-J/[2sin(pi/(2L))]                     for odd L.

Thus the odd half-filled doubling defect is zero. For even L,

    L(F_(2L)-2F_L)/J
       = L tan(pi/(4L)) sec(pi/(2L)) <= 4-2sqrt(2).

The last inequality follows since tan x/x and sec(2x) increase on 0<x<=pi/8. Equality holds at L=2. This sharper free constant supersedes the coarse Section 1 constant at Delta=0.

[Proved, inherited] At Delta=-1, every minimizing-flux sector has energy -JL/4 and the defect is zero. C129's bound for matched sectors is Q<=(Delta+1)L^2, with Q in units of M. On odd half-filled doubling, an equal-amplitude half-filled state on 2L and C129's sector lower bound give instead

    F_(2L)(L)-2F_L((L-1)/2)
       <= J(Delta+1)[L(L-1)/(2L-1)-1],

which is no greater than the even relaxed ceiling J(Delta+1)L/2. Thus Q<=(Delta+1)L^2 remains an upper bound for the present half-filled defect at both parities. This only controls a size-dependent band, not a fixed neighborhood of -1.

## 10. Why the other routes do not close the remaining interval

[Proved, obstruction to this argument] For 0<Delta<1, mathcal K is positive. Then sigma=a1-mathcal K*mu is not guaranteed positive on the whole line. The midpoint representation in a globally increasing counting coordinate and (4.1) cannot be used without a new estimate. The norm inequality ||mathcal K||_1<1 alone does not repair that missing input. This proof therefore does not establish the theorem on a compact interval crossing Delta=0.

[Proved, logical limitation] Even when D(L) is nonnegative and nonincreasing, a value D(L_0)<=A/L_0 gives D(L)<=A/L_0 for L>=L_0, not A/L. For a scalar illustration, E(L)=eL-b L^(-eta), 0<eta<1, is concave and has positive decreasing D(L)=b(2-2^(-eta))L^(-eta), but L D(L) diverges. This is a counterexample to the inference, not an XXZ counterexample. C125 provides a stronger spectral countermodel even with prescribed infrared behavior. Its positive spatial transfer construction is an even-ring statement for Delta>=0 and does not automatically apply to a fixed magnetization sector or to odd minimizing fluxes. Entrywise positivity of a six-vertex transfer matrix is not positive semidefiniteness of the crossed transfer operator.

[Proved, conditional endpoint extension] If one proves a common finite C for all Delta in an interval (Delta_*,1), uniformly in L, then the same C holds at Delta=1. For each fixed L, finite-matrix spectral continuity gives F_L(Delta)->F_L(1). Take this limit in the finite-L inequality. No monotonicity in Delta is needed. A bound C(Delta) diverging at 1, or constants on compact intervals whose upper endpoint stays below 1, is insufficient.

Indeed, Weyl's inequality only gives

    |F_L(Delta)-F_L(Delta')| <= JL |Delta-Delta'|/4,
    |L D_L(Delta)-L D_L(Delta')| <= J L^2 |Delta-Delta'|.

This does not control the supremum in L by continuity. Ground energy concavity in Delta also gives no sign for D_L', which is the difference of two longitudinal correlations at different lengths. Neither the requested monotonicity nor a separate Heisenberg proof has been established here.

[Assumed, unresolved] The open parts are 0<Delta<=1 apart from the isolated point Delta=1/2 treated next, and a constant uniform as Delta approaches -1 from above. The compact theorem on [-1+delta,0] cannot be patched to C129's crossover band without an additional estimate in the intervening scales. The numerical ceiling Q<=4 is not proved for the whole interval.


## 11. An isolated positive-anisotropy point from an exact twisted energy

[Proved using a published exact-energy identity] At Delta=1/2 in our convention, the bound holds for all L with

    L(F_(2L)-2F_L) <= (2pi^2/9)J.                            (11.1)

The input is the supersymmetric XXZ energy. [Batchelor, de Gier and Nienhuis, cond-mat/0101385](https://arxiv.org/pdf/cond-mat/0101385), equations (1)-(2) and the paragraph following (6), state the ground energy -3L/4 for their Hamiltonian -2 sum[XY-(1/2)ZZ], with periodic boundary conditions for odd L and twist 2pi/3 for even L. The half-filled or neighboring sectors are the ones involved. Dividing their Hamiltonian by two gives (2.1) at our Delta=+1/2 and energy -3JL/8. Their wavefunction enumeration claims are conjectures, but this energy identity is stated separately as known. [Yang and Fendley, cond-mat/0404682](https://arxiv.org/html/cond-mat/0404682), Section 2.2, explains the supersymmetric ground-state mechanism, with care about translation sectors and twists.

Thus, writing e_*=-3J/8,

    F_L=e_*L                                           (odd L),
    e_L^(L/2)(Phi_m+2pi/3)=e_*L                        (even L).

At even L, apply C129's twist bound to the actual minimizing-flux ground state:

    0 <= e_*L-F_L <= (JL/2)[1-cos(2pi/(3L))]
                         <= pi^2 J/(9L).

In particular F_(2L)<=2Le_* by flux minimization. Subtracting twice the lower bound proves (11.1) at even L. At odd L the defect is nonpositive. For L=2 the special-twist energy is also verified directly: it is -J/4-J cos(pi/3)=-3J/4. This removes any ambiguity from the oriented two-bond convention.

This is a consequence of the exact supersymmetric energy and the elementary twist estimate, not a claim of a newly proved supersymmetry theorem. It does not extend to a neighborhood of Delta=1/2 by finite-matrix continuity alone.
