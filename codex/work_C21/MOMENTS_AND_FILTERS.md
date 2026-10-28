# Local boundary moments and delayed Jacobi filters

## 1. A necessary locality condition

[proved] Let H be a nearest-neighbor open-cell Hamiltonian, let HG=0,
and let F_a act at its leftmost site. Its source moments are

    M_r[a,b]=G^dagger F_a H^r F_b G.

Induction using HG=0 gives

    H^r F_b G=ad_H^r(F_b)G,
    M_r[a,b]=G^dagger F_a ad_H^r(F_b)G.                    (1)

The induction step is H O G=[H,O]G. Each commutator can enlarge the
support by at most one site. Thus the physical operator compressed in
(1) is supported on the first r+1 sites. No frustration-free assumption
is used. For a nonzero code energy, replace H by H-E0.

[proved, restricted obstruction] Suppose the code subspace is a product
of local subspaces. Its isometry may include an arbitrary global change
of logical basis. Compressions of operators supported on disjoint sets
then commute: first compress using the product isometry, where the
operators act on disjoint tensor factors, and finally conjugate both
by the same logical unitary. Therefore left and right moments of orders
r and q commute entrywise whenever their (r+1)-site and (q+1)-site
neighborhoods are disjoint. This also applies to Hermitian combinations
of syndrome entries.

This is not a theorem about arbitrary entangled code subspaces.
Compressions of commuting physical operators need not commute on such
a subspace. In particular, the result does not rule out the history
encoding still sought by route A.

## 2. C19 violates this condition for a long product code

[proved using C19's weights] Write b_a=sum_e e X_ea for C19's first
filter moment. Its polynomial formula gives

    b_(j+1)=-Q_j(1)/sqrt(h_j)
            =-binom(j+2,j)/sqrt(h_j),                     (2)

which is nonzero for every j. Define the Hermitian logical operators

    R_u=sum_a b_a eta_(2a-2+u),
    L_w=-sum_(a,v) b_a J_a[w,v] gamma_(2a-2+v).

Up to fixed signs and the nonzero factor rho, these are Hermitian
combinations of the XZ and ZX entries of the two first source moments.
The Y-syndrome chemical-potential filter does not enter them.

The Majorana identities give

    [eta_u,gamma_v]=-2i delta_uv P,
    [R_u,L_w]=2i sum_a b_a^2 J_a[w,u] P.                  (3)

For C19's stated hopping target,
J_a[2,1]=(t_a/2)cos(theta_a)>0. Consequently [R_1,L_2] is nonzero.

In a physical cell of length ell>=4, the first two sites and last two
sites are disjoint. Equations (1) and (3) contradict a product-subspace
code with the original C19 source weights. This remains true after
any common logical basis rotation. The objection is stronger than
observing that a literal Jordan-Wigner string reaches a distant qubit.

[proved scope] This excludes the original weights with a long product
code and strictly boundary-supported sources. It does not exclude a
local model with different weights. The Clifford cursor in LOCAL_CURSOR.md
is an example with a product code whose early moments involve only
nearby logical operators, consistent with (1).

## 3. Delaying the filters

[proved] Choose positive integers n and s. Let k=n+s, E_e=e for
e=1,...,k, C_ef=1/(e+f), and p_e=1/k. Define

    Q_j^(s)(x)=sum_(l=0)^j binom(j+2s,l) binom(j+1,j-l)
                              x^l (x-1)^(j-l),
    h_j^(s)=(j+1)/[(2j+2s+2)(j+2s+1)],
    f_j^(s)(u)=x(1-x)^s Q_j^(s)(x)/sqrt(h_j^(s)),
    x=exp(-u),   j=0,...,n-1.

Expand f_j^(s)(u)=sum_(e=1)^k X_(e,j+1) exp(-e u). Then

    sum_e e^r X_ea=0       for r=0,...,s-1,
    X^T C X=I_n.                                           (4)

For completeness, the all-n proof uses the weight w_s=x(1-x)^(2s).
The displayed binomial polynomial satisfies

    Q_j^(s)=(-1)^j/[j! w_s] * d^j/dx^j
                       [x^(j+1)(1-x)^(j+2s)].

Its leading coefficient is binom(2j+2s+1,j). Integrating j times by
parts gives zero against all lower-degree polynomials. All endpoint
terms vanish because of the indicated powers. The norm is

    integral_0^1 w_s (Q_j^(s))^2 dx
      =binom(2j+2s+1,j) B(j+2,j+2s+1)
      =(j+1)/[(2j+2s+2)(j+2s+1)].

Changing variables x=exp(-u) yields (4)'s Cauchy orthogonality.
The zero of order s at u=0 gives its moment identities. The first
nonzero moment is

    sum_e e^s X_(e,j+1)
      =(-1)^s s! binom(j+2s,j)/sqrt(h_j^(s)).                (5)

This follows from x=1-u+O(u^2), Q_j^(s)(1)=binom(j+2s,j), and
the s-th derivative of the exponential expansion. At s=1 these are
C19's original filters and norms.

[proved optimality for this constraint] If all n independent filters
have their first s moments zero on k distinct energies, the corresponding
Vandermonde constraints have rank s, provided k>=s. Their nullspace
has dimension k-s. Thus k>=n+s is necessary. The family above saturates
that bound. It is a bound for simultaneous moment delay, not a universal
pole bound for every possible local encoding.

## 4. Fields, positive weights, and exact parity

[proved extension of C19] Use these X in C19's D_e^L and D_e^R. Their
sum vanishes and their Cauchy Gram matrix is unchanged, so the hopping
contraction remains rho^2 K_hop with arbitrary C19 phases and strengths.
The Clifford square-root and reflection arguments from C19 apply
without alteration once positivity is established.

To delay the chemical-potential filter as well, take

    d_e=(-1)^(e-1) binom(s,e-1) for 1<=e<=s+1,
    d_e=0 otherwise.

Its exponential sum is x(1-x)^s. Its overlap with p is

    g_s=p^T C d
       =(1/k) integral_0^1 sum_(f=1)^k x^f (1-x)^s dx
       =[B(2,s)-B(k+2,s)]/k >0.                            (6)

The last equality sums the geometric series and cancels one factor
of 1-x. Positivity also follows directly from the integral. At s=1,
(6) reduces to C19's 1/[2(k+2)].

For F=-sum_a mu_a Z_a/2, set

    sigma=-rho^2/g_s,
    M_e^R=p_e I+rho D_e^R,
    M_e^L=p_e I+rho D_e^L+sigma d_e |Y><Y| tensor F.

With r_max=max_e sqrt(2 sum_a X_ea^2), and C19's t_a<2, it suffices to take

    rho=min(1/(4k r_max),
            sqrt(g_s/(4k max_e|d_e| ||F||))).                (7)

Omit the second bound if F=0. The X,Z and Y blocks then each have
eigenvalues at least 3/(4k), and the weights sum to identity. C19's
entrywise contraction gives

    H2=a0 I+rho^2(K_hop+I tensor F),  a0=-3p^T C p.         (8)

The transpose on C19's left syndrome matrices is retained. Every
syndrome-odd coefficient multiplies a logical-odd operator, and F is
even, so exact product parity survives. The inherited finite-matrix
reflection has an exactly 2^n-dimensional zero-energy code and dimension
16(k+1)2^n. No locality is being claimed for that reflection.

## 5. What delay repairs, and what it does not

[proved] All moments below s are now scalar. At order s, the argument
of (3) applies with the nonzero coefficients (5). The delayed family
therefore still cannot be realized on a product-subspace code if the
two radius-s boundary neighborhoods are disjoint. In a cell of length
ell this gives the necessary condition

    ell < 2(s+1), or s >= ceil((ell-1)/2).                  (9)

Together with the simultaneous-delay bound, this requires
k>=n+ceil((ell-1)/2) for this particular product-code repair.
It does not say that these conditions suffice. The local spectral
weights, port geometry, and legality of a controller remain to be built.

[numerical and exact arithmetic] filter_checks.py verifies (4)-(6)
using Fraction at n=1,2,3,5,8 and s=1,2,4,8, twenty cases. Direct
matrix contractions at n=1,2,3 and s=1,2,3 verify (8), positivity,
normalization and parity. Absolute second-order matrix errors are
below 5.6e-16. Because rho^2 becomes small, the largest error divided
by rho^2 is 9.4e-7, at n=3,s=3. This conditioning is recorded rather
than hidden by reporting only absolute error. The exact arithmetic
and analytic identities establish the cancellations independently.
