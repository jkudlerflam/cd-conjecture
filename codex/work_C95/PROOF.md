# C95. Thermodynamic charge comparison

[Status] C92 equation (26) is not proved here. This note bounds the signed spatial sum of the canonical contact kernel uniformly in volume using exact free-fermion form factors, bounds the thermodynamic vacancy-number contribution, proves that the adjacent-sector crossing contribution to the first-order finite-volume residual is subextensive, and identifies a separate obstruction to a temperature-uniform Taylor-remainder argument. The tagged-contact spatial first absolute moment and the full relative thermal estimate remain open. No full mobile defect theorem for d=13 is asserted.

[Convention] A tag governs the following argument until the next tag. Proved denotes a derivation supplied here. Imported denotes a result retained with its original scope. Conditional denotes an implication whose additional hypothesis is stated. Numerical denotes a finite calculation. Open denotes an estimate not established here. All energies use t>0. All writes and Python execution stay in this folder.

## 1. The requested statement and its distinct inputs

[Imported, C92 equations (25)-(26)] In ordered material coordinates,

    H_sigma(J,s) = H_0 + s V + J X_sigma,
    H_0 = -t sum_j (E_j^+ E_(j+1)^- + E_j^- E_(j+1)^+),
    V = sum_j g_j,
    X_sigma = sum_j A_j h_sigma,j,
    A_j = |0><0|_j,
    Hmat_sigma = sum_j h_sigma,j,
    ||h_sigma,j|| <= M.

The frozen subtraction is

    R_sigma = Fmob_sigma - Fgap - Fmat_sigma(J a_0),
    a_0 = <A_j>_(J=0,T,s).

The target, uniformly in the length pattern and in a compact source interval, is

    R_sigma - R_good >= -J epsilon(J/t) W(sigma) - B_N,
    epsilon(x) -> 0,
    B_N/N -> 0 at every fixed T>0.                         (1)

The thermodynamic limit precedes T->0. W is the total material length in wrong cells plus the number of pinned selected cells. In the application b=15 and d=13. C92 already treats selected-cell pins, so the wrong-length estimate can be sought without those pins.

[Imported, C94 section 7] At fixed vacancy number, for zero-port product clock ground states,

    chi_N(r) = <0|(A_0-a) Q(H_0-e_0)^(-1)Q(A_r-a)|0>,
    chi_N(r) = k_N(r) - c_N/N,
    K_N = sum_r (1+dist_N(r,0)) |k_N(r)|,

and

    |B_sigma-B_good| <= 4 M^2 b W (K_N+|c_N|).             (2)

This is a conditional estimate on a ground-state second-order coefficient. Its premises do not include the thermal auxiliary trace or the nonzero port.

[Proved, distinction of limits] Even a proof of sup_N t(K_N+|c_N|)<infinity would establish only the coefficient estimate in (2). A finite-volume remainder O_(N,sigma)(J^3) need not be uniform in N, T, or the wrong-cell length. Moreover, the Taylor coefficient at J=0 and positive T uses normalized trace on the material space, while (2) uses a frozen ground vector. Section 6 gives an explicit example in which this difference prevents a uniform Taylor remainder after the frozen subtraction.

## 2. Exact free equation of state and a uniform density curvature

[Imported, C94 equation (7)] With v vacancies, L=N+v, and rho=N/L, the free ground energy and mean contact are

    E_N(v) = -2t sin(pi N/L)/sin(pi/L),
    a_N(v) = rho - sin^2(pi rho)/(N L sin^2(pi/L)).         (3)

[Proved] Put q=v/N=(1-rho)/rho. On every compact subinterval of 0<rho<1, (3) gives

    E_N(v) = N e(q) + O(t/N),
    a_N(v) = a(q) + O(N^(-2)),
    e(q) = -2t sin(pi rho)/(pi rho),
    a(q) = rho - sin^2(pi rho)/(pi^2 rho).                 (4)

The bounds, and any fixed number of derivatives of their smooth interpolations in q, are uniform on that compact interval. To obtain them, use

    z/sin z = 1 + z^2/6 + O(z^4),  z=pi rho/N,

in both formulas (3). In particular the finite-volume energy error is O(t/N), not an extensive error.

Since d rho/dq=-rho^2, differentiation yields

    e'(q) = (2t/pi)[pi rho cos(pi rho)-sin(pi rho)],
    e''(q) = 2 pi t rho^3 sin(pi rho) > 0,
    a_q = -rho^2 + (2rho/pi) sin(pi rho)cos(pi rho)
                         - sin^2(pi rho)/pi^2.            (5)

The equilibrium tilt is

    s = -e'(q)
      = (2t/pi)[sin(pi rho)-pi rho cos(pi rho)].           (6)

Its derivative with respect to rho is 2 pi t rho sin(pi rho)>0. Thus the free minimizer is unique at every interior source.

For t/2<=s<=3t/4, it lies strictly between rho=1/3 and rho=2/3. Indeed (6) at those two endpoints gives respectively

    s(1/3)/t = sqrt(3)/pi - 1/3 < 1/2,
    s(2/3)/t = sqrt(3)/pi + 2/3 > 3/4.

Throughout that larger density interval,

    e'' >= pi sqrt(3) t/27,
    |a_q| <= 4/9 + 2/(3pi) + 1/pi^2 < 1.                 (7)

Here the second bound uses |2 sin x cos x|<=1. Consequently

    0 <= c_dens(q) := a_q^2/(2 e'') < 3/t.                (8)

This is a uniform density-mode estimate, not a bound on the full contact kernel.

At half filling q=1,

    s=2t/pi,
    a=1/2-2/pi^2,
    e''=pi t/4,
    a_q=-(1/4+1/pi^2),
    c_dens = (2/(pi t))(1/4+1/pi^2)^2.                   (9)

The factor 1/2 in (8) is the one-sided ground-state resolvent convention. A full static Duhamel susceptibility has twice this coefficient.

## 3. What the thermodynamic vacancy-number average does

### 3.1 Exact finite-temperature decomposition

[Proved] Let Z_(N,v)(T) be the free trace in a fixed vacancy sector, and define normalized trace weights

    w_v = exp(-s v/T) Z_(N,v)(T) / Z_N(T,s).

Write a_v(T)=<A_j>_(N,v,T). Translation symmetry makes this independent of j. Since A_i preserves v, the connected Duhamel kernel in the tilted ensemble splits as

    Ktilt_ij = sum_v w_v Kcanon_ij(v)
               + (1/T)[sum_v w_v a_v(T)^2
                              -(sum_v w_v a_v(T))^2].     (10)

This follows by splitting the trace over v before integrating imaginary time and then adding and subtracting sum_v w_v a_v(T)^2. It uses no approximation. The last term is independent of i,j. At a finite-volume two-sector crossing it is the term isolated by C94.

### 3.2 Thermodynamic Hessian identity

[Proved] The following identity expresses the density part without needing the full canonical contact susceptibility. Let f(T,q,lambda) be the canonical free energy per material site when a uniform scalar source lambda sum_i A_i is added. At lambda=0,

    f_lambda(T,q,0) = a(T,q).

At an interior, stable minimizer define

    g(T,s,lambda) = min_q [f(T,q,lambda)+s q].

Implicit differentiation of f_q+s=0 gives

    q_lambda = - f_qlambda/f_qq,
    g_lambdalambda = f_lambdalambda - a_q(T,q)^2/f_qq.
                                                               (11)

Thus the full static response per material site gains

    a_q(T,q)^2/f_qq(T,q)                                   (12)

when vacancy number can relax. This is the thermodynamic version of the constant term in (10). The one-sided convention assigns half of (12) to c_dens(T,q).

For completeness, (11) requires differentiability only at lambda=0, at each fixed T>0. The uniform scalar contact is local in physical free-fermion coordinates, as shown in section 5. Its fixed-T two-point functions are finite sums of Wick contractions of thermal propagators. The Fermi function is analytic in a strip about real momentum, giving exponential spatial decay, while imaginary time ranges over the finite interval [0,1/T]. The scalar response coefficients at lambda=0 are therefore finite. No temperature-uniform convergence radius in lambda is asserted. Equation (11) can equivalently be read as the Schur complement of these thermodynamic response coefficients.

### 3.3 Evaluation at positive temperature

[Proved] The free physical-fermion dispersion is epsilon(k)=-2t cos k. Diagonalizing the hopping by Fourier transformation gives

    P(T,mu) = T integral_(-pi)^pi
                 log[1+exp((mu+2t cos k)/T)] dk/(2pi),
    n(k) = 1/[1+exp((-2t cos k-mu)/T)],
    rho = integral n(k) dk/(2pi),
    g_1 = integral cos k n(k) dk/(2pi),
    D = d rho/d mu = integral n(k)(1-n(k)) dk/(2pi T).
                                                               (13)

The physical canonical free energy is f_ff=mu rho-P. The free gap free energy per material site is

    f(T,q,0) = f_ff(T,rho)/rho.                            (14)

The root multiplicity and the boundary twist do not affect this thermodynamic identity. This is the free version of the root and bounded-seam transfer used in C92. Differentiating (14), rather than identifying finite-ring spectra, gives

    f_q = -P,
    f_qq = rho^3/D,
    s = P(T,mu).                                         (15)

The contact sum counts adjacent occupied pairs. Wick contraction gives

    a(T,rho) = [rho^2-|g_1|^2]/rho,
    a_rho = 1+g_1^2/rho^2 - (2g_1/rho)(g_(1,mu)/D).
                                                               (16)

Here g_1 is real. Combining (12), (15), and a_q=-rho^2 a_rho yields

    c_dens(T,s) = rho D a_rho^2/2.                        (17)

As T->0, uniformly over a compact interior interval of mu,

    D -> 1/[2pi t sin(pi rho)],
    g_1 -> sin(pi rho)/pi,
    g_(1,mu)/D -> cos(pi rho).

These limits can be obtained by changing variable from k to epsilon near the two Fermi points. The kernel n(1-n)/T has total integral one in energy, its first centered moment vanishes, and its second moment is O(T^2). The smooth density of states and its derivatives are bounded away from the band edges. Contributions outside fixed neighborhoods of the Fermi points are exponentially small. This also proves uniform convergence of the displayed derivatives.

The source interval in section 2 stays a positive distance from both band edges. Equations (7)-(8) and uniform convergence therefore imply the existence of T_0>0 such that

    sup_(t/2<=s<=3t/4, 0<T<T_0) t c_dens(T,s) <= 4.        (18)

No explicit optimal value of T_0 is needed for this statement.

[Scope] Equations (10)-(18) control the conserved-density part of the response. They do not identify an all-distance decomposition chi_N=k_N-c_N/N with bounded K_N. In particular, setting c_N equal to (17) by definition does not establish that the remainder k_N has a summable first spatial moment. That would be a separate theorem about the tagged contact.

[Conditional, quadratic comparison] In the zero-temperature product-clock comparison, if the canonical kernel has the decomposition and locality required by C94, density relaxation adds c_dens(q)/N to its one-sided response. The associated change is c_dens(q) E_sigma^2/N. C94 proves

    |E_sigma^2-E_good^2|/N <= 4 M^2 b W.

Thus the additional thermodynamic density contribution to the difference of quadratic coefficients is at most

    4 M^2 b W c_dens <= 12 M^2 b W/t                     (19)

in the source window of section 2. Equation (19) does not supply the missing canonical part. The finite-temperature bound (18) does not extend the product-clock argument to the full thermal material trace.

## 4. Adjacent finite-volume crossings contribute only O(J) at first order

[Proved] Let v/N remain in the compact interval of section 2. Taylor expansion of the smooth interpolations in (3)-(5) gives

    a_N(v)-a_N(v+1) = -a_q(v/N)/N + O(N^(-2)),
    s_N := E_N(v)-E_N(v+1)
         = -e'(v/N)-e''(v/N)/(2N)+O(t/N^2).               (20)

All constants are uniform on the compact interval. In particular

    |a_N(v)-a_N(v+1)| <= C_a/N                            (21)

for a fixed C_a and all sufficiently large N. The finitely many smaller N can be included by enlarging C_a.

At an isolated crossing of these two ground lines, C94's first-order residual is

    R_sigma(J,0,s_N)
       = -J |a_N(v)-a_N(v+1)| |E_sigma|/2 + O_N,sigma(J^2).
                                                               (22)

Since ||Hmat_sigma||<=MN, one has |E_sigma|<=MN. Therefore the magnitude of the displayed first-order term is at most J M C_a/2, uniformly over all material patterns. The difference of two such first-order terms is at most J M C_a. It is a global O(J) term, so it is subextensive.

This proves that the specific linear crossing mechanism from C94 cannot by itself produce an extensive violation of (1) along an interior-density sequence. It does not bound the remainder in (22). The assertion is about the explicitly displayed first-order term only.

At v=N, (20) gives the two useful limits

    s_N -> 2t/pi,
    N[a_N(N)-a_N(N+1)] -> 1/4+1/pi^2.                    (23)

[Scope] Taking T->0 at a fixed crossing still produces C94's divergent pointwise response. Taking N->infinity at fixed positive T instead gives the finite density response (12). These statements agree, but neither order of limits establishes the spatial first-moment bound or the full nonlinear comparison.

## 5. What the exact free-fermion contact calculation controls

### 5.1 A uniform bound on the signed spatial sum

[Proved] In the canonical gap ground state, let C_gap=sum_j A_j and

    B_C = <0|(C_gap-<C_gap>) Q(H_0-E_0)^(-1)
                                      Q(C_gap-<C_gap>)|0>.

Translation symmetry gives

    B_C = N sum_r chi_N(r).                              (24a)

There is an exact isometry between the cyclic-label-invariant gap subspace and the translation-invariant subspace of the physical hard-core hopping chain at L=N+v and N particles. Here is the normalization and the reason it intertwines the operators. A physical occupation set X has N choices of occupied root. A gap tuple g has L choices of physical root position. The resulting rooted objects are in bijection. If psi(g) is invariant under cyclic permutation, its value does not depend on which particle is chosen as root. Therefore

    sum_X |psi(g(X))|^2 = (L/N) sum_g |psi(g)|^2,
    (U psi)(X) = sqrt(N/L) psi(g(X)).

The image is translation invariant, and the inverse follows by reading the cyclic gaps. Moving a physical particle by one empty site transfers one vacancy between neighboring gaps. When the root changes, cyclic invariance removes the root choice. Hence U intertwines both hopping Hamiltonians. It also intertwines C_gap with

    C_ff = sum_x n_x n_(x+1).

The gap ground vector and C_gap applied to it are cyclic invariant. The resolvent in (24a) stays in this invariant subspace. Thus B_C is the physical free-fermion contact resolvent, with the parity boundary twist obtained from the hard-core chain. This statement is about the invariant subspace used by C_gap, not the full finite-ring spectrum for arbitrary tagged sources.

The filled sea has momenta k_j=2pi(j+theta)/L, theta=0 for odd N and theta=1/2 for even N, and dispersion epsilon_k=-2t cos k. Fourier transformation of C_ff gives, between the sea and an excitation with two holes k_1,k_2 and two particles p_1,p_2,

    <p_1 p_2; k_1 k_2|C_ff|FS>
       = (2/L)[cos(p_1-k_1)-cos(p_2-k_1)]                 (24b)

up to the basis sign, subject to p_1+p_2=k_1+k_2 modulo 2pi. This is obtained by antisymmetrizing the nearest-neighbor Fourier potential 2cos q. The magnitude is at most 4/L. There is no off-diagonal one-particle-hole term, because C_ff is translation invariant and such a term would require p=k. The remaining diagonal term is removed by Q.

Let F be the occupied sea and define

    S_L = L^(-2) sum_(k in F, p outside F)
                                      1/(epsilon_p-epsilon_k).

Sum the squared matrix elements (24b) divided by the excitation energies. Overcounting with ordered holes and dropping the requirement that the second particle lie outside F only increases the positive sum. Since each particle-hole energy is positive,

    1/(epsilon_(p1)+epsilon_(p2)-epsilon_(k1)-epsilon_(k2))
                    <= 1/(epsilon_(p1)-epsilon_(k1)).

For fixed k_1,p_1 there are at most L choices of k_2, and momentum conservation determines p_2. It follows that

    B_C/L <= 16 S_L.                                    (24c)

To bound S_L, write k and p for the absolute momenta in [0,pi]. If rho=N/L is between 1/3 and 2/3, then k<=k_F=pi rho<=p and

    epsilon_p-epsilon_k
      =4t sin((p+k)/2) sin((p-k)/2)
      >= (2t/pi)(p-k).                                  (24d)

The first sine is at least 1/2 throughout this domain. The second inequality uses sin x>=2x/pi on [0,pi/2]. The occupied and empty absolute momenta adjoining the Fermi surface are separated by 2pi/L. Index their distances from those edges by r,s>=0. Then p-k=2pi(r+s+1)/L, and each absolute pair has multiplicity at most four. Consequently

    S_L <= (1/(tL)) sum_(r=0)^R sum_(s=0)^S 1/(r+s+1)
        <= (R+S+1)/(tL) <= 1/t.                         (24e)

For the middle inequality, group terms by r+s. The number in each group is at most r+s+1. The last inequality follows because the two radial lists together contain at most L distinct radial levels. This bound includes endpoint levels of multiplicity one.

Combining (24a), (24c), and rho>=1/3 proves the finite-volume bound

    0 <= sum_r chi_N(r) = B_C/N <= 48/t.                 (24f)

No lower bound on the many-body bath gap was used. The constant is deliberately loose. This is a bound on the zero-momentum eigenvalue of the circulant susceptibility matrix. It does not control its other eigenvalues, absolute row sum, or first absolute spatial moment.

### 5.2 Why the tagged contact needs an additional estimate

[Proved] Let x_1<...<x_N be ordered material positions on a physical ring, and let n_x be the free-fermion occupation operator. On an occupation basis vector,

    A_j = 1_(x_(j+1)-x_j=1),
    sum_j A_j = sum_x n_x n_(x+1).                       (24)

The second identity follows by counting every adjacent occupied pair once. The first is a tagged-particle observable. Specifying j selects an ordered particle, not a fixed physical coordinate.

The distinction already appears at equal time. Put D_x=n_x n_(x+1), and for an oriented arc from x to y let

    Q_(x,y] = sum_(z in (x,y]) n_z.

For 1<=r<=N-1 the following diagonal operator identity holds:

    sum_j A_j A_(j+r)
      = sum_x sum_(y forward from x, y!=x)
                    D_x D_y 1_(Q_(x,y]=r).               (25)

If D_x D_y is nonzero, x and y are occupied and begin contacts. Their material-label separation is the number of occupied sites in (x,y], which is Q_(x,y]. This proves (25) configuration by configuration.

The counting projector can be written as

    1_(Q=r) = integral_0^(2pi) exp[i theta(Q-r)] dtheta/(2pi).
                                                               (26)

For a quasi-free state with one-particle covariance C, including the filled sea or a free grand-canonical thermal state, its counting generating function on a physical interval I is

    <exp(i theta Q_I)>
      = det[1-C+C exp(i theta P_I)].                      (27)

To see (27), write the quasi-free density operator as the normalized second quantization of exp(-h), use Tr Gamma(B)=det(1+B), and divide det(1+exp(-h)exp(i theta P_I)) by det(1+exp(-h)). Density insertions in (25) are obtained by differentiating additional diagonal counting fields.

Thus even the equal-time tagged-contact correlation involves counting projectors or determinants whose interval size grows with the physical separation. The imaginary-time tagged-contact correlation additionally follows material order through the evolution. Replacing it by the two-point function of D_x at a fixed physical separation discards this information.

[Open, first missing estimate] No uniform estimate for that imaginary-time tagged determinant or its form-factor representation is established here. In particular, this note does not prove

    sup_(N, interior v/N) t sum_r (1+dist_N(r,0))
                         |chi_N(r)+c_N/N| < infinity.    (28)

An exact free dispersion and Wick's rule for local physical observables do not, without this additional step, prove (28). The number-mode calculation in section 3 uses the sum in (24), so it avoids the tagged-observable problem. The required spatial estimate cannot avoid it in the same way.

[Imported, heuristic only] C89's continuum argument gives cancellation of the static smooth-density tail and then predicts faster descendant and harmonic tails. C94 explicitly labels this as power counting rather than an all-distance theorem. It is not promoted to a proof here.

## 6. Frozen subtraction does not imply a uniform Taylor remainder

[Proved, exact diagnostic model] Consider an independent two-level bath and a two-level material variable,

    H_0 = t tau_z,
    A = (1+tau_x)/2,
    h = diag(0,1),
    H(J) = H_0 + J A h.

This bath is gapped, A is a projector, and ||h||=1. In the uncoupled thermal state a_0=1/2 at every T. Define the same frozen residual as in (1),

    R(J,T) = F(H(J))-F(H_0)-F(Jh/2).

The h=0 block has energies +/-t. The h=1 block has energies

    J/2 +/- sqrt(t^2+J^2/4).

Consequently, with beta=1/T and d(J)=sqrt(t^2+J^2/4),

    R(J,T) = -T log {
       [1+exp(-beta J/2) cosh(beta d(J))/cosh(beta t)]
                     /[1+exp(-beta J/2)] }.              (29)

This is an exact full trace, with no vacancy-number crossing.

Expand at J=0 for fixed positive T. First,

    d(J)=t+J^2/(8t)+O(J^4),
    cosh(beta d(J))/cosh(beta t)
        =1+ beta tanh(beta t) J^2/(8t)+O_T(J^4),
    [1+exp(beta J/2)]^(-1)=1/2-beta J/8+O_T(J^3).

Substitution into (29) gives

    R(J,T) = -tanh(t/T) J^2/(16t)
               +tanh(t/T) J^3/(64tT)+O_T(J^4).          (30)

The cubic Taylor coefficient diverges as 1/T. Taking m independent copies multiplies R and every coefficient by m. This obstruction survives the thermodynamic limit of those copies.

For every fixed J>0, the h=0 block is the unique material ground sector, since

    J/2-sqrt(t^2+J^2/4) > -t.

Both the mobile and frozen traces therefore select h=0 as T->0, and

    lim_(T->0) R(J,T)=0.                                 (31)

In contrast, the J=0 thermal Taylor term tends to -J^2/(16t). If R2 denotes that quadratic term, then

    lim_(T->0) (R-R2)/|R2|=1                             (32)

for every fixed positive J. In particular there is no function eta(J/t)->0 such that |R-R2|<=eta(J/t)|R2| for every sufficiently small positive T. There is also no bound |R-R2|<=C J^3/t^2 with a temperature-independent C.

[Proved, compatible uniform integrated bound] The diagnostic does not violate an estimate of the form (1). For d>=t,

    0 <= T log[cosh(d/T)/cosh(t/T)] <= d-t,
    d(J)-t = (J^2/4)/(sqrt(t^2+J^2/4)+t) <= J^2/(8t).

The ratio in (29) is a normalized positive sum of 1 and the cosh ratio. It follows that

    -J^2/(8t) <= R(J,T) <= 0                             (33)

uniformly in T. For m copies the bound is -m J^2/(8t)<=R_m<=0. Thus a small integrated correction can coexist with a nonuniform Taylor series.

[Scope] This is not a counterexample to (1) for C88. It disproves a generic inference that a free-bath correlation bound and frozen subtraction automatically give a temperature-uniform Taylor remainder. The material thermal weights reorganize on the scale T~J. Removing finite-volume vacancy crossings does not remove that effect. A successful expansion for C88 may resum the frozen material dynamics and control a bath-induced interaction directly. Such an argument has not been supplied here.

The same diagnostic applies to a relative residual: choose a reference with h=0 on the same material space. Its residual vanishes identically. Replication then gives an extensive relative example, although it is still not the C88 chain.

## 7. Exact resummed formulation and what remains to be bounded

[Proved, exact identity] Set

    Href_sigma = H_0+sV + J a_0 Hmat_sigma,
    Y_sigma = sum_i (A_i-a_0) h_sigma,i.

The two factors in Href commute. Its partition function is Zgap Zmat_sigma(Ja_0). The interaction-picture identity gives

    exp(-R_sigma/T)
       = < T_tau exp[-J integral_0^(1/T) Y_sigma(tau) d tau]
                                                       >_ref. (34)

No scalar replacement of h_sigma,i is made. Its time dependence in (34) is generated by J a_0 Hmat_sigma. The one-point expectation of Y vanishes, because the free bath contact is a_0.

Equation (34) incorporates the material thermal weights at the actual value J/T. A connected expansion of its logarithm, if convergent, involves all multi-time bath correlations and all ordered material traces. A bound on the static two-point function proves neither their convergence nor a comparison proportional to W.

[Conditional, sufficient bound] One sufficient conclusion from a resummed argument would be

    |[R_sigma-R_good] - [R_sigma^(2,res)-R_good^(2,res)]|
       <= (J^2/t) eta(J/t) W + B_N,
    |R_sigma^(2,res)-R_good^(2,res)| <= C (J^2/t) W,
    eta(x) -> 0,
    B_N/N -> 0.                                         (35)

Here the retained quadratic expression is evaluated with the frozen material dynamics at Ja_0, not with normalized material trace at J=0. Equations (35) would imply (1) with epsilon(x)=x[C+eta(x)]. A weaker direct bound on the full logarithmic difference in (34) would also suffice. There is no requirement that its Taylor series at J=0 converge uniformly in T.

[Open, second missing estimate] No estimate of (34) or (35) uniform in wrong-cell length, material thermal state, nonzero port, N, and low T has been proved. Even a total-volume bound C J^2 N/t on each residual would not suffice. It would have to cancel on arbitrarily long selected runs, including the different marker offsets produced by wrong cells, leaving W and one global subextensive error. C94 performs that matching only for its product-state quadratic Gram matrix.

## 8. Consequences for d=13 and consistency checks

[Imported, established part of the application] At infinite grammar penalty, C92 proves a positive contact lower bound, activated removal of selected length-15 pins including the full legal auxiliary trace and nonzero port, and exclusion of entire marker-free and one-marker sectors at low temperature. Those results remain available. The finite-K relative correction is not required by the present task.

[Conditional] If (1) were established in the full legal trace, and its epsilon fit below the frozen selector margin, the C92 length-weighted activity argument would give

    0 <= p_full(T)-p_good(T) <= C T exp(-delta/T)=o(T^2)

in the prescribed thermodynamic order of limits and source window. Neither (18) nor (21) supplies (1), so this remains conditional.

[Proved consistency checks] The contact formula satisfies a->1 at full material occupancy. The interior estimates deliberately exclude that endpoint, where e'' vanishes. At half filling (5) reproduces s=2t/pi and the inherited contact 1/2-2/pi^2. All density-response coefficients have units 1/t. The finite-volume contact jump is O(1/N), and multiplication by an extensive material energy gives only O(J), as required for a global finite-size allowance. The exact two-level bound (33) is quadratic in J despite the divergent cubic Taylor coefficient, so there is no contradiction between that diagnostic and the desired weaker comparison.

[Numerical] Running python3 -B checks.py in this folder checks (17) by small thermal quadratures, (23) directly from the exact finite formulas, the form-factor calculation (24b) against three C94 gap-resolvent row sums, and the two-level trace (29). All checks pass. CHECKS.md records the values. No large computation, cluster job, network command, or interacting d=13 simulation was performed.

[Final status] The thermodynamic vacancy-number obstruction can be isolated and bounded. The free tagged-contact spatial estimate and the full quantum relative remainder are separate unresolved steps. The requested full mobile defect comparison cannot be concluded from the results established here.
