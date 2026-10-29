# C65. A sharp-window obstruction and a trace-compatible smooth regulator

## 1. Verdict and conventions

[proved, negative result for a specified construction] The finite fermion windows of C62 cannot be used, at every finite cutoff, in the absolute spatial moment norm of C59 TP1 with constants uniform in spatial volume. A second-order connected four-field kernel has a nonzero static momentum cusp at every positive temperature. Its first spatial moment already fails, hence so does the fifth-moment norm. This is a failure of an actual coefficient, not just of a propagator majorant.

[proved, repair of this obstruction] Sections 5 and 6 give a smooth fermion regulator of the same C62 Hamiltonian, prove convergence of its full Gibbs trace, and remove this cusp and the associated single-scale moment obstruction. The interaction remains quartic. This does not require identifying the scalar of BM's different Euclidean model with a trace.

[assumed, remaining theorem] The uniform all-depth joint comparison for the smooth construction is not proved here. In particular, this note does not prove TP or the fixed-interaction XXZ thermal formula. The negative result concerns uniform control before removal of the sharp regulator. It does not rule out first removing that regulator at fixed volume and temperature, as C62 does, and constructing an expansion afterwards.

[proved work convention] All momenta below are measured from a Fermi point. Write A>0 for C62's free chiral speed, g=2 pi B for its Hamiltonian interaction coefficient, and

    H = H_0(A) + V,
    V = (g/L) sum_p chi(p) rho_+(p)rho_-(-p).             (1)

Here chi is real, even and C-infinity, is one near zero, and vanishes for |p|>=kappa. Both densities are sea normal ordered. The parameters A,B and the dressed speed v=sqrt(A^2-B^2) are those already matched in C62. They are not identified with BM's bare Euclidean parameters. We do not recompute their conversion or the oscillator pressure.

[assumed, imported representation] We use C62 Section 3's filled-sea current representation and its form inequality. On a spatial circle, with antiperiodic fermions,

    |<psi,V psi>| <= b <psi,H_0 psi> + C_L ||psi||^2,
    b=|B|/A<1,
    C_L=sum_(m>0) (2 pi m/L)|B chi(2 pi m/L)|.           (2)

The finite-excitation space is a form core. The free Gibbs trace is finite for each fixed L,beta>0. These statements are proved or explicitly imported in C62 and are retained here. Constants C_L in the trace arguments need not be uniform in L. They are not RG norm estimates.

## 2. The finite-temperature sharp-window bubble

### 2.1 Exact coefficient

[proved] Take the thermodynamic version of a symmetric finite window |k|<Lambda. It is obtained from C62's finite windows along sequences with 2 pi N/L tending to Lambda. Write

    n(k) = 1/(1+exp(A k/T)),
    G_Lambda(omega,k)=1_(|k|<Lambda)/(-i omega+A k),
    omega=(2n+1)pi T.

Define the positive static density response by minus the fermion loop,

    Pi_Lambda^T(p)
       = -T sum_omega integral dk/(2 pi)
                    G_Lambda(omega,k)G_Lambda(omega,k+p).
                                                               (3)

For real p!=0, a partial fraction in the two denominators, followed by the free occupation sum, gives

    Pi_Lambda^T(p)
       = integral_(|k|<Lambda, |k+p|<Lambda) dk/(2 pi)
                           [n(k)-n(k+p)]/(A p).          (4)

One can also derive (4) without a conditionally convergent single-propagator sum. In a finite free Fock space integrate the connected imaginary-time density correlation. For two distinct energies E,E', its value is [n_F(E)-n_F(E')]/(E'-E). Take the spatial Riemann-sum limit. The product sum in (3) is absolutely convergent at fixed Lambda,T.

For 0<p<2Lambda, cancellation of the overlapping integration intervals gives

    integral_(-Lambda)^(Lambda-p) [n(k)-n(k+p)] dk
       = integral_(-Lambda)^(-Lambda+p) n(k)dk
                              -integral_(Lambda-p)^Lambda n(k)dk
       = p-2 integral_(Lambda-p)^Lambda n(k)dk.          (5)

The last step uses n(-k)=1-n(k). Changing k to k+p in (4) shows Pi(-p)=Pi(p). Thus for p near zero,

    Pi_Lambda^T(p)
       = [1-2n(Lambda)]/(2 pi A) + c_(Lambda,T)|p| + O(p^2),
    c_(Lambda,T)=n'(Lambda)/(2 pi A)
       = -n(Lambda)[1-n(Lambda)]/(2 pi T) < 0.           (6)

In particular, the right derivative minus the left derivative is 2c_(Lambda,T), not zero. For fixed Lambda, c_(Lambda,T) is of order -exp(-A Lambda/T)/(2 pi T) as T tends to zero. Its units are inverse energy times inverse momentum, as required for the coefficient of |p| in a response of units inverse speed.

At T=0 and |p|<Lambda, the numerator in (5) is p, so

    Pi_Lambda^0(p)=1/(2 pi A).                           (7)

Consequently Pi_Lambda^T-Pi_Lambda^0 has the same cusp (6). Smallness of its coefficient does not restore differentiability.

### 2.2 The cusp survives ultraviolet preparation

[proved] Choose C_<(omega,k)=chi_<(omega,k)/(-i omega+A k), where chi_< is C-infinity, equals one in a neighborhood of the origin, and has support inside |k|<Lambda/3 and a fixed bounded frequency region. Subtract it from the full window covariance,

    C_> = G_Lambda-C_<.                                 (8)

The connected loop of C_> is the loop of G_Lambda minus its two cross loops with C_<, plus the C_< loop. At fixed T the cross loops involve finitely many fermionic frequencies because C_< has compact frequency support. At sufficiently small |p|, the other line's sharp spatial cutoff is identically one on those supports. Every cross-loop integrand is then smooth in p on a fixed compact k domain, with denominators bounded away from zero by pi T. Differentiation under its integral is legitimate to every fixed order. The same holds for the C_< loop.

Thus the static loop Pi_>^T differs from Pi_Lambda^T by a C-infinity function near p=0. Its derivative jump is still 2c_(Lambda,T).

For completeness, the zero-temperature high loop is smooth near zero. Compare its sharp spatial cutoff with a smooth compact cutoff equal to one on a larger neighborhood of the low support, for example |k|<=Lambda/2. Their cross terms with C_< agree at sufficiently small |p|. Their full-loop difference is supported where k and k+p have the same sign. The convergent frequency integral of G(omega,k)G(omega,k+p) is zero there at T=0, by the divided occupation formula in (4), or by a residue calculation. Thus their high static loops coincide. The smooth high covariance has its infrared singularity removed, smooth compact spatial support and integrable large-frequency product tails, also after any fixed spatial derivative. Its loop is C-infinity near p=0 by differentiation under the integral. This proves the assertion for the sharp high loop at T=0 as well. Therefore the thermal replacement Pi_>^T-Pi_>^0 retains (6).

This argument also applies to the sum of any exact smooth shell decomposition of (8). At least its summed coefficient fails the required norm. Splitting the same coefficient into more scales cannot give a uniformly summable absolute moment majorant for it.

## 3. Failure of the effective-kernel norm, not only of a line bound

### 3.1 An actual generated quartic

[proved coefficient statement] Integrate the high covariance (8), keeping low external fields. In the connected expansion at order g^2, retain the four minus-chirality fields at the two interaction vertices and contract their four plus-chirality fields. The only connected contraction is the plus density loop. If

    R_-(t,x) = integral dy chi_check(x-y) rho_-^<(t,y),

the resulting nonconstant term is, up to the fixed overall sign convention for the effective action,

    (g^2/2) integral dz dz' Pi_(+,>)(z-z') R_-(z)R_-(z'). (9)

Its static transfer multiplier is g^2 chi(p)^2 Pi_(+,>)^T(p), apart from that constant prefactor. Since chi=1 near zero, its cusp is unchanged. Tadpole products are removed by connectedness. There is no other connected order-g^2 diagram with four external minus fields from two g2 vertices that could cancel this loop.

Fermionic antisymmetrization does not cancel the cusp. With all four external Matsubara frequencies equal, choose spatial momenta

    k1=a, k2=a-p, k3=b-p, k4=b,

inside a fixed low-momentum patch, with a-b a fixed small nonzero number. They obey k1-k2+k3-k4=0. The antisymmetric density-density vertex is proportional to

    F(p)-F(a-b),    F(p)=chi(p)^2 Pi_(+,>)^T(p),          (10)

with an inessential nonzero convention factor. The exchange term is smooth as p passes through zero and the first term has (6). Equivalently, in position space the two density pairs have distinct points; their quartic monomial does not vanish. Only the fully coincident local same-chirality quartic vanishes by the Pauli principle.

### 3.2 Exact failed estimate

[proved] Let a pinned kernel, including its possible delta measures, have finite first spatial absolute moment. Its Fourier transform is continuously differentiable in each independent spatial momentum. Indeed multiplication by the corresponding coordinate difference is integrable, so differentiation and continuity follow by dominated convergence. For a time-periodic kernel the same argument applies to its zero-frequency coefficient. In particular,

    integral_0^beta dt integral_R dx (1+q|x|)^5 |K(t,x)| < infinity
                                                               (11)

implies continuous first spatial differentiability of K_hat(0,p), for every fixed q>0. The multi-field pinned norm of TP1 implies the corresponding assertion in each independent external spatial momentum. Polynomial factors from a finite set of derivative-field labels preserve this differentiability.

Equations (6), (9) and (10) contradict this necessary consequence of (11). Therefore the order-g^2 prepared four-field kernel, and its thermal replacement, are not in the common fifth-moment Banach space at any finite Lambda and T>0, after the spatial thermodynamic limit. Subtracting local constants or affine quadratic coefficients cannot remove a |p| cusp. Nor can a smooth change of field amplitude.

The estimate that fails is already the ultraviolet coefficient estimate

    || [g^2] W_(4,>,Lambda)^T ||_(q,5) <= C,

and also its purported thermal version

    || [g^2](W_(4,>,Lambda)^T-Per_beta W_(4,>,Lambda)^0) ||_(q,5)
                                      <= C (T/q)^p.     (12)

Here [g^2] means the Taylor coefficient, Per_beta means periodization with the prescribed fermionic signs, and ||.||_(q,5) is the pinned coefficient norm with fifth spacetime moment at a fixed preparation scale q. Both left sides are infinite in the thermodynamic limit. Thus (12) fails for every p, including p>2+theta. A cutoff-independent analytic connected expansion in that norm would imply a finite coefficient bound by Cauchy's estimate. It is impossible for this regulator and this norm.

This also rules out a uniform finite-volume version. If the finite-circle fifth-moment bounds were uniform in L with Lambda fixed, their distributional infinite-volume limit would have a finite first weighted variation, by tightness and lower semicontinuity of weighted variation. The explicit Riemann-sum limit (4) does not. Thus some sequence of finite-volume norms diverges. The uniform-regulator requirement allows sequences with N proportional to L and is contradicted by this sequence.

[proved scope] This is not a counterexample to TP after regulator removal, to a weaker momentum norm, or to the thermal coefficient. For fixed p, the edge term vanishes as Lambda tends to infinity. C62 removes its window first at fixed L,T. That order of trace limits is valid and is not contradicted. What fails is using these finite windows as a uniformly bounded family in the absolute moment space before taking those limits. Exponentially small times an infinite moment remains an infinite moment.

### 3.3 A simpler propagator warning

[proved] A sharp endpoint also spoils the usual spatial line majorant. For a smooth symbol a(k) with a(Lambda) nonzero, integrating

    integral_(-Lambda)^Lambda e^(ikx)a(k) dk

by parts gives a boundary term [e^(iLambda x)a(Lambda)-e^(-iLambda x)a(-Lambda)]/(ix), followed by an O(x^-2) remainder if a has two integrable derivatives. This prevents arbitrary polynomial spatial decay. It motivates concern about the regulator but by itself would only refute that line estimate. Sections 2 and 3 instead show failure in the actual generated kernel, where cancellations have already been taken.

## 4. What BM's smoothing does and does not resolve

[assumed, imported primary definitions] BM 0907.2837, equations (29)-(33), uses a smooth frequency-momentum fermion cutoff and a spacetime interaction v0. Appendix D distinguishes its parameter convention from the Hamiltonian Luttinger one. Its Section 3 compares infrared running data and correlations. These definitions and the correlation-comparison scope were checked in the primary source. They are not a finite-temperature scalar theorem invoked here.

[proved] At a fixed smooth fermion regulator and T>0, the analogue of (3), with smooth compact spatial weights on each fermion line, is C-infinity in p. Every spatial derivative has a fixed compact integration region; the Matsubara product and its spatial derivatives have a summable large-frequency tail. Thus BM's smooth fermion cutoff removes the moving sharp-overlap boundary responsible for (6).

Smoothing only the interaction would not do this. In (9) a spacetime-smoothed interaction inserts v0_hat(0,p)^2 at the external transfer. If that multiplier is smooth and nonzero at p=0, the derivative jump is multiplied by v0_hat(0,0)^2 and persists. It is the smooth fermion regulator, not the word "retarded", that repairs this particular estimate.

[assumed, unresolved scalar identification] BM's spacetime interaction is not the instantaneous Hamiltonian (1). No identity between its full cylinder functional, with all vacuum constants, and the Gibbs trace of (1) follows from the stated correlation comparison. A positive Hamiltonian or transfer-operator realization of that precise Euclidean regulator has not been constructed here. An independent o(T^2) comparison of its ground-subtracted scalar would suffice, but has not been proved here either. No claim that every spacetime-smoothed reference lacks a Hamiltonian realization is made.

The smooth Hamiltonian construction below avoids this extra scalar-identification obligation while retaining the spatial smoothing that (12) requires.

## 5. A smooth quartic regulator with the same trace limit

### 5.1 Sea-preserving attenuation

[proved construction] Work at fixed L. Realize the filled sea as the vacuum of excitation operators d_j. Above the sea d is the original fermion annihilator; below it d is a hole annihilator. Including both chiralities,

    H_0 = sum_j E_j d_j^dagger d_j,    E_j=A|k_j|>0.     (13)

Choose an even C-infinity function r with 0<=r<=1, r(s)=1 for |s|<=1 and r(s)=0 for |s|>=2. Put r_Lambda(k)=r(k/Lambda).

Introduce one auxiliary copy of the excitation one-particle space. The map

    j_Lambda f = r_Lambda f direct_sum sqrt(1-r_Lambda^2) f

is an isometry. Its fermionic second quantization J_Lambda=Gamma(j_Lambda) is also an isometry, from the original excitation Fock space to the doubled one. Define, first on finite-mode polynomials or forms,

    Phi_Lambda(O)=J_Lambda^* (O tensor I_aux) J_Lambda.   (14)

This map is unital and preserves positivity. It preserves quadratic-form inequalities on the common energy form domain. To see its action on polynomials, extend j_Lambda to a two-copy one-particle unitary. A physical annihilator becomes r_j d_j+sqrt(1-r_j^2) a_j. Take the auxiliary vacuum matrix element. For an excitation-normal-ordered monomial, every auxiliary term vanishes unless no auxiliary operator occurs. Hence

    Phi_Lambda(:d_(j1)^# ... d_(jl)^#:)
       = product_(i=1)^l r_Lambda(k_ji)
                                     :d_(j1)^# ... d_(jl)^#: . (15)

For the original sea-normal-ordered density this gives

    rho_(omega,Lambda)(p)
       = sum_k r_Lambda(k)r_Lambda(k+p)
                              :c_(omega,k)^dagger c_(omega,k+p): . (16)

The two chiral densities have no cross contractions and are even. Their auxiliary vacuum expectation factorizes, so

    Phi_Lambda(rho_+(p)rho_-(-p))
                         =rho_(+,Lambda)(p)rho_(-,Lambda)(-p). (17)

In particular define

    V_Lambda=Phi_Lambda(V),
    H_Lambda^sm=H_0+V_Lambda.                            (18)

This is a quartic interaction in finitely many fermion modes, with the exact normal-order constants and quadratic terms retained. All other modes remain free. It is not multiplication of the whole many-body interaction by an occupation projector, which would generate arbitrary degrees.

### 5.2 Stability uniform in the smooth regulator

[proved] Apply (14) to the two form inequalities for V in (2). Directly from (13)-(15),

    Phi_Lambda(H_0)=sum_j E_j r_Lambda(k_j)^2 d_j^dagger d_j <= H_0.

Thus

    -b H_0-C_L <= V_Lambda <= b H_0+C_L,
    (1-b)H_0-C_L <= H_Lambda^sm <= (1+b)H_0+C_L.         (19)

The constants are the same as (2), independently of Lambda. For precision on unbounded forms, J_Lambda maps D(H_0^(1/2)) into the first-copy energy form domain because the pulled-back first-copy energy is Phi_Lambda(H_0). Apply (2) there and pull back. This proves (19) without extending a bounded-operator identity beyond its domain.

The relative bound b<1 defines closed forms on the common free form domain. At fixed Lambda the interaction is a bounded finite-mode polynomial. The lower bound and the compact embedding of the free form domain imply compact resolvent for each H_Lambda^sm. Min-max gives

    E_j(H_Lambda^sm) >= (1-b)E_j(H_0)-C_L.              (20)

Since the right side has a finite Gibbs sum, all the traces are finite and uniformly dominated at each fixed L,beta.

### 5.3 Trace convergence, including the scalar

[proved] For every finite-excitation vector phi, V_Lambda phi=V phi once Lambda is sufficiently large. Indeed a density transfer in (1) has |p|<kappa. Acting on a finite excitation, a sea-normal-ordered density either moves an existing excitation by at most kappa or creates a particle-hole pair across the sea within kappa of the Fermi point. A product of one density of each chirality still has only finitely many possible excitation modes at fixed L. If the original vector uses |k|<=M, the relevant modes are within M+2kappa. All their r factors are one for sufficiently large Lambda. This also holds for all matrix elements in a fixed finite-excitation subspace.

Let E_j and E_(j,Lambda) denote ordered eigenvalues of H and H_Lambda^sm. We prove

    E_(j,Lambda) -> E_j for every fixed j.              (21)

For the upper bound, approximate the first j+1 eigenvectors of H in its form norm by vectors in the finite-excitation core, using (2). On their common finite-dimensional trial space the forms of H_Lambda^sm and H agree for large Lambda. Min-max and then arbitrarily accurate core approximation give limsup E_(j,Lambda)<=E_j.

For the lower bound, take a subsequence along which E_(j,Lambda) approaches its liminf, and choose its first j+1 orthonormal eigenvectors. Equations (19) and the upper bound give uniform free-form norms. The free form embedding is compact, so a further subsequence of each of these finitely many vectors converges strongly in Hilbert norm and weakly in free form norm. Their limits remain orthonormal. Their eigenvalues have convergent subsequences as well.

For a finite-excitation test vector phi, H_Lambda^sm phi=H phi eventually. The eigenvector equation therefore passes to the limit,

    <H phi,psi_i>=e_i <phi,psi_i>.

The limiting psi_i is in the common form domain by weak compactness. By the boundedness of the H form in the free form norm and the form-core property, this equality extends to all form-domain test vectors. The representation theorem for the closed form then gives H psi_i=e_i psi_i. There are j+1 orthonormal eigenvectors with eigenvalues no greater than liminf E_(j,Lambda). Min-max yields E_j<=liminf E_(j,Lambda), proving (21).

Now use (20). The summable sequence

    exp(beta C_L) exp[-beta(1-b)E_j(H_0)]

dominates exp[-beta E_(j,Lambda)] independently of Lambda. Dominated convergence over j gives

    Tr exp(-beta H_Lambda^sm) -> Tr exp(-beta H),
    min spec H_Lambda^sm -> min spec H.                (22)

Monotonicity in Lambda is neither needed nor asserted.

At fixed Lambda the active-mode Hamiltonian is finite dimensional and the spectator modes are free. Its Duhamel formula is an exact Grassmann Wick expansion after multiplication by the full free trace Z_0. Each vertex leg has an r_Lambda factor, so equivalently contractions of smeared fields use

    G_Lambda^sm(omega,k)=r_Lambda(k)^2/(-i omega+A k).    (23)

Normal ordering uses the smeared sea contraction in (16), not a new temperature-dependent counterterm. The finite-dimensional identity fixes all scalar factors. Equation (22) therefore identifies the total normalized scalar after regulator removal with C62's Hamiltonian trace. Taking L to infinity next and then T to zero imports C62/C59's result

    f_ref(T)-e_ref=-pi T^2/(6v)+o(T^2).                (24)

This proves trace compatibility of a smooth quartic fermion regulator. It does not establish a uniform radius or coefficient norm for its connected logarithm.

## 6. Which analytic bounds the smooth regulator now satisfies

### 6.1 The bubble is smooth

[proved] For (23), write w_Lambda(k)=r_Lambda(k)^2. The static response is

    Pi_Lambda^sm,T(p)
       = integral dk/(2 pi) w_Lambda(k)w_Lambda(k+p)
                               [n(k)-n(k+p)]/(A p).    (25)

The divided difference extends smoothly to p=0, since

    [n(k)-n(k+p)]/p = -integral_0^1 n'(k+s p) ds.

The other factors are smooth and have compact support. Differentiating under the integral proves C-infinity regularity. It is even in p by a change of variable, hence its derivative at zero is zero. There is no analogue of (6). Smooth infrared subtraction retains this property.

The assertion is finite-regulator smoothness, not a bound on every derivative uniform as T tends to zero for the unsplit massless covariance. Such bounds must be stated shell by shell.

[proved stronger static check] The specific bad coefficient of Section 3 has a uniform repair. Fix a transfer support |p|<=p0 and Lambda>=Lambda0>2p0. The cutoff-free static loop is 1/(2 pi A) at every T: integrate the divided difference over all k, or integrate -n' and use its endpoint values. Thus

    Pi_Lambda^sm,T(p)-1/(2 pi A)
      = -1/(2 pi A) integral dk integral_0^1 ds
          [w_Lambda(k)w_Lambda(k+p)-1] n'(k+s p).     (25a)

On the support of the bracket, |k+s p|>=Lambda-2p0. For each fixed m, differentiating under this integral gives a bound

    sup_(|p|<=p0) |partial_p^m (Pi_Lambda^sm,T-1/(2 pi A))|
        <= C_m(1+T^(-m)) exp[-A(Lambda-2p0)/T],       (25b)

for 0<T<=T0, with constants independent of Lambda. To verify it, distribute the m derivatives between w and n'. A derivative of w costs at most C_j Lambda^(-j). For |k+s p| away from zero, a derivative of order j of n' is bounded by C_j T^(-j-1) exp[-A|k+s p|/T]. Its k integral over the two tails costs at most C T, leaving T^(-j). Sum the finitely many allocations and use Lambda>=Lambda0. This proves (25b).

At T=0 the divided difference is supported between 0 and -p, where both w factors are one, so the smooth static loop equals 1/(2 pi A) on this transfer support. Multiplying (25a) by the compact smooth chi(p)^2, with support inside |p|<=p0, and integrating by parts m+2 times in p gives a finite spatial absolute moment of order m for its inverse Fourier transform. Equation (25b) bounds that moment by C_(m,P,q) T^P for every fixed P, at a fixed preparation scale q and uniformly for Lambda>=Lambda0. This follows because every fixed negative power of T times exp[-A(Lambda0-2p0)/T] is bounded by C_P T^P.

The cross terms with a fixed low covariance are independent of Lambda once their momentum supports lie in r_Lambda=1. Hence this bound also controls the regulator-dependent part of the prepared static loop. It is a check of the troublesome order-g^2 coefficient. It is not a full-frequency, all-order preparation theorem.

### 6.2 Common single-scale bounds and thermal power

[proved] Let f_q(omega,k) be a fixed smooth rescaled annular cutoff, with |omega|+A|k| comparable to q. For 0<q<=c Lambda, define

    g_(q,Lambda)(omega,k)
         = f_q(omega,k) r_Lambda(k)^2/(-i omega+A k).    (26)

Constants can depend on c,A and the cutoff functions. With omega=q u, k=q z/A, every fixed derivative of the rescaled numerator is bounded uniformly in Lambda and q, because each derivative of r contributes q/(A Lambda)<=c/A. The rescaled denominator is bounded away from zero on the annulus. Repeated Fourier integration by parts gives, for every fixed integer M,

    |g_(q,Lambda)(t,x)| <= C_M q(1+q|t|+q|x|/A)^(-M),
    integral dt dx (1+q|t|+q|x|/A)^m |g_(q,Lambda)(t,x)|
                                                    <= C_m/q. (27)

To prove the second bound use the first with M>m+2 and rescale the two coordinates. The factors of A are absorbed in the constants. The absolute Fourier integral is at most Cq, so the square-root Fourier factorization supplies a Gram representation with product of vector norms at most Cq. Sampling allowed fermionic frequencies gives the same bound for T<=c' q: the number of sampled frequencies in the annulus is O(1+q/T), and their measure is T.

The antiperiodic periodization of (27) is the thermal shell. For |t|<=beta/2, the nonzero images have |t+m beta|>=(|m|-1/2)beta. Using (27) with enough extra powers gives any prescribed fixed power of T/q, with any fixed spatial moment remaining. Thus the covariance hypotheses of the collaborator's WINDING_ESTIMATE hold uniformly in Lambda for these shells. Its single-integration theorem, including its parameter derivatives and vacuum root, applies once the endpoint norms it requires are supplied.

[assumed, imported one-step result] The four-node THERMAL_LOCALIZATION projection has p=3 with the stated five moments. Its proof is not repeated. Hence the sharp-window obstruction to using p=3 is absent in this construction, and p>2+theta permits any 0<theta<1 at this single-step level. Neither an irrelevant defect gain nor a complete RG trajectory follows merely from (27).

[proved limit of the bound] Equation (27) covers shells with q comparable to or below the spatial fermion regulator. The Hamiltonian covariance has arbitrarily large frequencies even at finite Lambda. Its entire ultraviolet integration, including the frequency tail, must be prepared with an appropriate determinant construction. We have not replaced that task by a false two-dimensional Gram bound for an uncut covariance. Uniform bounds on the elementary infrared shells also do not bound the effective kernels produced after all higher shells have been integrated.

## 7. The remaining joint estimate after this repair

[assumed, exact unproved requirement] Use (18) for the trace-defined reference and BM's smooth regulator for the comparison reference. Remove regulators only after an appropriate common construction, or prove uniform preparation before removal. Beyond the single-scale statements above, the needed estimate is an absolutely summable normalized renormalized-tree expansion in which one reference difference and one thermal replacement are retained simultaneously. At root h, q_h=2^(-h), its full kernel norm must obey

    ||X_h^0-Y_h^0|| <= C q_h^theta,
    ||(X_h^T-Y_h^T)-Per_beta(X_h^0-Y_h^0)||
                                   <= C q_h^theta(T/q_h)^p, (28)

with compatible local covariance and field-normalization coordinates. The notation Per_beta in (28) includes the common finite-temperature localization identification; it is not a claim that actual thermal kernels equal periodized zero-temperature kernels. The pure reference shift is bounded by C(T/q_h)^p. The vacuum maps must retain their Gaussian factors and satisfy

    ||D E_(h,T)||+||D^2 E_(h,T)|| <= C q_h^2,
    ||D(E_(h,T)-E_(h,0))|| <= C q_h^2(T/q_h)^p.        (29)

A sufficient tree majorant, in the notation of C62 Section 5, retains

    C^n |g|^n 2^(-sigma max(h-k,h-l)) q_k^alpha (T/q_l)^p
                                                               (30)

after ordinary renormalized dimensions and endpoint weights are summed. Here k and l are defect and thermal attachment scales. A single fixed reserve of field weight and coordinate moments must pay all ordinary Taylor distances and the fixed marks, uniformly in tree depth. Equations (28)-(30) are assumed targets, not consequences of the trace argument.

[proved implication using prior work] C62's double scale sum handles (30) for theta<min(alpha,sigma). The collaborator's MOMENT_RESERVE Section 6 pays a fixed thermal mark once ordinary distance routing is controlled. C57/C59 then propagate the four histories and sum the scalars for p>2+theta. Those arguments are retained without rederivation. The present obstruction occurs before they can be applied to the sharp regulator. The smooth regulator removes that obstruction but does not supply the ordinary comparison, its critical parameter domains, or the normalized differentiated tree sum.

[assumed, ultraviolet subproblem for the smooth construction] One concrete earlier obligation is a common fixed-scale ultraviolet effective action U_Lambda^T, analytic in a cutoff-independent neighborhood of the matched parameters, with a bounded stronger field/moment norm and thermal derivatives of power p. It must include the exact scalar. For the sharp windows its four-field coefficient fails (12). For (18) it is no longer refuted by that coefficient; its all-order bound and its comparison to BM remain unproved here. C62's affine-critical-interpolation example still prevents inferring global derivative domains from two controlled endpoints.

[proved application boundary] The C48 massive matching and C58 ensemble transfer are separate conditional inputs for the encoded chain. This note neither repeats nor strengthens them. It verifies a regulator obstruction and a regulator repair for the single-channel reference. It verifies no new microscopic hypothesis of the multichannel or retarded-action application.

## 8. Checks and provenance

[numerical] checks.py compares (4) with the boundary identity (5), checks the nonzero one-sided slope (6), verifies a growing weighted spatial integral of an explicitly isolated cusp, and checks that the smooth version (25) has a vanishing first derivative at zero. It also checks sea-preserving attenuation on a small finite CAR algebra and the inherited two-sided operator bound. These are finite controls of the new lemmas, not a numerical test of TP.

[proved scope] No large computation, lattice RG simulation, subagents, shell network commands, messages to authors, or edits outside work_C65 were used. Primary sources were accessed with the browser. The trace and cutoff arguments have not received independent or human review. RESULT.md records the verdict and STATE.md is the local handoff because the project-level memory is outside the permitted write scope.
