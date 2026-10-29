# C78. Mixed density rings and the remaining fermionic comparison

## 1. Verdict and conventions

[assumed, unresolved target] Lemma TP has not been proved for interacting XXZ. No mixed-chirality coefficient violating TP has been established. The results below integrate the alternating two-density-loop subseries of both high chiralities, including its vacuum scalar, on a uniform analytic disk. They do not integrate the complete high action. Section 6 identifies the first additional mixed vacuum coefficients. Section 8 states the remaining infrared obligation without replacing it by a Gaussian density calculation.

[proved partial result] A single spatial smearing of a high chiral bubble belongs uniformly to weighted L1 and L2. Its thermal replacement has the same bounds with any fixed power of T and two speed derivatives. The L2 bound controls the shared time variable in the mixed order-g^2 vacuum coefficient. Together with convolution estimates it gives an analytic sum of all alternating bubble rings. These statements concern infinite spatial volume and a fixed preparation scale, not a new uniform finite-circle theorem.

[assumed, imported inputs] Use C62's trace-defined Hamiltonian and its matched parameters, and C65's smooth sea-preserving regulator. Write

    V = g integral dt dx dy v0(x-y) rho_+(t,x) rho_-(t,y),
    g=2 pi B,       A>|B|,       v=sqrt(A^2-B^2),
    Fourier(v0)=chi,     supp chi subset [-kappa,kappa].       (1)

Here chi is real, even, smooth, and one near zero. The two densities are sea normal ordered. C62 already supplies the conversion from the BM parameters to A,B and the complete trace evaluation. Those results are not rederived or identified with BM's differently regulated Euclidean functional.

[assumed, imported covariance estimates] Use C71's inversion-even compact low multiplier f, whose spatial support is within |k|<=K, and its splitting

    c_(+,Lambda)(omega,k)=w_Lambda(k)/(Ak-i omega),
    u_+(omega,k)=f(omega,k)/(Ak-i omega),
    a_(+,Lambda)=c_(+,Lambda)-u_+,
    a_-(t,x)=a_+(t,-x).                                  (2)

The minus low multiplier is reflected with the covariance. Take Lambda>=Lambda0>K+kappa. All speed derivatives below are through order two, with A in a compact positive interval and the cutoff profiles fixed. C71 proves the common bound

    |partial_A^j a_(omega,Lambda)(z)|
        <= C_N |z|^(-1)(1+|z|)^(-N),       j<=2,           (3)

for every fixed N. The same bound holds after cutoff removal. C71 also fixes the distributional bubble using the Hamiltonian prescription. No equal-time contact term is chosen anew here.

[proved convention] Let

    B_(omega,Lambda)(z)=-a_(omega,Lambda)(z)a_(omega,Lambda)(-z),
    K_(omega,Lambda)=v0 *_x B_(omega,Lambda).              (4)

The B product in the removed-cutoff limit is the C71 distribution. At finite Lambda it is fixed by the regulated Wick expansion. Spatial convolution in (4) has only one v0, whereas C71's four-field coefficient uses two. Its proof applies to one smearing as well. Below this extension and the stronger L2 estimate are proved explicitly. P_beta denotes bosonic periodization for these density kernels. Kernel convolutions without a subscript use both time and space.

## 2. The stronger norm needed to glue two bubbles

[proved lemma] For every fixed M and j<=2,

    ||(1+|t|+|x|)^M partial_A^j K_(omega,Lambda)^0||_1
      +||(1+|t|+|x|)^M partial_A^j K_(omega,Lambda)^0||_2
                                                       <= C_M. (5)

The constants are independent of Lambda>=Lambda0. The zero-temperature kernel is exactly independent of Lambda on this range. Each weighted L2 assertion in (5) is new relative to the L1 assertion used in C73.

[proved] For the cutoff independence, the full plus density response on the transfer support is C71's formula

    Pi_+^0(Omega,p)=p/[2 pi(Ap-i Omega)].                  (6)

In a cross term with u, one momentum lies within K and the other within K+kappa. All regulator factors are therefore one. The same is true in the u/u term. Multiplication by chi(p), rather than chi(p)^2, still restricts the transfer to |p|<=kappa. This proves the stated independence for a single smearing and all the displayed derivatives.

[proved] It remains to bound the limiting distribution. Put c(t,x)=1/[2 pi(At-i x)] and choose a smooth compact eta equal to one near the origin. The distributional identity imported from C71 is B_full=-(i/(2 pi)) partial_x c. It gives

    K_+^0 = -(i/(2 pi)) v0' *_x (eta c)
       +(i/(2 pi)) v0 *_x ((partial_x eta)c)
       +v0 *_x [eta(-2cu_+ + u_+^2)]
       +v0 *_x [(1-eta)a_+^2].                           (7)

Every term other than the first and the cu term is bounded with arbitrary spatial decay on a bounded time interval, and decays to arbitrary order for large time. This uses smoothness of u near the origin and (3) outside that neighborhood. For the two singular terms, if the support of eta is inside |x'|,|t|<=R, then, for 0<|t|<=R,

    integral_(|x'|<=R) |c(t,x')| dx'
                      <= C_R [1+|log |t||].              (8)

To see this, split |x'|<=|t|, where the integral is bounded, from |t|<|x'|<=R, where its majorant is C/|x'|. Schwartz decay of v0 and v0', with bounded x' in these terms, then implies for any fixed N

    |partial_A^j K_+^0(t,x)|
       <= C_N (1+|x|)^(-N)
          {1_(|t|<=1)[1+|log |t||]+(1+|t|)^(-N)}.       (9)

Changing the fixed support radius only changes the constant. A derivatives of c are bounded by C_j/|z|, so (8) applies to them too. Finally integral_0^1 |log t|^r dt is finite for r=1,2. Integrating (9), with N chosen larger than M+3, proves both parts of (5). Reflection proves the minus statement. This proof does not assert a uniform pointwise bound at t=0.

[proved qualification] Two independent spacetime L1 bounds do not control integral K_+(z)K_-(z) dz. The two L2 bounds do. Thus (5) supplies the missing gluing estimate for the particular two-bubble graph, without taking the absolute value of its unsmeared fermion loops.

## 3. Thermal replacement in the stronger norm

[proved lemma] For beta>=1, every fixed M,P>1 and j<=2, define

    Delta K_omega=K_(omega,Lambda)^T-P_beta K_(omega,Lambda)^0,
    w_beta(t,x)=1+d_beta(t,0)+|x|.

Then, on one time period times the spatial line,

    ||w_beta^M partial_A^j Delta K_omega||_1
       +||w_beta^M partial_A^j Delta K_omega||_2
                                                <= C_(M,P) beta^(-P). (10)

The corresponding norms of K^T and P_beta K^0 are bounded without beta^(-P). The constants are uniform in Lambda>=Lambda0. Fixed first and second variations of v0 in a sufficiently strong fixed Schwartz norm are also permitted.

[proved] At finite Lambda the exact difference before spatial convolution is

    B_beta(t,x)-P_beta B_0(t,x)
       =sum_(n!=m) (-1)^(n+m)
                  a_Lambda(t+n beta,x)a_Lambda(t+m beta,x). (11)

The equal-image terms have canceled, including their short-distance distribution. Work on |t|<=beta/2. When one image index is zero, the other is nonzero. For any fixed P,S, (3), with N increased once, gives

    sum_(m!=0) |partial_A^j a_Lambda(t+m beta,x)|
              <= C_(P,S) beta^(-P)(1+|t|+|x|)^(-S).       (12)

Indeed |t+m beta|>=(|m|-1/2)beta and |t+m beta|>=|t|. Distribute decay powers among these two inequalities and |x|, then sum the remaining summable power of |m|. The terms with both indices nonzero satisfy a nonsingular bound of the same form, with arbitrarily large S, after enlarging N again. The restriction n!=m only reduces this bound.

Multiplying (12) by the possible zero-image singularity in (3) shows that the magnitude of (11), including its derivatives, is bounded by

    C beta^(-P) [|z|^(-1)(1+|z|)^(-S)
                                    +(1+|z|)^(-S)].     (13)

Convolving in x with v0, and repeating the integral in (8) near z=0, gives the analogue of (9) with the prefactor beta^(-P). Away from |t|<=1, take S large enough to retain the desired temporal and spatial weights. This proves (10) for both L1 and L2. The derivative allocations have total order at most two and a fixed number of terms. Varying v0 changes only its fixed convolution constants.

For periodization itself, the principal image has (5); all other images of (9) have summably small weighted norms on the principal cell. Minkowski's inequality yields a uniform L2 bound and the triangle inequality gives L1. Equation (10) then bounds K^T. These estimates also justify removing Lambda by dominated convergence, including in L2. No bound on an entire high covariance determinant is used.

[proved check] Unlike the bare absolute bubble, the thermally subtracted loop has at most a single 1/|z| singularity at any point. The spatial smearing converts it into a logarithm in t, which is square integrable. There is no contradiction with C71's logarithmic divergence for the unsmeared absolute square of a propagator.

## 4. Coupled integration of the entire bubble-ring subseries

[proved definition and combinatorics] Retain, in the connected density cumulant expansion of each high chirality, only its two-density cumulant B_omega. This selects a well-defined subseries of the original fermionic vacuum expansion. It is not the assertion that the high density is Gaussian. At order 2m the connected pairings of the two species form one alternating ring. With spacetime volume divided out, their sum is

    R_0(g)=sum_(m>=1) g^(2m)/(2m)
                         (K_+^0*K_-^0)^{*m}(0),
    R_T(g)=sum_(m>=1) g^(2m)/(2m)
                         (K_+^T*beta K_-^T)^{*beta m}(0). (14)

Here R is the contribution to log Z per spacetime volume. Its free-energy contribution is -R. For the combinatorial factor, introduce formal commuting density variables with covariance diag(B_+,B_-) and bilinear coupling g v0 sigma_x. Gaussian algebra, equivalently the sum over their pair contractions, gives -Tr log(I+g diag(B_+,B_-)v0 sigma_x)/2. Odd traces vanish. Each even trace has two equal cyclic blocks, producing 1/(2m) in (14). This formal identity counts exactly the pair-cumulant sector of the fermionic expansion, regardless of the other density cumulants. It neither assumes positivity of the auxiliary covariance nor changes the normalization of Z.

[proved uniform analyticity] For any N>=2 and kernels k_1,...,k_N with the norms above,

    |(k_1*...*k_N)(z)|
                  <= ||k_1||_2 ||k_2||_2 product_(i>=3)||k_i||_1. (15)

For the first two factors use Cauchy-Schwarz at any translation. Convolution with the remaining factors preserves the L-infinity bound. The same proof works on the cylinder. Applying (5) and (10) gives a single C with

    |[g^(2m)] R_0|+|[g^(2m)] R_T| <= C^(2m)/(2m).       (16)

Thus both series are analytic for |g|<1/C, independently of Lambda and beta>=1. Speed or shape derivatives through order two cost at most a fixed polynomial in 2m. On any smaller disk the differentiated series still converge. The order-g^2 coefficient is the actual mixed high vacuum graph, including its factor 1/2, not just a positive upper bound for it.

[proved thermal estimate] For every fixed P>1, on a sufficiently small fixed coupling disk,

    |partial_A^j (R_T-R_0)| <= C_P |g|^2 T^P,
    j<=2.                                                (17)

First replace K_omega^T by P_beta K_omega^0 in each factor. A telescoping expansion has 2m terms. Bound the replaced factor with its L2 estimate in (10), one other factor with its uniform L2 bound, and all others in L1. This contribution is at most (2m) C^(2m) beta^(-P).

The remaining cylinder convolution of periodized zero-temperature kernels is

    sum_(ell in Z) (k_1*...*k_(2m))(ell beta,0).          (18)

This is the periodization identity for integrable convolution kernels. The continuous representative of the convolution is supplied by the two L2 factors. Its ell=0 term is the zero-temperature scalar. For ell!=0, the edge coordinates sum to (ell beta,0), so at least one edge has temporal length at least |ell|beta/(2m). Restrict that factor to this tail, use its weighted L2 norm, use any other factor in L2, and the rest in L1. Equation (15) bounds the resulting convolution at the specified point. Summing the possible long edges gives

    |(k_1*...*k_(2m))(ell beta,0)|
       <= C_P^(2m) (2m)^(P+1) beta^(-P)|ell|^(-P).      (19)

The weighted L2 norms of (5) justify this step. Summation over nonzero ell converges for P>1. The derivative allocations add only another polynomial in 2m. Finally sum (19) and the replacement bound against |g|^(2m)/(2m), reducing the disk once if needed. This proves (17). Moments are chosen in terms of P and derivative order once, not in terms of m.

[proved cutoff statement] The zero-temperature ring kernels are cutoff independent above Lambda0, while the thermal kernels converge in the norms used in (15). The same geometric majorant gives locally uniform convergence of R_T and its indicated derivatives as Lambda tends to infinity. This is a uniform scalar construction for this subseries. It does not supply the free Gaussian trace or the scalar contributions of the omitted mixed cumulants. The free trace remains the one fixed in C65.

## 5. What the mixed loops contribute before a high/low split

[proved coefficient calculation] The full, cutoff-removed free density responses are

    Pi_+(Omega,p)=p/[2 pi(Ap-i Omega)],
    Pi_-(Omega,p)=p/[2 pi(Ap+i Omega)],
    Pi_+ Pi_-=p^2/[4 pi^2(Omega^2+A^2 p^2)].             (20)

At bosonic thermal frequencies the same response follows by integrating the divided Fermi occupation difference over the full line. The integral is p, independently of T. The static and dynamic limits at the origin need not agree. The value at the isolated point is immaterial for the thermodynamic momentum integral.

The actual order-g^2 free-energy coefficient is consequently

    f^(2)(T)=-g^2/2 T sum_Omega integral dp/(2 pi)
                         chi(p)^2 Pi_+(Omega,p)Pi_-(Omega,p). (21)

Using T sum_Omega (Omega^2+E^2)^(-1)=coth(E/(2T))/(2E), subtracting T=0, and combining positive and negative p gives

    f^(2)(T)-f^(2)(0)
       =-g^2/(8 pi^3 A) integral_0^infinity
                        dp p chi(p)^2/(exp(Ap/T)-1)
       =-g^2 T^2/(48 pi A^3)+o(T^2).                    (22)

For the final step set p=Tx/A. The rescaled chi tends to one, and x/(exp x-1) is an integrable majorant. Its integral is pi^2/6 by expansion into positive exponentials. Thus the mixed coefficient is nonzero and is of order T^2. It cannot be discarded as an irrelevant scalar. It agrees with the order-g^2 expansion of the dressed velocity already fixed in C62.

[proved ring resummation] Applying the same ring combinatorics to (20) gives the full-density ring determinant

    f_ring(T)=1/2 T sum_Omega integral dp/(2 pi)
          log[(Omega^2+u(p)^2 p^2)/(Omega^2+A^2 p^2)],
    u(p)^2=A^2-B^2 chi(p)^2.                            (23)

The branch is the analytic one at B=0. Since |B|<A, the logarithmic series converges absolutely in coupling on compact subdisks. Its frequency integral is controlled by the first term, which decays as Omega^(-2) for fixed compactly supported p. The ground-subtracted answer is

    f_ring(T)-f_ring(0)
      =T/pi integral_0^infinity dp
           log[(1-exp(-u(p)p/T))/(1-exp(-Ap/T))].         (24)

One derivation of the frequency sum differentiates the logarithm with respect to u, uses the sum in (22), and integrates back from u=A. The integration constant is zero there. Dominated convergence, as in C59, gives

    f_ring(T)-f_ring(0)
       =-pi T^2/6 [1/sqrt(A^2-B^2)-1/A]+o(T^2).        (25)

[assumed imported identification, not a new off-shell theorem] C62's current representation says that the full unsplit Hamiltonian density sector has the oscillator trace. Equation (23) reproduces its interaction part. This is a consistency check on the mixed-loop normalization, not a second proof of the previously settled trace theorem. It gives no identity between the high-field determinant in (14) and the full effective fermionic action. In the full current theory the higher connected density cumulants vanish with the specified contact convention; C73 shows why that fact does not persist for the high covariance a.

[proved current response within the Gaussian density sector] With the interaction sign in (1), inversion of the two-by-two quadratic density form gives

    Pi_g=(Pi_0^(-1)+g chi sigma_x)^(-1),
    (Pi_g)_(+-)=-g chi Pi_+Pi_-/(1-g^2 chi^2 Pi_+Pi_-).  (26)

Its denominator has poles at the sound dispersion u(p)|p|. This displays how the two chiralities are coupled. The marginal coupling g is already present in (1); integrating a decoupled g=0 model does not generate it. Both (22) and (26) vanish in the corresponding interaction quantities at g=0. The determinant is even in g, whereas its off-diagonal response is odd, as required by the two-species quadratic form.

## 6. The first mixed coefficients omitted by a ring-only integration

[proved exact algebra] Let C_(omega,r) be the connected cumulant of r high density insertions at positions z_1,...,z_r. Normal ordering and inversion symmetry give C_1=0 and C_3=0. At separated points

    C_(omega,2)(z_i,z_j)=B_omega(z_i-z_j),
    C_(omega,4)(z_1,...,z_4)
       =-sum_(directed 4-cycles sigma)
                         product_(i=1)^4 a_omega(z_(sigma_i)-z_(sigma_(i+1))). (27)

The cycles are counted modulo rotation. The minus sign is the connected fermion-loop sign (-1)^(4-1). C73 proves that C_4 generally survives the high/low split. Its formula is a finite distribution after the stipulated smearings, not a new ultraviolet divergence.

Put x_i=(t_i,x_i^space), y_i=(t_i,y_i^space), and write B_omega(ij) for the appropriate pair kernel. Define

    P_omega(1,2,3,4)=B_omega(12)B_omega(34)
                    +B_omega(13)B_omega(24)
                    +B_omega(14)B_omega(23).             (28)

The part of [g^4] log Z per spacetime volume outside the ring sector is

    N_4=1/4! integral_pinned product_(i=1)^4
                                  [dt_i dx_i^space dy_i^space v0(x_i^space-y_i^space)]
          [ C_(+,4)(x_1,...,x_4) C_(-,4)(y_1,...,y_4)
             +C_(+,4)(x_1,...,x_4) P_-(y_1,...,y_4)
             +P_+(x_1,...,x_4) C_(-,4)(y_1,...,y_4) ].  (29)

The notation integral_pinned means removing the common time and spatial translation integral, equivalently dividing by spacetime volume. One may pin t_1 and x_1^space. Formula (29) is understood first at finite regulator. Its zero-temperature cutoff-removed value is fixed by C73's coefficient stabilization theorem.

[proved combinatorial proof] The fourth density moment of each species is C_4+P because C_1=C_3=0. Multiply these two moments, since the two free species are independent. The cumulant of the four interaction vertices removes disconnected unions of two double-vertex components. Among the nine products in P_+P_-, the three identical pair partitions are precisely these disconnected components. The six different pair partitions are connected alternating rings, already included in (14). Every term with a C_4 is connected across all four interaction vertices and survives. The exponential expansion contributes 1/4!, giving (29). This accounts for all vacuum contractions at this order.

[proved limitation of a proposed Gaussian shortcut] Dropping C_4 would change the high-field generating functional. C73 supplies an explicit allowed low-cutoff family with nonzero C_4, so the deletion is not an identity justified by the full-covariance loop cancellation. Formula (29) locates the first missing mixed vacuum terms in that shortcut. It does not show that their fully integrated sum is nonzero for every v0 or every cutoff, and it does not show a failure of TP. In particular, C73 proves that its zero-temperature value is finite at each fixed perturbative order after external projection where needed.

[assumed, unresolved bound] Equations (5)-(19) do not bound the complete sum of mixed graphs containing C_4,C_6 and higher even cumulants together with the open fermion chains. The required estimate is exponential growth in perturbative order in one fixed stronger field/moment norm, with the thermal replacement and scalar derivatives included. A bound on each isolated cumulant or on each fixed-order graph is insufficient. There are two families of fermion permutations to sum, and discarding their determinant structure can cost a factorial. No factorial growth of the actual coefficients is inferred from that poor bound.

## 7. A scale-by-scale check for Gaussian density determinants

[proved restricted comparison] There is no mixed-ring obstruction from replacing a spatial interaction smoothing by a smooth spacetime smoothing once the low-momentum stiffness data agree. The following calculation proves this statement only for the density determinants explicitly defined here.

Let chi_H(p)=chi(p), and let chi_E(Omega,p) be a smooth inversion-even multiplier with chi_E(0,0)=1. Fix the same stable A,B in both determinants. In a neighborhood of zero, chi_H=1 and

    chi_E(Omega,p)-1=O(Omega^2+p^2).                    (30)

The first derivatives at zero vanish by inversion symmetry; Taylor's formula with bounded second derivatives proves (30). Define

    D_j(Omega,p)=Omega^2+[A^2-B^2 chi_j(Omega,p)^2]p^2,
    F(Omega,p)=1/2 log[D_H(Omega,p)/D_E(Omega,p)].        (31)

Take B in a fixed disk small enough for both denominators on real shells to remain comparable to Omega^2+p^2, and choose the analytic logarithm at B=0. On a smooth annulus of radius q_h, (30) gives D_H-D_E=O(B^2 q_h^4), hence

    |F|<=C |B|^2 q_h^2.                                (32)

After rescaling (Omega,p)=q_h(s,r), every fixed derivative of the rescaled F has the same q_h^2 gain. This follows either from the logarithmic power series or differentiating its uniformly invertible denominators. Fixed derivatives with respect to A and shape parameters retain the gain. B derivatives have the reduced coupling power appropriate to differentiating B^2.

Let eta_h be a fixed smooth rescaled annular partition and

    delta e_h(T)=T sum_(Omega=2 pi n T) integral dp/(2 pi)
                                      eta_h(Omega,p) F(Omega,p). (33)

Then, for any fixed integer P>1 and T<=q_h,

    |delta e_h(0)|<=C |B|^2 q_h^4,
    |delta e_h(T)-delta e_h(0)|
                     <=C_P |B|^2 q_h^4 (T/q_h)^P.       (34)

The first bound is the two-dimensional shell measure q_h^2 times (32). For the second, bosonic Poisson summation writes the sum-integral difference as the sum of time Fourier transforms at ell/T, ell!=0. Integrating by parts P times in Omega gives C |B|^2 q_h^(4-P)(T/|ell|)^P. Summing |ell|^(-P) proves (34). The compact smooth annular support removes all boundary terms. This proof also establishes the indicated parameter derivatives.

[proved consequence and limitation] On these determinant shells the smoothing defect has an extra q_h^2. For any 0<theta<1 one may weaken q_h^4 to q_h^(2+theta), and P=3 then has P>2+theta. Thus the two powers requested in TP are compatible in this explicitly computed mixed sector. Equations (31)-(34) do not identify chi_E's determinant with BM's full reference functional or construct the interacting fermionic histories. They use fixed A,B, not scale-dependent histories. They also do not supply the terminal or ultraviolet scalar comparison of the actual models. Claiming TP from (34) would omit precisely those hypotheses.

## 8. Where the actual BM and XXZ comparison stops

[assumed, primary-source input] [BM 0907.2837](https://arxiv.org/html/0907.2837), equations (81), (89), and (93)-(96), provides the constant-history beta-function decay, preparation of its specified Euclidean reference, and matching of limiting marginal and field-normalization data. These passages were inspected through the browser. They concern the model and quantities defined there. No temperature-sensitive vacuum comparison for the C65 trace regulator is imported from them.

[proved logical distinction] C62's conversion of A,B from the matched v,K removes the wrong-speed error in the full reference. Equations (22)-(25) check a mixed coefficient that would otherwise survive as an erroneous T^2 term. Matching those numbers does not identify the complete nonconstant effective potentials. Open chains and high C_4 terms are still present. Nor does matching two zero-temperature critical trajectories give derivative bounds on their parameter rectangles, as C62's explicit recurrence counterexample already demonstrates.

[assumed, remaining constructive statement] Let U denote the complete two-chirality high integration, with all open chains, all mixed cumulants such as (29), and every vacuum contribution. Before an infrared theorem can use the trace reference in this scheme, establish on a fixed coupling disk

    sum_(n>=1) |g|^n ||U_n^0||_(H,M) < infinity,
    ||D_a^j(U^T-P_beta U^0)||_(H',m) <= C T^P,    j<=2, (35)

with its scalar per volume separately included and with compatible parameter directions. H>H' and M are fixed independently of perturbative order and later infrared depth. Appropriate vanishing coupling prefactors apply to undifferentiated terms. The present proof supplies (35) only for the ring vacuum subseries, in its stated thermodynamic setting. C73 supplies the one-chirality statement and complete fixed-order zero-temperature stabilization. Neither is the full estimate (35).

[assumed, remaining infrared statement] With the full prepared reference and XXZ actions in common normalized slow-field coordinates, TP still requires at q_h=2^(-h)

    ||X_h^0-Y_h^0|| <= C q_h^theta,
    ||Y_h^T-P_beta Y_h^0|| <= C(T/q_h)^P,
    ||(X_h^T-Y_h^T)-P_beta(X_h^0-Y_h^0)||
                                <= C q_h^theta(T/q_h)^P, (36)

and the complete vacuum-map bounds

    ||DE_(h,T)||+||D^2E_(h,T)|| <= C q_h^2,
    ||D(E_(h,T)-E_(h,0))|| <= C q_h^2(T/q_h)^P.          (37)

The common localization identification is included in P_beta. The physical counterterms remain the zero-temperature ones. A single fixed reserve of moments and field weights must handle both marks, history derivatives, and scalar roots, with P>2+theta. The ultraviolet and terminal estimates and appropriate thermodynamic-regulator control are part of TP, not optional additions.

[assumed, imported conditional implications] C62 already sums the shared-short-memory factors. The collaborator's MOMENT_RESERVE pays a fixed thermal mark conditional on ordinary routing, and C57/C59 close the four-history recursion and scalar sum conditional on (36)-(37). C65 fixes the trace normalization. Their combination would give the desired XXZ thermal coefficient if the missing constructive estimates hold. No microscopic or retarded-endpoint hypothesis is newly verified here.

[proved final status] The leading mixed vacuum coefficient and the complete alternating ring subseries are finite and thermally controlled at fixed high scale. They do not obstruct TP. The first non-ring mixed coefficients are displayed in (29), but no estimate on their all-order completion has been proved and no divergent coefficient has been found. The requested full two-chirality integration and actual infrared comparison remain unfinished.

## 9. Checks and work boundary

[numerical, checks specified] checks.py independently checks alternating-pairing combinatorics, the complete order-four connected cumulant decomposition using a finite Grassmann algebra, and the second-order thermal integral. It also checks the q_h^2 Gaussian determinant defect by direct evaluation and the bosonic Poisson shell difference for an explicit annulus. These controls do not test (35)-(37) for XXZ.

[numerical, completed] The final run passes. The connected alternating-pairing counts at orders 2,4,6,8 are 1,6,120,5040, giving the factors 1/(2m). Rational arithmetic verifies (29) against the logarithm of the finite density-source polynomial. The second-order thermal integral at A=1.3, g=0.2, T=0.05 divided by its limiting T^2 coefficient is 0.9999999999998982. The frequency-logarithm quadrature residual is at most 6.33e-17. The Gaussian shell comparisons are finite controls, not interval-certified estimates.

[numerical, failed check and correction] The first run imposed the unsupported test condition |delta e_h(T)-delta e_h(0)|/[q_h^4(T/q_h)^3]<0.01. At q_h=0.125 and T/q_h=0.2 the ratio was 0.01759624465840341, so that assertion failed. No analytical argument fixed the constant to 0.01. The revised script reports the normalized ratios, verifies their finiteness and decreasing absolute errors over the sampled thermal scales, and retains the other exact and convergence assertions. The largest reported ratio in the final sample is 0.026165663261909018. This correction does not prove a numerical uniform constant or alter (34).

[proved scope] All new files and all Python execution are inside work_C78. There are no shell network commands, large computations, subagents, author messages, or edits to parent files. The local STATE.md supplies the handoff under the explicit write restriction. These derivations have not received independent or human review.
