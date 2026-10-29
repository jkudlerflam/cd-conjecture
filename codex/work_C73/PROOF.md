# C73. Chiral preparation before coefficient norms

## 1. Result and scope

[proved partial result] There is a cutoff-independent analytic integration of one high chirality of the C65 Hamiltonian reference, including its closed density cycles and open fermion chains. Its thermal replacement has any fixed power P>1 in every fixed polynomial moment norm, with two speed derivatives. This statement is at a fixed preparation scale. The other high chirality has not been integrated by this theorem.

[proved partial result] After projecting external fields onto a fixed slow momentum support, every fixed-order zero-temperature coefficient of the full reference preparation, including its vacuum coefficient, becomes exactly independent of the sea cutoff once that cutoff exceeds a constant times the perturbative order. This is a coefficientwise statement. Its proof does not give an order-independent analytic radius for the full preparation.

[proved scope qualification] The new kernel estimates below are stated in infinite spatial volume, on the line or thermal time cylinder. No additional uniform finite-spatial-volume theorem is claimed. The exact finite-regulator algebra is used to fix the identities, while the norm bounds are proved in the displayed thermodynamic kernel spaces. This is another reason not to identify them with the complete regulator-uniform TP statement.

[assumed, unresolved target] The full analytic two-chirality preparation and the subsequent lattice/reference estimates of C59 Lemma TP remain unproved here. No coefficient counterexample to TP is established. In particular, a nonzero cutoff remainder in a Ward identity is not a divergent coefficient. Section 8 states the remaining estimate and why the partial theorems do not imply it.

[assumed, imported inputs] Use the reference, trace normalization and parameter conversion of C62 and C65. In their conventions g=2 pi B, A>|B|, and the dressed velocity is sqrt(A^2-B^2). The interaction is

    V = g integral dt dx dy v(x-y) rho_+(t,x) rho_-(t,y),
    v_hat = chi,     supp chi subset [-kappa,kappa].       (1)

The function chi is real, even and smooth, and is one near zero. Thus v has all polynomial spatial moments. The densities are normal ordered relative to the free sea. This is the spatially smoothed Hamiltonian reference, not BM's different spacetime-smoothed functional. C65 fixes its total scalar by the finite-regulator Duhamel identity and the Gibbs trace limit.

[assumed, imported covariance bounds] Use exactly C71's splitting

    c_Lambda(k0,k) = w_Lambda(k)/(A k-i k0),
    u(k0,k) = f(k0,k)/(A k-i k0),
    a_Lambda = c_Lambda-u,       a = c-u.                 (2)

Here w_Lambda=r(k/Lambda)^2, f is smooth, compactly supported and inversion even, f=1 near zero, and its spatial support is inside |k|<=K. The regulator is on a plateau throughout the low support. C71 proves, uniformly for Lambda>=Lambda0 and A in a compact positive interval,

    |partial_A^j a_Lambda(z)| <= b_N(z),
    b_N(z)=C_N |z|^(-1)(1+|z|)^(-N),   j=0,1,2.          (3)

The same bound holds for a. At every fixed finite Lambda these kernels are also bounded at short distance, with a bound depending on Lambda. The discontinuous equal-time value is specified by the Hamiltonian convention. Self-contractions are zero in the sea-normal-ordered, inversion-symmetric splitting. The opposite chirality is obtained by spatial reflection.

[proved construction] Set the preparation scale to one. For an l-field translation-invariant coefficient measure W_l define

    ||W_l||_m = integral |W_l(0,X_2,...,X_l)|
                            (1+diam(X))^m,
    ||W||_(H,m) = sum_(l even>=2) H^l ||W_l||_m.          (4)

The diameter uses spacetime distance, with circular temporal distance at temperature T. The first field is pinned. Instantaneous and coincident-density delta measures are allowed. Symmetrization and fermion antisymmetrization are normalized averages, so their operator norms in (4) are at most one. We give explicit ordered representations of coefficients and bound their variation before that final average. There is no norm on individual bare loop lines in (4).

Periodization means inserting antiperiodic fields into the line functional and integrating the common anchor time over one period. This also specifies its signs for nonlocal kernels. Density insertions are periodic. The free Gaussian partition factor is kept separately, as in C65. No scalar is adjusted by a choice of normalization in this note.

## 2. All long chiral cycles have uniform absolute moments

[proved lemma] Let r>=3. For a directed cycle with z_1=0, z_(r+1)=z_1 put

    L_(r,Lambda)(z_1,...,z_r)
                      = product_(i=1)^r a_Lambda(z_i-z_(i+1)).

For each fixed m and j<=2,

    ||partial_A^j L_(r,Lambda)||_m
                              <= C_m^r r^(m+3).          (5)

The power of r is deliberately loose. The bound is uniform in Lambda, and these kernels converge in the same norm as Lambda tends to infinity.

Proof. Let b dominate the derivatives in (3), taking N large enough for the chosen m. Put b_m(z)=(1+|z|)^m b(z). Both b and b_m belong to L^1 and L^(3/2). At the origin this follows from

    integral_0^1 R^(1-p) dR = 1/(2-p),   p<2.

At infinity choose N>m+3. Interpolation gives a bound C_m for their L^p norms simultaneously for every 1<=p<=3/2.

Set p_r=r/(r-1). It belongs to [1,3/2]. The generalized Young inequality gives

    (b_1*...*b_r)(0) <= product_i ||b_i||_(p_r),
    sum_i 1/p_r = r-1,                                (6)

for nonnegative b_i equal to b or b_m. One can obtain this form by iterating Young's inequality, with intermediate reciprocal exponents 1-j/r after j factors, and ending in L^infinity. The convolution has a continuous representative at this endpoint, or its value can be obtained by truncation and monotone convergence. The same bound holds at every translation of its last argument.

The integral of the unweighted cycle is the left side of (6) with all b_i=b. For the weight, every pair of vertices can be joined along the cycle, hence

    diam(z_1,...,z_r) <= sum_i |z_i-z_(i+1)|.

For m>=1 use

    (1+sum_i |e_i|)^m <= r^(m-1) sum_i (1+|e_i|)^m.

Thus each summand in the weighted integral has one b_m factor and r-1 b factors. For m=0 use (6) directly. At most r^2 allocations arise from two speed derivatives. Enlarging C_m proves (5).

Pointwise a_Lambda tends to a away from the origin. The product majorant used above is integrable in all relative coordinates. Dominated convergence proves convergence in (4), also after the two derivatives. This establishes the claim without treating a squared propagator as integrable. QED.

[proved distinction] For r=2 the Young exponent would be two, and (3) is not in L^2 locally. This is the exceptional cycle already treated in C71. Its spatially smeared coefficient, with the Hamiltonian distributional extension, has bounded moments. The trace-log cycle of length one is the normal-order tadpole and is zero here. There is therefore no new isolated-cycle ultraviolet singularity at length three or greater.

[proved path lemma] An open chain of r vertices has r-1 propagators. Its pinned integral is bounded by ||b||_1^(r-1), because its successive edge coordinates are independent. With a fixed diameter moment and two speed derivatives, the bound is C_m^r times a fixed polynomial in r. Dominated convergence again gives cutoff removal. This is the open-field counterpart of (5).

## 3. Thermal winding of an arbitrary cycle

[proved lemma] Let L_(r,Lambda,beta) be the product of r antiperiodized propagators, r>=3. Let Per_beta L_(r,Lambda) mean periodization of all r vertex times as a density kernel. For beta>=1, P>1, and j<=2,

    ||partial_A^j[L_(r,Lambda,beta)-Per_beta L_(r,Lambda)]||_(m,beta)
               <= C_(m,P)^r r^(m+P+5) beta^(-P).         (7)

No moment order depending on r is required.

Proof. Expand each edge propagator as its image sum. The r integer image labels can be changed into r-1 vertex lifts and one total winding ell around the cycle. Unfold the r-1 vertex lifts. The ell=0 contribution is precisely the periodized zero-temperature cycle. The remaining terms carry sign (-1)^ell and have ell!=0. After taking absolute values and pinning one vertex, they are bounded by convolutions with constraint

    e_1+...+e_r = ell beta e_t.                         (8)

Circular external time distances are at most the lengths of paths in the unfolded graph. Thus the same diameter-weight allocation as in Section 2 applies to the e_i. At least one edge satisfies |e_i|>=|ell| beta/r.

Choose N in (3) sufficiently large in terms of m,P only. Uniformly in 1<=p<=3/2 and R>=0,

    ||b_m 1_(|z|>=R)||_p <= C_(m,P)(1+R)^(-P).          (9)

For R<=1 use the full L^p bound. For R>=1 take out (1+R)^(-P) and use the L^p norm of b_(m+P), with N increased once. This proves (9) including nonintegral fixed P.

Allocate the diameter moment to edge j and the long-edge condition to edge i. Apply (6) to the r factors, using (9) on edge i, and using b_m there if i=j. Summing these choices, and then the speed derivative choices, costs at most the polynomial displayed in (7). Since

    (1+|ell| beta/r)^(-P) <= r^P beta^(-P)|ell|^(-P),

the remaining sum over ell!=0 is 2 zeta(P). This proves (7). The bounds justify the unfolding and all image sums by absolute convergence. They also allow cutoff removal. QED.

[proved path identity] An open chain has no independent winding. Its r-1 edge image labels are in bijection with r-1 vertex lifts. The product of edge signs telescopes to the signs of the two external fermionic endpoints. Consequently the thermal open-chain functional is exactly its multivariate periodization at fixed parameters. It has zero thermal defect of the kind in (7). This statement is about an unlocalized coefficient, not about the thermal RG history or its localization coordinates.

[assumed, imported two-cycle estimate] C71 Section 5 proves the r=2 replacement after the two spatial potential smearings, with any P>1 and j<=2. It deals with its nonintegrable equal-image square before taking norms. Combined with (7), this covers every closed chiral cycle. The normal-order tadpole is zero at both temperatures by inversion symmetry; it is not subtracted using a temperature-dependent counterterm.

## 4. Analytic preparation of one high chirality

[proved theorem] Keep all minus-chirality fields external for this operation, and integrate only the plus high field eta with covariance a_Lambda. Let psi denote the plus slow field and set

    R(t,x)=integral dy v(x-y) rho_-(t,y),    M=g R.

The resulting interaction functional F_(+,Lambda) is analytic in g on a disk independent of Lambda and beta>=1. For every fixed H,m and P>1 there is a smaller fixed disk on which

    ||F_(+,Lambda)^0||_(H,m) <= C |g|,
    ||partial_A^j(F_(+,Lambda)^T-Per_beta F_(+,Lambda)^0)||_(H,m)
                                      <= C |g|^2 T^P,
    j=0,1,2.                                             (10)

The cutoff limit exists in these norms, locally uniformly on that disk. Fixed first and second variations of v in a stronger fixed Schwartz norm are also allowed. Derivatives with respect to g itself have the corresponding reduced powers of g, rather than the prefactor written in (10).

Proof of the algebra. At finite Grassmann regulators, Gaussian completion gives

    integral dmu_a(eta)
      exp[-(eta^+ + psi^+) M (eta^- + psi^-)]
       = det(1+a M)
           exp[-psi^+ M(1+a M)^(-1) psi^-].              (11)

This is verified by completing the quadratic form a^(-1)+M. Multiplication of the Schur complement gives

    M-M(a^(-1)+M)^(-1)M = M(1+a M)^(-1).

The entries of M are even elements of the remaining Grassmann algebra and commute, so the ordinary matrix identities apply. Singular Gaussian regulators can be handled by a nonsingular approximation or by the polynomial Wick identity. Normal-order constants are those in (2), not additional terms in (11).

Taking minus the connected logarithm yields

    F_+ = psi^+ sum_(r>=1) (-1)^(r-1) M(a M)^(r-1) psi^-
                   - sum_(r>=1) (-1)^(r+1) Tr[(a M)^r]/r. (12)

Thus an open chain with r potential insertions has 2r+2 external fermion fields, and a cycle with r insertions has 2r external fields. The r=1 chain is the original quartic. The r=2 chain is the nonlocal six-field contribution highlighted in C59, now retained. The trace term with r=1 vanishes. The trace terms with odd r vanish as well, but this extra cancellation is not required for the convergence bound.

Proof of the norm estimate. For an open chain integrate its independent relative plus coordinates using the path lemma. Every smearing adds a factor bounded by

    V_m = integral dx (1+|x|)^m |v(x)|.

For a cycle, pin the first remaining minus density coordinate. The smearing at that vertex integrates the plus anchor, and the other smearings integrate their spatial offsets. The cycle lemma controls the remaining relative spacetime coordinates. The output density fields are duplicated at the same coordinates, represented by finite delta measures. The elementary inequality 1+|x+y|<=(1+|x|)(1+|y|) and the diameter allocation in Section 2 control the smearing moments. Hence, for r>=3, the coefficient bounds have the form

    chains:  H^(2r+2) |g|^r C_m^r r^d V_m^r,
    cycles:  H^(2r)   |g|^r C_m^r r^d V_m^r/r,           (13)

with d fixed in terms of m and the derivative count. The two-cycle is instead bounded using C71, with phi=v*v and one spatial derivative of phi in its proof. A fixed stronger norm of v pays that derivative. The r=1 local chain is bounded directly.

There is no factorial count hidden in (13). Formula (12) is the trace-log identity, already summing the permutations of a closed connected fermion contraction. Equivalently, a symmetric r-density kernel has (r-1)! cyclic contractions and the exponential contributes 1/r!, leaving 1/r. Taking an unnormalized sum of all those cycles after dropping 1/r! would destroy this accounting. The normalized coefficient convention in (4) preserves it.

Choose |g| C_m H^2 max(1,V_m)<rho<1. The sums of r^d rho^r converge. This proves the first estimate in (10). Equations (7) and the path identity prove the second after one fixed reduction of the g disk to pay r^(m+P+5). The two-cycle uses the imported estimate. All thermal corrections start at r=2. Derivatives in A or v distribute over at most r factors twice, so cost only a fixed polynomial. This establishes their bounds without a loss at each value of r.

Dominated convergence from Section 2, the path lemma, and C71 for the bubble proves the cutoff limit term by term. The same geometric majorant permits summation and proves locally uniform analytic convergence. Equation (11) identifies this convergent series with the finite-regulator analytic germ; it retains every term of that germ for the integration performed. QED.

[proved scalar qualification] No interaction vacuum scalar is produced by this first integration alone: every term in (12) still contains at least two minus fields. A cycle is a scalar only with respect to the integrated plus variables. It becomes a physical vacuum term only after the remaining minus fields are contracted. The full free trace factor is still required. Equation (10) is not a bound on the final vacuum map in TP3.

[proved norm qualification] The minus fields in (10) are external algebraic variables, comprising both their eventual high and low parts. A coefficient norm for these variables does not justify integrating them against the unbounded continuum high covariance. The finite-scale Grassmann theorem requires its own covariance/determinant bound. C71's logarithm prevents inserting the entire remaining high covariance into that theorem with the earlier absolute majorant.

## 5. Ward identities after the high/low split

[proved algebraic identity] Write D(k)=A k_1-i k_0 and a(k)=(1-f(k))/D(k), after removal of the sea cutoff. Since D(k+q)-D(k)=D(q),

    D(q) a(k)a(k+q)
      = a(k)-a(k+q)+f(k)a(k+q)-f(k+q)a(k).              (14)

Multiply both sides by D(k)D(k+q) to check the identity. The last two terms are fixed-low-cutoff insertions. They persist when Lambda tends to infinity. At finite Lambda the more general formula is

    D(q) a_Lambda(k)a_Lambda(k+q)
          = F_Lambda(k+q)a_Lambda(k)
                              -F_Lambda(k)a_Lambda(k+q),
    F_Lambda=w_Lambda-f.                               (15)

Replacing the right side of (14) by its first two terms is therefore not an exact Ward identity of the prepared action. The missing terms are not small merely because the sea cutoff is large. They can be finite kernels, as the next calculation demonstrates.

[proved full-cycle identity] For distinct complex coordinates Z_i=A t_i-i x_i and r>=3,

    sum_(oriented cycles sigma on {1,...,r})
           product_(i=1)^r 1/(Z_(sigma_i)-Z_(sigma_(i+1))) = 0. (16)

Cycles are counted modulo cyclic rotation, not reversal. For a proof, view the sum as a rational function of Z_1, keeping the other Z_i distinct. A possible pole at Z_1=Z_j comes only from cycles containing the adjacent pair 1,j. Pair the cycles obtained by inserting 1 immediately before or immediately after j into each cycle on the other r-1 vertices. Their residues are opposite and their remaining factors coincide at Z_1=Z_j. Hence all residues vanish. Each term has two factors involving Z_1 and tends to zero at infinity. A rational function without poles and vanishing at infinity is zero. For r=3 this is also immediate by reversal. This proves (16) as a separated-point identity. We do not use this argument to choose contact terms of singular distributions; the two-cycle contact term remains the Hamiltonian one in C71.

[proved odd-cycle identity] Every inversion-odd covariance, including a_Lambda and its antiperiodization, has zero symmetrized cycle of odd length r>=3. Reversing a directed cycle multiplies its product by (-1)^r. The reversed cycles have the same fermion-loop permutation sign. Pairing them proves the cancellation before taking any norm.

[proved nonzero four-cycle example] The even high cycles need not obey (16). For A=1 take first the rapidly decreasing low multiplier f_G(k)=exp[-(k_0^2+k_1^2)]. Its high covariance is

    a_G(t,x)=exp[-(t^2+x^2)/4]/[2 pi(t-i x)].             (17)

To verify (17), use

    (1-exp[-|k|^2])/(k_1-i k_0)
            = (k_1+i k_0) integral_0^1 exp[-s|k|^2] ds.

The inverse Fourier heat kernel is (4 pi s)^(-1) exp[-|z|^2/(4s)]. Applying the numerator's derivatives and performing the s integral gives (17), with the Fourier convention of C71.

At z_j=(j-1,0), j=1,...,4, put a_d=exp[-d^2/4]/(2 pi d). Summing the six oriented cycles gives

    S_4(a_G)=2[-a_1^3 a_3+a_1^2 a_2^2+a_1 a_2^2 a_3]
      = 2/(2 pi)^4 [-exp(-3)/3+exp(-5/2)/4+exp(-9/2)/12]
      >0.                                              (18)

The first two displayed contributions already have positive sum, since exp(1/2)>4/3. With a_d=1/(2 pi d), the same expression vanishes, agreeing with (16).

This example can be placed inside the allowed compact plateau cutoff class. Choose smooth radial f_eps, equal to one in a shrinking disk, equal to f_G on an expanding annulus, and zero outside a larger disk. It can be chosen between zero and one. Then

    integral d^2k |f_eps-f_G|/|k| -> 0.

Near zero, the difference is O(|k|^2); its radial integral tends to zero. At infinity use the Gaussian tail. Fourier inversion shows u_eps-u_G tends uniformly to zero. Thus a_eps at the finitely many nonzero separations in (18) tends to a_G, and its four-cycle remains nonzero for sufficiently small eps. These f_eps satisfy the compactness, plateau and inversion assumptions in (2).

[proved interpretation] Length four is the first length greater than two at which the symmetrized high cycle can survive in this class. It enters (12) at g^4, with eight remaining minus fields. The example proves that a general preparation cannot impose the full-covariance zero identity on all its high cycles. It does not show a norm obstruction. Section 2 bounds this very coefficient, and Section 3 bounds its thermal replacement. Spatial smearing is still applied when it is inserted into (1).

[assumed, infrared Ward requirement] BM's vanishing beta function involves localized sums, cutoff correction identities, and matching, rather than deleting each high cycle separately. The new estimates do not establish its differentiated four-history vacuum extension. In particular, (14) is an elementary line identity, not a substitute for the complete interacting Ward system.

## 6. Exact cutoff stabilization at each fixed order

[proved theorem] Fix the preparation covariance and project every final external slow field with a smooth spatial multiplier s(k), equal to one on the low covariance support and supported in |k|<=K_e. This leaves the functional on slow fields unchanged. Let U_(n,Lambda)^0 denote the coefficient of g^n of the full high integration, in this representative, including its scalar per spacetime volume. For

    R_n=max(K,K_e)+(n+1)kappa+1,

one has

    U_(n,Lambda)^0=U_(n,R_n)^0,       Lambda>=R_n.         (19)

The statement concerns the thermodynamic zero-temperature coefficient. It does not assert equality of finite-temperature functions at a finite sea cutoff. Every coefficient in (19) has finite fixed polynomial moments, and its first two A derivatives do as well. No bound on their growth with n is asserted.

Proof, graph organization. At fixed perturbative order there are finitely many connected Wick graphs. Each vertex of (1) has one incoming and one outgoing field in each chirality. For each chirality, its internal contractions form disjoint open chains and closed cycles. An open chain ends on final external slow fields. A closed one-cycle is the zero normal-order tadpole. At each vertex the change of fermion spatial momentum has magnitude at most kappa, since it is the interaction transfer. This remains true whether the field at that vertex is high or low.

On an open chain, every line momentum differs from its external endpoint momentum by a sum of at most n such transfers. Thus every high line on that chain lies in |k|<=K_e+n kappa. All its w_Lambda factors are one for Lambda>=R_n.

On a closed cycle expand each high line c_Lambda-u. If at least one line is u, its momentum lies in |k|<=K. Every other line then lies in |k|<=K+n kappa. Again all c_Lambda weights are one, and the term is independent of Lambda. This expansion is finite at fixed n.

It remains to treat a cycle of length r>=2 consisting only of c_Lambda lines. Fix its interaction transfer momenta and frequencies. If its successive spatial shifts are s_j and frequency shifts are nu_j, its internal frequency integral has the form

    I_r(k)=integral d omega/(2 pi)
             product_(j=1)^r [A(k+s_j)-i(omega+nu_j)]^(-1). (20)

The product decays as |omega|^(-r), so its contour can be closed for r>=2. If all k+s_j have the same sign, all poles lie in the same half-plane. Closing in the opposite half-plane gives I_r(k)=0. This argument also holds for repeated poles. It therefore has support, as a function or distribution in the other transfers, only where

    -max_j s_j <= k <= -min_j s_j.                     (21)

Throughout this interval every |k+s_j| is at most max s-min s<=r kappa<=n kappa. All the regulator weights, independent of omega, equal one there. Thus the pure full cycle also agrees for every Lambda>=R_n. Unlike (16), this step does not require summing different orientations of the loop.

Loop-frequency integration can be performed at fixed transfers independently for each closed chiral component. Transfers and base loop frequencies are independent graph momentum variables. At a finite spatial regulator these frequency integrals are convergent for cycles of length at least two. Mixed u cycles have the compact-frequency factor, and open chains have no independent closed-loop frequency. This justifies the operation before further transfer integration; it is also obtained by first using the time-ordered finite-regulator Wick formula. No shift of a divergent one-line frequency integral is used. The same-sign zero in (20) persists under A differentiation while A stays positive.

One can avoid any ambiguity at exceptional external transfers by applying this argument directly to the difference of two plateau cutoffs. Its nonzero spatial weights require at least one line momentum outside R_n. Since the spread is at most n kappa, all line momenta then have the same sign and are bounded away from zero. The frequency contour vanishes there without an infrared singularity. Thus the difference is zero also after testing against smooth external transfer functions, not merely at generic nonzero momenta. In particular no additional contact distribution is left to choose.

We have exhausted the possible chiral components. Their finite sum proves (19) in external Fourier variables, and hence as coefficient distributions. For the scalar there are no open chains, so the same proof applies. The external multiplier s is essential to the assertion about open chains. Without choosing a slow-field representative, an arbitrary off-shell external momentum could itself be outside every fixed cutoff plateau.

Proof of the finite-moment assertion. Use the finite cutoff R_n on the right side of (19), without expanding c-u again. Its high covariance is bounded and has arbitrary polynomial decay away from the origin. In a connected graph choose a spanning tree on the n interaction vertices using n-1 of its fermion contractions. Use coordinates (t_i,x_i,y_i), and put d_i=y_i-x_i. Each spatial potential contributes |v(d_i)|. With the d_i fixed, integrate the independent relative t_i,x_i along the spanning tree, using the weighted L^1 norm of a_(R_n). Every remaining propagator is bounded by its weighted L^infinity norm. Integrate the d_i with the corresponding weighted norm of v. One vertex time and spatial coordinate are pinned; for a vacuum graph their integral is instead the spacetime volume divided out of its scalar.

Any fixed output diameter moment is controlled by distances along this tree and the n spatial offsets. Increasing their moment order once by a fixed amount pays it. External spatial projections are convolutions by a Schwartz function and preserve these bounds. The finitely many speed-derivative allocations use the same finite-cutoff bounds. Thus every graph at this fixed n has finite moments, and so does their finite sum. This completes the proof. QED.

[proved limitation] Equation (19) does not prove cutoff-independent analyticity. Its replacement cutoff R_n grows with n, and the finite-cutoff L^infinity and determinant constants used in the last paragraph grow with R_n. There is no established bound

    ||U_(n,R_n)^0||_(H,m) <= C^n                        (22)

with a single C. A sequence of finite stabilized coefficients may still have no positive radius of norm convergence. Conversely, a poor upper bound growing faster than exponentially is not evidence that the actual coefficients grow that way.

[proved consequence for a coefficient search] Within this fixed reference, projection convention, and zero-temperature ultraviolet question, no finite perturbative order can give a new log Lambda or power-of-Lambda divergence of the integrated coefficient norm. The norm at that order equals a finite-cutoff norm by (19). This conclusion does not exclude failure of an infrared depth estimate, a uniform sum over orders, or a thermal estimate for the fully integrated coefficients. It also does not cover arbitrary retarded endpoints or a local unsmoothed continuum interaction with unbounded transfer momentum.

## 7. Why the two chiralities still need a constructive bound

[proved diagnostic] Even after integrating each separate density loop before taking its norm, contracting their resulting kernels can be ill-controlled by their individual L^1 norms. For example the order-g^2 vacuum graph is proportional to

    integral dt dx dy phi(x-y) B_+(t,x) B_-(t,y),
    phi=v*v,    B_omega=-a_omega(z)a_omega(-z).           (23)

The two bubbles share their time difference t. A separate spacetime L^1 bound on each smeared bubble is not a bound on this shared-time product. Algebraic cancellation, smoothing and time estimates have to be kept through the gluing. This is a failure of an inference, not a claim that (23) diverges. Equation (19) proves its integrated zero-temperature coefficient is finite.

More generally, (12) has arbitrary powers of the remaining density. They are not arbitrary bounded commuting sources once those fermions are integrated. The contraction count of two ungrouped chiral Wick permutations can have factorial growth. The exact determinant identities of the original Grassmann integral prevent one from identifying that count with the growth of the actual coefficients. They also explain why replacing the second integration by a positive sum of all such contractions is not a completed proof.

[heuristic next organization] A full construction could retain the source-resolvent formula (11) inside a second determinant expansion, or use a multiscale expansion that groups closed chiral cycles before its positive bounds. Both must handle the sharing of vertex times between the two chiralities and the fixed-low-cutoff insertions in (14). The analytic bounds in Section 4 provide reusable endpoints for that construction, but do not supply its determinant or connected-graph estimate.

## 8. The exact remaining statement for XXZ

[assumed, unresolved ultraviolet estimate] For the full reference high integration, prove a connected expansion on a fixed complex disk in g such that, in one fixed stronger moment/field norm,

    sum_(n>=1) |g|^n ||U_n^0||_(H,M) < infinity,
    ||D_a^j(U^T-Per_beta U^0)||_(H',m) <= C |g| T^P,
    j<=2,   H>H',   M fixed,                            (24)

and prove the corresponding full interaction scalar bounds per spacetime volume, retaining the free trace. Coupling derivatives have their own prefactors as in Section 4. The reference covariance and interaction parameter variations must belong to the actual domains used in the later comparison. Section 4 proves (24) only for the first chirality; Section 6 proves only the zero-temperature coefficientwise cutoff removal for the complete action.

[assumed, unresolved infrared estimate] For weak XXZ off half filling, after zero-temperature counterterm and finite field-amplitude matching, C59 TP still requires, at q_h=2^(-h),

    ||X_h^0-Y_h^0|| <= C q_h^theta,
    ||Y_h^T-Per_beta Y_h^0|| <= C(T/q_h)^P,
    ||(X_h^T-Y_h^T)-Per_beta(X_h^0-Y_h^0)||
                              <= C q_h^theta(T/q_h)^P,  (25)

with P>2+theta, together with

    ||D E_(h,T)||+||D^2 E_(h,T)|| <= C q_h^2,
    ||D(E_(h,T)-E_(h,0))|| <= C q_h^2(T/q_h)^P.          (26)

The periodization here includes the common finite-temperature localization identification. A norm of source cycles alone cannot replace the complete effective-kernel/history norm. One must control the actual parameter rectangles, ordinary localization distance routing, endpoint degrees, normalization ratios and the vacuum roots. The moment and field reserves must be fixed independently of h.

[assumed, imported conditional reductions] The collaborator's MOMENT_RESERVE Section 6 already allocates one fixed thermal mark once ordinary routing is established. C62 performs the shared-short-memory double scale sum. C57/C59 give the four-history closure and the scalar summation. C65 fixes the reference trace, whose thermal coefficient C59/C62 evaluate. Those results are used with their stated hypotheses; none proves (24)-(26). They imply the desired XXZ thermal coefficient if these missing estimates and terminal control are supplied.

[proved final status] This run establishes a uniform analytic one-chirality preparation and a complete zero-temperature fixed-order cutoff stabilization theorem for the reference. It neither proves TP nor locates a coefficient that refutes it. The surviving g^4 high cycle is a counterexample only to discarding all cycles above length two. It obeys the prepared norm and thermal estimate. No microscopic hypothesis of the multichannel encoded-chain application is newly verified here.

## 9. Checks and source boundary

[numerical] checks.py independently checks the contour formula (20) by frequency quadrature, its finite-cutoff plateau for several cycles, the high-cutoff line identity (14), the nonzero four-cycle (18), and finite-image cycle unfolding. These are small controls of algebra and coefficients, not an interacting XXZ flow or a test of (24)-(26). checks.json records the output.

[numerical, truncation correction] The first finite-image test used vertex lifts from -5 to 5. At cycle length three and beta=2 its omitted image tail gave residual 1.30e-11, failing the 2e-12 tolerance. Expanding the lifts to -8 through 8 and the winding range to -16 through 16, with the tolerance unchanged, gives maximum residual 1.08e-16 across the final cases. The contour quadrature residual is 2.22e-16 and the plateau residual is 6.08e-18. These errors are numerical controls, not certified interval bounds.

[proved analytical checks] The isolated-cycle Young exponent is below two exactly for r>=3, identifying the exceptional bubble. All open-chain thermal images are accounted for by external lifts, leaving no unassigned winding. The two-cycle uses its fixed Hamiltonian contact convention. The fixed-order plateau bound depends on n and therefore does not imply analyticity. At g=0 the one-chirality effective interaction and its thermal replacement vanish. At fixed nonzero g no thermal XXZ conclusion is drawn. Speed derivatives are controlled only on a compact positive interval, consistent with excluding a band edge.

[proved work boundary] All new files are inside work_C73. The parent project memory is not modified under the explicit folder restriction. No subagents, shell network commands, large computations or author correspondence were used. The arguments have not received independent or human review. SOURCE_AUDIT.md distinguishes the primary-source context from the new proofs above.
