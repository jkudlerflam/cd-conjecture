# C50. Constructive RG interface and partial extensions

[proved, status] This note does not prove S(b). It proves a finite-scale integration theorem for general even kernels, a retardation localization estimate, and conditional estimates that specify the remaining all-scale problem. No microscopic power-law exponent or gapless-mode count is claimed as an established result.

[assumed, inputs] Fix n and put N=n+1. Import C39/C42's finite spectral data, C44's rooted length/payload representation, C45's normalization W, and C46 Section 5's exact charges and momentum classification. Do not assume C45's RG estimates. Assume a C48/C49 reduction supplies the logical theory. At the time of the initial inspection those runs supplied no completed action theorem. The results below say what a usable output must contain.

[proved, conventions] A proved tag refers to the stated mathematical implication, including any explicit assumptions. An assumed tag marks an input not established here. A read tag refers to SOURCE_AUDIT.md. Numerical checks are supplementary. No probability estimates are used. All writes stay in work_C50.

## 1. The proposed implication and its quantifiers

[assumed, desired theorem] Let p_alpha, alpha=1,...,2N, be distinct simple free Fermi points with signed velocities u_alpha. There are N positive and N negative velocities. Let delta>0 bound every symmetry-allowed non-density quadratic or quartic mismatch. For the stronger unrestricted version, delta is the full pair-sum margin of C46 (24). Require the complete interaction, not only its quartic part, to be small in an explicitly stated coefficient norm. The intended conclusion is a fixed-coupling infinite-volume, zero-temperature Luttinger theory with 2N chiral branches, power-law correlation asymptotics, and summable marginal beta functions.

[proved, quantifiers] The required smallness has the form

    for every fixed reference R_n, there exists epsilon_0(R_n)>0,
    for every interaction in the specified epsilon_0 ball,
    estimates hold uniformly as L,beta tend to infinity.             (1.1)

It is not a single epsilon_0 valid as n tends to infinity or min|u_alpha| tends to zero. One chooses couplings before L,beta. This is the ordering in C45. The auxiliary parameters may enter the conversion to the normalized kernel norm, but should not enter the subsequent RG constants if that norm and the normalized reference are fixed.

[read] SOURCE_AUDIT.md identifies the applicable statements in MP and PSS. In particular the many-channel reference beta-function theorem already exists. The genuinely new task is an extension of the lattice-to-reference comparison to the assumed kernel family.

[assumed, sufficient format of input] In retained units use a Grassmann action

    V(psi) = sum_(p even, p>=2) (1/p!)
             sum_labels integral K_p(X) psi_1 ... psi_p.           (1.2)

Kernels can be finite complex measures in relative times, so instantaneous delta functions are allowed. Translation invariance factors off one space-time volume. Sum absolute coefficients over all orbital, charge and creation/annihilation labels. Kernels must preserve the exact charges. Source kernels must be given separately if source correlations are required. A scalar kernel is retained in the generating functional but does not enter the normalized correlation argument.

[proved, precise meaning of "general even"] Here it means every even field degree compatible with the exact charges. Parity alone is not the assumed symmetry. The inherited payload U(1) and rooted length-number symmetry are the reason anomalous number-changing terms are excluded.

## 2. Normalizing retarded kernels without losing the auxiliary hierarchy

[proved, lemma 1] Suppose the physical coefficient measure, with one time fixed at zero, has moments

    sum_labels integral (1+Delta T)^s (1+ell(X))^s |K_p|
        <= epsilon W a_p,                                         (2.1)

where T=sum_(i>=2)|t_i|, ell is a spatial connected length, and s>=2. Define retained time t'=Wt, without rescaling the lattice fermion. For an ordinary density kernel,

    K'_p(t'_1,...,t'_p)=W^(-p) K_p(t'_1/W,...,t'_p/W).              (2.2)

For a measure this means the corresponding pushforward, with the overall factor W^(-1) after anchoring. Then, for 0<=j<=s,

    sum_labels integral (T')^j (1+ell(X))^s |K'_p|
        <= epsilon a_p (W/Delta)^j.                               (2.3)

Proof. The p time differentials in the action each contribute W^(-1). Integrating the p-1 relative times contributes W^(p-1), leaving W^(-1). A jth time moment contributes W^j. Inequality (2.1) bounds the physical jth moment by epsilon W a_p Delta^(-j). Multiplying proves (2.3). Approximation by measures, or the defining pushforward identity, gives the same result for delta functions. No covariance estimate is used.

[proved, consequence] If W/Delta<=1, the normalized zero, first and second time moments are bounded by epsilon a_p independently of Delta. The spatial weights are unchanged. Thus a bound already supplied in the strong physical norm (2.1) removes explicit inverse auxiliary-gap factors from the subsequent endpoint estimates. A claim of exponential decay without its prefactor, degree sum and moments is insufficient to apply this lemma.

## 3. Separated velocities at the propagator level

[proved, lemma 2] Let f be a smooth function supported in a fixed annulus c<=sqrt(q0^2+q1^2)<=C, and define the continuum single-scale covariance

    g_(u,E)(t,x) = integral d k0 d k/(2pi)^2
                  exp[-i(k0 t+k x)] f(k0/E,u k/E)/(i k0+u k),
    E>0, u!=0.                                                   (3.1)

For any nonnegative integers a,b,M there are constants depending on f,a,b,M, but not on E or u, such that

    |partial_t^a partial_x^b g_(u,E)(t,x)|
      <= C E^(1+a+b)/|u|^(1+b)
         [1+E|t|+E|x|/|u|]^(-M).                                 (3.2)

Moreover

    integral dt dx |g_(u,E)| <= C E^(-1),
    integral dt dx (|t|+|x|)^s |g_(u,E)|
      <= C_s E^(-1-s) max(1,|u|)^s.                              (3.3)

Proof. Substitute k0=Eq0 and k=Eq1/u. The Fourier multiplier becomes E^(-1) f(q)/(iq0+q1), and the measure contributes E^2/|u|. Its inverse Fourier transform is Schwartz because the annulus avoids zero. Differentiation supplies E^a(E/|u|)^b. Integration by parts gives (3.2). For (3.3), use y0=Et and y1=Ex/|u|, whose Jacobian is |u|E^(-2), and integrate the Schwartz majorant. This proves every displayed bound.

[proved, determinant bound] Write the multiplier in (3.1) as m(k). Factor it into |m|^(1/2) and its phase times |m|^(1/2). Translates of these factors give a Gram representation of the covariance. Each squared Hilbert norm equals integral |m| <= C E/|u|. Therefore an r by r covariance determinant is bounded by

    (C E/v_min)^r,       v_min=min_alpha |u_alpha|>0,              (3.4)

for finitely many branches, using their Hilbert direct sum. The same estimate holds after multiplication by a positive interpolation matrix t_ij with t_ii=1. Factor t_ij as inner products of unit vectors and tensor these vectors into the Gram representation. The Gram determinant inequality then gives (3.4).

[proved, interpretation] No denominator |u_alpha-u_beta| occurs. Degenerate speeds cause no problem here. Small absolute speeds do enter the field/determinant bounds. Spatial derivative estimates also cost inverse speeds. This proves a propagator statement, not an all-scale interacting theorem.

[proved, finite-volume qualification] The continuous estimate (3.2) is not asserted at arbitrarily small energy on a finite spatial ring. A spatial zero mode can remain when E<<|u|/L. In Section 5 we instead use a fixed positive matching energy and a smooth covariance, for which periodization and sufficiently large L,beta give uniform constants. The final infrared shell at arbitrary aspect ratio is part of the unproved full extension.

## 4. Taylor localization for retarded interactions

[proved, lemma 3] Anchor one coordinate of a translation-invariant p-field kernel. Put

    R(X)=sum_(i>=2)(|t_i|+|x_i|),
    M_j(K)=sum_labels integral R(X)^j |K(X)|,  j=0,1,2.            (4.1)

If M_j is finite, the Fourier transform satisfies

    |Khat(q)-Khat(0)| <= Q M_1(K),
    |Khat(q)-Khat(0)-D Khat(0)[q]| <= (Q^2/2) M_2(K),          (4.2)

where Q=max_i max(|q_i,0|,|q_i,1|). Here q denotes independent external deviations after imposing total momentum conservation. Constants arising from selecting those independent variables can instead be absorbed into Q.

Proof. Write Khat(q)=integral exp(i sum_i q_i dot X_i) K(dX). The scalar phase z obeys |z|<=Q R. Use |exp(iz)-1|<=|z| and |exp(iz)-1-iz|<=z^2/2 for real z, then integrate the total variation. Dominated convergence supplies the derivative and proves the second formula, including singular measures.

[proved, power-counting consequence] In 1+1 dimensions a fermion has engineering dimension 1/2. A local p-field term has coupling dimension 2-p/2. Subtracting a constant from a quartic kernel gives one external momentum and hence dimension -1. Subtracting constant and linear terms from a quadratic kernel gives two external momenta and again dimension -1. Every local p>=6 term has dimension at most -1 before anomalous dimensions. On an energy-E patch, Q<=C E max(1,v_min^(-1)), so the gain has an explicit velocity cost. This is why two moments are enough for these Taylor remainders. It does not establish an interacting anomalous-dimension bound.

[proved, inherited momentum step] C46 Section 6 already proves support exclusion for forbidden quadratic and quartic branch assignments after shrinking patches below the mismatch margin. Apply it with the dressed roots only once their displacement is controlled. This is an exact support statement, valid at every perturbative order. Higher-degree vertices may contract into permitted density quartics and diagonal quadratic terms. Their irrelevance as uncontracted vertices does not make those induced terms zero.

[proved, local marginal form] A surviving quartic has equal incoming and outgoing branch multisets. Anticommuting the fields writes it as a multiple of rho_alpha rho_beta, rho_alpha=psi_alpha^+ psi_alpha^-. If alpha=beta its strictly local Grassmann monomial vanishes. Nonzero point-split same-branch terms belong to derivative terms or quadratic renormalizations after contact contractions. Diagonal quadratic constants and first derivatives must be included among the running parameters. No microscopic U(1) identity sets them to zero.

## 5. A finite-scale theorem for arbitrary-degree kernels

[proved, definition of the norm] Work first on a finite space-time torus with finitely many Fourier modes. Let Vol=beta L. For a list X of space-time points, let T(X) be the infimum of lengths of connected networks containing X, in the periodic metric d((t,x),(s,y))=|t-s|_beta+|x-y|_L. Auxiliary network vertices are allowed. For a singleton T=0. Put w_s(X)=(1+T(X))^s, s>=0. The properties used below are

    T(Y)<=T(X) if Y is a subset of X,
    T(union_i X_i)<=sum_i T(X_i)+sum_(ij in tree) d(x_ij,x_ji).
                                                                    (5.1)

In the second inequality x_ij is a selected point of X_i. Both follow by taking unions of connected networks, then taking infima. Since 1+sum a_i<=product(1+a_i), the weight of the union is bounded by the product of vertex weights and line weights.

For an expansion (1.2), define

    ||V||_(rho,s)= sum_p rho^p/(p! Vol)
                   sum_labels integral w_s(X)|K_p(X)|.             (5.2)

A scalar term c has norm |c|/Vol. The integral over all coordinates in (5.2) divided by Vol equals the anchored integral by translation invariance. One may use any explicit coefficient representation of the polynomial. The bound does not rely on finding the minimal representation or on cancellations between equivalent representations. Antisymmetrization by averaging does not increase this norm.

[proved, theorem 4] Let C be a translation-invariant, number-conserving Grassmann covariance with the following two bounds, including all flavour sums in K_s.

    sup_(X,label) sum_label' integral dY (1+d(X,Y))^s |C(X,Y)| <= K_s,
    |det C_t(X_i,Y_j)| <= D^r for every r by r interpolated minor.  (5.3)

Here the interpolation matrix is any positive Gram matrix of unit vectors associated with the clusters. Set

    R=max(1,rho,sqrt(D)),  H=8 exp(2) R,  K=max(1,K_s).

If V is even and ||V||_(H,s)<=epsilon with K epsilon<1, the connected effective action

    F_C(V)(psi)=log E_C exp[V(psi+zeta)]                          (5.4)

has an absolutely convergent connected expansion, modulo the convention that its scalar is divided by Vol. In the coefficient norm,

    ||F_C(V)||_(rho,s) <= epsilon/(1-K epsilon).                   (5.5)

The constants are independent of Vol and of the number of Fourier modes whenever (5.3) and the input norm are uniform. This is a finite-integration theorem. It does not assert preservation of the same field radius under infinitely many integrations.

[proved, proof of theorem 4] Expand the logarithm in truncated expectations,

    F_C(V)=sum_(r>=1) (1/r!) E_C^T(V;...;V).                      (5.6)

The fermionic tree interpolation identity expands the rth truncated expectation over trees on r labelled interaction vertices. Each of its r-1 tree edges contracts one field at each endpoint. The fields not contracted by tree edges and not retained as external fields form a determinant with the interpolated covariance in (5.3). This is the finite-algebra identity used in PSS (4.109); its hypotheses are Gaussian Grassmann integration and even vertices, not a two-body interaction. We use this identity as an established algebraic tool, not as an assumed RG bound.

Fix vertex degrees p_1,...,p_r and P=sum_i p_i. For a fixed tree with degree d_i at vertex i, assigning fields to its incident edges costs at most product_i p_i^(d_i). Selecting external fields costs at most 2^P. Allow another factor 2^P for orientations and signs before taking absolute values. If q fields remain external, the determinant has (P-2(r-1)-q)/2 rows. Its absolute value times the external radius is bounded by

    D^((P-2(r-1)-q)/2) rho^q <= R^P.                            (5.7)

Terms with an inadmissible parity or negative row number are zero and need not be counted.

The weight of the external coordinates is bounded using (5.1). Integrate each kernel in its coordinates relative to a chosen vertex center. With those offsets fixed, integrate all centers but one along the tree, starting at its leaves. Each line integration costs K_s. This remains true with arbitrary offsets because a translation of its relative argument leaves the line integral unchanged. The remaining common center contributes Vol. The endpoint weights are the weights in (5.2). Thus the coordinate bound is

    Vol K_s^(r-1) product_i [anchored weighted kernel norm at i].
                                                                    (5.8)

The same proof applies to q=0, for which division by Vol is precisely the scalar convention. Measures with time contacts are handled by total variation and Fubini; no product of singular measures at coincident arguments is introduced, since the contracted covariance is bounded at this finite scale.

For r>=2, the labelled-tree sum obeys the weighted Cayley identity

    sum_trees product_i p_i^(d_i)=(product_i p_i) P^(r-2).        (5.9)

For r=2 it is immediate. For r>=3 the Prufer sequence has r-2 entries, each entry i supplies p_i, and every vertex supplies one additional p_i. This proves (5.9). Since product_i p_i<=exp(P), and P^(r-2)/r!<=exp(P) for P>=1, the tree sum divided by r! is at most exp(2P). The latter inequality follows from P^r/r!<=sum_j P^j/j!=exp(P).

Combining (5.7)-(5.9) absorbs all degree factors into (4 exp(2)R)^P, which is smaller than H^P. The 1/p_i! factors in (1.2) are retained in the endpoint norms. Equivalently, grouping the output into antisymmetric kernels with 1/q! normalization gives the same external-subset count. Summing independently over endpoint degrees and labels yields

    ||(1/r!) E_C^T(V;...;V)||_(rho,s) <= K^(r-1) epsilon^r.      (5.10)

For r=1, direct Wick contraction gives the bound epsilon, also covered by (5.10). Sum the geometric series to get (5.5). At finite cutoff the connected expansion is the Taylor expansion of the logarithm at V=0. Absolute convergence supplies its analytic branch throughout the stated ball. If all bounds persist along a cutoff limit, the resulting coefficient series converges in the same norm. That last sentence requires the uniform bounds and convergence of the input covariance and kernels, not just finite-cutoff existence. This proves the theorem.

[proved, corollary 5, finite crossover] Suppose the assumed action is already defined below a fixed temporal cutoff Lambda in retained units. Choose a positive matching energy E_d below all root-patch, band-isolation and quartic-detuning scales. Split off a smooth covariance C_> supported between E_d and Lambda, including all bands in that window. At fixed reference data, its multiplier is smooth and supported away from every free pole. Fourier integration by parts gives finite weighted L1 constants K_s; factoring the multiplier gives a finite Gram constant D. The matrix-valued version follows from polar decomposition and the finite number of orbitals. Periodizing the smooth kernel preserves weighted L1 bounds; for sufficiently large L,beta its sampled compact Fourier multiplier has a uniform Gram bound. Theorem 4 therefore integrates this entire crossover window at once, uniformly in L,beta, provided epsilon is small in its specified H-radius norm.

Proof of the volume assertion. Periodization expresses the kernel as a sum of translated infinite-volume kernels with possible fermionic signs. Taking absolute values removes the signs, and summing a fundamental domain is bounded by the infinite-volume integral/sum. Periodic distances are at most their lifted distances. The discrete Gram norm is a Riemann sum of a smooth compactly supported absolute multiplier, bounded uniformly once both meshes are sufficiently fine. There is no infrared singularity because E_d>0. This proves the corollary.

[proved, what corollary 5 resolves] Extremely small E_d changes K_s,D and the permitted epsilon. It does not make them infinite at fixed n. No hypothesis of comparable velocities is used. General even, time-nonlocal vertices are permitted with the degree-weighted norm stated above. This gives a rigorous finite crossover result, not merely a count of the number of crossover shells.

[assumed, limitations for the proposed application] The supplied logical output must meet that norm at radius H and must already have a controlled ultraviolet representation. A radius-2 norm alone does not automatically satisfy a larger-radius estimate. A local Hamiltonian ultraviolet limit requires its own determinant/time-ordering estimates; the elementary compact-frequency Gram argument does not prove them uniformly as Lambda tends to infinity. These requirements have not been verified for an actual C48 action. Nor does (5.5) iterate indefinitely with a fixed radius, since it trades H for rho. An infrared proof needs dimensional rescaling, localization and the reference cancellation together.

## 6. Higher-degree vertices change marginal initial data at first order

[proved, lemma 6] There is no general extension with a beta-function error necessarily quadratic in the norm of every added irrelevant kernel. A single six-field vertex can induce a quartic at first order.

Proof. Split fields into retained psi and integrated zeta and consider an endpoint containing

    epsilon rho_1(x) rho_2(x) zeta_3^+(y) zeta_3^-(z).             (6.1)

Gaussian expectation gives

    epsilon c(y,z) rho_1(x) rho_2(x),
    c(y,z)=E_C[zeta_3^+(y) zeta_3^-(z)].                         (6.2)

Choose a finite relative displacement for which c is nonzero and translate the entire vertex together. A conjugate endpoint can be added for the physical reality condition. The resulting action preserves number and translation. Smooth relative kernels can replace the contacts. The first term of the connected logarithm is (6.2), so its size is O(epsilon), not O(epsilon^2). Such an endpoint is a term in the expansion of an ordinary six-field interaction after splitting its fields into psi+zeta. This proves the claim. It does not claim a gap or a relevant instability.

[proved, consequence] The extended comparison estimate must allow an O(epsilon_irr) shift of the matched quartic and quadratic parameters. The inherited C45 drift budget already allows linear irrelevant contributions. Reusing the O(lambda^2) coefficient of a special quartic microscopic expansion without rematching would be incorrect. A forward-scattering beta function can tend to zero while its total accumulated shift is nonzero.

## 7. A sufficient all-scale comparison estimate, and its conditional closure

[assumed, the missing estimate] After the finite matching, let g_h collect marginal forward couplings and velocities, with h<=0 decreasing into the infrared. Use a fixed collective norm, including flavour sums and relative velocity weights. Let B_h=||g_h-g_(h+1)||. A sufficient comparison estimate has the following form for constants A,b,kappa,theta>0 with theta<kappa:

    B_h <= A epsilon 2^(theta h)
       + b epsilon sum_(k=h+1)^0 2^(-kappa(k-h))
                                  sum_(j=h+1)^k B_j.             (7.1)

The first term includes curvature, retarded remainders and arbitrary-degree endpoints, after their local induced terms have been matched. The second measures the dependence of the reference beta function on a nonconstant history. The reference constant-history cancellation is the published input. The full bound (7.1) for the supplied kernel class is not proved here or asserted to be a stated theorem of either paper.

[proved, theorem 7] If (7.1) holds on every finite range h_min<=h<=0 with the same constants, then for

    b epsilon C_(kappa,theta)<=1/2,
    C_(kappa,theta)=2^(-(kappa-theta))/
       [(1-2^(-kappa))(1-2^(-(kappa-theta)))],                     (7.2)

one has uniformly in h_min

    B_h<=2A epsilon 2^(theta h).                                 (7.3)

Consequently any compatible infinite-scale trajectory has a limit g_* and

    ||g_h-g_*|| <= [2A epsilon/(1-2^(-theta))] 2^(theta h).        (7.4)

Proof. At h=0 the double sum is empty. Assume (7.3) at all j>h. Extending the nonnegative sums in (7.1) to infinity bounds their contribution by

    2Ab epsilon^2 2^(theta h)
       sum_(d>=1) 2^(-kappa d) sum_(r=1)^d 2^(theta r).

Reverse the order of summation. The double sum equals C_(kappa,theta) in (7.2). Equation (7.2) bounds the contribution by A epsilon 2^(theta h), closing induction. Sum (7.3) over j<=h to get (7.4), with a harmless larger endpoint convention. Finite-volume families still need convergence or compactness and consistency to produce a trajectory; the theorem is a uniform bound, not a substitute for constructing the RG map. This proves every stated implication.

[proved, why this is the central gap] A bound B_h<=C epsilon^2 without a positive scale power has an infinite sum as h tends to minus infinity. Smallness at one scale and power counting do not repair it. Conversely (7.1) and the inherited root/stiffness reserves would bound the full marginal drift. Merely writing a differential majorant x'<=Cx^2 gives a finite RG-time guarantee, not an infrared theorem at fixed coupling.

[proved, meaning of vanishing beta function] The appropriate conclusion is summable decay such as (7.3) and a limiting fixed point. It is not the statement that every lattice Wilsonian beta coefficient is identically zero at every scale. Lemma 6 already exhibits a finite-scale change that the latter statement would incorrectly forbid.

## 8. Fixed physical parameters and the Fermi-point equation

[read] PSS uses added-and-subtracted quadratic terms. The original Hamiltonian stays fixed, and the roots shift. A requirement to tune independent microscopic chemical potentials is not a hypothesis to import from that procedure.

[proved, lemma 8] Let e_alpha be fixed bare dispersions on disjoint intervals I_alpha=[p_alpha-r,p_alpha+r], with e_alpha(p_alpha)=mu and |e'_alpha|>=v_min/2 throughout. Let nu_alpha(p') be real continuous functions on the product of these intervals, with

    max_alpha |nu_alpha(p')|<=C_nu epsilon,
    2C_nu epsilon/v_min<r.                                    (8.1)

Then there exists p' in the product satisfying

    e_alpha(p'_alpha)-mu+nu_alpha(p')=0                         (8.2)

for every alpha. Every such solution satisfies

    max_alpha |p'_alpha-p_alpha|<=2C_nu epsilon/v_min.          (8.3)

If nu is Lipschitz in the sup norm with constant L_nu<v_min/2, this solution is unique.

Proof. Each e_alpha is strictly monotone with inverse Lipschitz constant at most 2/v_min. Its image contains [mu-v_min r/2,mu+v_min r/2]. Define F_alpha(p')=e_alpha^(-1)(mu-nu_alpha(p')). Equation (8.1) ensures F maps the compact convex product interval into itself. Continuity and Brouwer give a fixed point. The inverse Lipschitz estimate gives (8.3). If the stronger Lipschitz condition holds, F is a contraction with constant 2L_nu/v_min, proving uniqueness. The bare e_alpha and mu have not changed.

[proved, margin and finite-ring consequences] If the right-hand side of (8.3) is less than delta/8, every four-root mismatch initially at least delta remains at least delta/2 by the triangle inequality. Rounding the roots to the ring mesh costs at most pi/L per root and hence at most 4pi/L in a quartic mismatch. Once L is sufficiently large at fixed delta, this loses only a further reserved fraction. Equation (8.2) then has a residual bounded by sup|e'_alpha| pi/L if the counterterm is held fixed during rounding. This is a finite-mesh statement, not an estimate on interacting energy levels.

[assumed, missing physical input] For the supplied general action, the existence, reality, uniform smallness and continuity of the counterterm map nu(p') must come from the constructive RG. Lemma 8 removes a purely finite-dimensional inverse problem once those facts hold. It does not construct that map. A frequency-dependent critical self-energy must not be assumed C1 at the origin without proof.

## 9. Source normalization and correlation remainders

[proved, lemma 9] Suppose a source or field normalization obeys

    log Z_(h-1)-log Z_h=eta log 2+r_h,
    |r_h|<=C epsilon 2^(theta h),  h<=0, theta>0.              (9.1)

Then for a nonzero finite constant Z_bar,

    Z_h=Z_bar 2^(-eta h)[1+O(epsilon 2^(theta h))], h->-infinity.
                                                                    (9.2)

Proof. Sum (9.1) from h through zero, with the index shifted at most once according to the initial convention. The sum of r_j converges absolutely. Absorb its full sum and the initial Z into log Z_bar. The missing tail has absolute value bounded by a geometric series C'epsilon 2^(theta h). Exponentiating proves (9.2). This proof also specifies the sign convention for eta; no sign is inferred from a source's notation.

[proved, lemma 10] Suppose the difference between a matched correlation and its reference correlation has a shell expansion R(x)=sum_(h<=0) R_h(x) satisfying

    |R_h(x)|<=C 2^(a h) [1+2^h |x|]^(-M),
    a>0, M>a.                                                  (9.3)

Then for |x|>=1,

    |R(x)|<=C' |x|^(-a).                                       (9.4)

Proof. Choose h_x with 1/2<2^(h_x)|x|<=1. For h<=h_x, bound the denominator below by 1 and sum the geometric series in 2^(a h). For h>h_x, bound it below by (2^h|x|)^M and sum 2^((a-M)h). Both sums are bounded by constants times |x|^(-a). This proves (9.4).

[proved, conditional power-law implication] If the reference two-point function for branch alpha has a nonzero leading amplitude of degree -(1+eta_alpha), and if the matched observable has nonzero overlap with that field, then a comparison (9.3) with a=1+eta_alpha+theta proves that the logical correlation has the same leading exponent. The angle dependence can contain several velocities; a single Lorentz-invariant denominator is not assumed. The free limit has eta_alpha=0. The anomalous exponents depend on the matched forward matrix, not just on n or the free velocities.

[assumed, missing source result] No bounds (9.1) or (9.3), no matched forward matrix, and no nonzero source overlap for the full supplied retarded theory have been established here. Theorem 4 applies to even interactions at a finite integration step; inserting and controlling fermion or composite sources through the entire infrared iteration is additional work. Current power laws alone would not determine every charged-vertex exponent. C51 separately addresses what a complete Gaussian limit would imply for modes and central charge.

## 10. What has and has not been extended

[proved] Sections 2-6 establish normalization uniformity under a stated massive-kernel moment bound, unequal-velocity propagator estimates, Taylor localization for measure-valued retarded kernels, finite-scale determinant integration for arbitrary even degree in a specified field-radius norm, and the necessity of linear rematching of irrelevant endpoints. Sections 7-9 prove conditional numerical-sequence, root-inversion and source-summation implications. These are genuine bounded statements and do not use S(b).

[assumed, remaining RG theorem] To finish the route one must establish all of the following for one explicitly specified input norm.

1. A uniform ultraviolet initialization, if the input is a Hamiltonian or an action without a temporal cutoff, and a summed degree bound sufficient for the determinant estimates.
2. A renormalized multiscale expansion with arbitrary-degree retarded endpoints that preserves a usable norm after localization and field rescaling.
3. A reference comparison estimate of the strength of (7.1), including endpoint degree sums, velocity dependence and all irrelevant insertions.
4. A real continuous small counterterm map for fixed bare parameters, and compatible finite-volume limits.
5. Source comparison bounds strong enough to identify nonzero leading power-law amplitudes and their reference exponents.

[proved, separation from massive matching] These are requested under the assumption that C48 or C49 has supplied its logical output. Their absence is not repaired by assuming S(a). Conversely, the failure to prove them here is not evidence that the logical model is gapped. No counterexample to the proposed phase has been found.

[proved, final mathematical status] The requested many-flavour Luttinger theorem, with n+1 gapless modes at fixed coupling, remains unproved in C50. The finite-scale theorem and the primary-source audit narrow the missing work, but replacing items 2-5 by the phrase "the usual RG argument" would assume the central part of the question.
