# C55. Single-channel retarded initialization and an all-scale linear endpoint theorem

[proved, status] The nonlinear Luttinger theorem is not proved here. Sections 2–5 give an off-half-filled initialization and first infrared step, using the collaborator's ultraviolet theorem and C50's finite-integration theorem. Sections 6–8 prove an all-scale result for the term linear in an arbitrary-degree initial action, including a continuous first-order Fermi counterterm. Section 9 states the remaining nonlinear lemma. Section 10 proves a counterterm existence implication that requires only continuity, with its hypotheses explicit. Section 11 explains the additional source conclusion needed for the exponents.

[convention] A proved tag applies to the stated implication, including its named hypotheses. Imported results are identified separately. No statement below identifies the central charge, the thermal coefficient, or the selected physical type sector. All coordinates are logical cell coordinates and all energies are in hopping units.

## 1. Input and choice of model

[imported] Use the nearest-neighbor Jp/Jz tuning in [SINGLE_CHANNEL_RG.md](../../contrib/decisive_attack_20261028/demonstrator/SINGLE_CHANNEL_RG.md), Section 1. Its reference hopping and longitudinal coefficient are

    t = Jp Jz kap^2/[4D(D^2-kap^2)],
    u = 3 Jp^2 kap^4/[32D(D^2-kap^2)(4D^2-kap^2)].                 (1)

Thus 2u/|t| can be made small at fixed nonzero t by choosing the ratio first. The physical bandwidth unit is W=|t| lambda^2. Choose lambda only after that ratio, the cell gap, the decay weights and all crossover constants. The matched correction has norm O(lambda/|t|). C53's diagonal counterterm gives another tunable reference, but its BPPB term has physical range four. It is not the nearest-neighbor option used here.

[assumed, matching input] Take the zero-temperature finite-cell action and its symmetry-preserving source construction as supplied by the finite-cell matching theorem, with the strengthened connecting-path norm in the collaborator's Section 5. The input is

    S = S_free + V_XXZ + R,
    ||V_XXZ|| <= C delta,       ||R|| <= C lambda/|t|,              (2)

in a specified, sufficiently large field radius and a positive exponential connecting-path weight. Kernels can be measures in relative times. This includes instantaneous vertices. The action is invariant under translation, logical U(1), spatial reflection and the reality symmetry inherited from the Hermitian real microscopic Hamiltonian. These are hypotheses on the full matched action, not just its quartic truncation. Off half filling there is no particle-hole symmetry at fixed field.

[proved, field implementation] The finite-cell conserved charge Q of the collaborator's Section 3 commutes with its cell Hamiltonian and acts on the code as 1/2-n_f. Adding a uniform field of order W times Q therefore gives the desired logical field and preserves the code. Since Q is bounded on the finite cell, the gap changes by at most twice the field magnitude times ||Q||. It stays positive for sufficiently small lambda. This field is a sum of on-site charge terms in the physical cell pattern, so it does not require the BPPB counterterm. Finite matching constants can be chosen uniformly on a fixed compact interval of dimensionless logical fields by reserving half the original gap.

[definition] Fix an interval I compactly contained in (0,pi), not meeting pi/2. For p in I put

    e_p(k)=2(cos p-cos k),       v_p=2 sin p,
    v_min=min_(p in I) v_p>0,
    d_* = min_(p in I){dist(2p,2pi Z),dist(4p,2pi Z)}>0.         (3)

A staggered fermion rotation fixes the hopping sign if necessary. The two free roots are +p and -p. All smallness constants below may depend on I and the finite cell. No estimate is claimed uniformly near a band edge or half filling.

## 2. Off-half-filled ultraviolet integration

[imported] The collaborator's Sections 5.1–5.3 and [UV_AUDIT.md](../../contrib/decisive_attack_20261028/exact/UV_AUDIT.md) prove a determinant and decay initialization for general small retarded endpoints. C50, PROOF Section 5, proves the finite-scale connected-integration estimate. We use those results rather than infer a new ultraviolet theorem from infrared power counting.

[proved extension] Their off-half-filled constants can be made uniform on I. Set

    chi_G(k0,k)=exp[-k0^2-e_p(k)^2],
    g_UV=(1-chi_G)/(-ik0+e_p(k)).                              (4)

The denominator cancels in the identity

    g_UV=(ik0+e_p(k)) integral_0^1 exp[-s(k0^2+e_p(k)^2)] ds.    (5)

For |Im k|<=a, |e_p(k)|<=E_a=2+2 cosh a. Substitution in the collaborator's heat-kernel estimate gives finite exponential row integrals at any fixed smaller spatial rate and fixed temporal rate, uniformly in p. Nothing in that estimate uses cos p=0.

For completeness the convolution filter is bounded by the elementary operator-series estimate. Write e_p=2cos p-S-S^(-1), where S is the shift. Its convolution row norm is at most A=4. Submultiplicativity gives

    ||exp[-e_p^2]||_row <= exp(A^2) <= exp(16).                (6)

The time heat filter has total variation at most one, also after antiperiodic periodization. The audited row-filter determinant lemma therefore gives the safe constant

    delta_UV = 2 sqrt(1+exp(16)).                             (7)

The improved half-filled constant is unnecessary. The constant in (7) also controls interpolation minors. Spatial and temporal periodization preserve the decay estimates. Thus the published arbitrary-vertex integration theorem, with the collaborator's weight reserve, integrates (4) uniformly in volume and inverse temperature. It produces a small all-degree action, with the same U(1), reflection and reality symmetries. It does not set its quadratic constant to zero.

[proved, compact-patch crossover] Choose smooth, inversion-even branch cutoffs chi_omega(k0,k-omega p), with disjoint supports, equal to one sufficiently close to each pole, and with deviation widths a_* satisfying 4a_*<d_*. They can be chosen uniformly on I. The desired retained covariance is

    g_patch = sum_omega chi_omega/(-ik0+e_p(k)).               (8)

The remaining covariance C=chi_G g-g_patch has an integrable Fourier multiplier and a smooth space-time kernel with all fixed polynomial moments. Near a pole chi_omega=1, and (chi_G-1)g is smooth by (5). Away from the poles the denominator is separated from zero on compact sets; at large frequency chi_G g has Gaussian decay. These observations, followed by integration by parts, give the claimed moments. Square-root factorization of the absolute Fourier multiplier gives a finite Gram bound. The bounds are uniform in p in I. Periodization gives the corresponding uniform finite-regulator bounds once the meshes resolve this fixed crossover.

Apply C50's finite-integration theorem to C. Alternatively split C into the gapped compact-patch complement of chi_G g and the regular correction chi_patch(chi_G-1)g and apply it twice. This explains why a Gaussian cutoff alone must not be used for exact momentum-support exclusions. A finite number of field-radius losses is harmless, since the matching norm can be required at any prescribed finite radius before choosing lambda.

[proved, norm used below] Following the collaborator, use

    w_q(X)=min_(terminal trees T) product_(ij in T)(1+r_ij)^q,
    r_ij=|x_i-x_j|+|tau_i-tau_j|.                            (9)

Exponential connecting-path control implies this norm at any fixed q after a geometric loss in field radius. The connected-network proof of their Section 5.3 gives output w_q from input w_(2q) and covariance moment 2q. Hence the compact crossover produces such a tree norm. It is not necessary, or correct, to assert exponential locality for a compact smooth Fourier cutoff.

We will also use the anchored radius R_X=sum_(i>1)r_i1. For p external fields and a tree minimizing (9), R_X<=p sum_(e in T)r_e. Consequently

    (1+R_X)^s <= p^s product_(e in T)(1+r_e)^s <= p^s w_q(X)
    whenever q>=s.                                         (10)

A fixed loss of field radius absorbs p^s. This proves the moment conversion needed below. It does not prove that a fixed q pays for every repeated localization in an infinite nonlinear tree.

## 3. Degree summation and the microscopic hierarchy

[definition] Write an even, number-conserving action in antisymmetric coefficient form, with 1/p! included at degree p. Put

    ||V||_(H,s)=sum_(p even>=2) H^p/p!
                 sum_labels integral (1+R_X)^s |K_p|(0,dX'). (11)

Use an explicit coefficient representation, not the infimum over representations. A scalar is excluded from this norm and can be tracked separately. At a finite integration a sufficiently large H pays for field choices and determinant powers. The same convention will be used for output kernels. In Gaussian integration formulas V denotes the signed potential in exp[V], so V=-(V_XXZ+R) for (2). Inverse-action coefficients are extracted from minus the connected logarithm. All norm estimates are independent of that sign.

[imported] For the actual finite-cell residual, [DEGREE_REPARAMETERIZATION.md](../../contrib/direct_attack_20261028/vacancy/DEGREE_REPARAMETERIZATION.md) supplies the support/order mechanism. If physical order l has degree at most D l and starts at l=3, its normalized degree-p kernel admits the bound C^p epsilon^max(2,p/2-1), after choosing a possibly very small epsilon and a fixed artificial analytic parameter. The static weak XXZ term must also be included in that smallness budget. The hierarchy survives a finite integration because contractions only remove fields.

[proved qualification] This is additional microscopic information. A bound ||K_p||<=eta H^(-p) alone does not imply powers eta^(p/2-1). Fix p>4 and send eta to zero to see the failure. The linear theorem below requires only (11); the prospective nonlinear lemma may use the stronger microscopic hierarchy. Time rescaling does not imply Banach analyticity in physical lambda through lambda=0. We hold physical lambda fixed while using the artificial parameter.

## 4. Localization and the first infrared integration

[proved] Translation invariance imposes total momentum conservation modulo 2pi. For two external legs the off-diagonal branch assignment transfers 2p; for four legs, a non-forward assignment transfers 2p or 4p. The total deviation is bounded by respectively 2a_* or 4a_*. Equation (3) therefore excludes these assignments from the compact patches. This argument concerns the external legs of the effective vertex, irrespective of its internal contractions. No Diophantine condition on higher multiples of p is needed for this quadratic/quartic statement.

After number conservation and reflection, the localized action consists of

    sum_omega psi_omega^+ [a + b(-ik0) + c omega q] psi_omega^-
         + g psi_+^+ psi_-^+ psi_-^- psi_+^-.               (12)

The reality symmetry makes a,b,c,g real in this convention. The time coefficient is absorbed into field normalization. The spatial coefficient changes the common speed. The constant a is relevant and must be tuned or absorbed into a shifted p. Off half filling it is not removed by the half-filled symmetry argument in the collaborator's note.

There is only one local quartic because two identical creation or annihilation fields at a point give zero. For degree at least six the engineering dimension 2-p/2 is at most -1. This statement does not require all high-order momentum resonances to be excluded.

[proved] If M_j(K)=integral R_X^j |K|, elementary Taylor remainder bounds give

    |Khat(q)-Khat(0)| <= Q M_1(K),
    |Khat(q)-Khat(0)-D Khat(0)[q]| <= Q^2 M_2(K)/2,          (13)

where Q bounds the independent external momentum/frequency deviations in the dual norm. The inequalities follow by applying |exp(iy)-1|<=|y| and |exp(iy)-1-iy|<=y^2/2 under the kernel integral. They also hold for finite measures. With Q=O(E), a quartic remainder has a gain E, and a quadratic remainder a gain E^2. Relative to their engineering dimensions both are irrelevant by one power.

[proved using C50] Integrate any fixed first shell of (8), including all contractions of every initial degree. Its row moment and Gram bound are finite. C50 gives an analytic connected action with norm at most epsilon/(1-K epsilon), after a fixed radius loss, if K epsilon<1. Equations (10) and (13) then localize that output as in (12). Thus arbitrary-degree retardation is absorbed by the ultraviolet integration, compact crossover and first infrared step. No truncation of initial degree is used.

[limitation] K grows as the shell energy decreases, and this estimate loses field radius. Repeating the same finite-step inequality without a renormalized dimensional bound does not prove an infinite expansion. Sections 6–8 avoid that iteration only for the coefficient linear in the initial action.

## 5. Free shell estimates with a local cancellation

[definition] Work first on the infinite lattice at zero temperature. Use a smooth cutoff chi(q0,q) even under (q0,q)->(-q0,-q), supported sufficiently near zero, and a shell partition f_h=chi(gamma^(-h)q)-chi(gamma^(-h+1)q), h<=0, gamma>1. Fixed velocity factors can be included in chi. Rescale the fixed crossover energy to one. Put E_h=gamma^h and

    D_omega(q)=-iq0+omega v_p q,
    e_p(omega p+q)=omega v_p q + r_omega,p(q),
    |r_omega,p(q)|<=C q^2,
    g_(h,omega)(z)=integral e^(iq.z) f_h(q)
                       /[D_omega(q)+r_omega,p(q)] d^2q/(2pi)^2.   (14)

The deviation q is within one nonwrapping branch chart. Constants are uniform in p in I. The cutoff can be chosen small enough that |r|<=|D|/2 on every shell. For spatial estimates the Fourier integral on the chart is restricted back to integer x.

[proved, shell lemma] For each fixed number of derivatives and moments,

    ||g_h||_infinity <= C E_h,
    ||g_h||_L1 <= C E_h^(-1),
    ||partial_z g_h||_infinity <= C E_h^2,
    |g_h(z)| <= C E_h^2(1+|z|).                            (15)

The L1 estimate includes the sum over lattice x. The first three bounds follow by q=E_h u, with |D+r| bounded below by c E_h on the annulus; integration by parts gives a decaying envelope C_M E_h(1+E_h|z|)^(-M). For the last bound, introduce g_h^lin with denominator D alone. Since f_h is inversion-even and 1/D is inversion-odd,

    g_h^lin(0)=0.

Moreover

    1/(D+r)-1/D = -r/[D(D+r)].                             (16)

The absolute multiplier on the right is O(1) on a shell of measure O(E_h^2). Hence |g_h(0)|<=C E_h^2. Integrating the derivative bound along a straight line from 0 to z proves the last estimate in (15). A continuous Fourier interpolation in x justifies this line argument before restricting to integer x. The physical Fermi phase exp(i omega p x) has modulus one, so it does not change this bound.

[proved, cumulative Gram bound] Let C_h=sum_(k=h+1)^0 g_k. Its Fourier absolute integral is bounded uniformly in h, since each shell costs C E_k and sum_(k<=0)E_k<infinity. The same holds for

    C_h+s g_h,       0<=s<=1.                              (17)

Square-root factorization gives Gram minors bounded by D0^r with a fixed D0. This uses an absolute Fourier integral of the cumulative covariance, not its L1 norm in position space. The latter is not uniformly bounded as h tends to minus infinity.

[scope] The strengthened local cancellation in (15) is proved in the thermodynamic zero-temperature shell representation. Arbitrary finite grids need not be inversion-symmetric about p. We do not assert this cancellation uniformly for their last shell. At fixed h the thermodynamic/zero-temperature limit of the smooth shell exists, so the result below applies to the coefficientwise bulk limit followed by h->-infinity. Uniform nonlinear finite-volume estimates remain part of Section 9.

## 6. The all-scale one-endpoint theorem

[definition] For the compact-patch action V define

    F_h(tV)=log E_(C_h) exp[t V(psi+zeta)],
    A_h V = (d/dt) F_h(tV)|_(t=0) = E_(C_h) V(psi+zeta).     (18)

Here the parameter t multiplies the entire action. It is unrelated to the hopping t in (1). A_h is the derivative at the Gaussian theory, not the derivative about an interacting XXZ background. Every initial degree and every possible number of contractions inside that one vertex is retained.

[proved, theorem L] Fix s>=0 and an output field radius rho. There is a finite input radius H and a constant C, independent of h<=0, such that

    ||A_h V||_(rho,s) <= C ||V||_(H,s),
    ||(A_(h-1)-A_h)V||_(rho,s)
                  <= C E_h^2 ||V||_(H,s+1).              (19)

Consequently A_h V converges in the weaker norm, and

    ||(A_(-infinity)-A_h)V||_(rho,s)
             <= C ||V||_(H,s+1) E_h^2/(1-gamma^(-2)).     (20)

The same estimates apply after projecting to a fixed degree, and to its Fourier derivatives of order at most s. The constants include the finite branch and creation/annihilation label sums. This theorem requires a sufficiently large but finite field radius once, not once per shell.

Proof. Start with finite degree and compact relative-coordinate support. In the Wick average of a degree-p monomial choose its retained legs. Number conservation makes the contracted fields a determinant of size r<=p/2, rather than a factorial sum over pairings. There are at most 2^p external-subset choices; harmless factors for labels and ordering can be bounded by another fixed constant to the pth power. The remaining determinant is bounded by D0^r from Section 5. With the 1/p! coefficient convention of (11), retaining these subsets and then regrouping into antisymmetric output kernels introduces no additional p! factor. It is also possible to leave each subset as an explicit output coefficient representation.

This bounds the first line of (19) by a sum of

    p^s [C0 max(1,rho,sqrt(D0))]^p/p!
                                      integral (1+R_X)^s |K_p|.  (21)

The polynomial p^s allows for changing the anchor to a retained leg. Indeed each new distance is at most the sum of two old anchored distances, so the new sum is at most p R_X. Choose H with a fixed geometric reserve over the bracket in (21). The degree sum is then bounded by the input norm. No covariance position integral is taken in this estimate: all contracted coordinates are already integrated inside the initial kernel.

To bound a shell increment, interpolate its contracted matrix as M(s)=M(C_h+s g_h). The exact determinant identity is

    det M(1)-det M(0)
        = integral_0^1 sum_(i,j=1)^r (g_h)_ij Cof_ij M(s) ds.     (22)

Each cofactor is an (r-1)-row Gram minor bounded by D0^(r-1). There are at most p^2 terms. The marked shell line has endpoints among the original vertex coordinates. Equation (15) bounds it by

    C E_h^2(1+r_ij) <= C E_h^2(1+R_X).                       (23)

For r=0 the increment vanishes. For r>=1, multiply (21) by Cp^2 E_h^2(1+R_X), increasing its geometric reserve if necessary. Sum over p. This proves the second line of (19). Equation (17) makes the cofactor estimate valid throughout the interpolation. The sum over h converges geometrically and proves (20).

Finally remove the degree and coordinate truncations. Equations (21)–(23) give an absolutely summable dominating majorant in the input norm, so the truncated operators converge uniformly and the estimates pass to the limit. Finite measures are allowed because the compact-frequency covariance and all functions used here are bounded and continuous at coincidences. This proves the theorem.

[proved, why this is stronger than finite-scale integration] The result sums the infinitely many contractions of one high-degree endpoint at all infrared scales without multiplying the finite-step radius losses. It also proves an extra E_h^2 gain, rather than just the E_h pointwise shell bound. It does not sum graphs with two or more initial interaction vertices.

## 7. Linearized flow and counterterms

[proved] Apply theorem L with s=2. Extract from -A_h V the inverse-action coefficients in (12), denoting them a_h,b_h,c_h,g_h^loc. Then

    |a_h-a_infty|+|b_h-b_infty|+|c_h-c_infty|
                      +|g_h^loc-g_infty^loc| <= C epsilon E_h^2, (24)

where epsilon=||V||_(H,3). Every limit is O(epsilon). Choosing the relevant bare counterterm -a_infty cancels the limiting constant at this order. Its remaining dimensionless coefficient obeys

    |nu_h^lin| = |a_h-a_infty|/E_h <= C epsilon E_h.          (25)

The linear marginal increments are summable and bounded by C epsilon E_h^2. The time and velocity shifts are finite. Their linear tail estimates use one coordinate moment from the Fourier derivative and one from the marked shell line; the norm in (24) is more than sufficient.

The same theorem with s=2 bounds the Taylor remainder of the output quadratic kernel uniformly, and with s=1 bounds that of the quartic. Hence (13) continues to hold uniformly in h at this linearized level. Higher-degree coefficients obey the uniform field-radius bound in (19). Thus the free dimensional classification is stable under the entire one-endpoint contraction sum.

[proved, parameter continuity for the localized coefficients] Consider kernels fixed in physical coordinates, or an initial family continuous in an appropriately dominated pinned-kernel norm. In a one-endpoint contraction, internal lines conserve branch labels. A diagonal quadratic or forward quartic output is branch neutral, so only branch-neutral input assignments contribute. Their Fermi phases depend on relative coordinates after pinning the output anchor, with no absolute-anchor phase. Differences of those phases are bounded by a constant times |p-p'| R_X. Dominated convergence therefore applies with one additional relative-coordinate moment. For each fixed h, the shell integrals and determinants in (18) are continuous in p. The convergence estimate (20) is uniform on I, so a_infty(p), b_infty(p), c_infty(p) and g_infty^loc(p) are continuous. Reality and reflection make the quadratic constant real and the same at the two roots. This does not assert continuity in an anchor-supremum norm for a branch-charged initial endpoint. The collaborator's FERMI_SHIFT_EXISTENCE.md, Section 6, gives the corresponding neutral-root argument for finite multi-vertex trees.

[proved, first-order Fermi root at fixed physical field] Let p0 be the free root at the chosen physical field. Decompose its free inverse propagator as

    e_p0(k)=e_p(k)+d(p),       d(p)=2(cos p0-cos p).          (26)

At linear order the root equation is d(p)+a_infty(p)=0. If |a_infty(p)|<=C epsilon and the interval [p0-r,p0+r] is contained in I, set r=2C epsilon/v_min. Since d'(p)=2 sin p>=v_min, d(p0-r)<=-2C epsilon and d(p0+r)>=2C epsilon. The continuous function in the root equation changes sign. It has a zero p_* and

    |p_*-p0| <= 2C epsilon/v_min.                           (27)

This constructs the counterterm in the linearized theory at fixed physical field. Inserting that O(epsilon) counterterm inside nonlinear graphs changes terms of higher order, which have not been bounded by theorem L. Equation (27) is not a proof of the full physical Fermi point.

## 8. What the linear theorem does not imply

[proved limitation of this argument] For one endpoint every marked covariance line joins two coordinates controlled by the same initial kernel. This is the reason (23) costs only one fixed kernel moment. With two vertices, their relative center is instead integrated by contraction lines. A bound E_h^2(1+|z|) on a line at separation |z| of order E_h^(-1) is only of order E_h. It cannot be integrated by the initial kernel moment. Renormalized tree power counting and localization are needed to recover the gain in that situation.

For example, naively integrating an infrared shell line costs ||g_h||_1=O(E_h^(-1)). Pairing this with only its O(E_h) Gram scale gives no small factor for a marginal cluster. Repeating the estimate over scales gives an unsummable constant. The cancellation in the reference beta function is therefore indispensable. A proof that every microscopic correction is irrelevant cannot be substituted for it.

[proved] A high-degree endpoint can shift a quartic at first order in its coefficient by contracting one pair. This remains true in one channel when its legs are at distinct space-time points. Consequently the matched forward coupling can differ by O(||R||), even if a pure quartic microscopic model has a higher-order matching correction. Theorem L includes this finite shift. It establishes vanishing linear shell increments, not a zero finite matching shift.

## 9. The exact remaining nonlinear lemma

[assumed, lemma N, not proved] The following would close the requested RG extension. Initial kernels are those obtained in Sections 2–4, in a fixed polynomial tree norm w_q and sufficiently large fixed field radius, with the microscopic degree hierarchy of Section 3. Include sources for the 2p density bilinear and a nearest-neighbor pair. Constants may depend on I, q and the cell, but not on infrared depth or the regulators.

(N1) At every finite depth the localized Gallavotti–Nicolò expansion is absolutely convergent whenever the prescribed running couplings are sufficiently small. The single-step wavefunction normalization is solved continuously on that history domain and obeys |log(Z_(h-1)/Z_h)|<=C epsilon, so its ratios stay in a fixed neighborhood of one without assuming convergence of Z_h. Its sum and marked sums have regulator-independent bounds. In particular, no extra field-radius or endpoint-moment loss accumulates with the number of scales. The final finite-volume shell is included, or a precise cofinal thermodynamic limiting prescription is supplied.

(N2) Trees containing an initial irrelevant endpoint, or a lattice-versus-linear covariance remainder, have normalized marginal contribution at scale h bounded by A epsilon gamma^(theta h), for some theta>0. Here epsilon bounds all initial perturbations, so the bound permits contributions linear in R. A history insertion at scale k>h has an additional factor C epsilon gamma^(-kappa(k-h)), with kappa>theta. Equivalently, after using the reference constant-history cancellation, the actual marginal beta obeys

    |b_h^mu| <= A epsilon gamma^(theta h)
       + B epsilon sum_(k=h+1)^0 gamma^(-kappa(k-h)) |mu_k-mu_h|
       + B epsilon sum_(k=h)^0 gamma^(-kappa(k-h)) |nu_k|.     (28)

The relevant beta, in nu_(h-1)=gamma nu_h+b_h^nu, obeys

    |b_h^nu| <= A epsilon gamma^(theta h)
       + B epsilon sum_(k=h)^0 gamma^(-kappa(k-h)) |nu_k|.     (29)

The constants are uniform on the indicated history domain. The beta functions are continuous there. Their finite-depth definitions agree on common higher-scale histories. For the Fermi-point map one needs either a continuous choice of the counterterm solution or the coupled Brouwer construction of Section 10. No differentiability in p is required by that construction.

(N3) The same renormalized bounds with one and two marked density or pair endpoints give a summable normalized difference from a matched one-channel reference theory. In a normalization with matched marginal limit, require

    |log(Z_(O,h-1)/Z_(O,h))
                   -log(Z_(O,h-1)^ref/Z_(O,h)^ref)|
                                <= C epsilon gamma^(alpha h)   (30)

for O=fermion, 2p density, pair, with some alpha>0. Nonzero source matching at the initial scale must persist. The following is one explicit sufficient two-insertion bound. After factoring the density oscillation and matching the limiting source amplitude, write the correlation difference as a finite-crossover term plus sum_(h<=0) R_(O,h)(z), and require

    |R_(O,h)(z)| <= C epsilon E_h^(2X_O+alpha)
                                      /[1+(E_h |z|)^M],
    M>max_O(2X_O+alpha),                                  (30a)

where X_O is the matched reference dimension and the crossover term decays faster than |z|^(-M). This suffices for an O(epsilon |z|^(-2X_O-alpha)) remainder. An equivalent marked-tree estimate implying (30a) would do. The coefficient bound, its scale sum and the nonzero limiting source overlap are all required; (30) alone does not control a position-space correlation.

[definition] mu in (28) comprises the forward coupling and the velocity parameters in a scheme where their linear-reference constant histories have no unsuppressed beta. Wavefunction normalization itself is not a marginal coupling with a summable beta: it has an anomalous power and is treated through its normalized ratios and (30). A formulation matching the limiting speed first and making the residual velocity tend to zero is equivalent, but must establish that matching as part of its history construction.

[proved, reduction to a concrete endpoint problem] Pure localized forward trees and their reference Ward identities are unchanged by the origin of the action. The extra trees in N2 must touch a general initial endpoint or a covariance correction. Thus the new proof obligation is the uniform renormalized bound for those trees and their marks. A sufficient technical implementation is a localization routing expansion with at most C^p terms per p-leg endpoint and at most b coordinate factors on each edge of a chosen endpoint tree, with b independent of depth and p. Then w_q with q>=b pays the endpoint factors by the collaborator's Section 5.3 inequality. The radius reserve pays C^p. The missing point is proving that routing and its compatibility with nested subtractions, determinant derivatives, history marks and source marks. The inequality for a tree with bounded loads does not itself establish that the loads remain bounded in the RG construction.

[scope] N1–N3 are not asserted as an abstract theorem in BM, C50 or the collaborator's note. They are the remaining lemma for this action. Theorem L verifies the one-initial-vertex, no-other-interaction part of N2, with the stronger E_h^2 shell estimate in the bulk limit. Sections 2–4 verify its finite initialization. They do not verify the multi-vertex part of N1–N3.

## 10. A continuity/Brouwer closure, conditional on N1–N2

[imported context] The collaborator's [COUNTERTERM_INDUCTION.md](../../contrib/direct_attack_20261028/vacancy/COUNTERTERM_INDUCTION.md), [FERMI_SHIFT_EXISTENCE.md](../../contrib/direct_attack_20261028/vacancy/FERMI_SHIFT_EXISTENCE.md), and [GN_MARK_CALCULUS.md](../../contrib/direct_attack_20261028/vacancy/GN_MARK_CALCULUS.md) already prove a conditional contraction construction, finite-cutoff physical matching and diagonal compactness. They remain the fuller account of history-difference estimates. The variant below replaces contraction by a joint finite-dimensional self-map. It proves existence under (28)–(29) and continuity, without a history Lipschitz estimate or a unique counterterm branch. It does not remove the absolute renormalized bounds needed to define that map.

[proved, abstract sequence theorem] This section separates the existence argument from the missing RG estimate. Assume continuous beta functions satisfying (28)–(29), causal in the history from h to 0. Work in a domain with |mu_0|<=c0 epsilon, and let 0<alpha<min(theta,kappa). Define

    H_alpha = 1/[1-gamma^(-(kappa-alpha))],
    D_alpha = gamma-gamma^(-alpha),
    C_alpha = 1/[(1-gamma^(-kappa))(gamma^(kappa-alpha)-1)].   (31)

The last constant sums the double history series. Use coordinates

    d_h=mu_(h-1)-mu_h,
    |d_h|<=R_mu epsilon gamma^(alpha h),
    |nu_h|<=R_nu epsilon gamma^(alpha h).                    (32)

For finite indices -N+1<=h<=0 set nu_(-N)=0. Construct mu_h from mu_0 and the d's at strictly higher scales. Define a continuous map by

    (T d)_h = b_h^mu(mu,nu),
    (T nu)_h = -sum_(k=-N+1)^h gamma^(k-h-1) b_k^nu(mu,nu).  (33)

The nu formula is the solution of its relevant recurrence with the stated terminal boundary, for the beta values evaluated on the proposed history.

From (32), for k=h+n,

    |mu_k-mu_h| <= R_mu epsilon gamma^(alpha h)
                                      sum_(r=1)^n gamma^(alpha r).

Multiply by gamma^(-kappa n) and sum n>=1. Interchanging the nonnegative sums gives

    sum_(n>=1) gamma^(-kappa n) sum_(r=1)^n gamma^(alpha r)
                 = C_alpha.                               (34)

Similarly the nu history sum is at most H_alpha R_nu epsilon gamma^(alpha h). Equations (28)–(29) imply

    |(Td)_h|/gamma^(alpha h)
       <= A epsilon+B epsilon^2(C_alpha R_mu+H_alpha R_nu),
    |(Tnu)_h|/gamma^(alpha h)
       <= [A epsilon+B H_alpha R_nu epsilon^2]/D_alpha.    (35)

The second denominator follows by summing

    sum_(k<=h) gamma^(k-h-1) gamma^(alpha k)
                         = gamma^(alpha h)/D_alpha.        (36)

Take R_mu=2A and R_nu=2A/D_alpha. Choose epsilon small enough that

    B epsilon(C_alpha R_mu+H_alpha R_nu)<=A,
    B epsilon H_alpha R_nu<=A.                             (37)

Also require that c0 epsilon+R_mu epsilon/(1-gamma^(-alpha)) and R_nu epsilon lie inside the small-history domain of N1. The product of balls (32) is compact and convex in finite dimensions, and (35) maps it into itself. Brouwer gives a fixed point. At that point both actual recurrences hold, with the prescribed mu_0 and terminal nu boundary.

Let N tend to infinity. Successively extract subsequences for every fixed index h. Bounds (32) give a diagonal subsequence. Each fixed beta function depends only on finitely many indices h,...,0, so its continuity passes the recurrences to the limit. The relevant formula (33) also passes to its infinite sum, since its omitted tail is bounded by a geometric series using (29), (32). Thus

    nu_h=-sum_(k=-infinity)^h gamma^(k-h-1)b_k^nu,
    |nu_h|<=R_nu epsilon gamma^(alpha h),
    |mu_h-mu_infty|<=R_mu epsilon gamma^(alpha h)
                                           /(1-gamma^(-alpha)). (38)

This proves existence and summable marginal flow under N1–N2. It does not prove uniqueness or a continuous selection of the fixed point as a function of p.

[proved, coupled physical-root construction] A continuous selection is unnecessary for existence at fixed physical field. Include p as one more Brouwer coordinate, with |p-p0|<=r. Require the exact physical matching hypothesis, as in the collaborator's H3,

    nu_0^phys(p)=A(p)d(p)+u(p),
    s(p,nu_0)=A(p)^(-1)[u(p)-nu_0],                         (38a)

where A,u are continuous, A and its inverse are uniformly bounded, and |u|<=C_u epsilon. At each p all other initial physical kernels are held fixed while only the localized relevant coefficient is varied. This retains their full dependence on d(p) and avoids modifying the physical Hamiltonian. Here d(p) is (26), and the relevant ball bounds |s|<=C_s epsilon. The equation d(p)+s(p,nu_0)=0 is precisely nu_0^phys(p)=nu_0. Replace the p component of the map by

    T_p = d^(-1)(-s(p,nu_0)).                              (39)

Here d is strictly increasing on [p0-r,p0+r], d(p0)=0. Choose r=2C_s epsilon/v_min and keep that interval in I. Its image contains [-2C_s epsilon,2C_s epsilon], so (39) is defined and stays inside the interval. Apply Brouwer to (33),(39) jointly, allowing mu_0 to be the continuous physical initial marginal vector mu_0(p). The same subsequence argument retains a limiting p and the matching equation because (38a) uses only p and the initial counterterm coordinate. The root displacement is at most r. Verifying (38a) and the continuity of the actual finite-depth beta functions is part of the application, not a consequence of the inequalities alone.

This proves that continuity and uniform beta bounds suffice. It does not assume that the set of fixed points has a continuous branch. It also shows why a finite-order self-energy, or a finite-volume fixed point without uniform estimates, is insufficient.

## 11. Density and pairing exponents, conditional on N3

[read, published input] BM, [arXiv:0907.2837](https://arxiv.org/html/0907.2837), Section 3, equations (63)–(65), starts its infrared construction from all-degree ultraviolet kernels. Equations (76)–(82) separate reference and lattice beta functions. Equations (83)–(95) perform source matching. Theorem 1.1 and equations (52)–(53) give the density/pair dimension relation. This is a theorem for the paper's microscopic interaction, not a stated arbitrary-retarded-kernel theorem. Its reference theory and identities are imported here, but the extension of the marked-tree bounds is N3.

[proved implication] If N1–N3 hold, match the limiting forward coupling and speed to the same reference model. The analytic reference matching has nonzero derivative at the free point and hence is locally invertible. Summing (30) gives a convergent product for each source ratio. If its initial ratio is nonzero, the limit is nonzero, since the sum of absolute logarithms of successive ratios is finite. The two-insertion remainder in N3 then transfers the reference power law without changing its exponent.

To check the last step explicitly, put beta_O=2X_O+alpha>0 and split the sum in (30a) at E_h |z|=1. The low-energy part is bounded by C epsilon sum_(E_h<=1/|z|)E_h^beta_O<=C' epsilon |z|^(-beta_O). The high-energy part is bounded by C epsilon |z|^(-M) sum_(E_h>1/|z|)E_h^(beta_O-M). Because M>beta_O, this sum is controlled by its lowest-energy term and has the same bound. The finite-crossover term is smaller. This proves the stated improved decay, uniformly in direction in the anisotropic space-time norm.

In the reference normalization put tau=g_ref/(4pi v). Its dimensions are

    X_density=(1-tau)/(1+tau)=K,
    X_pair=(1+tau)/(1-tau)=1/K.                            (40)

Thus the oscillating density correlation has a term proportional to cos(2p_* x) r^(-2K), and a local spinless pair correlator has leading power r^(-2/K), where r^2=x^2+v^2 tau_E^2. The exponents K and 1/K are scaling dimensions; the two-point decay powers contain the factor two. The smooth density component has its separate dimension-one current contribution. Formula (40) uses the matched reference coupling, not the bare XXZ coefficient or the physical port lambda.

[proved, source overlap check at the free point] Write f_x approximately as exp(ipx)psi_+ + exp(-ipx)psi_-. The density f_x^dagger f_x has cross-branch bilinears with unit coefficients at momenta +/-2p. For the nearest-neighbor pair,

    f_x f_(x+1) = [exp(-ip)-exp(ip)]psi_+ psi_- + derivatives
                = -2i sin p psi_+ psi_- + derivatives.      (41)

The coefficient is bounded away from zero on I. Finite analytic UV matching preserves both overlaps at sufficiently small perturbation. Their survival through infinitely many scales is the source-product assertion of N3, not a consequence of (41) alone. For dressed microscopic observables, the source action produced by massive matching must also satisfy N3 and have the corresponding nonzero overlap.

[limitation] A generic retarded U(1)-invariant action need not have BM's bare instantaneous lattice current identity. Local gauge variation can insert time-nonlocal current terms. The exponent implication above uses matching to the reference source problem, not an unjustified Ward identity for the bare retarded action. No susceptibility or Drude formula for an unmatched physical density is inferred here.

## 12. Outcome and checks

[proved] Retarded all-degree kernels can be initialized and localized off half filling. Their complete Gaussian linear response has a uniform all-scale contraction sum, summable marginal shifts and a continuous first-order counterterm. The nonlinear extension remains N1–N3, specifically the uniform renormalized endpoint and marked-endpoint estimate with no cumulative norm loss. The Brouwer argument is proved once those estimates and their continuous initialization are supplied.

[checks] The free limit of (24)–(27) has zero shift, and (40) gives K=1. Equation (41) is nonzero precisely away from the band edges covered by (3). The marked-line estimate uses the correct shell amplitude E_h in two dimensions and improves it to E_h^2 at fixed separation by inversion symmetry; it does not assert this gain at separations of order E_h^(-1). The counterterm is divided by one energy power in (25), consistent with its relevant dimension. The half-filled and band-edge limits are outside the uniform constants rather than claimed as consequences.

[numerical] `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B checks.py` passes. The smooth-shell quadrature tests two fillings, three separations and ten energies; |g_h(z)|/E_h^2 approaches finite values and the odd linear origin residual is below 1e-18. Determinant interpolation residual is 5.56e-17. The history-sum identity agrees to 8.89e-16 and the finite relevant recurrence to 8.68e-19. A synthetic continuous root lies within the proved displacement bracket, with residual 1.59e-17. These are controls of the displayed lemmas, not a computation of the microscopic RG flow. Details are in [checks.py](checks.py) and [checks.json](checks.json).

[scope] All writes are confined to work_C55. Parent project memory is unchanged under the task's explicit write restriction. This note has no independent referee review. It does not establish the first rigorous Luttinger theorem through ports.
