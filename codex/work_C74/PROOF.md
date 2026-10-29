# C74. Taylor projectors across contractions

## 1. Status and inputs

[assumed, unresolved target] This is a partial result. I do not prove C60 (17) for a general GN tree, do not complete corrected Lemma N, and do not exhibit a contraction pattern that forces unbounded counters under every reorganization. The first merger below defeats commuting the old projector through the contraction unchanged. It admits a bounded repair, so it is not the requested obstruction.

[proved below] There are three results. Section 2 gives an all-depth field-derivative induction through Gaussian contractions, without a bound on coefficient moments. Section 3 gives the exact extra terms when a Taylor projector is moved through a contraction. Sections 4 through 7 prove the simultaneous spatial and scale-cut bound for arbitrary-depth quartic pair ladders. The ladder proof retains the BBF determinants and gives C60 (17) on that family.

[assumed, imported] Use C60's fixed real cutoffs, free lattice dispersion, exact quadratic absorption and order of regulator limits. Use C61's endpoint norm, source subtraction and oriented incidence convention. The fixed-order covariance and Gram estimates are those of MARKED_COVARIANCE. TREE_NORM_EXTENSION supplies the terminal-tree norm. These are inputs, not consequences of the present calculations. No matching, reference-cancellation, source or thermal assertion is reproved.

[definition] Larger scale indices mean shorter distances. Write E_s=2^s, with initial scale zero. R_m subtracts the total Taylor jet of degree less than m at zero deviation momenta on the conservation hyperplane. The orders are m=2 for quadratic and m=1 for quartic kernels. R_0 is the identity. Engineering derivatives are accounted for separately, as in C60.

[definition] All displayed covariance norm estimates use the normalized shell covariances. Explicit products of field-normalization ratios are factored out, as in the localization-factor question C60 (17). Their separate bounds are imported from MARKED_COVARIANCE, not included in the bare factor 2^h below.

## 2. Field derivatives stay bounded through general contractions

[definition] At a fixed finite tree and fixed infrared cutoff, retain each term as

    P(p_ext) F(p_ext),
    P = product_f product_mu p_(f,mu)^(a_(f,mu)).                 (1)

The momenta in P keep their original field labels. Conservation is imposed on F and P together. Do not replace every inherited pivot momentum in P by a sum of the other momenta. Let a_f=sum_mu a_(f,mu) and a=sum_f a_f. F includes all coefficient integrals, Taylor parameters and Gaussian contractions below that vertex.

[proved, contraction-compatible derivative induction] Starting with a_f=0 for localization derivatives, there is an exact expansion at every finite depth with

    a_f <= 2                                                   (2)

on every original field until its first contraction or the root. A tree covariance has at most four localization derivatives, two from either end. A determinant row or column has at most two. This assertion does not bound moments or derivatives of the original coefficient measures hidden in F.

Proof. Consider one fixed BBF contraction and one monomial term from each child. Separate the factors P into those on contracted and uncontracted fields. Factors on uncontracted fields pull outside the Gaussian expectation. Factors on contracted fields differentiate their covariance at the indicated field position. A tree field occurs on its selected tree line and is absent from the residual determinant. A remaining contracted field occurs in one determinant row or column. Thus this operation does not increase a_f and does not expand the determinant into pairings.

Call the output C and the surviving monomial P_ext, of total degree a. If the output has more than four fields, no Taylor operation is added. Otherwise homogeneity gives

    R_m(P_ext C) = P_ext R_(m-a)_+ C.                           (3)

If a>=m, nothing is added. If a<m, put r=m-a. Taylor's formula in independent external components gives

    R_r C(q) = sum_(|gamma|=r) (r/gamma!) q^gamma
                 integral_0^1 (1-t)^(r-1) D^gamma C(tq) dt.     (4)

In this case the total degree after the operation is a+r=m<=2. Consequently every field degree is at most two. If a>=m, all field degrees retain their previous bounds. This proves (2) inductively over contractions and localizations. Once contracted, an original field cannot acquire a new derivative at a later vertex.

There are at most six independent scalar components at a quartic vertex in two spacetime dimensions, and two at a quadratic vertex. The coefficient mass in (4), including the Taylor integral, is d^r/r!, with d<=6 and r<=2. It is uniformly bounded. Unary repetitions with no intervening contraction or merge reduce by idempotency. Each remaining merge is charged to a branch vertex, and each unary contraction removes fields. Their number is O(n+sum_i p_i). Thus the numerical coefficient mass of this particular factoring and Taylor procedure is at most C^(n+sum_i p_i). This count leaves the differentiated coefficients D^gamma C unexpanded. It is not the coefficient-mass bound for their eventual spatial routing.

[proved conditional on the covariance input] For a BBF matrix written as

    G_ij = <u_i tensor A_(x_i), u_j tensor B_(y_j)>,

derivatives on field i or j replace the corresponding Fourier vector by its momentum-multiplied vector. With row orders a_i<=2 and column orders b_j<=2, Gram's inequality gives

    |det G_decorated| <= C^N E_s^(N+sum_i a_i+sum_j b_j).        (5)

The intersection of a differentiated row and column has both derivatives. No cofactor expansion is needed for these field derivatives. The tree-line estimate is

    ||partial^(a+b) g_s||_1 <= C E_s^(a+b-1), a,b<=2.           (6)

Affine evaluations of positions at fixed Taylor parameters do not change the Fourier vector norms. These observations prove the stated Gaussian derivative bounds.

[proved, limitation] Equation (4) differentiates C, which can already contain derivatives of an endpoint from a lower Taylor operation. The induction counts exposed momenta on fields, not those coefficient derivatives. Taking (2) as a bound on endpoint moments would repeat the mistake identified in C72 Section 5. The proof also makes no assertion that coordinate moments created by D^gamma C are already on actual shifted tree lines.

[assumed, regularity boundary] The algebra can first be performed with compactly supported endpoint relative coordinates and fixed regulators. For each finite tree all operations are defined. Extending the resulting representation uniformly to the prescribed endpoint norm requires the missing moment bound. In particular, no arbitrary-order moment assumption is imposed on the fixed conversion legs to bypass that issue.

## 3. The extra jets produced by a contraction

[definition] Consider a contraction in independent variables c and q,

    (C_W F)(q) = integral dc W(c,q) F(c,q).                      (7)

c denotes integrated variables, and q denotes the surviving external variables. W includes the other child coefficients, covariance factors and the unexpanded residual determinant. The following identities hold whenever differentiation through this finite regulated integral is justified. They extend to several c variables without change. They are not assertions that every momentum routing has the spatial anchor property.

[proved, contraction formula] For any m>=1,

    R_m^q C_W F
      = C_W R_m^q F
        + sum_(|alpha|<m) q^alpha/alpha!
            R_(m-|alpha|)^q
              integral dc W(c,q) D_q^alpha F(c,0).             (8)

To prove it, decompose F=R_m^q F+T_(m-1)^q F. All q derivatives of the first summand of order less than m vanish at q=0. Leibniz's rule shows that multiplication by W and integration preserve these vanishing jets. R_m therefore acts as the identity on C_W R_m^q F. In the finite jet sum, multiplication by q^alpha lowers the remaining subtraction order by |alpha|, which gives (8).

[proved, scalar cases] For one external component, write F_j(c)=partial_q^j F(c,0), and let C_W(q) act on functions of c. Then (8) says

    R_1 C_W F = C_W R_1^q F + [C_W(q)-C_W(0)]F_0,             (9)

    R_2 C_W F = C_W R_2^q F
       + [C_W(q)-C_W(0)-q C_W'(0)]F_0
       + q[C_W(q)-C_W(0)]F_1.                              (10)

For m=1 the additional term vanishes only if the contraction is independent of q on F_0, or if that particular jet is annihilated. Neither condition follows from the existence of a Gram bound.

[proved, quartic child followed by quartic parent] Put F=R_1^(c,q) K=K(c,q)-K(0,0). Equation (9) becomes

    R_1^q C_W R_1^(c,q) K
      = integral dc W(c,q)[K(c,q)-K(c,0)]
        + integral dc [W(c,q)-W(c,0)][K(c,0)-K(0,0)].          (11)

This is a contraction version of the survival split. The first term keeps a surviving q factor. The second uses a contracted c factor and a new q difference of the contraction. The old projector cannot in general be collapsed across C_W. Section 5 realizes both terms with one spatial representation and separate cut labels.

[proved, compatibility with incoming derivatives at one event] Applying at most two q derivatives to (8) and using D^beta R_m=R_(m-|beta|)_+ D^beta gives a finite exact formula. For m<=2 it uses at most four derivatives of F and W at this one event. This bound is relative to the current F and W. Repeating it on coefficients already differentiated by earlier events gives no depth-independent order on the original endpoints. In addition, expanding derivatives of W may differentiate another child coefficient. Equation (8) therefore adds data to C72's retained flag instead of proving that its state closes unchanged.

## 4. A legal first merger with a nonzero regenerated contact

[definition] Use the two branches R,L of the single-channel action. Take the retarded quartic endpoint A and a local forward vertex B with fields

    A: psi_R^+(a+u), psi_L^+(a+u), psi_L^-(a+r), psi_R^-(a),
    B: psi_R^+(b),   psi_L^+(b),   psi_L^-(b),   psi_R^-(b).    (12)

Here r and u may be taken along Euclidean time. The coefficient is a finite measure kappa(dr,du), with finite exponential terminal-tree norm. Its amplitude can be multiplied by the required epsilon^2. This is a subclass of the admitted quartic endpoint space, not an assertion that the microscopic matching produces a specified coefficient in it. Charge is conserved, and each coincident pair has distinct branch labels, so Pauli antisymmetry does not kill (12).

Contract the two annihilation fields of A with the two creation fields of B at scale s. Retain the two creation fields of A and the two annihilation fields of B. There are four external fields because 4+4-2(2)=4. Choose the R contraction as the BBF tree edge. The residual determinant is the one-by-one matrix

    det[lambda g_(L,s)(a+r-b)],  0<=lambda<=1.                 (13)

The tree line is g_(R,s)(a-b). Its BBF parameter integral is 1/2 times the product of the two lines. The other tree choice gives the other half with the same fixed overall Grassmann sign. Work with the stated labeled tree choice and absorb its factor 1/2 into C_s. No determinant has been replaced by an absolute Wick sum.

[definition] Set z=a-b and define

    K(c,Q) = integral kappa(dr,du) exp(i c dot r-i Q dot u),
    W_s(c,Q) = (1/2) ghat_(R,s)(Q-c) ghat_(L,s)(c).           (14)

Q is the total deviation momentum of the surviving creation pair, with the opposite momentum on the annihilation pair. The direct Fourier transform of g_R(z)g_L(z+r), with external displacement z+u, gives (14). The factors (2 pi)^(-2) are included in dc. Child contact localization subtracts K(0,0). Parent contact localization is Q=0, since the output depends only on Q. Thus its twice-regularized amplitude is the left side of (11).

[proved, nonzero contact before the parent subtraction] Choose u=0 and the temporal measure with masses A/2 at r and -r, where A>0 and r!=0. Then

    K(c,Q)=A cos(r c_0),    K(0,0)=A.

Use the free lattice covariance, fixed reflection-related branch cutoffs and a reflection-symmetric dispersion. In deviation coordinates,

    e_R(-c_1)=e_L(c_1),    f_(R,s)(-c)=f_(L,s)(c).

At Q=0,

    (C_s R_1 K)(0)
      = (A/2) integral dc
          f_(L,s)(c)^2 [cos(r c_0)-1]
                    /[c_0^2+e_L(c_1)^2] < 0.                  (15)

The denominator is positive on the nonzero shell. The numerator is nonpositive and strictly negative on a set of positive Lebesgue measure, since a nontrivial smooth shell contains an open set and r!=0. All integrals are finite. The coefficient measure has finite exponential norm for fixed r, and an arbitrary small prefactor puts it in N0. A symmetry completion can be considered separately; the nonzero assertion is for this labeled component of the absolute tree problem.

[proved, interpretation] Equation (15) disproves the shortcut that a child quartic remainder stays in the parent contact ideal after its carrier ports are contracted. It is the first two-child quartic merger needed to see this particular failure. It does not force a counter to grow. Equation (11) repairs it, and the next section pays the repair with bounded loads. No claim is made about a nonzero projected amplitude for every arbitrarily long external scale span.

## 5. Spatial realization of both terms and their two cuts

[definition] For this calculation write

    A_+(x)=psi_R^+(x) psi_L^+(x),
    B_-(b)=psi_L^-(b) psi_R^-(b).

These are ordered Grassmann products used to display the calculation. A derivative of A_+ is always expanded by the ordinary even Leibniz rule into its two field derivatives. Its use does not change the coefficient convention or introduce a new fundamental field.

[proved, exact position formula] Apart from the fixed BBF factor and Grassmann sign, (11) paired with external fields is the sum

    I_surv = integral kappa(dr,du) dz
       g_R(z) g_L(z+r)
       [A_+(b+z+u)-A_+(b+z)] B_-(b),                         (16)

    I_reset = integral kappa(dr,du) dz
       g_R(z) [g_L(z+r)-g_L(z)]
       [A_+(b+z)-A_+(b)] B_-(b).                            (17)

One can verify this without Fourier manipulations. Before parent subtraction the child remainder is

    g_R(z){g_L(z+r) A_+(b+z+u)-g_L(z) A_+(b+z)}B_-(b).

Parent contact evaluation replaces each A_+ by A_+(b). Subtract it and add and subtract g_R(z)g_L(z+r)A_+(b+z). The result is (16)-(17). This proof keeps the action on all coordinate-dependent exterior fields.

[proved, interpolation and load bounds] In (16), use

    A_+(x+u)-A_+(x)
       = integral_0^1 u dot grad A_+(x+theta u) dtheta.         (18)

There is one endpoint moment, no line moment, no covariance derivative, and one derivative on a root field. In (17), use

    g_L(z+r)-g_L(z)
       = integral_0^1 r dot grad g_L(z+eta r) deta,
    A_+(b+z)-A_+(b)
       = integral_0^1 z dot grad A_+(b+theta z) dtheta.        (19)

There is one endpoint moment r, one moment z on the actual R tree line, one derivative in the L determinant and one derivative on a root field. Both parameter measures have mass one. The spacetime component sums and the two-field Leibniz rule cost fixed constants. No derivative acts only on a covariance while leaving a preattached weight at an old argument.

[proved conditional on the fixed-order bounds] Normalize the root external fields at scale h<s. The derivative in (18) or the second formula in (19) costs E_h. The remaining factors obey

    ||g_R||_1 sup|g_L| <= C E_s^(-1) E_s = C,
    ||z g_R||_1 sup|grad g_L| <= C E_s^(-2) E_s^2 = C.       (20)

The suprema are uniform in z+eta r. They can equivalently be obtained from the one-row derivative Gram estimate for (13). For the fixed metric defining the endpoint norm, |r|+|u| is bounded by a constant times the terminal-tree length of (12). Since t<=exp(t) for t>=0, the endpoint exponential norm pays the one moment in each term with a fixed constant. Thus (16) and (17) are each bounded by

    C ||kappa||_exp E_h                                    (21)

in the quartic coefficient estimate, with its local coupling factors retained. Endpoint scaling is E_0=1.

[definition] The root estimate here and in (26) is the tested, derivative-decorated coefficient estimate used in (17). One can implement it by testing each root field against functions with sup norm at most one and first derivatives bounded by C E_h, after pinning the root anchor. It is not a claim that the untested position-space total variation of F_L tends to zero with an independently chosen h.

[proved, anchor compatibility] Pin b. The actual tree argument is z=a-b, so the a integral is its unshifted anchor integral with unit Jacobian. Only the L determinant is shifted by eta r, and its Gram bound is uniform in that shift. In (16) and (17) the root test field is evaluated at an affine combination of positions within the already formed cluster. It has not yet been contracted. On the lattice the z sum is an integer anchor difference; no fractional substitution of lattice variables is used.

[proved, separate cut records] The endpoint is at zero, the first contraction at s and the root sink at h. In (16) assign the endpoint moment u to a path 0 to h. It crosses both required quartic intervals. In (17) assign r to the derivative in the L determinant at s, and z to the external derivative at h. These are the paths

    0 to s,    s to h.                                     (22)

Their scale factors are respectively E_s/E_0 and E_h/E_s. Each required cut is crossed once, and the product is E_h/E_0, in agreement with (21). All line moments, derivatives, determinant bounds and anchor integration belong to the same displayed terms. This supplies the scale record absent from projector collapse alone.

## 6. Arbitrary-depth quartic pair ladders

[definition] Start with (12), and successively attach local forward vertices B_1,...,B_L at scales

    0 > s_1 > s_2 > ... > s_L > h.

At every attachment contract the current two annihilation fields with the creation fields of the new vertex. Keep the original creation pair and the new annihilation pair. Choose one of the two cross lines as the BBF tree line at each attachment and retain the other as its one-by-one determinant. Apply the quartic contact remainder to each output. Running local vertices carry their given small amplitudes; they are not initial irrelevant endpoints.

[proved, exact depth induction] After the first attachment and its remainder, let the coefficient be F_1(Q). Then F_1(0)=0. Every later attachment acts by multiplication in the conserved pair momentum,

    C_j F(Q) = B_j(Q) F(Q),
    B_j(Q) = (1/2) integral dc
                      ghat_(R,s_j)(Q-c) ghat_(L,s_j)(c).      (23)

Indeed the two fields removed from the previous cluster carry total momentum Q and transfer that same total momentum to the two new external fields. Equivalently the position kernel is convolved with g_(R,s_j)(z)g_(L,s_j)(z)/2. No separate individual external momentum enters a local B_j endpoint. Therefore

    R_1[B_j F_(j-1)] = B_j F_(j-1),
    F_L(Q) = F_1(Q) product_(j=2)^L B_j(Q).                    (24)

The first equality follows from F_(j-1)(0)=0, and proves the assertion by induction. An inherited factor Q can be kept outside every later contraction. This is an induction through genuine Gaussian contractions, not a chain on one unchanged kernel. It works because these particular contractions are multiplication operators on the preserved momentum ideal.

[proved, all-depth spatial bound] Use (16)-(19) only at the first attachment. Put their derivative on the surviving creation pair, and keep it there until scale h. Subsequent attachments remove annihilation fields only. They add no moment or derivative, and no further nontrivial Taylor interpolation. The derivative of the creation pair remains a sum of two single-field derivatives.

Let b_j be the position of B_j. The later tree arguments are z_j=b_(j-1)-b_j. Holding r,u and the first Taylor parameters fixed, the transformation from the anchors to z_1=a-b_1 and z_j, j>=2, is triangular with determinant of modulus one. Each new factor is bounded by

    ||g_(R,s_j)||_1 sup|g_(L,s_j)| <= C.                      (25)

The determinant representation is retained at every step. The root creation fields keep the affine positions introduced at the first step; no lower ladder covariance is incident on them. Consequently the first-contraction anchor condition is preserved. Combining (21) and (25) gives

    ||F_L||_root <= C^L ||kappa||_exp
                        product_j |lambda_j| E_h,            (26)

where lambda_j denotes the prescribed local amplitude, distinct from the BBF interpolation variables. This bound has no constant per unit scale span.

[proved, cut assignment at every depth] For the first term keep the path 0 to h. For the second term keep 0 to s_1 and s_1 to h. The latter path crosses every lower ladder cut. Hence every required spine cut has F_v=1, including cuts whose contact subtraction became redundant in (24). The total gain is

    2^[-(0-s_1)] product_(j=1)^(L-1) 2^[-(s_j-s_(j+1))]
                               2^[-(s_L-h)] = 2^h.           (27)

Equations (26)-(27) prove C60 (17) on this family. The at most two BBF tree choices per attachment, component sums and first Taylor choices cost C^L, which is within C^(sum_i p_i+n). The maximum endpoint moment order is one, tree-line moment order one, covariance derivative order one, and field derivative order one. Increasing L changes none of them.

[proved, scope] The theorem allows arbitrarily large independent spine scale spans and arbitrary finite L. It does not require nonvanishing of the final kernel at all choices of scales, since it is an upper bound. It covers the specified pair-conserving ladder, not all quartic trees. A merger that exchanges different port subsets need not act by multiplication in one conserved pair variable, and (24) then has no stated justification.

## 7. Checks and the unresolved induction

[proved, consistency checks] If kappa is concentrated at r=u=0, both (16) and (17) vanish, as required for a localized initial quartic. If r=0 only the surviving term remains. If u=0 only the reset term remains. The latter is the nonzero-contact example (15). Setting a scale gap to zero in (22) or (27) gives factor one on that gap and leaves the identities unchanged. For L=1, (26)-(27) reduce to (20)-(22). The quartic output degree and U(1) charge are preserved at every ladder attachment.

[proved, what has closed] Section 2 closes the exposed field-derivative counter through arbitrary finite GN contractions. Sections 5-6 close all counters and the incidence record for a contraction family where the surviving momentum ideal is preserved after the first reset. These two statements are logically different, and neither proves the general spatial induction by itself.

[assumed, remaining general assertion] At a general merger, (8) retains jet coefficients for several independent port variables. A subsequent contraction can act on these jets and on their differentiated sibling coefficients. To conclude (17), one still needs a representation in which incoming moments are processed before the necessary cross-component sums are separated, with bounded moments on original endpoint edges and actual tree-line arguments, no uncontrolled coordinate weights on determinant entries, and an oriented record covering every cut. The coefficient bound must apply after those operations, not merely to (4) before its differentiated coefficient is opened.

[proved, why no obstruction follows here] The contact (15) is not a counter-growth theorem: (16)-(22) repair it with fixed counters. The arbitrary-depth ladder shows that repeating this sort of Gaussian attachment need not produce growth. A nonzero contraction commutator alone therefore cannot answer the task's alternative. I have not proved that any remaining pattern defeats every bounded representation.

[assumed, corrected Lemma N boundary] Even a general source-free proof of (17) would still need C61's N-comparison and N-source in the same scheme. In particular the interaction-prefactored source bound concerns D_O(u)-D_O(0), while D_O(0) retains its free lattice descendants. The ladder and derivative theorems prove neither of those comparison hypotheses.

[numerical, exact-arithmetic controls] checks.py verifies (8)-(11) on finite polynomials with an explicit integration functional, checks (16)-(17) before interpolation, tests the derivative counter transition including contraction and repartition of fields, verifies repeated multiplier contractions and all ladder cut counts, and checks the two-child BBF interpolation identity. These finite checks do not verify a general GN routing or an interacting RG.

[record, primary context] The field labels, weighted-line differentiation and paired scale organization are discussed in Sections 3.3-3.4 and 3.9-3.11 of [Benfatto and Mastropietro](https://arxiv.org/html/cond-mat/0004226). I checked those passages through the browser. Their construction is not invoked as an arbitrary-retarded-endpoint theorem or as a substitute for the unresolved induction above.

[record] The four requested work folders were read in full. All new files are in work_C74. Local STATE.md supplies the handoff under the prohibition on parent edits. No subagents, large computations, shell network commands or independent review were used.
