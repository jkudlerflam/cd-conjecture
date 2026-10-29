# C72. Commuting moments with a survival-adapted Taylor decomposition

## 1. Status and scope

[assumed, unresolved target] I have not proved the simultaneous tree induction requested in C72. I have also not constructed a legal tree whose loads are unbounded under every reorganization. C60 (17) and the source-free part of corrected Lemma N remain unresolved here.

[proved below] There is an explicit finite-order formula that processes incoming moments before exposing surviving momentum factors. It includes the derivatives of the Taylor kernels and works with arbitrary partitions of the independent port momenta. For one application it has kernel differentiation order at most two, bounded polynomial degree on each field, bounded coefficient mass, and affine relocations with permitted anchors. Sections 2 through 5 prove these assertions. They do not assert that the input representation is closed under a GN contraction.

[proved below] A complementary obstruction explains why that last qualification matters. Keeping a particular surviving factor explicit after two derivatives in other port directions can require a third derivative of the underlying kernel. No formula using only derivatives through order two and uniformly bounded coefficient mass can avoid this while retaining that factor. Section 6 proves this statement for an explicitly defined class of formulas. Dropping the factor removes the obstruction. It is not an obstruction to every reorganization, nor a counterexample to the strong endpoint norm.

[assumed, imported inputs] Retain C60's fixed real cutoffs, exact quadratic absorption, infinite-volume limit at fixed infrared cutoff, and separate treatment of engineering derivatives. Retain C61's corrected source subtraction and its oriented path convention. Endpoint measures and fixed conversion legs have the norms stated in C61 and the collaborator's TREE_NORM_EXTENSION. Covariance and Gram bounds are inputs from MARKED_COVARIANCE. The matching, normalization, reference-comparison and source results are not rederived here.

[record] All files in C61, C64 and C70 were read. The S and S2 backgrounds, C60, the current endpoint-majorant argument and the obligations ledger were consulted. This investigation uses their settled local identities. It does not replace a missing induction with their conditional consequences.

## 2. Partial Taylor projectors

[definition] Eliminate one pivot momentum using conservation. Write the remaining momentum variables as q=(c,s), with c the port components to be contracted before a specified ancestor and s the components surviving there. A component can be a time or spatial deviation momentum. The formulas below are identities in independent variables on this conservation hyperplane.

Let T^s_{r-1} mean Taylor expansion in s about zero, through total s-degree r-1, with c held fixed. Set

    R^s_r = 1-T^s_{r-1},       R^s_0=1.

Define R^c_r in the same way. Write R_r for total Taylor subtraction in all components. Evaluations such as H(c,0) retain all c dependence. A multi-index factorial is the product of its component factorials.

[proved, flag decomposition] For every integer r>=1,

    R_r H(c,s)
      = R^s_r H(c,s)
        + sum_(|alpha|<r) s^alpha/alpha!
             R^c_(r-|alpha|) [D_s^alpha H](c,0).          (1)

To prove this, first write H=R^s_r H+T^s_{r-1}H. The first summand has no total Taylor coefficient of degree less than r, since all of its s derivatives of order less than r vanish at s=0. Its total R_r is therefore itself. In the second summand, multiplication by s^alpha raises total degree by |alpha|. Subtracting the remaining total jet is consequently the same as applying R^c_(r-|alpha|) to its coefficient. This proves (1) using only finite jets. Analyticity is unnecessary.

[proved, integral realization] For a block x with d components and r>=1,

    R^x_r H(x)
      = sum_(|gamma|=r) (r/gamma!) x^gamma
          integral_0^1 (1-t)^(r-1) D_x^gamma H(tx) dt.     (2)

Apply the one-variable Taylor theorem to h(t)=H(tx). The chain rule gives

    h^(r)(t)=sum_(|gamma|=r) (r!/gamma!) x^gamma
                                      D_x^gamma H(tx).

The integral remainder has factor 1/(r-1)!, which gives (2). Spectator variables remain fixed. The measure (1-t)^(r-1)dt has mass 1/r, so the sum of the absolute numerical coefficients, including this integral mass, is

    sum_(|gamma|=r) 1/gamma! = d^r/r!.                    (3)

[proved, cases used here] In one c and one s variable the first two orders are

    R_1 H = R^s_1 H + R^c_1 H(c,0),                       (4)

    R_2 H = R^s_2 H + R^c_2 H(c,0)
                              + s R^c_1 H_s(c,0).       (5)

The c-only coefficients in (4) and (5) are independent of s. This independence must be retained when incoming moments are processed. Replacing these formulas by a radial integral in all variables loses that visible structure, even though it leaves the function unchanged.

## 3. Incoming moments must precede the split

[proved, partial commutation rule] Let E_s H(c,s)=H(c,0). For beta=(beta_c,beta_s),

    D^beta (1-E_s)H = (1-E_s)D^beta H,   |beta_s|=0,
    D^beta (1-E_s)H = D^beta H,          |beta_s|>0.        (6)

In the first case c differentiation commutes with evaluation in s. In the second case at least one s derivative annihilates E_s H. In particular, keeping (1-E_s) after a derivative in s would introduce an erroneous subtraction E_s D^beta H.

[proved, full commutation rule] For b=|beta|,

    D^beta R_m K = R_(m-b)_+ D^beta K,                    (7)

where (m-b)_+=max(0,m-b). Differentiating each term of the finite total Taylor polynomial proves the identity. For b>=m the differentiated Taylor polynomial vanishes.

[proved, explicit moment-first interpolation] Set r=(m-b)_+ and H=D^beta K. If r=0, use H without another split or interpolation. If r>0, use

    D^beta R_m K
      = R^s_r H(c,s)
        + sum_(|alpha|<r) s^alpha/alpha!
             R^c_(r-|alpha|) [D_s^alpha H](c,0),          (8)

and expand only the displayed partial remainders using (2). The resulting kernels are evaluated at (c,t s) or (t c,0). Their differentiation orders on the original K are all

    b+r = max(b,m).                                      (9)

For the first term this follows from b+r derivatives. For a jet term the count is b+|alpha|+(r-|alpha|). Thus m,b<=2 require only derivatives of K through order two. No division by t or by a momentum component occurs.

Equation (8) is also an explicit prescription for recombining derivatives of an already expanded Taylor integral. Restore its partial remainder by (2), add all terms belonging to the original R_m, apply (7), and then use (8). These are equalities before absolute values. This prescription requires access to the entire indicated sum. It does not authorize using (8) on one separated coefficient from that sum.

[proved, example with a differentiated survival kernel] For m=2 and beta=partial_c, differentiating (5) gives

    R^s_2 K_c + R^c_1 K_c(c,0)
                                      + s K_sc(c,0).     (10)

The first and third terms combine into R^s_1 K_c. Indeed

    R^s_2 K_c + s K_sc(c,0)
       = K_c(c,s)-K_c(c,0).

The moment-first representation is therefore

    R^s_1 K_c + R^c_1 K_c(c,0).                          (11)

Expanding the first term of (10) would have used K_css. Formula (11) needs only K_cs and K_cc. The cancellation crosses two components of the survival decomposition. It is not just an integration by parts within the first integral.

[proved, two incoming moments] Differentiating (5) twice in c gives

    R^s_2 K_cc + K_cc(c,0) + s K_scc(c,0)
       = K_cc(c,s).                                     (12)

Directly expanding the first term on the left would introduce fourth derivatives. The right side has order two. This illustrates why a second application of a formula to individual previously separated kernels does not follow from the bound (9).

## 4. Existing homogeneous factors

[definition] Let P_a(q) be homogeneous of degree a. For field-load statements take P_a to be one monomial in the labeled field momentum components; general homogeneous polynomials are handled term by term with their coefficient norm. Fix beta with b=|beta|<=2. For delta<=beta set

    d=|delta|,
    P_delta=D^delta P_a,
    H_delta=D^(beta-delta) F,
    r_delta=max(0,m-b-a+d).                              (13)

Zero polynomials P_delta are omitted.

[proved, explicit formula with factors] C61's normal form, followed by (1), gives

    D^beta R_m(P_a F)
      = sum_(delta<=beta) binom(beta,delta) P_delta
                                         Phi_(r_delta) H_delta,    (14)

where Phi_0 H=H, and for r>0,

    Phi_r H
      = R^s_r H
        + sum_(|alpha|<r) s^alpha/alpha!
                R^c_(r-|alpha|) D_s^alpha H(c,0).         (15)

For completeness, (7) first moves D^beta across R_m. Leibniz distributes the derivatives between P_a and F. Every surviving P_delta is homogeneous of degree a-d, so

    R_u(P_delta G)=P_delta R_(u-a+d)_+ G.

This proves (14). Equation (15) is (1). Thus (14) includes derivatives hitting both the exposed factor and its coefficient, with the requisite recombination before the final integral expansion.

[proved, local counters] For each nonzero delta term define

    B=b-d+r_delta,
    A=a-d+r_delta,
    J=a+b-d+r_delta.                                     (16)

B is the number of derivatives on F after (2). A is the total degree of the remaining momentum monomial after that expansion. J is the algebraic cut count of C61, including the b incoming paths. One has

    B<=2,       A<=max(a,2),       J>=m.                  (17)

If r_delta=0, B=b-d<=2, A=a-d<=a and J=a+b-d>=m. If r_delta>0, substitute its definition to obtain

    B=m-a<=2,       A=m-b<=2,       J=m.

These cases prove (17). In particular, if the incoming monomial has degree at most two on each field, the output does too. When r_delta=0 no field degree increases. When r_delta>0 the total output degree is at most two.

[proved, local coefficient mass] At a two- or four-field remainder there are at most six independent scalar spacetime momentum components. The block sizes in (3) are therefore bounded. There are finitely many alpha and gamma indices with total degree at most two. If each of at most four fields has incoming degree at most two, then a<=8, and Leibniz differentiation of P_a through order two has bounded coefficient mass as well. Hence (14)-(15), with their parameter integrals included, cost a fixed constant C_0 independent of scale spans and of tree depth. This statement assumes the bounded input field degrees. Without that assumption the coefficients can grow with a.

[proved, boundary of the counter] B counts derivatives on F as supplied to this application. It does not count derivatives of an original endpoint hidden in a previously constructed F. For example, F may already be an integral of K_ss. Applying two derivatives to F is then an operation on fourth derivatives of K. Equations (10)-(12) show how a full local sum can repair this. They do not show that the full sum is still available after intervening contractions and different ancestor operations.

## 5. What the local formula preserves in position space

[proved, measure interpretation] Adopt the Fourier convention

    Khat(q)=integral exp(-i q dot r) K(dr).

D^nu Khat is the transform of (-i r)^nu K(dr). An evaluation Khat(c,t s) is the transform of the pushforward under

    r_c -> r_c,       r_s -> t r_s.

An evaluation at s=0 sends those relative coordinates to zero. This description remains valid at t=0 and never changes integration variables by dividing by t. Derivatives and evaluations in (8) therefore have an exact interpretation for finite coefficient measures with the indicated moments. Equations (14)-(15) have the same interpretation with the factor P_delta acting on the external test fields.

[proved, local anchor compatibility] Let y be the child anchor and write its original ports as y+r_i. The pushforwards just described move a port to y+t r_i or to y, leaving the anchor coefficient one. With the pivot as origin, these are permitted affine relocations within the same child. A later BBF line between children w,w' consequently retains the form

    y_w-y_w' + rho(relative coordinates,parameters).

Conditional on the first-contraction pivot property of C61, the reduced incidence determinant remains of modulus one. The spatial lattice uses the shifted lattice bound, as in C61. This verifies anchor compatibility for this local formula. It does not verify the pivot property for a globally regrouped forest.

[proved, local endpoint payment] When K itself is an endpoint coefficient, (9) places at most two coordinate factors on its relative coordinates. Expand each relative coordinate along the fixed endpoint tree. Every edge then has load at most two. With metric lengths rho_e,

    product_e rho_e^(b_e)
       <= product_e b_e! exp(rho_e)
       <= 2^(p-1) exp(sum_e rho_e).                      (18)

This is the endpoint payment already established by C61. It applies to (8) directly when the input is the endpoint. If the input is an integrated subtree kernel F, its total second moment does not identify the scales of the lines paying it. Applying (18) to that subtree as though it were a scale-zero endpoint would discard the issue in C60 (17).

[proved, conditional path accounting] Suppose the b incoming moments in (14) have admissible paths in the actual spatial representation, and P_a has labeled derivative paths. In the delta term, concatenate the d consumed pairs as in C61. The remaining b-d incoming moments and r_delta newly created moments account for B. Their b outer sinks, a-d surviving old sinks and r_delta new sinks give J. Thus the local count is compatible with the required cut demand m. A path assignment is part of the hypothesis of this assertion. Algebraic momentum degree alone does not specify it.

## 6. A precise obstruction to keeping every surviving factor

[definition] For scalar s define the divided-difference coefficient

    Q_s K(c,s)=[K(c,s)-K(c,0)]/s,       s!=0,
    Q_s K(c,0)=K_s(c,0).                                (19)

This is the coefficient of the explicit surviving factor s in a quartic survival term. Smoothness follows from its integral formula. The question here is whether partial_c^2 Q_s K can be represented using only derivatives of K through order two, with uniformly bounded coefficients as s approaches zero.

[proved, restricted impossibility theorem] No identity of the form

    partial_c^2 Q_s K(c,s)
      = sum_(|nu|<=2) integral D^nu K(A_nu(c,s,t))
                                                   d mu_nu(c,s;t)  (20)

can hold for every smooth K if the total variation summed over the measures is bounded by a constant C independent of K and of (c,s) near (0,0). Here the evaluation points are real, and may be arbitrary functions, not just affine ones. The measures may depend on (c,s). This is a statement about this bounded two-derivative representation class.

Proof. For an integer N>=1 take

    K_N(c,s)=N^(-2) exp(iN(c+s)).                         (21)

Every derivative of total order at most two has magnitude at most one at every real point. The right side of (20) therefore has magnitude at most C. At (c,s)=(0,0), however,

    partial_c^2 Q_s K_N(0,0)
       = K_N,ccs(0,0) = -iN.                            (22)

Taking N>C is a contradiction. If (20) is only stated for s!=0 but has a uniform variation bound near zero, continuity of the left side gives the same contradiction by s approaching zero. This proves the theorem.

[proved, the obstruction disappears on recombination] The whole survival component obeys

    partial_c^2 [s Q_s K]
       = K_cc(c,s)-K_cc(c,0).                            (23)

This uses two evaluations of second derivatives, with coefficient mass two. Dividing (23) by s returns the obstruction. Alternatively, adding the differentiated contracted component from (4) gives

    partial_c^2 [R_1 K]
       = [K_cc(c,s)-K_cc(c,0)] + K_cc(c,0)
       = K_cc(c,s).                                     (24)

Thus a formula with bounded order is available, but it need not keep the previously selected surviving factor explicit. This is a concrete reason to allow both the carrier labels and the component grouping to change in a successful induction.

[proved, qualification concerning the endpoint norm] Each K_N is the Fourier transform of a point measure with weight N^(-2) and relative displacement of size N. Its exponential terminal-tree norm grows with N. Multiplying (21) by a compensating exponential would put it in a bounded endpoint class and also suppress (22). Therefore (22) does not contradict the endpoint norm, C60 (17), or a construction allowing a fixed third moment. It rules out only the two-derivative, fixed-surviving-factor representation (20). No universal no-go statement follows.

## 7. A depth-independent bound for a fixed nested flag

[definition] Fix a single smooth kernel K in D independent scalar momentum components. Let A_1,...,A_L be nonempty coordinate subsets totally ordered by inclusion, with repetitions allowed. At step j apply a partial remainder R^{A_j}_{m_j}, with m_j in {1,2}, at the same origin. No contraction, change of independent variables or multiplication by another kernel is performed between these operations. This is a restricted operator chain, not a GN tree with merging vertices.

[proved, depth reduction] The product of these projectors reduces to at most two factors. For each order k present, retain only the smallest coordinate set B_k on which order k occurs. If both orders are present and B_2 is contained in B_1, the order-one projector is redundant. Otherwise the reduced product is

    R^{B_2}_2 R^{B_1}_1,       B_1 properly contained in B_2.       (26)

To prove this, partial Taylor projectors on coordinate subsets commute. One can verify this by expanding the finite jets in their common coordinates: both compositions retain exactly the multi-indices satisfying both total-degree cutoffs. The coefficients in coordinates outside both sets are spectators. The corresponding remainders also commute. If A is contained in B and p>=q, then

    R^A_p R^B_q = R^A_p.                                  (27)

Indeed the B-jet of order less than q is killed by R^A_p, and the B-jet of R^A_p K of that order vanishes. Equation (27) removes every larger set of a given order. It also removes the order-one factor if B_2 is contained in B_1. This proves the assertion. As an induction in L, retain at most the smallest set for each of the two orders after every insertion and apply (27). The number of retained projectors is at most two at every depth.

[proved, explicit interpolation for the remaining pair] Put B=B_1 and C=B_2 minus B_1. Regard coordinates outside B_2 as spectators. Applying the flag formula (1) to R^{B_2}_2 and then R^B_1 gives

    R^{B_2}_2 R^B_1 K
      = R^B_2 K
        + sum_(j in B) q_j R^C_1 [partial_j K](q_C,0_B).   (28)

The B-degree-zero term of (1) is killed. Its B-degree-one terms survive, and its B-remainder is unchanged. Formula (2) now represents (28) using derivatives of K of total order two only. Its numerical coefficient mass is bounded by

    |B|^2/2 + |B| |C| <= D^2.                            (29)

The same fixed-order conclusion holds for a single retained projector. Increasing L does not change this bound.

[proved, incoming moments] For any beta,

    D^beta R^A_m = R^A_(m-|beta_A|)_+ D^beta.              (30)

This follows by differentiating the finite partial jet. Move D^beta through the entire product using (30), omit the resulting order-zero identities, and repeat the preceding reduction. The residual orders are still at most two and their coordinate sets remain nested. For b=|beta|<=2, the final interpolation therefore uses derivatives of K through order at most b+2<=4, with coefficient mass bounded independently of L. In cases where the surviving projectors see all derivative directions this improves to the order-two count of Section 3. Four is the uniform bound for arbitrary derivative directions in this restricted chain.

[proved, interpretation] Arbitrarily many nested coordinate-subset projectors on one fixed kernel do not by themselves force unbounded moment order. The reduction is an actual induction on the chain length. It does not handle intervening Gaussian contractions or the spatial paths needed to charge distinct scale cuts. In particular, reducing many projectors to one does not justify deleting their requested scale gains.

## 8. Why this is not yet an induction over the GN tree

[proved, distinction between two operations] After a child kernel has been split into s Q_s K and K(c,0)-K(0,0), a later contraction can treat these components differently. The first has an exposed surviving momentum factor; the second does not. C70 uses that distinction to omit an ancestor quartic contact subtraction in the first component. Equations (10)-(12) and (24) require reuniting components before assigning bounded kernel loads. That reunion must also restore or otherwise account for the different ancestor operations.

[assumed, missing closure assertion] The needed inductive object is consequently larger than a list of monomials P_a times independently bounded Taylor kernels. It must keep enough of the unsplit operator expression to perform (8) or (14) while retaining all spatial source/sink paths and contraction data. The following are not proved here:

1. After applying (14) at a vertex, all descendant coefficient sums required for its next application can be recovered with total coefficient mass at most C^(sum_i p_i+n).
2. Changing or removing surviving factors as in (23)-(24) admits a simultaneous reassignment of their spatial paths, with cut coverage on all affected ancestor edges.
3. Those regroupings preserve the first-contraction pivot condition, including when several ancestor test fields share child ports.
4. The locally bounded derivatives on F in (17) imply bounded derivatives and moments on original endpoint measures, actual shifted tree lines and Gram rows after arbitrary depth.

[proved, why the local counter cannot establish item 4] The order bound in (17) can be applied to an already differentiated coefficient and then bound differentiation order relative to that coefficient. This does not bound total differentiation order relative to its original endpoint. Equation (12) supplies an explicit order-four intermediate term with a correct order-two local recombination. An induction must prove availability of that recombination after every merge. Repeating the numerical inequality B<=2 without such an object is insufficient.

[proved, implication if closure is supplied] If the four assertions hold in one representation, C61's path identity gives

    sum_f a_f s(f)-sum_l b_l s(l)
       = -sum_v F_v ell_v <= -sum_v z_v ell_v.            (25)

The coefficient-mass bound then gives C60 (17). Endpoint payment follows from (18), with the corresponding fixed-order line and Gram estimates. This is a conditional implication, not a new proof of the four assertions.

[assumed, additional boundary for corrected Lemma N] Even a completed source-free (17) would leave C61's reference-comparison and source hypotheses to be verified in the same scheme. In particular, the interaction-prefactored source difference is D_O(u)-D_O(0). The free lattice-to-reference difference D_O(0) remains separate. Source normalization and the two-source estimate do not follow just from a source-free spatial routing theorem.

## 9. Checks and provenance

[proved, consistency] With no s variables, (8) is the ordinary Taylor formula after moment commutation. With no c variables it is the partial remainder itself. If b>=m, (8) has no Taylor integral and the original contact polynomial is annihilated. If a>=m and b=0, (14) adds no new moment or derivative. The examples (10)-(12) agree both at s=0 and at c=0 with direct differentiation. These are identities, not asymptotic tests.

[numerical, finite exact-arithmetic checks] checks.py tests the flag decomposition, its integral coefficients, the moment-first formula and the formula with homogeneous factors on finite polynomials. It also checks the local counters and the cross-component cancellations (10)-(12), verifies the derivative magnitudes in (21)-(22) by exact integer arithmetic, and tests (26)-(30) on nested flag chains through depth 40. All checks pass. Counts are in checks.json and PROGRESS.md. These checks do not construct a GN tree or test an interacting RG.

[record, primary-source context] Benfatto and Mastropietro keep weighted propagators under the differentiation operators in Section 3.9 and state the paired scale estimate in Sections 3.10-3.11. Their discussion uses a prescribed localization-label scheme. I checked those passages and Sections 3.3-3.4 in the [primary paper](https://arxiv.org/html/cond-mat/0004226). They are not being invoked as an abstract theorem for the current endpoint class.

[record] No subagents, large computations or network shell commands were used. All writes are in work_C72. The local STATE.md supplies the handoff because the task prohibits editing parent memory. No independent review is claimed.
