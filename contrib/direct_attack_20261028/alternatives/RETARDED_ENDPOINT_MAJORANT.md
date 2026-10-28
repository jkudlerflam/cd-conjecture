# Arbitrary-degree retarded endpoints in a normalized GN forest

Status: **[proved endpoint-majorant lemma under the stated normalized covariance bounds; independently reviewed by the vacancy worker and root coordinator]**. The normalized covariance, derivative and Gram estimates are inputs from [MARKED_COVARIANCE.md](../infrared/MARKED_COVARIANCE.md). This note supplies the complementary endpoint, factorial, field-label and localization argument. It does not assume a marked beta-function estimate as an input, and does not prove a reference Ward identity, a current/stress limit, or a CFT. The independent reviews were mathematical AI reviews, not formal verification or human refereeing.

Author: alternatives worker, 2026-10-28. AI-assisted mathematical work; no human review or formal verification.

## 1. Endpoint space and precise claim

Work on the line, or first at a fixed infrared cutoff in the volume regime of MARKED_COVARIANCE. There are finitely many branches. A degree-p initial endpoint, p even, is an antisymmetric coefficient measure K_p on its relative spacetime coordinates; all species and creation/annihilation labels are included in its total variation norm. Delta-local temporal factors are permitted. Coefficients use one consistent ordered-monomial convention throughout.

For the external locations P define T(P) to be their minimum spanning-tree length in the positive metric

    d((x,t),(y,s))=a|x-y|+kappa|t-s|,   a,kappa>0.

The input is

    H^p sup_anchor integral exp[T(P)] d|K_p|
        <= A epsilon^d(p),
    d(p)=max(2,p/2-1),        eta=sqrt(epsilon).          (1)

A and the metric constants are fixed. H is chosen large enough for the finite branch conversion and the desired output field weight. The intended free quadratic term is excluded from (1). At the starting scale, extract the quadratic constant and linear coefficients and the resonant quartic constants into the relevant/marginal coordinates; the initial remainder has R_2 K_2 and R_1 K_4 in the kernel/field sense used below. Local forward and relevant running endpoints have their own magnitude at most eta. At scale zero their prescribed smallness is as in COUNTERTERM_INDUCTION. An unextracted constant quadratic or quartic endpoint is not assigned an irrelevant short-memory gain.

The static finite-degree perturbation is included: DEGREE_REPARAMETERIZATION (7a)--(7b) gives (1) after choosing its independent small parameter before the final microscopic coupling. The order-three matched residual satisfies the same hierarchy. [TREE_NORM_EXTENSION.md](../infrared/TREE_NORM_EXTENSION.md) supplies the stronger spatial/time tree weight directly from the massive polymer, with no degree-dependent conversion from a diameter norm.

The following equivalent decorated-endpoint input is also allowed and is useful for band projection. Start with a kernel satisfying (1) in the original finite orbital labels. Attach one fixed linear conversion kernel L_f to each field leg, with its species sum and weighted L1/total-variation moments through order four bounded by a fixed constant. The actual branch endpoint is the convolution of K_p with these p kernels. Its bounding representation, rather than its already-convolved coefficient, is retained in the forest. This extension is proved in Section 4. Smooth compact branch cutoffs give such finite moments but generally do not give exponential Fourier decay; no claim that they preserve (1) literally is needed.

**Claim.** Given the normalized covariance bounds, a fixed reduced GN tree and its scale labels obey the usual renormalized dimensional bound, with a factor C^p times the norm in (1) for every initial degree-p endpoint. All field labels, determinant constants, anchored tree choices, localization choices and one or two specified marks are summable with one fixed field-radius reserve. The constants do not grow with the number of scales. In particular, the sum of trees with at least one initial endpoint has a factor O(epsilon^2) and an initial-scale short-memory gain after normalizing its quadratic/quartic root to beta-function units.

The marks here are endpoint selections/replacements, choices of one or two field labels, or the one-coordinate/one-common covariance and normalization marks already bounded in MARKED_COVARIANCE. A replacement coefficient must satisfy the corresponding norm (1), with its own prefactor. Derivatives of moving Fermi-center phases are not asserted in this lemma.

## 2. The anchored-tree factorial is exponential in fields

At a truncated expectation with s>=2 children, let q_i be the number of fields from child i contracted at this vertex, and Q=sum_i q_i. A connected anchored tree uses d_i distinct fields at child i, where d_i is its degree in the abstract child tree. Ignoring signs and creation/annihilation restrictions only enlarges the count.

For fixed d_i>=1, sum_i d_i=2s-2, the number of abstract labeled trees is

    (s-2)! / product_i (d_i-1)!.

There are at most (q_i)_(d_i) ways to assign distinct field labels to the incident edges. Including the factor 1/s! from the truncated expectation, their total is

    [1/(s(s-1))] [z^(s-2)]
        product_i {q_i(1+z)^(q_i-1)}
      = [product_i q_i/(s(s-1))] binom(Q-s,s-2)
      <= 4^Q.                                            (2)

Here (q)_d is a falling factorial. If a q_i vanishes the connected contribution vanishes. Edge orientations give another harmless exponential factor. For s=1 there is no anchored-tree sum; all contractions remain in the determinant. This proves (2) without replacing a determinant by its factorially many Wick pairings.

The Gram bound for the remaining matrix costs C to the number of its rows, hence C^Q. Across the GN forest,

    sum_v Q_v = sum_initial/leaves p_i - p_root,           (3)

because each field is contracted once. Thus these constants telescope to C^(sum p_i), rather than being paid once per scale for every surviving field.

## 3. Localization zeros: the bounded-load construction

The following kernel-field formulation is important when K_p is a measure. Taylor localization is performed on external fields, or equivalently on the output kernel as a distribution; it never differentiates the initial measure.

Use remainder order m=2 for a quadratic cluster and m=1 for a quartic cluster. Let R_m=1-T_(m-1), with T the Taylor polynomial at the localized momenta. Multiplication of the position-space kernel by a relative-coordinate monomial of degree b is momentum differentiation. Consequently, for b<=m,

    x^b R_m K = R_(m-b)(x^b K),                         (4)

with the appropriate tensor/multi-index version and R_0=1. For b>=m, the subtracted polynomial is killed. For a quartic contact subtraction, any external-coordinate difference kills the contact term. Formula (4) is an exact distributional identity, not an estimate.

Fix localization points using the field's next contraction toward the root. For a quadratic cluster choose its difference field to be a field contracted at the next genuine merging/contraction vertex, and use the other field as pivot; if both contract together, fix either choice deterministically. For a quartic cluster, if an ancestor retains exactly one of its four fields, use that field as pivot. Any two such singleton choices coincide by nesting. Otherwise retain a deterministic choice, with the common-pivot convention in the mixed case below. Use the unique paths in the following augmented tree: the BBF tree joins child clusters, and each initial endpoint is equipped with a minimum spanning tree of its field locations. An endpoint-tree edge is only a coordinate connection; it is not an additional propagator or contraction.

Process clusters from ancestors to descendants. A coordinate zero joining two external positions is written as the sum of edge displacements on its tree path. When the path enters a child, its entry and exit are external field positions of that child, because the parent connects children only through their uncontracted fields. Before expanding the zero further inside the child, apply (4). Each zero already entering that child reduces its required subtraction order. A quadratic child therefore requires at most two incoming-plus-new zeros; a quartic child at most one new zero and at most two total; a higher-degree child produces no new subtraction zero. If two incoming zeros are present they are kept as a pair, rather than re-expanding a new second-order remainder.

This gives an induction with the following coordinate invariant: each recursive kernel carries at most two unallocated zeros; after distributing them into its child paths, every line or endpoint-tree edge carries coordinate degree at most two. The derivative-carrier invariant needs a separate check: a zero routed outside a child does not erase a difference field that still survives in that child. The following cases supply that check.

- Two quadratic regularizations cannot repeatedly load the same surviving field across a genuine merge: the chosen difference field is contracted at the next such vertex. Equal-valence unary repetitions are removed by idempotency.
- For two quartic regularizations, the localization of the earlier difference field at the later contact point vanishes. The additional contact subtraction is therefore zero, rather than another derivative on that field.
- The remaining case is a quadratic cluster above a quartic cluster with both quadratic fields surviving from the quartic cluster. Choose their pivot consistently. The quadratic relative-coordinate zero kills the quartic contact subtraction by (4). Remove the duplicate quartic difference operation. Its dimensional gain is retained because that zero can be evaluated at the shorter, descendant-cluster scale; the resulting scale ratio supplies precisely the gain of the erased quartic remainder.

The pivot rule also protects the pivot from being independently differentiated along the same nested chain. A quartic/ quartic change can relocate a pivot, but cannot additionally differentiate the corresponding already-loaded difference field, since its contact subtraction has just vanished. These cases exhaust the possible nontrivial two/four-field remainders. Hence a field needs at most a second-order difference/derivative, not an order growing with tree depth. A covariance line can receive derivatives at both ends, a fixed total order covered by MARKED_COVARIANCE. None of this case analysis restricts the degree of an initial endpoint.

Unary repetitions do not introduce extra choices: if no fields are removed between two equal-valence localized clusters, R_m^2=R_m. Along a reduced unary segment, a nontrivial remainder can first occur at degree four and then at degree two, at most twice. A segment that contracts fields is charged to those removed fields in (3). This is the reason the construction does not spend a norm radius at every RG scale.

The routing count is local and does not require injective reconstruction from final edge loads. At a vertex with s children, a zero path has at most s-1 pieces on its BBF edges and at most s pieces inside child clusters. At most two zeros therefore give at most (2s-1)^2 choices. For s=1 this is one; for s>=2 it is bounded by C^(s-1). Since sum_v(s_v-1)=n-1 in a rooted tree with n endpoints, these choices cost C^n in total, regardless of depth. At a degree-p initial endpoint, each of its at most two incoming zeros chooses one of at most p-1 edges in the fixed endpoint tree, for a factor at most (p-1)^2<=C^p. Coordinate directions, pivot and remainder labels cost a fixed factor per effective low-valence remainder vertex, whose number is O(n+sum p_i). The resulting total is C^(sum p_i+n). This counts different routing histories even if they lead to the same final loads; no unproved uniqueness claim is needed.

The Taylor parameters are integration variables, not additional field labels. Each retained first- or second-order difference uses one parameter in [0,1], with weight 1 or 1-t of total mass at most one. The derivative-carrier construction bounds their number by a constant times the original number of fields. An interpolated position keeps its original field label and its prescribed pivot and parameter; it is never recounted as a new endpoint when a later BBF graph is chosen. Thus continuous interpolation creates neither a degree factor nor an independent path-record multiplicity.

These are the algebraic ingredients of the bounded-zero construction in [Benfatto--Mastropietro, Sections 3.3, 3.6--3.9](https://arxiv.org/pdf/cond-mat/0004226). The added endpoint tree changes only where a zero can be routed at a terminal kernel; the induction depends on the two/four-field cluster subtractions and the external-field nesting, not on a quartic bound for the terminal degree. Its bounded load is explicit above rather than inferred from power counting.

## 4. Paying endpoint moments without factorials

For a chosen endpoint spanning tree write its nonnegative weighted edge lengths as r_e, so sum_e r_e=T(P). Every routed coordinate factor has bounded degree b_e<=2. Therefore

    product_e r_e^b_e <= product_e (b_e! exp(r_e))
                      <= 2^(p-1) exp[T(P)].              (5)

Changing from a physical coordinate component to the metric length introduces only fixed powers of a^-1 and kappa^-1 per loaded edge, hence C^p. Two additional coordinate-moment marks, if required, can be assigned to the same tree with b_e<=4, replacing 2^(p-1) by 24^(p-1) and another C^p routing count. Ordinary history marks generate no such extra coordinate moments.

Integrating (5) against the endpoint total variation gives C^p times (1). No moment of order growing with p is bounded by a single diameter moment of that order; that shortcut would produce a factorial. The exponential tree weight pays each bounded edge load separately.

For a decorated endpoint, adjoin its p conversion legs to the original coordinate tree. Each attached leg has a bounded moment of the required order, and there are only p extra edges. The routed load on it is bounded by the same two-zero induction (or four including the optional extra two coordinate marks). Integrating the leg kernels therefore costs another C^p, preserving the degree hierarchy. A smooth band eigenvector on a fixed isolated patch, multiplied by a smooth extension equal to one on the field support, has a Fourier kernel with these moments by finite-order integration by parts; temporal delta factors are harmless. Species sums are finite. This handles branch projection without a false exponential-decay claim for a smooth compact cutoff.

For a Gaussian tree line, the same bounded moment/derivative load is paid by its fixed-order scale-dependent L1 bounds. The following operation-order argument is necessary; moving an arbitrary integration variable toward a different cluster's anchor would instead produce a possibly singular Jacobian.

At cluster v, first take its truncated Gaussian expectation and contract its internal fields. R_v then acts only on the remaining external fields. A Taylor-interpolated field at v cannot be an endpoint of any Gaussian line already contracted at v or in a descendant. Suppose it first contracts at an ancestor u. The whole cluster v, including its localization pivot, belongs to a single child w of u. Both the original point and every interpolation pivot therefore lie in w, even when they belonged to different immediate children of v. Ancestor-first processing above refers to the allocation of zeros in this already-constructed forest; it does not reverse this Gaussian-contraction order.

Write every original coordinate in w as y_w+r_f, with y_w its common anchor. Its interpolated field coordinate has the form

    (1-t)(y_w+r_f)+t(y_w+r_pivot)
        =y_w+[(1-t)r_f+t r_pivot].                       (5a)

Repeated permitted relocations remain affine combinations whose coefficients sum to one. For the BBF tree T_u between children of u, a line between w and w' thus has argument

    y_w-y_w' + rho_line(relative coordinates,parameters).

Pin one child anchor. The map from the other child anchors to the signed line differences has exactly the reduced incidence matrix of T_u, with determinant +/-1 in each continuous coordinate. Offsets may depend on internal relative coordinates and the Taylor parameters; they do not change that determinant. Integrate the outer child anchors holding these offsets fixed, then proceed recursively to the child relative variables. The complete change of continuous variables is block triangular with these incidence blocks, so no factor 1/(1-t) occurs.

For the discrete spatial coordinate, fractional interpolation offsets need not map the integer lattice onto itself. No false discrete change of variables is used: the covariance bound holds uniformly for real offsets, and its shifted lattice L1 sum has the same scale bound. The incidence matrix is unimodular, so the unshifted integer anchor differences still give independent lattice sums. At finite temperature the same statement uses periodized kernels and circular distance.

Initial coefficient measures are disintegrated into a common anchor and relative coordinates, or represented by the equivalent translation-covariant pushforward measure. Formula (5a) changes only the arguments of fields that later become propagators; it does not differentiate or rescale the initial relative-coordinate measure. Scalar products and bilinear kernel-measure contractions are therefore bounded by total variation and the same uniform shifted covariance estimates. The Gram estimate treats non-tree lines and their fixed derivative load. This proves the anchor integration step for singular as well as ordinary retarded endpoints.

## 5. Field survival labels and the scale sums

Before the small normalization correction, the positive renormalized dimension of an internal cluster of even field degree p>=2 is

    d_R(p)=p/2-2+z(p),
    z(2)=2, z(4)=1, z(p)=0 for p>=6.

It satisfies both d_R(p)>=1 and d_R(p)>=p/6. The fixed independent ratio domain in MARKED_COVARIANCE is chosen so its growth, proportional to p times the edge scale length, uses only a fixed fraction of these dimensions. Split the remaining gain into three fixed positive reserves: a field-label reserve, a scale-summation reserve, and a marked-path reserve.

For the field-label reserve use, for example, sigma^(sum_v |P_v|) with a fixed 0<sigma<1 furnished by a small fraction of p/6. Each original field follows one ancestral path until it is contracted or remains external at the root. Ignoring constraints only enlarges the label sum. Summing its possible stop depths gives

    sum_(ell>=0) sigma^ell = (1-sigma)^-1.

Multiplying over the original fields bounds all compatible surviving-field sets by C^(sum p_i). This is a single global allocation from the dimensional gains, not a succession of field-radius reductions.

Remove trivial unary vertices from the combinatorial skeleton; their scale spans remain on its edges. Branch vertices number at most n-1. Unary contraction vertices remove at least two fields and hence number at most half the total initial degree; effective remainder vertices add only O(n+sum p_i). The geometric scale sums from the fixed residual dimension cost C per skeleton edge, again at most C^(sum p_i+n). Tree topologies and finite branch labels cost another such exponential. The marked-path reserve gives 2^(-alpha(k-h)) for a designated dependency at scale k, with a fixed alpha>0. A tree touching the initial scale k=0 gains 2^(alpha h).

This supplies the renormalized tree majorant with the same root engineering factor 2^(h(2-p_root/2)) as usual. Quadratic relevant roots divide by their canonical 2^h, and a quadratic derivative root pays one inverse scale; both leave dimensionless beta estimates. Quartic roots already have zero engineering dimension.

## 6. Degree sums and one or two marks

All preceding constants have the form C^(sum p_i) once the covariance parameters and the metric have been fixed. Fixed powers of the total field count from one or two determinant/history marks are absorbed by increasing C once. An endpoint-field selection contributes at most p or p^2 and is handled in the same manner.

The initial hierarchy then gives the sharp useful summation

    sum_(p even>=2) C^p epsilon^d(p) <= C' epsilon^2.     (6)

Indeed p=2,4,6 are a finite set with d(p)=2. For p>=8,

    d(p)-2=p/2-3 >= p/8,

so their tail is bounded by epsilon^2 times a geometric series in C epsilon^(1/8). The same argument works with any fixed polynomial p^r multiplying the summand. Choose epsilon after C is fixed. A strict fixed reserve between the initial H and the output field weight pays the external-field factors and any remaining finite conversion constants.

Every running local endpoint has degree two or four and costs C eta. Trees with n total endpoints and at least one initial endpoint are consequently summed by a geometric endpoint-count series with one designated O(epsilon^2) factor. Designating one or two initial endpoints, or differentiating scalar parameters multiplying them, introduces at most n^2, still summable. Replacing a designated endpoint by a direction satisfying (1) replaces its prefactor A by the norm of that direction; multilinearity gives the identical estimate.

For history differentiation the initial data stay fixed. A once-marked analytic history dependence of an initial-endpoint tree can alternatively be bounded on an eta-sized complex disk: its epsilon^2=eta^4 budget becomes eta^3. A second common-history mark leaves eta^2. These are stronger than the eta and 1 required by GN_MARK_CALCULUS. The marked covariance estimates preserve the same scale dimensions and short-memory path reserve; their fixed degree polynomials have already been paid above.

The constants may be very poor as the finite branch count, code or memory metric changes. They are nevertheless chosen once, before epsilon and the microscopic coupling, and do not depend on the infrared cutoff or on field degree.

## 7. Scope and independent review

This construction supplies the endpoint side of the ordinary normalized GN bound, including arbitrary-degree temporally nonlocal coefficient measures. It uses the geometric tree norm already supplied by massive matching and the explicit covariance lemmas, rather than postulating G1/G2.

The bounded-zero routing argument in Section 3 was the most delicate step. The vacancy worker independently checked its carrier-load bookkeeping, mixed two/four-field cancellation, and local routing count against the primary construction. The root coordinator challenged the interpolation Jacobian; the resulting first-contraction-vertex and incidence-matrix argument was then checked independently by both reviewers. Both accepted the revised argument. Replacing these steps by an informal appeal to exponential moments would not suffice.

The normalization resolvent and history calculus use this lemma together with the separately audited covariance input. Reference comparison must still be made in its proper normalized cutoff scheme. Source-dependent comparison and a stress/pressure or conformal-limit argument are separate from this endpoint lemma. A further peer-chat audit may supplement the reviews recorded here; no such additional verdict is claimed yet.
