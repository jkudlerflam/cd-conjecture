# C61. Scale allocation with a retained incidence record

## 1. Scope and input ledger

[proved below] A cut criterion implies C60 (17). The elementary Taylor identity with an incoming coordinate moment is also proved, including every product-rule term. These results do not by themselves prove that the proposed spatial routing has the required cut criterion. Section 6 states the missing simultaneous construction.

[assumed, corrected Lemma N] Section 2 states the corrected nonlinear target. Its source assertion is a separate hypothesis. A bare source is not assigned the smallness of an interaction endpoint. No interacting Luttinger theorem is asserted here.

[imported] C55 initializes the off-half-filled single-channel action and controls one initial endpoint under Gaussian integration. C60 Section 3 proves the degree hierarchy after finite initialization. FINITE_CELL_MATCHING and TREE_NORM_EXTENSION supply massive matching and the exponential terminal-tree norm. DEGREE_REPARAMETERIZATION gives the stronger residual hierarchy under its schedule. These results are not rederived.

[read] CURRENT_OBLIGATIONS still calls for review of simultaneous zero and derivative routing. RETARDED_ENDPOINT_MAJORANT Section 3 claims bounded loads; Section 4 proves anchor incidence for permitted interpolations. MARKED_COVARIANCE supplies fixed-order line and Gram bounds. NORMALIZATION_RESOLVENT and GN_MARK_CALCULUS are implications from raw tree bounds. CONSTANT_HISTORY_COMPARISON retains a separate reference-scheme hypothesis. None is used to assume the missing simultaneous construction.

[convention] Work at a fixed infrared cutoff, remove volume and inverse-temperature regulators there, then send the cutoff down. No assertion about a spatial zero mode below the finite-volume spacing is needed. Use C60 Section 2's fixed real branch cutoffs and exact absorption convention. Scalars have a separate density norm. Writes stay in work_C61.

## 2. Corrected statement of Lemma N

[definition] Fix a compact Fermi-point interval away from band edges and half filling. Set E_h=2^h for h<=0. Use one antisymmetric coefficient convention with the degree factorial included. The source-free action has prescribed small local quadratic and forward-quartic coordinates and initial irrelevant kernels. One sufficient strong endpoint input is

    H^p sup_anchor integral exp[T(P)] d|K_p|
       <= A epsilon^max(2,p/2-1), p even, p>=2.        (N0)

T is the terminal minimum-spanning-tree length in a fixed positive space/time metric. Fixed band-conversion legs can instead be retained in the bounding representation, with fixed moments. A smooth compact Fourier cutoff is not assumed to preserve exponential decay. Keep the bare free dispersion in the covariance; its curvature is not interaction-proportional.

[assumed, N-action] The source-free normalized renormalized expansion and its specified one- and two-history marks converge uniformly in infrared depth on a fixed small history domain. Solve the normalization equations there, with single-step ratios in a fixed neighborhood of one. Spend one fixed field-radius reserve. In each routed tree, bound the coordinate moments, derivatives, anchor integrals and scale gains in the same representation. A sufficient scale bound is C60 (17), followed by

    d_R(p)=p/2-2+z(p),
    z(2)=2, z(4)=1, z(p)=0 for p>=6.                  (N1)

Thus d_R(2)=d_R(4)=1 and d_R(p)>=p/6. Keep explicit engineering derivatives in the dimensional count. They cannot be omitted from the localization weight and also used to pay z(p).

[assumed, N-comparison] In the same cutoff, absorption and history convention, initial-irrelevant-endpoint and covariance-difference trees have short memory. A sufficient beta-function formulation is

    |b_h^mu| <= A delta 2^(theta h)
       + B eta sum_(k>h) 2^(-kappa(k-h)) |mu_k-mu_h|
       + B eta sum_(k>=h) 2^(-kappa(k-h)) |nu_k|,

    |b_h^nu| <= A delta 2^(theta h)
       + B eta sum_(k>=h) 2^(-kappa(k-h)) |nu_k|,       (N2)

with 0<theta<kappa. Here eta bounds local interactions and delta bounds ultraviolet forcing. The reference's ultraviolet kernels can be O(eta), so delta need not equal the stronger physical residual budget O(epsilon^2). Reference constant-history cancellation is a separate signed-sum input. An absolute bound on each reference tree is not that cancellation. Include the history-difference bounds and continuity required by the cited conditional sequence theorems.

[definition, free subtraction] For specified small action data u, write the connected lattice and matched reference correlations as C_O(u;z) and C_O^ref(u;z), including the reference amplitude and oscillation convention. Define

    D_O(u;z)=C_O(u;z)-C_O^ref(u;z),
    F_O^free(z)=D_O(0;z),
    D_O^int(u;z)=D_O(u;z)-D_O(0;z).                    (N3a)

The zero datum retains the same bare source and free lattice dispersion. If u changes the Fermi point or free speed, specify the Gaussian comparator at those same parameters and subtract its difference there. This defines the subtraction; it does not prove a uniform estimate.

[assumed, N-source] The one-source normalized ratio difference after its own Gaussian subtraction is summable, and its matched initial overlap remains nonzero. Separately require

    D_O^int(u;z)=F_O^cross,int(u;z)+sum_(h<=0) R_O,h^int(u;z),

    |R_O,h^int(u;z)| <=
       C ||u|| E_h^(2X_O(u)+alpha)/[1+(E_h |z|)^M],    (N3b)

with alpha>0 and M above the exponents throughout a fixed small domain. The crossover is generated by finitely many fixed smooth integrations. Constants can depend on fixed bare-source moments. A proof comparing anomalous powers may need a smaller alpha to absorb logarithms. Subtraction at u=0 alone does not prove (N3b).

[proved, necessary normalization check] At u=0 the interaction remainder vanishes but F_O^free need not. C60 gives, for the nearest-neighbor pair at p=pi/3 and x=3m,

    F_pair^free(x)=3/[4 pi^2 x^2(x^2-1)].              (N4)

This infrared algebraic descendant is retained, not renamed a finite smooth crossover. IR_SOURCE_REVIEW likewise separates fixed bare-source moments from interaction-generated source terms. A source mark is not an interaction endpoint of size epsilon^2.

[proved, conditional consequence] If (N3b) holds, split its scale sum at E_h |z|=1. Below that scale, sum E_h^(2X+alpha). Above it, sum |z|^(-M) E_h^(2X+alpha-M). Both geometric sums are bounded by C ||u|| |z|^(-2X-alpha). This controls D_O^int. Compare the separate free term with the leading reference decay before deducing an exponent statement.

## 3. The scale-cut identity

[definition] Let tau be a finite reduced rooted scale tree. The root has scale h, initial endpoints scale zero, and each child v has h_v>h_parent(v). Put ell_v=h_v-h_parent(v). Each coordinate factor on a covariance line of scale t is a source of one unit at that scale. An endpoint moment is a source at zero. Each field derivative is a sink at its actual first contraction scale s, or at h for a root external field. Count localization derivatives only in this definition.

[definition] An oriented incidence record pairs sources with sinks at ancestor nodes. If a product-rule contraction consumes an intervening derivative and zero, concatenate their paths, retaining the original source and final sink. Let F_v count the paths crossing edge (v,parent(v)), with multiplicity.

[proved, cut identity] Any such record satisfies

    sum_f a_f s(f)-sum_l b_l s(l)
       = -sum_(v!=root) F_v ell_v.                    (1)

For one pair, the traversed edges have total length t-s, and its contribution on the left is s-t. Sum over pairs and interchange the two finite sums. Endpoint moments contribute zero on the left but are included among sources.

[proved, sufficient criterion] Suppose one exact routing expansion obeys

    F_v>=z_v on each required remainder edge,
    sum_A |c_A| <= C^(sum_i p_i+n),                    (2)

where coefficient mass includes Taylor-parameter integration. Apply (1) termwise, multiply by |c_A|, and sum. The result is C60 (17). A separate favorable pairing chosen after the spatial estimate would not suffice.

[proved, necessity with independent spans] If a nonzero term has F_w<z_w, coefficient independent of ell_w, all other spans fixed, and remains legal for arbitrarily large ell_w, its ratio to the target grows as 2^((z_w-F_w)ell_w). No fixed C^(sum p_i+n) repairs it. The legality and nonvanishing qualifications matter: shell support can kill a formal low-load tree.

[proved, consumed-pair bookkeeping] For t>=u>=s, concatenate a path t to u with one u to s. The scale exponent is preserved because

    (u-t)+(s-u)=s-t.                                  (3)

Every cut count is also preserved. Erasing both records would discard a real gain. This identity does not show that every proposed spatial routing permits such concatenation.

## 4. Exact local normal form

[definition] Let R_m=1-T_(m-1), with T total Taylor expansion at localized external momenta and R_0=1. Let P_a be homogeneous of degree a in external momenta, representing exposed external field derivatives. Let beta be a multi-index of order b. Work with sufficiently differentiable Fourier kernels; the distributional identities follow by duality when the required moments exist.

[proved] For every m>=0,

    partial^beta R_m(P_a F)
      = sum_(delta<=beta) binom(beta,delta)
          (partial^delta P_a)
          R_[max(0,m-b-a+|delta|)]
                         (partial^(beta-delta) F).     (4)

Omit zero terms. First,

    partial^beta R_m G=R_[max(0,m-b)] partial^beta G,

since differentiating the Taylor polynomial lowers its cutoff by b. Apply Leibniz. Each partial^delta P_a is homogeneous of degree a-|delta|. For a homogeneous polynomial P_c,

    R_r(P_c H)=P_c R_[max(0,r-c)] H,

as follows by multiplying the homogeneous Taylor components of H by P_c. These two facts give (4), componentwise also on the momentum-conservation hyperplane. No infinite-order analyticity is needed.

[proved, local accounting] For a nonzero term put d=|delta| and

    r=max(0,m-b-a+d),
    J=a+b-d+r=max(m,a+b-d),
    B=b-d+r.                                          (5)

If b<=2 and m<=2, then

    J>=m,  0<=B<=2.                                  (6)

Nonvanishing implies a>=d. If r=0, B=b-d<=2. If r>0, B=m-a<=2. The two cases also prove the equality for J.

[conditional interpretation] Suppose the algebra admits an incidence realization in the same spatial tree. There are b incoming moments whose derivative sinks lie outside the cluster. Of the a old derivative paths, d are concatenated with incoming paths. The remaining outgoing counts are a-d old derivatives, b outer sinks, and r new derivatives, totaling J. There are B incoming or newly created coordinate factors still to place in F. Equation (6) is the local algebra needed for a cut of demand m with at most two unallocated zeros.

It is not yet the geometric statement F_v=J. That requires actual paths through the same interpolated BBF tree. Homogeneous degree and the number B alone do not establish it.

[proved, cases] For a quartic child, m=b=1,a=0 gives r=0,J=1: its contact subtraction disappears, but the incoming path must still cross its edge. For a quadratic child, m=2,b=1,a=0 gives r=1,J=2, while b=2 gives r=0,J=2. For m=2,b=1,a=1, the d=0 term has r=0,J=2 and the d=1 term has r=1,J=2. The latter product-rule term cannot be discarded.

[proved, order of operations matters] In one dimension,

    R_1 K(q)=q integral_0^1 K'(tq) dt,
    partial_q^2 R_1 K(q)=K''(q).                       (7)

Differentiating the integral expression first gives

    2 integral_0^1 t K''(tq) dt
      + q integral_0^1 t^2 K'''(tq) dt.                (8)

The second expression contains a third coordinate moment before cancellation, though (7) needs only a second moment. Integrating the derivative of t^2 K''(tq) proves equality. Taking absolute values before reorganizing loses that cancellation. Applying the operator identities before estimating Taylor integrals is therefore essential. The proposed note instructs this for incoming zeros; a simultaneous construction must enforce it also while exposing inherited derivative carriers.

## 5. Geometry that can already be justified

[proved, conditional endpoint payment] Given endpoint-tree loads b_e<=2 and metric edge lengths r_e,

    product_e r_e^b_e <= product_e b_e! exp(r_e)
       <= 2^(p-1) exp(sum_e r_e).                     (9)

The endpoint norm (N0) pays this with a fixed exponential-in-p loss. With two further coordinate marks, b_e<=4 gives 24^(p-1). The calculation does not prove the load bound.

[proved, permitted anchor integration] Suppose a field first contracts at u and all its earlier pivots and relocations belong to one child w of u. Write original positions as y_w+r_f. Iterated affine interpolations preserve the coefficient one of y_w, so actual BBF line arguments between children are

    X_l=y_w-y_w'+rho_l.

Hold relative coordinates and parameters fixed, pin one child anchor, and transform the other anchors to the tree-line arguments. The linear part is the reduced incidence matrix of a tree. Remove a leaf and expand in its row to prove inductively that its determinant has modulus one. The continuous Jacobian is one. Parent-anchor integration followed by child-relative integration makes the full map triangular.

On a spatial lattice, apply unimodularity to unshifted integer anchor differences and use the covariance's uniform real-offset lattice L1 bound. Fractional shifts are not an integer change of variables. Coefficient measures are handled by total variation and translation-covariant disintegration. This proves the integration assertion for interpolations with the stated pivot property.

[proved, remaining moment issue] If a line originally joins x,y but its actual argument is X=(1-t)x+t p-y, then

    x-y=X+t(x-p).                                    (10)

Only X is paid by the shifted covariance moment at that line's scale. The offset must be routed into the child. Unit anchor Jacobian does not bound a moment of x-y by a moment of X uniformly in x-p; set X=0 to see the issue. This is not a counterexample to (17). It identifies a necessary step in realizing it.

[read] The proposed endpoint note combines original endpoint trees, BBF connections, ancestor-first allocation, x^b R_m=R_(m-b)x^b, and first-contraction pivots. These are compatible ingredients. It does not give the complete paired incidence record for (10), its repeated versions, and all product-rule concatenations. The assertion that an erased quartic gain survives needs that record, not only a comparison of coordinate length scales.

### 5.1 A local routing correction that does survive

[proved, local bound] Once the exact operator identities have put a line into weighted-kernel form, keep that form intact under its remaining differential and evaluation operators. For one coordinate component set W_s(X)=X^b g_s(X), with b bounded. The full product rule is

    partial_X^a W_s(X)
      = sum_(d=0)^min(a,b) binom(a,d) (b)_d
                   X^(b-d) partial_X^(a-d) g_s(X).    (11)

Every summand has the same scale cost. The supplied moment bound gives

    integral |X|^(b-d) |partial^(a-d) g_s(X)| dX
       <= C_(a,b) 2^[s(a-b-1)].                       (12)

For a pointwise bound the corresponding exponent is s(1+a-b), with a fixed weaker decay power. The derivative and moment cancellations occur at the same scale, so both exponents are independent of d. The tensor version follows by the multi-index Leibniz formula. Since a,b are fixed, coefficient losses are fixed.

[proved] Relocating the whole W_s evaluates both its polynomial and covariance at the same argument X(t). Differentiating only g_s and retaining the old coordinate polynomial is a different operation and produces the extra offset in (10). Thus keeping the weighted kernel intact avoids that local mismatch in a representation where the operator acts on W_s. This is not permission to change an operator that originally acted only on g_s: the preceding exact identities must first justify the weighted representation. Section 5's anchor argument applies to its actual argument X(t). This does not prove a uniform global bound on a or b.

[read, primary-source check] Benfatto--Mastropietro, [Sections 3.10--3.11 and 3.15](https://arxiv.org/pdf/cond-mat/0004226), retain derivatives of weighted propagators in (3.101), prove the line estimate in (3.116), and state the paired scale inequality in (3.83). Their (3.85) telescopes paired scale lengths. The double-derivative case is described without its full case analysis in Section 3.11. These passages support keeping the weights and derivatives together; they are not invoked as a replacement-endpoint theorem for the present arbitrary-degree retarded action.

## 6. Remaining construction and weakest criterion established here

[assumed, sufficient routing conditions] A completion would provide one exact expansion of a fixed GN/BBF tree with all of the following properties.

1. Route moments through actual external port coordinates. In (10), process the internal offset inside its child before expanding that child's Taylor remainder. Do not pay an original displacement as an actual shifted-line moment.
2. Retain oriented source/sink records for all normal-form product-rule terms, concatenate consumed pairs by (3), and verify F_v>=z_v on each remainder edge.
3. On these same terms, prove fixed bounds on endpoint/line coordinate and derivative loads, total coefficient mass C^(sum p_i+n), and the first-contraction pivot property used in Section 5.
4. Compress unary idempotent repetitions before counting. Keep engineering derivatives in their assigned dimension or in the incidence record, never both. If shell support removes a putative long surviving engineering-derivative path, supply that support argument.

[proved implication] These conditions give (17) by (1)--(2), endpoint moments by (9), and anchor integration by Section 5, with constants independent of depth. Reference cancellation and corrected source estimates remain separate. Equation (4) shows the compatible local algebra and (6) its order-two bound before geometric realization. Neither proves that realization.

[proved, scale criterion] The termwise criterion with no extra constant is

    sum_v (F_v-z_v) ell_v >=0.

One may weaken its right side to -c(sum_i p_i+n), with fixed c, and absorb the loss into C^(sum_i p_i+n). Together with the coefficient-mass bound this still implies (17). The cutwise form F_v>=z_v is sufficient; an isolated deficient cut with an independently unbounded span is excluded under the nonvanishing assumptions in Section 3. These are scale criteria, not a completed corrected routing algorithm for the action.

[unresolved] Conditions 1--3 have not been established simultaneously for every nested two/four-field pattern in this run. No legal, nonvanishing single-channel tree violating C60 (17) has been found either. Thus the requested proof-or-counterexample alternative remains unresolved.

## 7. Checks and handoff

[numerical, finite exact-arithmetic controls] checks.py compares (4) on polynomial kernels, checks (5)--(6), verifies (1) on branched trees with explicit paired paths, and checks (7)--(8). These finite algebra controls do not test arbitrary GN trees or interacting RG.

[proved, consistency] With z=0, an oriented record gives a nonpositive exponent. A single pair gives 2^(s-t), unchanged by concatenation. Endpoint sources at scale zero give the elementary derivative scale gain. The free source remainder (N4) persists at zero interaction while D_O^int vanishes by definition. None of these facts identifies a central charge.

[record] No subagents or large computations were used. Parent STATE.md is outside the task's permitted write scope; local STATE.md supplies the handoff. No independent review is claimed.
