# C70. A quartic child acting on the full exterior field

## 1. Status and conventions

[assumed, unresolved target] This note does not establish C60 (17) for the full single-channel routing, and does not exhibit a legal nonvanishing tree violating it. The missing assertion is a global bounded-load representation. A local deficient-cut counterexample cannot be obtained merely by including the previously omitted exterior-field derivative. That derivative has a favorable scale factor, as proved below.

[proved below] The spatial transposition formula extends to a quartic child with several surviving ports, including derivatives of affine ancestor test fields, weighted covariances, and determinant rows. Every term preserves existing source/sink paths and creates a path covering the child cut, or concatenates it with an existing path. This statement concerns one substitution with its incoming representation fixed. Its constants depend on the number of factors sharing the child coordinates and on their existing loads.

[proved below] A survival-adapted decomposition of a quartic remainder factors out the momenta of its surviving nonpivot ports. It proves that a later quartic contact subtraction is the identity on those components. For a later quadratic remainder, one surviving momentum lowers the remaining Taylor order from two to one. The statement holds before incoming coordinate moments are assigned. Section 6 explains why these two results are not a completed simultaneous construction.

[assumed, imported input] Use C60's fixed-cutoff convention and C61's source-free endpoint class. Covariance moment and Gram estimates are imported from MARKED_COVARIANCE. Endpoint coefficient measures have the terminal-tree norm of TREE_NORM_EXTENSION. C61's normal-form identity and cut identity, and C64's weighted-line substitution, are retained. No matching, source, reference-comparison or thermal assertion is added.

[definition] Larger scale indices mean shorter distances. A child moment has source scale t. Let s be the first ancestor scale outside the child. A covariance or exterior field reached by the child operator has contraction scale k<=s, with k=h for a root field. An initial endpoint moment has t=0. Scale records below are oriented toward the root.

[definition] Work first in one continuous coordinate, which can be Euclidean time. The extension to two spacetime directions uses directional derivatives. Spatial interpolation on the lattice has the meaning already specified in C64. All Taylor and BBF parameters are held fixed when differentiating coordinate arguments. Initial measures are paired with the exterior function and are never differentiated.

## 2. Full exterior transposition

[definition] A quartic child has coordinates x_1,x_2,x_3,x_4, with pivot a=x_4 and relative coordinates r_i=x_i-a for i=1,2,3. Let H be the entire exterior function of these ports, including all weights and ancestor test fields. Other coordinates and parameters are suppressed. Set

    x_i(u)=a+u r_i,  i=1,2,3,
    x_4(u)=a.

[proved] The transpose of contact subtraction is

    H(x_1,x_2,x_3,a)-H(a,a,a,a)
      = sum_(i=1)^3 r_i integral_0^1
                     partial_(x_i) H(x_1(u),x_2(u),x_3(u),a) du.       (1)

This is the fundamental theorem of calculus applied to H(x(u)). In particular, it is valid when H is Grassmann valued, component by component. Integrating (1) against a compactly supported coefficient measure is legitimate. For an exponentially weighted finite measure, polynomial growth bounds on the required derivatives allow the same conclusion by dominated convergence.

[proved, scope of the derivative] If a surviving field in H depends on x_i, its derivative occurs in (1). Neither survival of the original port nor assignment of other moments to a line makes that factor independent of x_i. Formula (1) acts on the whole exterior function.

[definition] A useful finite exterior representation has the form

    H(x)= product_e [Z_e(x)^b_e g_e(Z_e(x))]
          det G(x) product_(j=1)^q [partial^a_j phi_j](Y_j(x)).       (2)

Keep the Grassmann product in its original order. The Z_e are affine line arguments and the Y_j are affine field arguments. The coefficients of the affine maps may depend on fixed Taylor parameters. The matrix G is a BBF Gram matrix with affine field positions in its rows and columns. Coordinate-independent factors can be multiplied into (2).

[definition] Write

    L_ei=partial_(x_i) Z_e,
    A_ji=partial_(x_i) Y_j.

[proved] For each i the derivative in (1) is the sum of the following terms, with every factor evaluated at the same x(u).

1. A covariance derivative on e replaces its weighted factor by

       L_ei Z_e^b_e g'_e(Z_e).

2. A consumed weight on e replaces its weighted factor by

       L_ei b_e Z_e^(b_e-1) g_e(Z_e),

   and is absent when b_e=0.

3. A derivative on exterior field j replaces that factor by

       A_ji [partial^(a_j+1) phi_j](Y_j).

4. A determinant derivative replaces one affected row or column by its coordinate derivative, with the affine chain-rule coefficient of that position.

The proof is the product rule, the chain rule, and determinant multilinearity. Coordinate derivatives are even, so differentiating the ordered Grassmann product introduces no permutation signs. A coordinate appearing in both a row and a column contributes both replacements. Their intersection entry has the sum of its two coordinate derivatives. No expansion into Wick pairings is used.

[proved, ancestor Taylor fields are included] For example, suppose an ancestor factor is

    Phi_v(x)=partial^a phi(Y_v(x)),
    Y_v(x)=(1-v) U(x)+v V(x),

where U and V are affine functions of the surviving child ports. Then

    partial_(x_i) Phi_v
      = [(1-v) partial_(x_i)U+v partial_(x_i)V]
                     partial^(a+1)phi(Y_v).                         (3)

The ancestor integration in v is retained. Its coefficient is bounded whenever the affine coefficients of U,V are bounded, including at v=0 and v=1. Compositions of convex affine relocations remain convex affine relocations. Thus (3) needs no division by an interpolation parameter. It also covers a shared pivot, rather than declaring that pivot constant without checking its dependence.

[proved, local coefficient bound] At fixed parameters the sum of the scalar coefficients introduced by partial_(x_i), in absolute value, is bounded by

    sum_e |L_ei|(1+b_e)
       + sum_j |A_ji|
       + sum_rows |partial_(x_i) Y_row|
       + sum_columns |partial_(x_i) Y_column|.                       (4)

For b_e<=2, every line contribution is at most 3|L_ei|. Summing over the three relative coordinates and two spacetime directions gives another fixed factor. Taylor integration has mass one. Formula (4) is a bound for one operation, not a depth-independent bound on its product over the forest.

## 3. Cut records for all exterior derivative terms

[assumed, hypothesis of this local theorem] Before the substitution, every existing coordinate weight is accompanied by its oriented incidence path from its actual source scale to its derivative sink. Each new child displacement r_i has been assigned an admissible source at t. A covariance reached by the operation is outside the child, at scale k<=s. Every changed argument satisfies the first-contraction pivot condition of C61. No weights are assigned to determinant entries. These are properties of the incoming representation; this section does not construct it for an arbitrary forest.

[proved, full exterior cut theorem] Every term in Section 2 has at least one path through the child edge. Existing cut counts are retained. If a weight is consumed, concatenation retains its old cut counts. If an exterior field is differentiated, the new path ends at that field's actual later contraction scale, not at s.

Proof. First consider a covariance derivative at k. Its new moment and derivative have scale factor 2^(k-t). Add the path t to k. This crosses the child edge and possibly further ancestor edges. Existing weights and their paths are unchanged. A determinant row or column derivative has the same factor by the derivative Gram bound, so the same construction applies.

For an exterior field j first contracting at k_j<=s, (3) adds one derivative there. Its factor is 2^(k_j-t), so add a path t to k_j. If it is a root field, k_j=h. This is a legitimate path toward the root, with at least the child gain 2^(s-t). A later sink gives an additional favorable factor 2^(k_j-s).

Now suppose the derivative consumes a weight on a line at k. Label the b_e coordinate factors individually, even when they have the same coordinate component. The factor b_e in the product rule sums the choices of the consumed label. Let that label's old path run from k to an ancestor sink h_e. Add a temporary path t to k, concatenate it with k to h_e, and retain t to h_e. The remaining b_e-1 weight paths remain. The deletion was at a common scale k, so

    (k-t)+(h_e-k)=h_e-t.                                    (5)

The scale exponent and every preexisting cut count are preserved. The added child path crosses the child edge. This exhausts the terms in Section 2 and proves the theorem.

[proved, scalar specialization] With one line of scale s carrying b moments, all of their sinks at h, consider

    H(x)=Z(x)^b g_s(Z(x)) partial^b phi(Y(x)),
    partial_x Z=1,  partial_x Y=alpha.

One child coordinate r produces the three terms

    r integral_0^1 {
       b Z_u^(b-1) g_s(Z_u) partial^b phi(Y_u)
       + Z_u^b g'_s(Z_u) partial^b phi(Y_u)
       + alpha Z_u^b g_s(Z_u) partial^(b+1) phi(Y_u)
    } du.                                                    (6)

The first two have additional scale factor

    2^[(s-t)+b(h-s)].                                         (7)

Their records are the ones established in C64. The exterior derivative term has factor

    2^[(h-t)+b(h-s)]
       = 2^(h-s) 2^[(s-t)+b(h-s)].                            (8)

Its child count is one and its parent count is b+1. Hence it is stronger in scales than the requested counts one and b. At b=2 it contains a third exterior derivative. Its favorable scale factor does not prove a bound of two on the derivative load.

[proved, extension to different sinks] If the b weights have different ancestor sinks h_j, replace b(h-s) by sum_j(h_j-s). Consuming the j-th weight concatenates t to s to h_j. Differentiating an exterior factor with sink k adds t to k. The proof does not require equal sink scales or a chain with only two edges.

[proved, consistency] With no exterior dependence, alpha=0, (6) reduces to C64. With b=0 only a covariance derivative and an exterior derivative remain. If h=s, (7) and (8) coincide. Zero-length spans give factor one. These checks agree with the source/sink orientation.

## 4. Local spatial bounds and their limits

[proved conditional on the incoming geometry] Each residual polynomial in Section 2 is evaluated at the actual argument of its covariance. Thus the fixed-order moment bound applies to the weighted covariance itself. For instance,

    integral |Z|^(b-d) |partial^(c-d) g_k(Z)| dZ
        <= C_(b,c) 2^[k(c-b-1)],                              (9)

for c=1 and d=0,1 when the corresponding term is present. The exterior derivative term uses the unaltered line moment bound and an additional derivative at its exterior sink. There is no unmatched original displacement attached to a relocated covariance.

[proved conditional on the incoming geometry] A permitted affine interpolation at the first contraction of a field keeps coefficient one on that field's child anchor. For a BBF line joining children w,w', its argument is

    y_w-y_w'+rho(relative coordinates,parameters).

The child transposition in (1) preserves this property if its pivot and moved ports lie in the prescribed child of that first contraction. The reduced incidence determinant is then one in absolute value, by C61's argument. The discrete spatial sum uses its uniform real-offset bound. Differentiating an exterior test field changes neither the coefficient measures nor the anchor Jacobian.

[proved, one-event load statement] One substitution of a quartic remainder adds only one coordinate moment. It either adds one derivative to one covariance endpoint, determinant row/column, or exterior field, or consumes one existing coordinate weight. For fixed incoming loads and a fixed number of coordinate-dependent factors, all required moments, derivative orders and coefficient masses are bounded independently of scale spans.

[assumed, limitation] Iterating this statement with its maximal incoming derivative order does not give a uniform load bound. A bound that grows by one at every event is insufficient. Similarly, multiplying (4) without controlling the number of factors sharing each active port would not prove C^(sum_i p_i+n). Neither point is repaired by the favorable exponent in (8).

## 5. A decomposition adapted to surviving ports

[definition] Choose the quartic pivot to be a longest-surviving port, with a deterministic tie rule. Field survival sets are nested toward the root. Hence this pivot belongs to every nonempty surviving subset of the child's four ports. Remove its momentum using conservation. Separate the remaining independent child momenta into q_C for ports contracted before a specified ancestor and q_S for ports surviving there. The pivot is in the surviving set whenever that set is nonempty.

[proved] For a smooth child kernel K,

    R_1 K(q_C,q_S)
      = K(q_C,q_S)-K(q_C,0)
        + K(q_C,0)-K(0,0).                                 (10)

The first summand has the exact representation

    sum_(j in S) q_j integral_0^1
                     partial_(q_j) K(q_C,u q_S) du.          (11)

The second summand has the analogous representation using q_C and K(u q_C,0). Equation (10) is an addition and subtraction, and (11) follows by differentiating K(q_C,u q_S) with respect to u. It uses at most three momentum components per spacetime direction, since there are only four child ports.

[proved, effect at a quartic ancestor] Retain the factor q_j in each term of (11) until after the Gaussian contractions between this child and the ancestor. The port j is uncontracted there. Consequently its external momentum factor can be taken outside the contraction integral. The resulting ancestor kernel has the form q_j F_j(q_ext). Its contact localization is zero, since all external deviation momenta vanish at that contact. Thus

    R_1[q_j F_j]=q_j F_j.                                  (12)

This does not use a signed cancellation between Wick pairings. It holds for each contraction expression in which the factor is kept outside. The coefficient F_j can be nonlocal and can retain Taylor parameters.

[proved, effect at a quadratic ancestor] If both quadratic external ports survive from the quartic child, the independent external momentum is q. A surviving homogeneous factor P_a(q) satisfies

    R_2[P_a F]=P_a R_max(0,2-a) F.                         (13)

For a=1, only one more Taylor order is required. For a>=2, no additional Taylor order is required. This is the b=0 instance of C61's identity. If only one child port survives, it is the pivot and S is empty, so (11) supplies no incorrectly surviving derivative.

[proved, restricted consequence] Equations (10)-(13) avoid repeatedly differentiating the same surviving field through quartic ancestors in this decomposition before incoming moments are assigned. They explain how a term that would be written as a child derivative acting on an already expanded ancestor Taylor field can instead be organized with the ancestor contact subtraction already removed. For the mixed quartic/quadratic pattern, (13) gives the corresponding reduction of Taylor order.

[assumed, restriction] These formulas alone do not order the normalizations relative to the spatial moment allocation. In particular, the coordinate polynomial of an ancestor acts by momentum differentiation on (10)-(11), including their explicit momentum factors and their parameter-dependent kernels. The next section records this issue without assuming that (12) resolves it.

## 6. The remaining simultaneous step

[proved, imported normal form] With notation from C61, an incoming moment order b, inherited homogeneous derivative degree a, and d consumed derivatives give

    r=max(0,m-b-a+d),
    J=a+b-d+r>=m,
    B=b-d+r<=2,  when m,b<=2.                              (14)

For a quartic child m=1. This is an operator identity before the remaining Taylor integral is expanded. Sections 2-3 show how derivatives on the full exterior function can carry its paths. They do not identify every incoming moment in an iterated affine representation with the unexpanded operator on the left of (14).

[proved, explicit cancellation needed before absolute values] Even in one variable, with R_1 K(q)=q integral_0^1 K'(u q)du,

    partial_q^2 R_1 K(q)=K''(q),                            (15)

whereas direct differentiation of the Taylor representation gives

    2 integral_0^1 u K''(u q)du
      + q integral_0^1 u^2 K'''(u q)du.                    (16)

The third moment in the second term is removed only by combining the terms. Indeed the integrand sum is partial_u[u^2 K''(u q)], so its integral is K''(q). Thus applying a correct local product rule to previously expanded fields does not by itself enforce the B<=2 bound in (14).

[proved, multi-variable version] If R_1 K(q)=K(q)-K(0), any two momentum derivatives give partial_i partial_j K(q). Expanding the radial Taylor formula first gives terms with up to three derivatives of K. Differentiating the finite Taylor formula and integrating its parameter derivative recovers the same second derivative. The one-variable identity is therefore an actual instance of the multiport issue, not an effect specific to the choice of spacetime direction.

[assumed, unresolved construction] A proof of C60 (17) still has to supply one terminating rewriting procedure with the following joint properties. It must expose surviving factors as in (11), apply the normal form (14) before taking absolute values, route the resulting moments on actual shifted BBF lines and internal endpoint edges, and retain the source/sink labels through every cancellation such as (15)-(16). The procedure must bound the loads and its total coefficient mass, and preserve the first-contraction pivot property. The local extension here establishes cut coverage once such an incoming representation has been supplied. It does not establish that representation inductively.

[proved, why a counterexample has not been obtained] The extra term in (6) has an additional favorable scale gain. A third derivative in that expression is therefore not a deficient-cut counterexample. Nor is the third derivative in (16) such a counterexample, since (15) removes it before estimates. To refute C60 (17), one would need a legal nonvanishing projected tree with a scale deficit that cannot be removed by an exact bounded-load representation, and an independently unbounded span. No such tree is supplied here.

[assumed, final status] The requested proof-or-counterexample alternative remains unresolved. The proved advance is the full local exterior-field cut calculation and the survival-adapted quartic decomposition. The remaining issue is global compatibility of that decomposition with moment routing, not a missing local scale factor for the exterior derivative.

## 7. Checks and provenance

[numerical, finite exact-arithmetic controls] checks.py tests the full affine product rule with up to three child relative coordinates and up to four exterior factors, verifies the transposed integral identity, checks the different-sink cut records, and checks the cancellation (15)-(16). It also verifies (10)-(13) on finite polynomial kernels. These are algebra controls, not a search over legal GN trees or a numerical RG computation. The polynomial tests are scalar tests of operator identities and are not claimed as nonvanishing projected Grassmann trees.

[assumed, primary-source context] Benfatto and Mastropietro discuss weighted-line differentiation in Section 3.9, the pivot issue in Sections 3.3-3.4, and paired scales in Sections 3.10-3.11 of [their primary paper](https://arxiv.org/html/cond-mat/0004226). Those passages motivated checking the survival-adapted split. They are not used as a theorem for arbitrary replacement endpoints.

[record] C60, C61 and C64 were read, together with the task backgrounds, the current endpoint-majorant note, the covariance input and the obligations ledger. No subagents, cluster computations, shell network commands, or writes outside work_C70 were used. Local STATE.md provides the handoff under the task's write restriction. No independent review is claimed.
