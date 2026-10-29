# C64. Weighted-line transposition and two-cut coverage

## 1. Status and conventions

[assumed, unresolved target] The requested global routing theorem is not established here. No legal nonvanishing tree with a deficient cut and an independently unbounded scale span is exhibited. In particular, this note does not prove or disprove C60 (17).

[proved below] A local spatial substitution realizes the concatenation rule of C61, including every Leibniz term. It applies when a child Taylor operator reaches a BBF line already carrying ancestor moments. It preserves the ancestor and child cut counts simultaneously. Section 5 embeds the substitution in a quartic-to-quadratic GN configuration, with an actual BBF tree line and its unexpanded determinant. This is a restricted spatial realization, not a construction for every nested pattern.

[proved below] Transposing a remainder onto the covariance alone, after attaching an ancestor moment to that covariance, is generally a different operation. Section 3 gives the nonzero correction and the exact local repair. This is an operator-order error in that proposed shortcut. It is not a counterexample to the carefully ordered routing or to the physical action.

[assumed, imported inputs] Use the source-free endpoint class and fixed-cutoff convention of C60 Section 2 and C61 Section 2. The covariance moment and Gram estimates are those of MARKED_COVARIANCE, used at fixed infrared cutoff after removing the volume regulators. Initial coefficient measures have the exponential terminal-tree norm supplied by TREE_NORM_EXTENSION. Finite initialization and the degree hierarchy are imported from C55 and C60. None of the matching, normalization, reference comparison, source, or thermal estimates is proved anew here.

[definition] A larger scale index denotes a shorter distance. Write t >= s >= h for a child moment source, the first contraction of a child field, and an outer derivative sink. The cut between t and s is called the child cut, and that between s and h the parent cut. In the actual GN tree each interval can contain compressed unary vertices. The cut counts below apply on every edge of the corresponding intervals.

[definition] All displayed one-dimensional spatial identities can be used in the Euclidean time coordinate. This avoids any question about differentiating lattice positions. Tensor products and directional derivatives give their spacetime versions. Fixed Fourier interpolation supplies the spatial derivative version when needed. Physical branch phases have already been factored out, as in the supplied infrared convention.

[definition] R_m subtracts the total Taylor polynomial of degree less than m at zero external deviation momentum. R_0 is the identity. Its transpose on test fields subtracts the Taylor jet at a chosen anchor. Quadratic and quartic orders are respectively two and one. Engineering derivatives belonging to a prescribed local endpoint remain in the engineering dimension. They are not included among the localization derivatives in this note.

## 2. Exact transposition before interpolation

[proved, scalar identity] Let K be a compactly supported finite signed or complex measure of a relative coordinate r, and let H be a smooth test function. Define the distribution R_m K by

    <R_m K,H>
      = integral K(dr) [H(r)-sum_(j<m) r^j H^(j)(0)/j!].       (1)

With the Fourier convention Khat(q)=integral exp(-iqr)K(dr), this subtracts the Taylor polynomial of Khat at q=0. Equation (1) is also the definition by duality, so no derivative of K is taken. Taylor's integral theorem gives, for m>=1,

    <R_m K,H>
      = integral K(dr) r^m integral_0^1
            (1-u)^(m-1)/(m-1)! H^(m)(ur) du.                  (2)

For the endpoint measures used below the integrations are finite. For general exponentially weighted endpoint measures the same identity follows first with a compact relative-coordinate cutoff and then by dominated convergence, provided H and its needed derivatives have polynomial growth.

[proved, consequence for a weighted line] Suppose a coordinate factor has already been assigned to a tree line in the exterior test function H. It is part of H in (1), and therefore belongs inside the derivative and evaluation in (2). For example, with an absolute anchor a, a covariance endpoint y, and

    r=x-a,    Z=x-y,    Z_u=a+ur-y,
    W(Z)=Z^b g_s(Z),

the correct first-order substitution is

    W(x-y)-W(a-y)
      = r integral_0^1 [partial_Z W](Z_u) du.                 (3)

Both the monomial and the covariance are evaluated at Z_u. This follows from (1), rather than from an assumed freedom to relocate a coordinate factor after taking an absolute value.

[proved, why this applies to a BBF line] A field selected for a BBF tree edge occurs at one endpoint of that edge and is removed from the residual determinant. Thus its coordinate occurs in that covariance, in any weight assigned to that covariance, and in the external coefficient measures. In the setting where no other exterior factor depends on that coordinate, (1) acts on the weighted covariance alone. The coefficient measure is paired with the test function and is not differentiated. Section 5 verifies the stated independence in a concrete nested remainder configuration.

[proved, determinant alternative] If a field is instead retained in the determinant, its coordinate labels one row or one column. Differentiating that field replaces that row or column by the differentiated row or column. Multilinearity gives this identity without expanding into Wick permutations. The derivative Gram bound follows by replacing the corresponding Fourier vector by its momentum-multiplied vector. If a row and a column are differentiated, their intersection entry has both derivatives, as required. No coordinate moment is assigned to a determinant entry in the local substitution studied below.

[assumed, limitation of the local identity] If an exterior factor depends on the same child coordinate, its derivatives must also be included. In particular, an ancestor Taylor field that depends on a surviving child port cannot be silently treated as constant. The restricted theorem below states its independence hypothesis explicitly. Removing this restriction is part of the global problem.

## 3. A nonzero error in covariance-only relocation

[proved] Compare (3) with the operation that keeps the weight at its original argument and applies the child difference only to g_s. Their difference is

    Z^b [g_s(Z)-g_s(Z_0)]
       - [Z^b g_s(Z)-Z_0^b g_s(Z_0)]
      = (Z_0^b-Z^b) g_s(Z_0),
    Z_0=a-y.                                                (4)

For b=1 this is -r g_s(Z_0). For b=2 it is

    -r(2 Z_0+r) g_s(Z_0).                                    (5)

These are not identically zero. For a scalar check, choose Z_0=1, r=1 and g(Z)=1+Z. At b=2, the covariance-only expression is 4, the whole-weight difference is 10, and their difference is -6, in agreement with (5).

[proved, fixed-order repair] The correct expression at b=2 is

    r integral_0^1 [2 Z_u g_s(Z_u)+Z_u^2 g_s'(Z_u)] du.       (6)

It needs at most a second moment of g_s and a first derivative of g_s, together with one child moment r. No third moment is introduced. For each u the child anchor still has coefficient one in Z_u. Hence the allowed shifted anchor integral applies to the same argument as the moment bound.

[proved, interpretation] Formula (4) does not say that the full original field difference was wrong. It distinguishes two different placements of the ancestor weight relative to the child remainder. Once the ancestor weight has been distributed on the BBF line in the kernel-distribution expression, duality requires (3). Keeping the weight outside the child operator changes that routed contribution. Other route terms cannot be dropped on the supposition that these two expressions agree.

## 4. The local two-cut theorem

[assumed, hypotheses of the theorem] Consider one BBF line at scale s carrying b ancestor moments. Their b derivative sinks are outside the child, at scale h. A child remainder of order c introduces c coordinate moments with sources at scale t and c derivatives on one endpoint of this BBF line. Assume

    0<=b<=2,    1<=c<=2,    t>=s>=h,                         (7)

and assume the line is the only exterior factor depending on that endpoint coordinate. Its pivot is in the same child. The derivative is transposed before the child Taylor integral is expanded, as in (2). The child moments can already have been distributed to specified descendant lines or endpoint edges; that distribution is held fixed here.

[proved, theorem] Every nonzero Leibniz term of this substitution has child cut count c and parent cut count b. It admits a paired incidence record with these counts, preserves the scale exponent, and has bounded moment, derivative and coefficient cost independent of t-s and s-h.

Proof. The exact Leibniz formula is

    partial_Z^c [Z^b g_s(Z)]
      = sum_(d=0)^min(c,b) binom(c,d) (b)_d
                     Z^(b-d) partial_Z^(c-d) g_s(Z),         (8)

where (b)_d=b!/(b-d)!. Each term has the same line scale factor. Relative to the unweighted L1 line norm, the covariance estimates give

    || Z^(b-d) partial^(c-d) g_s ||_1
       <= C_(b,c) 2^[s(c-b)] ||g_s||_(1,scale),              (9)

where ||g_s||_(1,scale) denotes the reference scale factor 2^(-s), not a lower bound on the actual norm. Equivalently the absolute bound is C_(b,c)2^[s(c-b-1)]. The pointwise exponent is s(1+c-b). Both statements are uniform in allowed real offsets.

The c child moments supply 2^(-ct), and the b outer derivatives supply 2^(bh). Thus the total extra scale factor is

    2^[c(s-t)+b(h-s)],                                      (10)

for every d. At an initial endpoint t=0, the corresponding moments are paid by its fixed terminal-tree norm rather than by a covariance moment.

Before consuming a zero, label the c child paths t to s and the b parent paths s to h. When the term in (8) differentiates d powers of Z, concatenate d child paths with d parent paths. The surviving paths are

    c-d paths t to s,
    d paths t to h,
    b-d paths s to h.                                       (11)

Every child cut is crossed (c-d)+d=c times. Every parent cut is crossed d+(b-d)=b times. Their total scale exponent is

    (c-d)(s-t)+d(h-t)+(b-d)(h-s)
       = c(s-t)+b(h-s),                                    (12)

which agrees with (10). The c child sources and b final sinks were never erased. The deleted source and sink were both at scale s.

For c,b<=2 the sum of the nonnegative Leibniz coefficients is at most seven. The Taylor measure in (2) has mass 1/c!<=1. Tensor components in two spacetime directions multiply this by a fixed constant. Residual line moments are at most two, and residual line derivatives from this endpoint are at most two. Thus the local coefficient mass and loads are bounded. This proves the theorem.

[proved, example c=1,b=2] The term Z_u^2 g_s'(Z_u) in (6) has one t-to-s path and two s-to-h paths. The term 2 Z_u g_s(Z_u) has one t-to-h path and one s-to-h path. Both give

    F_child=1,    F_parent=2,
    extra scale factor = 2^[-(t-s)-2(s-h)].                 (13)

The second term has no derivative on g_s. Deleting its child path on that account would give the false count F_child=0. Its source is still the factor r outside the Taylor integral, paired now with an outer derivative.

[proved, example c=b=2] The terms in partial^2(Z^2g) have d=0,1,2 and coefficients 1,4,2. At d=2 both paths run directly from t to h. Thus both cuts are crossed twice even though the covariance has neither a residual derivative nor a residual moment. The gain remains in the two child moments and two outer derivatives.

[proved, multiple independent substitutions] For any finite collection of substitutions satisfying the hypotheses and using distinct field coordinates, apply (11) at each substitution. Path concatenation preserves every cut count outside the substitution as well as its two adjacent intervals. The coefficient cost is at most 7^M in one coordinate for M substitutions. This follows by finite induction and does not require choosing a new pairing after taking spatial bounds.

[assumed, boundary of this corollary] It does not show that an arbitrary GN forest can be decomposed into such substitutions, that its coordinate dependencies are independent, or that its number of effective substitutions is bounded by C(sum_i p_i+n). Those are separate requirements for the global coefficient-mass bound.

## 5. A spatial quartic-to-quadratic configuration

[definition] Write a spacetime point as its coordinate label and suppress the fixed branch phases. Take two number-conserving quartic vertices with fields

    A: psi_R^+(a), psi_R^-(x), psi_L^+(c), psi_L^-(d),
    B: psi_R^+(y), psi_R^-(b), psi_L^+(e), psi_L^-(f).       (14)

Retain psi_R^+(a) and psi_R^-(b). At scale s contract the other six fields. Choose the BBF tree edge from psi_R^-(x) to psi_R^+(y). The remaining four L fields give a two-by-two determinant. There is no Pauli or charge obstruction. The output has degree two because 4+4-2(3)=2.

[definition] A is an initial retarded endpoint at scale t=0, with its quartic contact subtraction. B can be a prescribed local forward endpoint, so y=e=f=b. One allowed coefficient component of A has x=a+r, c=d=a with a fixed nonzero temporal displacement r. Its coefficient measure has finite terminal-tree norm, bounded by a fixed multiple of exp(kappa|r|). Multiplying by an arbitrarily small coefficient enforces the initial degree bound. Conjugate and reflected components can be added to the action. The computation here is for the stated labeled component in the absolute tree expansion, not for a signed sum over its symmetry partners. It is not asserted that the particular matched microscopic action has a nonzero coefficient in this direction.

[proved, child localization] Localize A at a. With c=d=a the quartic remainder is represented by replacing psi_R^-(x) with

    psi_R^-(x)-psi_R^-(a).                                 (15)

For general c,d, this is one term in the telescoping decomposition of the full quartic contact subtraction, with the L fields in that term left at their original coordinates. Thus the same configuration is available for nonsingular endpoint kernels as well.

[proved, BBF determinant retained] For c=d=a and e=f=b, write z=b-a. Apart from an overall Grassmann sign, the determinant is

    D_s(z,lambda) = det [ g_L(0)          lambda g_L(-z)
                         lambda g_L(z)  g_L(0)        ]
                  = g_L(0)^2-lambda^2 g_L(-z)g_L(z),        (16)

with lambda in [0,1] the two-child BBF interpolation variable. The standard two-child forest interpolation integrates dlambda. The R tree edge supplies g_R(x-y). In particular, D_s is independent of x and y as independent field coordinates before the endpoint measures are paired with the test function. The R fields at those two positions were selected for the tree and are absent from the determinant. No Wick-pairing sum was used.

[proved, parent moment routing] Apply the quadratic remainder to the two external fields at a,b, using a as pivot. Its second-order test-field remainder contains (b-a)^2 and two derivatives on the retained B field. Before transposing (15), use the exact path decomposition

    b-a = (x-a)+(y-x)+(b-y).                               (17)

This is a path in the augmented spatial tree: an internal connection in A, the selected BBF edge, and an internal connection in B. For local B the last summand is zero. In the square of (17), retain the route term (y-x)^2. It has coefficient one. The other route terms remain part of the complete expansion and are not estimated or discarded in this computation.

[proved, exact spatial realization] In this selected route term the exterior function of x is now

    W(x-y)=(x-y)^2 g_R(x-y).                               (18)

The root Taylor field depends only on a,b and its own parameter. The determinant (16) is independent of x. The child subtraction (15), transposed onto this exterior function, therefore gives exactly

    r integral_0^1 [2 Z_u g_R(Z_u)+Z_u^2 g_R'(Z_u)] du,
    Z_u=a+ur-y.                                           (19)

The parent Taylor parameter is independent of u. Its field remains evaluated at a+v(b-a), with weight 1-v and with two external derivatives. Differentiating with respect to x does not differentiate that field. This verifies the independence hypothesis in Section 4, rather than assuming it for a general nested cluster.

[proved, simultaneous cut counts] The child endpoint supplies the source r at t=0. Both parent moments were placed on the line at s. The child derivative is first contracted at s, and the two parent derivatives have sink scale h. Equations (11)-(13) apply with c=1 and b=2. Thus the required quartic child cut has F=1 and the quadratic parent cut has F=2 on both terms of (19). This statement is uniform in the two scale spans for every nonzero contribution. It does not rely on assuming that the graph survives every possible support restriction.

[proved, anchor and moment compatibility] Write y=b for local B and pin a. For fixed r,u, the actual tree-line argument is Z_u=ur-b. The continuous time Jacobian for b to Z_u has modulus one. The moment in (19) is a power of that same Z_u. On the spatial lattice a shift is handled by the uniform shifted L1 estimate, without treating fractional points as integer variables. The determinant is bounded by its Gram estimate, independent of the anchor coordinates. The endpoint moment r is paid by the terminal-tree norm of A. Consequently the weighted-line bound, determinant bound, anchor integral and cut record all concern the same selected spatial representation.

[proved, nonzero coordinate integrand] This configuration is not a formal Grassmann monomial forced to vanish. For example, on the time axis a nontrivial shell covariance g is a decaying smooth function. The function W(z)=z^2g(z) is nonzero and cannot have a nonzero period. Hence W(r-z)-W(-z) is not identically zero for r!=0. With a compact smooth shell multiplier these functions are real analytic on the time axis. The determinant integrated over lambda is

    g_L(0)^2-(1/3)g_L(-z)g_L(z).                           (20)

If g_L(0)!=0, (20) tends to a nonzero constant as |z| tends to infinity. If g_L(0)=0 and the time-axis covariance is nontrivial, (20) is a nonzero analytic function. The product of two nonzero analytic functions on the line is not identically zero. This proves nonvanishing of the coordinate factor in (19) after its u and lambda integrals.

[proved, covariance qualification for this check] Choose a nonnegative shell cutoff even in q0, as permitted by the fixed-cutoff convention, and real free dispersion. The time-axis covariance is nontrivial. Indeed

    partial_tau g_s(0,0)
       = - integral d q0 d q1/(2pi)^2
                    q0^2 f_s(q0,q1)/(q0^2+e(q1)^2) < 0,   (21)

where the imaginary term is odd in q0. The shell is nonzero on a set with q0!=0. Compact frequency support supplies the analyticity used above. Small normalized perturbations are unnecessary for this nonvanishing test; the free covariance is an allowed input to the prescribed-history expansion.

[assumed, essential qualification] Nonzero coordinate integrand is weaker than a nonzero fully projected low-energy kernel at every root scale. A counterexample based on a deficient cut would additionally require survival of the shell support for an independently unbounded span and a nonzero final contribution. This note makes no such counterexample claim. The counts in (19) satisfy the desired inequalities anyway. The point of this configuration is to realize the derivative-consuming term and its repair inside a legal spatial contraction, not to refute (17).

## 6. What remains for the full induction

[proved, restricted conclusion] Sections 2-5 remove the local ambiguity about whether a child derivative should act on an ancestor line weight. Once that weight has been routed in the kernel-distribution expression, it must. The two-cut calculation then gives the cut counts without borrowing a gain from a different expansion.

[assumed, unresolved global invariant] To finish C64 one must construct, for every reduced GN/BBF tree, an exact representation with all of the following properties at once.

1. Incoming moments at a child are external-port coordinate differences to which C61's normal form applies before its Taylor integral is expanded.
2. Every transposed child operator either acts on a weighted covariance or is accounted for in the normal form of the remaining external test fields. All shared-coordinate terms are included.
3. The sources of consumed moments and sinks of consumed derivatives are concatenated at their actual common scale. The resulting record crosses every required remainder edge at least z_v times.
4. Moment and derivative loads stay bounded, all pivots satisfy the first-contraction anchor condition, and the total coefficient mass is at most C^(sum_i p_i+n).

[proved, why the present theorem is insufficient] A quartic child can retain two or more fields in an ancestor. In that case an ancestor Taylor test field can depend on several of its ports, and differentiating a port need not act only on one BBF line. Equation (19) verified the opposite situation by explicit field labels. Repeating its line-only rule in the multi-port situation would omit derivatives of the exterior test function. C61's polynomial normal form supplies algebraic counts for those terms, but an assignment of their spatial moments, relocated arguments and path records in one bounded-load representation is still required.

[proved, operation order to retain] One can avoid the local offset in C61 (10) by routing before transposition, as (17)-(19) demonstrate. This does not license replacing an already expanded shifted-line expression by the derivative of a weighted line. The equality must be established at the level of (1). Conversely, rerouting after a Taylor parameter has entered a displacement can create parameter-dependent child offsets. Applying the parameter-independent normal form to those offsets without another operator identity is not justified by (1)-(2).

[assumed, not a no-go statement] These observations do not exclude an exact induction combining distributional transposition, the polynomial normal form and suitable pivot choices. They identify the remaining interface. I did not construct that induction or find a legal counterexample to it.

[proved, conditional implication already established by C61] If the four global properties hold, the retained incidence identity gives

    sum_f a_f s(f)-sum_l b_l s(l)
       = -sum_v F_v ell_v <= -sum_v z_v ell_v.             (22)

Multiplying by the coefficients and summing proves C60 (17). The local theorem supplies a class of valid substitutions in such a proof, but is not a substitute for its global hypotheses.

## 7. Checks, provenance and scope

[numerical, finite exact-arithmetic checks] checks.py verifies the transposition identity, the full polynomial Leibniz formula, the two-cut records for all b,c in the specified range, and the error (4). All polynomial integrations are performed with rational arithmetic. The checks also integrate the two-by-two determinant polynomial in its BBF parameter. They are controls of the displayed identities, not a numerical RG or an exhaustive GN routing search.

[proved, consistency checks] Setting b=0 leaves only the child gain 2^[c(s-t)]. Setting c=0 in the algebraic formulas leaves only the parent gain 2^[b(h-s)]. At d=min(b,c), no nonexistent derivative or moment is retained. When s=t or h=s the corresponding scale factor is one and the same path accounting still holds. With c=1,b=2, the apparently derivative-free term retains the quartic cut gain through its concatenated path, as shown in (13).

[assumed, source boundary] The free lattice-source subtraction in C61 remains necessary. The present source-free calculation proves no density or pair exponent, history comparison, scalar thermal coefficient or central charge.

[read, primary-source context] Benfatto and Mastropietro apply a child operator to a weighted propagator in (3.79), and organize the scale pairing in (3.83)-(3.85). Their discussion of double derivatives does not spell out the full case analysis. This is consistent with the local operation above and is not invoked as an arbitrary-endpoint global theorem. See [the primary paper, Sections 3.9-3.11](https://arxiv.org/pdf/cond-mat/0004226).

[record] C55, C60 and C61 were read in full. Their settled finite-initialization results, source correction and reduced-tree identity are retained. The current endpoint note, covariance note and proof-obligation ledger were consulted. All writes stay in work_C64. Local STATE.md replaces a parent-memory edit under the task's write restriction. No delegated agents, independent review or large computations were used.
