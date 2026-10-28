# From normalized marked trees to the required history bounds

Status: **[proved transfer lemmas; the remaining covariance/tree input is
explicit]**. This note does not label the endpoint-only RG extension as
proved. It supplies the nonlinear calculus and exponent allocation that
would otherwise be hidden inside a claim of “continuity of GN trees.”

Author: OpenAI GPT-6-family Codex vacancy worker, 2026-10-28. AI-assisted;
no human review. The alternatives worker independently identified the
coefficient-difference term below and checked the published first-order
cancellations.

## 1. History coordinates and a sufficient marked estimate

Fix a root scale h. Write the history in coordinates

    u=g_h,        d_k=g_k-g_h,       n_k=nu_k,
    h<=k<=0,     d_h=0,
    w_(h,k)=2^[-a(k-h)],             a>0.                 (1)

At fixed physical Fermi centers, write the normalized beta function as

    B_h(g,nu)=F_h(u)+Psi_h(u,d,n),
    Psi_h(u,0,0)=0.                                      (2)

All norms below include the finite branch labels. Assume the following
bounds on a convex complex neighborhood of the provisional histories:

    |F_h(g_0)|<=C epsilon 2^(b h),
    ||D F_h(u)||<=C 2^(b h),                             (3)

and, for z_k either d_k or n_k,

    ||partial_(z_k) Psi_h||<=C eta w_(h,k),
    ||partial_u partial_(z_k) Psi_h||<=C w_(h,k).         (4)

Here eta=sqrt(epsilon), and constants are uniform in the number of
scales. For the second estimate, the u derivative moves the common
marginal history simultaneously, not one independently selected scale.
Choose a>=4theta and b>=2theta.

The analytic history domain is first extended formally so u and d are
independent, with the initial irrelevant endpoint kernels held fixed.
The physical fixed-g_0 relation `d_0=g_0-u` is imposed only afterwards.
Thus the common-history derivative does not covertly differentiate
the microscopic initial endpoint data.

**Lemma 1 [proved].** Equations (2)--(4) imply (G1) and the corrected
(G2) of COUNTERTERM_INDUCTION.

Proof. Integrate from (d,n)=(0,0) to the desired history. The first
bound of (4) gives

    |Psi_h(u,d,n)|<=C eta sum_k w_(h,k)(|d_k|+|n_k|).

Equation (3) bounds F_h(u) by
`C 2^(b h)[epsilon+|u-g_0|]`. This proves (G1).
For two histories, first change d,n while holding u fixed; this gives
the eta-weighted differences in (G2). Then change u with d,n fixed.
Since Psi_h(u,0,0)=0, the second bound of (4), integrated along the
segment from zero to d,n, gives

    ||partial_u Psi_h(u,d,n)||
                 <=C sum_k w_(h,k)(|d_k|+|n_k|).

The u difference therefore costs
`C|delta g_h|[Q_h(g,nu)+Q_h(tilde g,tilde nu)]`, with harmless enlargement
to include both histories. Finally use (3) on the F term. This is
precisely the corrected (G2). No Fermi-center derivative occurs. QED.

## 2. A more permissive ultraviolet history bound also closes

Some covariance conventions leave a first-order, scale-decaying
dependence on the marginal history. It need not be forced into (4).
One may replace the first bound in (4) for d_k by

    ||partial_(d_k) Psi_h||
                     <=C[eta+2^(b h)] w_(h,k),           (5)

while retaining `C eta w_(h,k)` for n_k. Make the corresponding
replacement in (G1)--(G2). The additional term in the absolute marginal
bound is

    C 2^(b h) sum_k w_(h,k)|g_k-g_h|
       <=C' 2^(b h) sum_(j>h)|b_j|.                     (6)

The analogous difference is bounded by the right side with
`|delta b_j|`. Since b>=2theta, both are absorbed by increasing the
constant of the strictly triangular K0 in COUNTERTERM_INDUCTION.
The factorial resolvent proof applies without any new smallness
requirement other than reducing epsilon after the new constant is fixed.

The relevant beta receives the same additional marginal-history term;
after the causal marginal estimate it is O(epsilon 2^(theta h)), and
its difference is O(eta ||delta nu||_theta 2^(theta h)). Thus the same
relevant contraction closes. The n_k estimate must still retain eta:
an arbitrary unsmall coefficient multiplying nu_h would not be removed
by a triangular marginal resolvent.

This variant means that exact confinement of every ultraviolet
first-order term to g_h is useful but unnecessary. Ultraviolet decay
and its marked version suffice.

## 3. Where the extra small factor comes from

Consider a convergent renormalized tree expansion on a complex domain
where all local quartic and relevant amplitudes have magnitude at most
2eta, and velocities stay a fixed distance from zero. Keep the physical
initial irrelevant endpoints fixed. Suppose these endpoints have total
size at most C epsilon^2 with the degree hierarchy established in
DEGREE_REPARAMETERIZATION. Recall epsilon^2=eta^4.

Separate the terms into:

1. The reference contribution at constant marginal history.
2. One-local-endpoint ultraviolet or curvature contributions.
3. Contributions with at least two local interaction endpoints.
4. Contributions with at least one initial irrelevant endpoint.

For a pure-local history-dependent tree in class 3, marking a local
amplitude removes at most one power of eta, leaving at least eta.
A second mark can remove that remaining power, leaving a bounded
coefficient. A velocity mark does not remove an interaction endpoint.
For class 4, Cauchy on a history disk of radius proportional to eta
costs at most eta^-1 per derivative: one and two marks leave respectively
eta^3 and eta^2, both better than the eta and 1 needed in (4).

These statements remain valid for sums, not just monomials. For example,
if the endpoint-count majorant is `sum_(n>=2) C^n eta^n`, then one and
two amplitude marks are bounded by

    sum_(n>=2) n C^n eta^(n-1)<=C' eta,
    sum_(n>=2) n(n-1) C^n eta^(n-2)<=C''.                (7)

Field-degree multiplicities from a mark are handled by the same field
radius reserve as in DEGREE_REPARAMETERIZATION. A second mark spends
another fixed part of that reserve once, not once per scale.

For the scale weight, restrict to a tree with the marked dependency at
scale k. If every renormalized non-root edge has gain at least delta>0,
reserve a<delta along the root-to-mark path. Its product is
`2^[-a(k-h)]`. The unspent dimensions still sum the tree scales.
The second, common-history mark can be estimated by a Cauchy derivative
in a single complex parameter shifting every marginal coordinate u.
Uniform convergence on the corresponding history polydisc bounds that
derivative independently of the number of scales. It therefore need
not consume a second root-to-k short-memory factor.

**Precise limitation.** This is a proof of (4) once the normalized tree
representation retains the stated endpoint counting and analytic
covariance bounds under those one- and two-mark operations. An absolute
unmarked tree bound alone does not establish those facts. In particular,
the running normalization ratios are functions of the history, not
independent arbitrary constants.

## 4. One-endpoint terms and the lattice/reference comparison

The reference covariance is odd in momentum, so its local quartic
tadpole vanishes. A one-local-quartic self-energy has no external
momentum dependence, hence contributes no velocity or wave-function
derivative. After the canonical relevant propagation and the carried
local quartic term are removed from the definition of the beta function,
the purely local reference remainder has no degree-one beta term.

There is also a scale-local observation. A one-local-quartic tadpole
cannot propagate through an intermediate renormalization operator:
R_2 annihilates its constant quadratic polynomial. R_4 annihilates an
uncontracted local quartic. Thus a genuine one-local-endpoint tree with
an intervening renormalized cluster vanishes. On the lattice a surviving
current-scale tadpole can be linear in the quartic coupling, but its
difference from the odd reference tadpole has a curvature/UV gain.
Such a term belongs to F_h or the ultraviolet allowance (5); it must
not be assigned the eta factor in (4) without this separation.

For the initial irrelevant kernels, one endpoint can self-contract to
a local quadratic or quartic term. This contribution is not zero, but
it belongs to class 4, is O(epsilon^2), and reaches the initial scale.
The root-to-initial path gives the ultraviolet decay. Keeping these
terms is necessary for the correct counterterm and marginal matching.

## 5. Polynomial losses along anomalous normalization paths

A separate elementary estimate explains why path lengths introduced
by a fixed number of marks need not destroy short memory. Suppose a
path has total scale length L>=0 and unmarked dimensional gain
`2^(-delta L)`. If r fixed marks produce at most C_r(1+L)^r times that
bound, then for every 0<a<delta,

    (1+L)^r 2^(-delta L)<=C_(r,delta-a) 2^(-a L).        (8)

This follows by maximizing `(1+L)^r exp[-(delta-a)(log 2)L]`.
Only r=1,2 is required here.

For example, a product of single-step normalization ratios with
`|log rho_j|<=C eta^2` has magnitude at most exp(C eta^2 L).
If the first and second marked logarithms are uniformly bounded by
C eta and C respectively, direct differentiation of the product costs
at most C eta L for one mark and C(L+eta^2 L^2) for two. After reducing
eta so C eta^2/log 2<delta-a, (8) absorbs those losses. This uses the
ratios and never sets Z_h or its anomalous flow equal to one.

The example is not a substitute for proving the marked-logarithm
bounds of the actual normalized covariance. It identifies exactly
what that proof must retain: the number of marks is fixed, the losses
are polynomial in the scale span, and the anomalous exponent is small.

## 6. Primary-source inputs and the remaining bridge

[Porta--Scola--Singh, Section 4.5](https://arxiv.org/html/2509.08665v1)
give determinant and covariance bounds in (4.111), renormalized
dimensions in (4.116)--(4.117), and single-mark short-memory/continuity
statements in Remarks 4.12--4.13. They retain normalization ratios in
the tree values. These estimates are stated after their running-flow
proposition, for the endpoints produced by their microscopic model.

[Mastropietro--Porta, Section 8 and Theorem 8.5](https://arxiv.org/html/2108.13708v3)
provide the reference constant-history beta decay. Equations (8.37),
(8.42), (8.44), (8.52), and (8.56) exhibit the normalized recursion,
odd covariance and quadratic reference corrections. Those primary
formulas justify the first-order separation above; they do not state
the full replacement-endpoint theorem.

To finish the ordinary GN bridge it is enough to verify, in the actual
normalized recursion, the following concrete statement: arbitrary
initial endpoints with the established hierarchy admit the same
renormalized tree and Gram bounds on a common complex history domain,
and one scale-specific history mark plus one common-history mark
preserve the path estimate and endpoint powers just described. Initial
retardation changes endpoint kernels, not the covariance, so it creates
no new covariance denominator. It still has to pass through the
localization and normalized marked-tree bookkeeping.

Once that statement and fixed-cutoff continuity are established,
Lemma 1 (or the ultraviolet variant), COUNTERTERM_INDUCTION and
FERMI_SHIFT_EXISTENCE give the physical RG trajectory. No additional
Fermi-point derivative theorem, uniqueness theorem, or microscopic
branchwise chemical-potential family is needed. Source/stress and CFT
conclusions remain separate.
