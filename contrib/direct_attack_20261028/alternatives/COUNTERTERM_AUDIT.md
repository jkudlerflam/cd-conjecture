# Independent counterterm and Fermi-shift audit

**Verdict: the revised abstract recursion theorem and conditional Fermi-shift existence theorem pass this independent audit.** Their application to arbitrary retarded initial kernels still requires the normalized marked-tree estimate identified below. This is an AI mathematical review, not formal verification or independent human refereeing.

Reviewed documents:

- [COUNTERTERM_INDUCTION.md](../vacancy/COUNTERTERM_INDUCTION.md), including the corrected coefficient-variation term in G2 and the additional operator K2.
- [FERMI_SHIFT_EXISTENCE.md](../vacancy/FERMI_SHIFT_EXISTENCE.md), including exact initial matching H3 and the strengthened continuity criterion.
- [DEGREE_REPARAMETERIZATION.md](../vacancy/DEGREE_REPARAMETERIZATION.md) and the initialization estimates it invokes.

## 1. Sequence estimates and indexing

With the convention `g_(h-1)=g_h+b_h`, the causal identity is exactly

    g_h=g_0+sum_(j=h+1)^0 b_j.

Thus its history difference contains increments indexed by `h+1,...,k`, not by `h,...,k`. Interchanging that finite sum with the short-memory kernel gives the stated K1 up to a fixed geometric constant. In the theta-weighted norm its operator norm is O(eta).

The normalized kernel of K0 is `C0*2^(theta h)*2^(theta j)` for `j>h`. An n-step path has strictly increasing indices. Replacing its intermediate `2^(2 theta j)` factors by `2^(theta j)` and summing ordered products proves

    ||K0^n||_theta <= (C0 S_theta)^n/n!,
    S_theta=1/(1-2^(-theta)).

This justifies the positive resolvent even when `||K0||` exceeds one. One must first use the finite triangular matrices and their uniform bounds, or the convergent Neumann series just proved; a small-norm assertion for K0 would be incorrect. The note uses the valid argument.

The relevant tail formula is

    nu_h=-sum_(j<=h) 2^(j-h-1) r_j.

Its norm is precisely `1/[2(1-2^(-(1+theta)))]`. Direct substitution gives `nu_(h-1)=2 nu_h+r_h`, with the signs in the note. Summing the marginal tail gives an O(epsilon*2^(theta h)) approach to its limit. These constants and off-by-one choices are consistent.

## 2. A necessary correction to the original G2

The first draft omitted the last term now present in G2. That stronger difference estimate is not a generic consequence of analytic beta functions with short memory.

An explicit scalar example is

    B_h(g)=g_h*(g_(h+1)-g_h), h<0;  B_0=0.

It satisfies the absolute history estimate on `|g_k|<eta`, after enlarging its constant by `2^(4 theta)`. Take `h=-N`, let `g_h=eta/4`, `g_(h+1)=eta/2`, and set `g_0=0`. Shift every `g_k` with `h<=k<0` by the same small number d, keeping `g_0` fixed. All recent history differences are unchanged, but

    delta B_h=d*eta/4.

The original G2 right side was at most a constant times

    |d|*[2^(-2 theta N)+eta*2^(-4 theta N)],

which is insufficient as N grows. No claim is made that this isolated monomial is a physical beta function; it is a counterexample to inferring the original G2 from analyticity and the absolute bound alone.

The corrected term is

    C |delta g_h| [Q_h(g,nu)+Q_h(tilde g,tilde nu)],
    Q_h=sum_(k=h)^0 2^(-4 theta(k-h))
                         (|g_k-g_h|+|nu_k|).

The absolute bootstrap is established before using G2 and gives `Q_h<=C epsilon*2^(theta h)` on the causal trajectories in the relevant ball. The additional difference operator is consequently

    (K2 f)_h=C epsilon*2^(theta h)*sum_(j>h) f_j,
    ||K2||_theta<=C epsilon S_theta.

After the K0 resolvent, K1 and K2 are jointly small. The corrected G2 therefore retains the O(eta) Lipschitz map and the relevant contraction. The author implemented this repair and acknowledged its provenance. The revised abstract theorem is valid conditional on its revised hypotheses.

## 3. Brouwer and removal of the infrared cutoff

The finite-scale boundary `nu_(-N-1)=0` gives exactly the truncated tail map written in FERMI_SHIFT_EXISTENCE. Its contraction factor is uniform in N and p by hypothesis. At each fixed N, continuous Picard iterates converge uniformly on the compact p-domain, so their fixed point is continuous without differentiating it.

The physical initial relation

    nu_0^phys(p)=A(p)[V p+Q(p)]+u(p)

is substantive. Keeping the other physical ultraviolet data fixed while varying only the local relevant coordinate ensures that equality of this coordinate really recovers the given action. The note does not silently introduce tunable chemical potentials.

For `|c_N|<=C_c epsilon`, the radius `R=2 C_c/v_*` and the condition `C_Q R^2 epsilon<=C_c` make `V^-1[c_N(p)-Q(p)]` a self-map of the stated cube. Brouwer requires continuity and a compact convex domain; it requires neither uniqueness nor a small derivative.

The diagonal argument passes to limits of p and of each fixed-scale history. Every beta function uses only finitely many coordinates at that scale, so H2 passes the recursion to the limit. Its backward boundary term is

    2^(-M-h-1) nu_(-M-1),

which tends to zero under the uniform weighted bound. The same bound dominates the tail series. Coordinatewise convergence is enough; convergence in the full weighted sequence norm at the moving terminal scale is unnecessary. Initial physical matching passes to the limit by H3. No continuous infinite-cutoff counterterm map or unique limiting center is used.

## 4. Fixed-cutoff continuity: a second correction

Pointwise continuity of an initial kernel and a uniform exponential L1 bound alone do not justify dominated convergence. For example, a smooth unit-mass bump of width p^2 centered at p can tend pointwise to zero as p decreases to zero while retaining its integral. This was a gap in the first proposed sufficient criterion, not in the abstract H2 assumption.

The revised note instead assumes pinned weighted total-variation/L1 continuity at each fixed degree, or common integrable domination with almost-everywhere continuous densities. This suffices, including for singular time-local kernels. A fixed coefficient measure multiplied by a continuously varying Fermi phase meets the domination condition.

At a branch-neutral beta root, pinning the root and contracting only equal branch labels makes the total endpoint branch charge zero. The phases then depend on relative positions. Fixed-cutoff propagator moment bounds control each tree; the uniform tree and degree majorants control their tails. This provides continuity without controlling a charged endpoint in a supremum over all absolute anchors. The author implemented this correction as well.

## 5. What the published normalized formulas actually supply

I opened the primary sources independently. In [Mastropietro–Porta, arXiv:2108.13708v3](https://arxiv.org/html/2108.13708v3), the discussion before (8.29) and (8.52) cancel the reference relevant term by momentum parity. Equations (8.37), (8.42), and (8.44) retain field normalization; (8.56) gives quadratic bounds for its infrared corrections and quartic beta function. Theorem 8.5 bounds the beta function on constant histories. Section 10 invokes tree comparison for the paper's lattice-generated kernels. It does not state a replacement-data theorem for every action considered here.

In [Antinucci–Mastropietro–Porta, arXiv:1708.08517v4](https://arxiv.org/html/1708.08517v4), (9.55) supplies field-degree endpoint powers, Remarks 9.1–9.2 give short memory and tree continuity, and (9.72) separates trees reaching the initial scale. These ingredients explain the intended extension, but do not by themselves state the reduced history derivative estimate below with arbitrary replacement endpoints and anomalous normalization retained.

A local quartic vertex cannot produce an external-momentum derivative through its one-vertex tadpole, and the symmetric linear-reference tadpole itself vanishes by oddness. The exact linear relevant propagation is already excluded from the beta function. Therefore first-order terms that survive in the replacement action must be assigned an initial-endpoint or curvature label and retain its scale gain. This assignment needs to remain explicit after normalization; setting the wave-function ratios to one is not allowed.

There is also an exact one-endpoint tree observation. A Gaussian expectation of a single local quartic yields only a local quartic, a constant quadratic tadpole, and a scalar. At an intermediate vertex the remainder operation removes the first two, while the scalar belongs to the vacuum sector. Consequently no such one-endpoint tree contributes to a marginal/relevant beta function through a genuine intermediate scale interval. At the adjacent root, its uncontracted quartic is canonical transport and its tadpole has no external-momentum derivative. A local quadratic endpoint similarly has only canonical transport and a scalar. These statements do not require oddness of the full lattice propagator; oddness is used only to cancel the adjacent reference tadpole. A surviving lattice-curvature tadpole is a decaying direct contribution.

The normalized lattice covariance still requires scrutiny: its curvature coefficient can contain accumulated Z ratios, so it cannot simply be declared a function only of the current velocity. Such history dependence is interaction-generated. A curvature gain `2^((1-O(eta))h)` can pay `2^(-4 theta(k-h))` for any `h<=k<=0` when `4 theta<1-O(eta)`. This gives room for the marked normalization estimate, but does not substitute for deriving that estimate from the solved ratio recursion.

## 6. A precise sufficient remaining lemma

Here is a concrete target, rather than an assertion that generic smallness closes the RG. Fix the initial data and let `u=g_h`, `d_k=g_k-u`. For either reduced marginal or relevant beta function, establish a uniformly convergent normalized GN representation of the form

    B_h = F_h(u) + Psi_h(u,{d_k},{nu_k}),
    Psi_h(u,0,0)=0,
    w_(h,k)=2^(-4 theta(k-h)),

with the following estimates on a fixed complex neighborhood of the provisional histories:

    |F_h(u)| <= C 2^(2 theta h) [epsilon+|u-g_0|],
    ||D_u F_h|| <= C 2^(2 theta h),

    ||D_(d_k) Psi_h|| + ||D_(nu_k) Psi_h||
                      <= C eta w_(h,k),
    ||D_u D_(d_k) Psi_h|| + ||D_u D_(nu_k) Psi_h||
                      <= C w_(h,k).

The formal history coordinates may be extended independently while deriving these estimates, then restricted to `d_0=g_0-u`. The initial kernels are held fixed in this differentiation. Integrating these derivatives gives G1 and the **corrected** G2: first vary deviations at fixed u, then vary u at fixed deviations. The mixed derivative and `Psi(u,0,0)=0` yield the Q_h term exactly. Thus this lemma would suffice for both audited existence theorems; a Fermi-point derivative estimate is not needed.

The estimate must include three concrete operations.

1. Mark one running localized endpoint, propagator, or solved normalization ratio at scale k. After removing the exact linear recursion and the vanishing single-local-vertex graphs, a remaining interaction factor must supply eta. A second mark in the common history value may remove this factor but must retain w_(h,k). Finite differences can be used instead of derivatives if the same marked bounds hold.
2. For every tree with a marked initial arbitrary-degree endpoint, or a curvature correction, retain the initial-scale gain after localization and either history mark. Endpoint degree reserves absorb field-label and mark multiplicities once; a fresh fixed field-weight loss at every RG scale would not suffice.
3. Solve and differentiate the single-step normalization identities inside the same majorant. The relevant quantities are ratios of Z's, whose logarithms have small increments; the accumulated Z itself need not converge. The denominators `1+z` stay away from zero, but that fact alone does not prove the history-weighted marked estimate.

The reference constant-history Ward estimate addresses F_h's reference part. Its analytic dependence provides the first derivative after choosing a smaller fixed neighborhood. The arbitrary initial endpoint norm and short-memory bounds address the other part only after the marked normalized tree sum is proved. The extra Q_h term shows why the unmarked bounds cannot simply be relabeled as the old G2.

No new microscopic chiral Ward identity has been identified as necessary. The unresolved issue is this specific extension of the normalized marked-tree majorant, including its allocation of the extra small factor, to the replacement endpoint class. The sequence fixed point, physical Fermi-center existence, and their continuity needs are now separately closed under explicit hypotheses. Current-source comparison, stress/pressure control, and a conformal scaling limit remain further tasks.

Reviewed 2026-10-28 by the alternatives worker. The author and root coordinator received both corrections while the proofs were being revised. No commits or shared-status edits were made by this worker.
