# C113. Smooth modulation, the anomaly, and the missing lattice hypotheses

2026-10-29

## 1. Status and definitions

[Proved] The smooth-modulation lemma in the supplied specification is false, even if one first chooses the representative attaining the gauge-invariant scale M. Section 7 gives an exact gapped counterexample for unrestricted representatives. Section 9 gives an exact variational counterexample at the ferromagnetic z = 2 point. Section 9.1 embeds these mechanisms into d = 4 models with a unique minimum-spread representative and M = 8. All use positive smooth periodic weights bounded away from zero, with no physical boundaries or cutoff singularities.

[Proved, conditional on stated matching] The continuum stress-tensor coefficient is

    kappa = 1/(24 pi).

If a chosen lattice representative satisfies pure stress-tensor matching for its weighted ground energy, the asymptotic lemma is equivalent to A <= 24 pi C M f(d). Without this additional matching, a weak-modulation lattice response bound is still sufficient under spectral stress-tensor matching, but a bound on A does not imply the lattice lemma. The distinction is substantive.

[Assumed, conventions] Lattice spacing and hbar are one. H_L = sum_i h_i is periodic, with h_i acting on sites i,i+1. Write e_L = E0(H_L)/L and e0 = lim e_L. A scalar shift of h is immaterial. M is the infimum of the half-spread after adding alpha I + o tensor I - I tensor o. All velocities below are positive speed magnitudes. A bond representative is fixed before its weights are applied.

[Proved] The finite-ring nonnegative deficit is

    D_L(w) = e_L sum_i w_i - E0(sum_i w_i h_i) >= 0.             (1.1)

Use a translation-averaged ground-state density matrix if necessary. Its expectation of every bond is e_L, so the variational principle proves (1.1). Replacing e_L by e0 on a finite ring adds a Casimir term even when w is constant. It is not an innocuous replacement at order 1/L.

[Definition] Fix a positive smooth function u of period one, with inf u > 0, and set w_i = u((i+1/2)/R), repeated on a ring of L = nR sites. Define the deficit per modulation period by

    D_R[u] = lim_(n -> infinity) D_(nR)(w)/n,
    G_R[u] = sum_(i=0)^(R-1) (w_i-w_(i-1))^2/w_i,
    I[u] = integral_0^1 u'(s)^2/u(s) ds.                       (1.2)

[Assumed] The thermodynamic limit in (1.2) is taken for each fixed periodic interaction. This limit exists for finite-range finite-dimensional chains. If a selected vacuum sector is imposed instead of the full ground energy, it must be specified separately.

[Proved] For these profiles,

    G_R[u] = I[u]/R + o(1/R).                                  (1.3)

Taylor expansion gives w_i-w_(i-1) = u'(i/R)/R + O(R^-3); the Riemann sum supplies (1.3). Shifting which neighboring weight appears in the denominator changes only the remainder. For nonconstant u, I[u] > 0. Uniformly multiplying all weights by b > 0 multiplies D_R and G_R by b.

[Definition] The pointwise asymptotic lemma for a specified class of representatives is

    limsup_(R -> infinity) D_R[u]/G_R[u] <= C M f(d)            (1.4)

for every fixed bond in that class and every fixed nonconstant positive smooth u. It does not assert a uniform wavelength threshold over bonds. A remainder merely o(1) is insufficient, because G_R itself is O(1/R). The remainder must be o(G_R).

## 2. Exact continuum calculation for all CFT sectors

[Assumed, continuum setting] Index independent chiral stress tensors by alpha. Their central charges are c_alpha >= 0, their speeds are v_alpha > 0, and their commuting positive-energy Virasoro representations contain the common vacuum. Use physical chiral energy densities t_alpha normalized by

    U(f) t_alpha(x) U(f)^-1
       = f'(x)^2 t_alpha(f(x)) - c_alpha/(24 pi) {f,x},
    {f,x} = f'''/f' - (3/2)(f''/f')^2.                        (2.1)

The physical Hamiltonian is sum_alpha v_alpha integral t_alpha dx. Set

    A = (1/2) sum_alpha c_alpha v_alpha.                       (2.2)

For a nonchiral sector with c_R = c_L = c, (2.2) gives c(v_R+v_L)/2. This is the velocity convention of the background. No charge current or global symmetry enters (2.1).

[Proved within the continuum assumptions] On a circle of circumference L with positive smooth w, define

    ell_w = integral_0^L dx/w(x),
    vbar_w = L/ell_w,
    f'(x) = vbar_w/w(x),    f(x+L) = f(x)+L.                  (2.3)

The symbol vbar_w is a dimensionless weight, not a sector velocity. Equation (2.1), followed by y = f(x), gives for each chirality

    U(f) [v_alpha integral w t_alpha dx] U(f)^-1
      = v_alpha vbar_w integral t_alpha(y) dy
        - c_alpha v_alpha/(24 pi) integral w {f,x} dx.        (2.4)

Indeed w f'^2 dx = vbar_w dy. Direct differentiation gives

    f''/f' = -w'/w,
    f'''/f' = -w''/w + 2w'^2/w^2,
    {f,x} = -w''/w + w'^2/(2w^2),
    integral_0^L w {f,x} dx = (1/2) integral_0^L w'^2/w dx.   (2.5)

The integral of w'' vanishes by periodicity. The uniform vacuum energy per chirality is -pi c_alpha v_alpha/(12L). Consequently

    E_CFT(w) = -pi A/(6 ell_w)
               - A/(24 pi) integral_0^L w'^2/w dx.           (2.6)

This is an exact unitary conjugacy calculation in the continuum theory, not a derivative expansion. Different velocities multiply the individual conjugacies and add in (2.2). The same diffeomorphism acts independently in the commuting chiral factors. Nontrivial modular invariants do not require a tensor factorization of the whole Hilbert space, provided the common vacuum and the commuting stress representations are present.

[Proved within the same assumptions] Relative to the uniform circle vacuum, rather than the infinite-line energy density, the deficit is

    D_CFT,L(w) = A/(24 pi) integral w'^2/w
                + pi A/6 [1/ell_w - W/L^2],
    W = integral_0^L w dx.                                   (2.7)

The second term is generally negative. For constant w, both terms sum to zero. For w = 1 + epsilon cos(2pi x/L), the two terms cancel exactly. The perturbation is a global conformal generator and leaves the circle vacuum a ground state. A formula omitting this term fails even at leading gradient order on a single circle.

[Proved] For n repetitions of a profile of wavelength R, the first term of (2.7), divided by n, is A I[u]/(24 pi R). The second term divided by n is O(1/(n^2 R)). Taking n -> infinity first yields

    D_CFT,R[u] = A I[u]/(24 pi R).                            (2.8)

This thermodynamic prescription removes the global finite-size ambiguity without any boundary condition at zeros of w.

[Proved in the cited literature] The sharp chiral quantum energy inequality is

    inf_psi integral G(x) <t_alpha(x)>_psi dx
      = -c_alpha/(12 pi) integral ((sqrt G)')^2 dx
      = -c_alpha/(48 pi) integral G'^2/G dx                  (2.9)

for the nonnegative test functions and positive-energy diffeomorphism-covariant CFTs covered by Fewster and Hollands, Theorem 4.1. Their theorem also establishes simultaneous sharpness for the two commuting chiral components. It independently fixes the line coefficient in (2.8). See https://arxiv.org/html/math-ph/0412028 .

[Proved, interpretation] Charged, neutral and unequal-speed sectors all contribute positively to (2.2). A difference of left and right central charges affects the gravitational anomaly and chiral transport, but does not replace the sum c_alpha v_alpha in this static energy bound. A stand-alone finite-dimensional one-dimensional lattice regulator imposes its own anomaly-consistency restrictions on which total chiral spectra are realizable. The calculation does not assume that every formal anomalous spectrum has such a regulator.

## 3. Curved background, local energy, and soft edges

[Proved within the continuum setting] For a nonchiral sector, the local velocity v w corresponds to the static metric

    ds^2 = -v^2 w(x)^2 dt^2 + dx^2
         = w(x)^2[-v^2 dt^2 + dy^2],    dy = dx/w.            (3.1)

The anomalous transformation of the stress tensor under this coordinate and Weyl change is the Schwarzian in (2.5). Thus the curved-background conformal anomaly and the conjugacy calculation are two descriptions of the same contribution. This is a Hamiltonian definition of the source coupling, which fixes the stress convention. With unequal speeds one uses the chiral conjugacies directly; a common Lorentz metric for all sectors is unnecessary. Moosavi's Lemma 5.3 gives the corresponding inhomogeneous-CFT conjugacy, https://arxiv.org/html/1912.04821 .

[Proved within the continuum setting] In the thermodynamic weighted vacuum, either (2.1) or functional differentiation of (2.8) gives

    <h(x)>_w - e0 = A/(24 pi) [2w''/w - w'^2/w^2],
    e0 - <h(x)>_w = A/(24 pi) [w'^2/w^2 - 2w''/w].            (3.2)

For the differentiation check,

    delta/ delta w integral w'^2/w
        = w'^2/w^2 - 2w''/w.

Multiplying (3.2) by w and integrating recovers (2.8), because integral w'' = 0. A local deficit need not be positive at every point even though the integrated deficit is nonnegative.

[Proved within a regulated continuum calculation] For w = x on a half-line, (3.2) gives A/(24 pi x^2). For w = x(ell-x)/ell on an interval,

    w'^2 - 2w w'' = 1,
    e0 - <h(x)>_w = A ell^2/[24 pi x^2(ell-x)^2],
    integral_epsilon^(ell-epsilon) w'^2/w dx
        = 2 log(ell/epsilon) - 4 + o(1).                     (3.3)

Thus the two soft edges give A log ell/(12 pi). These checks reproduce the stress contribution recorded by the modular work.

[Assumed, necessary limitation] In a lattice chain the cutoff neighborhoods of soft edges are not slowly varying relative to the weight. Boundary energy, oscillating bond observables, and gradient-induced fields require separate matching. C105 already shows why a raw local bond profile cannot be identified with (3.2) without smoothing. Neither (3.3) nor that smoothing removes a smooth nonuniversal gradient term. The bounded positive profiles in (1.2) avoid these endpoint issues altogether.

## 4. What infrared irrelevance does and does not prove

[Assumed, pure matching hypothesis P] For a fixed representative h, the renormalized coupling of w to the lattice energy has, to the needed accuracy, only the continuum stress source. It has no gradient source for an operator that contributes at order R^-1 or larger, and no extra local two-gradient identity term. In this convention the weighted ground-energy remainder obeys

    D_R[u] - D_CFT,R[u] = o(1/R).                            (4.1)

Uniform unweighted finite-size CFT scaling alone does not imply P. It is an additional lattice-to-continuum statement about a spatial source.

[Proved by infrared scaling, under regular conformal perturbation theory] A renormalized irrelevant operator O_j of scaling dimension Delta_j > 2, with coupling g_j, contributes to a one-period energy at first order on scale R as

    g_j R integral_0^1 ds R^(-Delta_j) F_j[u](s)
       = O(g_j R^(1-Delta_j)).                               (4.2)

The ratio to the anomaly is O(g_j R^(2-Delta_j)), which vanishes for each fixed bond. At higher order the dimensionless running coupling is g_j R^(2-Delta_j), so each further infrared insertion supplies another suppressed factor, assuming finite renormalized integrals and no dangerously irrelevant coupling. Intersector irrelevant couplings obey the same reasoning after keeping fixed nonzero velocities. Marginally irrelevant terms can instead give R^-1 times a vanishing inverse logarithm. They are subleading, but not a higher integer derivative order. Exactly marginal fixed-point parameters must be included in the sector data.

[Proved, qualification to (4.2)] This argument controls the long-distance part of the renormalized calculation. It does not control short-distance contact terms in products of O_j and energy-density insertions. Integrating out microscopic modes can leave an identity contribution to the source-dependent energy,

    E_extra[w] = -B_h integral w'^2/w dx.                    (4.3)

This term has the same gradient order as the anomaly. It vanishes for uniform w and is invisible to the uniform Casimir coefficient A. A lattice source has a preferred time slicing and is not automatically a generally covariant coupling to a relativistic metric. In particular, one cannot discard (4.3) by invoking the topological nature of the two-dimensional Einstein-Hilbert action.

[Proved, conditional on a local derivative expansion] The form in (4.3) is also forced for a local two-gradient term with one scalar source by E[bw] = bE[w]. After integrating a term proportional to w'' by parts, write the term as integral F(w)w'^2. Homogeneity gives F(bw) = F(w)/b, hence F(w) = const/w. Its coefficient is microscopic information, not determined by c_alpha and v_alpha. Section 7 realizes such a term without any continuum approximation.

[Proved] Total derivatives cannot be declared irrelevant without examining the differentiated operator. A lattice divergence produces a source proportional to w' for its underlying local operator. If O has CFT dimension Delta_O, the perturbation is -integral w' O. For Delta_O = 1 it can produce an additional two-gradient susceptibility term. If Delta_O < 1 its long-distance response can be more singular than the stress response.

[Heuristic, controlled only with the stated scaling assumptions] For a primary with Delta_O < 1 and a nonzero relevant source response, a gradient source of magnitude R^-1 has correlation length R^(1/(2-Delta_O)), smaller than R. Local hyperscaling predicts an energy gain per period of order

    R * (R^-1)^(2/(2-Delta_O))
      = R^(-Delta_O/(2-Delta_O)).                            (4.4)

This decays more slowly than R^-1. The formula requires a nonzero equation-of-state amplitude and control near source zeros. It is an explanatory scaling mechanism, not a general lattice theorem. The exact ferromagnetic counterexample below does not rely on (4.4).

[Proved, conclusion] The assertion requested in the task, that every correction from microscopic or non-CFT degrees of freedom is higher order in gradients, is false without additional source-matching assumptions. Infrared irrelevance of renormalized operators and absence of an ultraviolet analytic response are different assertions.

## 5. The exact implication and its limited converse

[Proved, conditional theorem] Suppose P holds for a specified representative, and suppose the asymptotic lattice lemma (1.4) holds. Equations (1.3), (2.8) and (4.1) give

    lim_(R -> infinity) D_R[u]/G_R[u] = A/(24 pi),
    A <= 24 pi C M f(d) = (C/kappa) M f(d).                   (5.1)

One nonconstant profile suffices for this implication. There is no symmetry assumption in the continuum calculation.

[Proved, conditional converse] Under the same P, if A <= C' M f(d), then for every fixed positive smooth nonconstant profile,

    D_R[u] <= [C'/(24 pi)] M f(d) G_R[u] + o(G_R[u]).         (5.2)

This is the precise sense in which the two statements are equivalent. A uniform remainder over bonds, over profiles whose minimum tends to zero, or over velocities tending to zero does not follow. Even when the coefficient is uniformly bounded, the wavelength needed for matching can depend on the bond.

[Proved] If the correct expansion is instead

    D_R[u] = [A/(24 pi) + B_h] G_R[u] + o(G_R[u]),             (5.3)

then a bound on A alone says nothing about B_h. For B_h >= 0 the lattice lemma is still sufficient for (5.1), but the converse fails. If gradient-induced relevant fields spoil even (5.3), the stronger lemma can fail for a single fixed bond. The examples in Sections 7 and 9 establish both types of failure.

[Proved, logical scope] A statement about A for conformal infrared limits cannot by itself imply a lemma for every nonconformal or gapped bond. Those bonds have additional response data which are absent from A. Declaring A = 0 in a gapped system does not make its weighted deficit zero.

## 6. A weaker sufficient reduction using positive spectral response

[Proved, finite ring with a simple ground eigenvalue] Let

    h_q = L^(-1/2) sum_j exp(-iqj)(h_j-e_L),
    S_L(q) = sum_(m>0) |<m|h_q|0>|^2/(E_m-E_0),
    Sbar_L(q) = [S_L(q)+S_L(-q)]/2.                          (6.1)

For w_j = 1 + epsilon cos(qj), with q neither 0 nor pi, momentum orthogonality and second-order perturbation theory yield

    D_L(w) = epsilon^2 L Sbar_L(q)/2 + O(epsilon^3),
    G_L(w) = epsilon^2 L [4 sin^2(q/2)]/2 + O(epsilon^3).
                                                                    (6.2)

Take epsilon -> 0 first at fixed finite L, then the thermodynamic limit at fixed q, and only then q -> 0. Thermodynamic differentiability and spectral convergence are extra assumptions when taking these limits in a critical system.

[Assumed, spectral stress matching S] The low-energy spectral measure of h_q contains the vacuum stress-tensor contribution with the normalization in (2.1), with errors negligible for its q^2 inverse-energy moment. Other components are orthogonal to this contribution in the scaling limit. This is weaker than matching the complete weighted ground energy. In particular it allows extra positive high-energy response.

[Proved within the CFT] The normalization in the next step can be checked without a Fourier-transform convention ambiguity. On a circle, write a chiral stress density as

    t_alpha(x) = (2pi/L^2)
       [sum_n L_n exp(2pi i n x/L) - c_alpha/24].

For the Fourier sign that creates L_(-n)|0>, n >= 2, the normalized physical density mode is (2pi v_alpha/L^(3/2)) L_(-n). The Virasoro commutator and L_n|0> = L_0|0> = 0 give

    ||L_(-n)|0>||^2 = <0|[L_n,L_(-n)]|0>
                    = (c_alpha/12)n(n^2-1),
    ||h_q|0>||^2 = (4pi^2 v_alpha^2/L^3)
                    (c_alpha/12)n(n^2-1).

At fixed |q| = 2pi n/L as L -> infinity this is c_alpha v_alpha^2 |q|^3/(24pi), and the excitation energy is 2pi v_alpha n/L = v_alpha |q|. The opposite chirality creates states at the opposite sign of momentum.

[Proved under S] Thus the chiral stress normalization gives a vacuum matrix-element weight for the L^(-1/2)-normalized density mode

    c_alpha v_alpha^2 |q|^3/(24 pi)

in the momentum direction created by that chirality, at energy v_alpha |q|. Dividing by its excitation energy and averaging q and -q gives

    Sbar_IR(q) = q^2 sum_alpha c_alpha v_alpha/(48 pi)
                = A q^2/(24 pi).                           (6.3)

Equivalently, expand the exact energy formula (2.8) to second order in epsilon. The spectral measure in (6.1) is nonnegative, so restricting it to these low-energy stress states proves

    liminf_(q -> 0) Sbar(q)/[4 sin^2(q/2)] >= A/(24 pi).       (6.4)

The restriction must be justified by S. A signed subtraction of overlapping operator amplitudes would not suffice. Currents that are orthogonal CFT primaries add response rather than cancel the stress contribution.

[Proved, conditional sufficient lattice inequality] An upper bound

    limsup_(q -> 0) Sbar(q)/[4 sin^2(q/2)] <= C M f(d)         (6.5)

therefore implies A <= 24 pi C M f(d) under S, even without P. It is enough to prove a lattice inequality for infinitesimal modulations, with remainders controlled in the order of limits just stated. A finite-amplitude lemma with only an unspecified pointwise o(G) remainder cannot silently be differentiated to obtain (6.5).

[Proved, exact continuity identity] Define j_j = i[h_(j-1),h_j]. Nearest-neighbor locality gives

    i[H,h_j] = j_j-j_(j+1),
    i[H,h_q] = (1-exp(iq)) j_q.                              (6.6)

For every excited eigenstate, (6.6) implies

    S_L(q) = 4 sin^2(q/2)
             sum_(m>0) |<m|j_q|0>|^2/(E_m-E_0)^3.           (6.7)

Thus (6.5) is a bound on a positive inverse-cube energy-current moment. It is not a norm bound on the bond gradient of A. Equation (6.7) also identifies the gapped contribution that the Casimir coefficient cannot see.

## 7. Exact gapped example and a conformal stack

[Proved] Let X,Z be Pauli matrices, g > 0, and lambda real. Take d = 2 and

    h_(g,lambda) = -(g/2)(Z tensor I + I tensor Z)
                  + lambda(X tensor I - I tensor X).        (7.1)

On every periodic chain,

    H = -g sum_i Z_i,    e0 = -g,    gap = 2g.                (7.2)

The unique ground state is the all-up product state. The uniform Hamiltonian and its gap do not depend on lambda.

[Proved] The gauge-optimized scale is M = g. Removing the divergence gives a bond with spectrum -g,0,0,g, hence M <= g. In any gauge and after any scalar shift, the expectations in |up up> and |down down> differ by 2g. The divergence has zero expectation in each of these states because its two single-site expectations agree. The spectral spread is at least their difference, so M >= g.

[Proved] For arbitrary positive periodic weights, put

    a_i = (w_i+w_(i-1))/2,    b_i = w_i-w_(i-1).

Then H_w is still a sum of independent single-spin Hamiltonians,

    H_w = sum_i [-g a_i Z_i + lambda b_i X_i],
    D_L(w) = sum_i [sqrt(g^2 a_i^2 + lambda^2 b_i^2)-g a_i]
           = sum_i lambda^2 b_i^2/
               [sqrt(g^2 a_i^2 + lambda^2 b_i^2)+g a_i].     (7.3)

No low-energy approximation enters this identity. Since a_i > 0,

    0 <= D_L(w) <= lambda^2/(2g) sum_i b_i^2/a_i.             (7.4)

For fixed g,lambda and the profiles of (1.2), Taylor expansion of the square root, with b_i/a_i = O(1/R), gives per period

    D_R[u] = lambda^2/(2g) G_R[u] + o(1/R).                   (7.5)

The quadratic coefficient is lambda^2/gap. It is produced entirely by gapped degrees of freedom and is not higher order in gradients.

[Proved, counterexample] M = g and d = 2 are fixed while lambda is arbitrary. Therefore no finite universal C f(2) can bound (7.5) in terms of the gauge-invariant M for all representatives. For any proposed constant, choose lambda first and then R sufficiently large. A remainder negligible as R -> infinity cannot repair the failure. This is a gapped failure with a unique ground state and a fixed positive gap.

[Proved, conformal extension] Tensor (7.1) with the usual XX bond on a separate spin at each site. The site dimension is d = 4. The XX conventions in the background give A = 1 and M_XX = 1/2; the gapped factor adds no CFT sector. The deficits add exactly because the two weighted Hamiltonians commute on separate tensor factors. Consequently the gapped factor adds lambda^2 G_R/(2g) while leaving A unchanged.

[Proved] In this stack M = 1/2+g. The upper bound follows by separately removing the divergence and adding the centered representatives. For the lower bound use the two extremal XX Bell states tensored with all-up or all-down states of the gapped bond, choosing their signs consistently. Their energies differ by 2(1/2+g), and each two-site state has equal left and right marginals, so every allowed divergence has zero expectation in both. This proves the lower bound. Thus even a chain with a nontrivial unitary CFT infrared does not have its complete leading deficit fixed by A under the stated assumptions.

[Proved, scope] These examples do not refute a lemma that first selects the minimum-spread representative. In (7.1), the representative lambda = 0 has zero deficit. They do refute applying such a lemma to an arbitrary representative with its separately optimized M, and they refute the requested converse from A alone as written. Section 9.1 supplies a separate example that also defeats the minimum-spread prescription.

## 8. A gapped bound with explicit gap and correlation dependence

[Proved, finite-dimensional spectral identity] Let V = sum_i (w_i-1)h_i and H_s = H+sV, for 0 <= s <= 1. Suppose the ground eigenvalue is simple throughout this path. Twice differentiating it gives

    E''(s) = -2 sum_(m>0) |<m_s|V|0_s>|^2/omega_m(s),
    D_L(w) = 2 integral_0^1 (1-s)
               sum_(m>0) |<m_s|V|0_s>|^2/omega_m(s) ds.      (8.1)

The first derivative at zero is e_L sum_i(w_i-1), which proves the second identity with the baseline in (1.1).

[Proved] Write u_i = w_i-1 and k_i = i[h_i,h_(i+1)]. Direct expansion of the commutator gives

    i[H_s,V] = sum_i (u_(i+1)-u_i) k_i =: K_u.               (8.2)

The s-dependent pieces cancel, since (1+s u_i)u_(i+1) - (1+s u_(i+1))u_i = u_(i+1)-u_i. Matrix elements of (8.2) replace each factor of V in (8.1) by K_u/omega. If the path gap is at least Delta_* > 0, then

    D_L(w) <= Delta_*^-3 sup_s Var_(0_s)(K_u).                (8.3)

The factor 2 integral_0^1(1-s) ds is one.

[Assumed, quantitative clustering data] Suppose the connected current covariances along the path obey the Schur bound

    sup_(s,i) sum_j |<k_i k_j>_(0_s,c)| <= K_*.               (8.4)

For example, a uniform bound C_j exp(-|i-j|/xi_*) gives

    K_* <= C_j [1+2/(exp(1/xi_*)-1)].                         (8.5)

A uniformly gapped local path with bounded interactions has exponential clustering under the usual unique-ground-state assumptions. Here (8.4) is the precise quantitative input, avoiding an unproved numerical value for its prefactor.

[Proved under these gap and clustering inputs] Apply the Schur bound to the quadratic form in (8.3):

    Var(K_u) <= K_* sum_i (w_(i+1)-w_i)^2,
    D_L(w) <= K_*/Delta_*^3 sum_i (w_(i+1)-w_i)^2
            <= w_max K_*/Delta_*^3 G_L(w).                   (8.6)

For a scalar-centered representative with ||h|| = m, ||k_i|| <= 2m^2. This provides the local size entering a clustering estimate, but m must not be replaced by M unless that representative attains M.

[Proved, limitation] The bound (8.6) controls the deficit by a gap and current correlations. It is not the universal desired bound. Neither the uniform gap of H alone nor the datum A = 0 supplies Delta_* and K_* along an arbitrary finite-amplitude path. Degenerate ground spaces require control of the effective perturbation inside the ground space before using the inverse-energy formula. Section 9 illustrates why omitting that condition is dangerous.

[Proved] If h_i-e0 I >= 0 for every bond and a common state annihilates all these positive operators, then every positive weighted sum has the same property and D_L(w) = 0. This frustration-free observation covers many gapped models and some gapless models exactly, without gap estimates.

## 9. XXZ and a rigorous z = 2 violation for an allowed representative

[Proved] Consider the spin-1/2 ferromagnetic bond

    h_i^F = -J S_i dot S_(i+1),    J > 0.                    (9.1)

Its triplet energy is -J/4 and its singlet energy is 3J/4. The fully polarized states minimize every bond. Therefore e0 = -J/4 and every positive weighted chain has D(w) = 0 in this representative. A one-magnon state has energy J(1-cos k) above the polarized ground state, so its dispersion is J k^2/2 + O(k^4). This is a nonconformal z = 2 point, not a CFT with a positive velocity.

[Proved] Now use the equally allowed nearest-neighbor bond

    h_i^(F,lambda) = -J S_i dot S_(i+1)
                    + lambda(S_i^z-S_(i+1)^z).               (9.2)

The uniform periodic Hamiltonian is still (9.1). M = J/2 in every gauge class. To prove the lower bound, compare a singlet and any triplet state with equal single-site marginals, for example the zero-magnetization triplet. Every divergence has zero expectation in each, while their energies differ by J. The original bond attains this half-spread.

[Proved, product-state witness] Take R > 2, 0 < epsilon < 1, lambda > 0 and

    w_i = 1+epsilon cos[k(i+1/2)],    k = 2pi/R.

Use the product of spin coherent states with unit vectors

    n_i = (cos(ki), 0, sin(ki)),    <S_i> = n_i/2.            (9.3)

Then n_i dot n_(i+1) = cos k, sum_i w_i = R, and

    w_i-w_(i-1) = -2epsilon sin(ki) sin(k/2),
    sum_i sin^2(ki) = R/2.

The exchange cost relative to e0 W is J R(1-cos k)/4. The telescoping term has weighted expectation

    (lambda/2) sum_i (w_i-w_(i-1)) sin(ki)
       = -lambda epsilon R sin(k/2)/2.

The variational principle therefore gives the rigorous lower bound

    D_R[u] >= lambda epsilon R sin(pi/R)/2
               - J R[1-cos(2pi/R)]/4.                      (9.4)

The same product witness can be repeated on any ring of nR sites, so (9.4) holds for the thermodynamic deficit per period. Taking R -> infinity yields

    liminf D_R[u] >= lambda epsilon pi/2 > 0,
    G_R[u] = I[u]/R + o(1/R) -> 0.                          (9.5)

Thus D_R/G_R diverges for one fixed d = 2 bond with fixed finite couplings. No finite C f(2), and no remainder o(G_R), can make the proposed lemma hold for (9.2). Its raw half-spread is also fixed, so replacing M by the raw half-spread does not fix this example.

[Proved, relation to the requested XXZ edge] For h = S^x S^x + S^y S^y + Delta S^z S^z, the point Delta = -1 is mapped to (9.1), with J = 1, by a pi rotation about z on alternating sites. The longitudinal divergence in (9.2) is unchanged. Even-length rings and positive weights are preserved. Thus (9.4) is an edge case at the ferromagnetic XXZ point itself.

[Assumed, standard XXZ infrared dictionary] For -1 < Delta < 1 write Delta = cos gamma, 0 < gamma < pi. The conventional zero-field Luttinger parameters are

    v = pi sin gamma/(2 gamma),
    K = pi/[2(pi-gamma)],
    chi = K/(pi v).                                         (9.6)

Sirker and Bortz give the velocity and susceptibility dictionary; their boson-radius convention uses a Luttinger parameter twice the K in (9.6). See https://arxiv.org/html/cond-mat/0511272 .

[Proved within that Gaussian infrared description] Add lambda(S_i^z-S_(i+1)^z) to the critical XXZ bond. The smooth density rho has static energy rho^2/(2chi). Its source from the weighted divergence is lambda w' rho. Completing the square gives

    integral [w rho^2/(2chi) + lambda w' rho] dx
      = integral w/(2chi) [rho+chi lambda w'/w]^2 dx
          - chi lambda^2/2 integral w'^2/w dx.               (9.7)

The minimizing density has zero total integral on a period because integral w'/w = 0. Hence the construction is also compatible with zero total magnetization. It supplies an additional coefficient chi lambda^2/2, besides v/(24 pi).

[Proved within (9.6)] Put delta = pi-gamma. As Delta -> -1 from above,

    v = delta/2 + O(delta^2),
    K = pi/(2delta),
    chi = delta^-2 [1+O(delta)],
    1+Delta = delta^2/2 + O(delta^4).                         (9.8)

The stress contribution tends to zero, while the improvement coefficient behaves as lambda^2/[4(1+Delta)]. Take R -> infinity at fixed Delta before approaching -1, so that the induced magnetization remains in the Luttinger regime. This is a nonuniform crossover, not a divergence of A. At Delta = -1 the rigorous coherent-state witness replaces the Gaussian response calculation.

[Proved, scope] The ordinary XXZ representative at Delta = -1 has zero deficit, and its M-minimizing representative need not retain the divergence. These two-state examples refute the literal all-representative lemma. The four-state construction next removes the possible repair by choosing a minimum-spread representative.

## 9.1 A four-state counterexample in the unique minimizing gauge

[Proved, stronger obstruction] The previous examples might suggest that minimizing the bond spread removes every dangerous gradient source. It does not. The following construction has a unique minimizing representative, up to scalar shifts, and retains the failure (9.5).

[Definition] A site has Hilbert space C^2_A direct-sum C^2_B, so d = 4. The bond preserves the A or B label of each site. In the four label blocks define

    h_AA = -S_1 dot S_2 + (S_1^z-S_2^z),
    h_BB = 16 S_1 dot S_2 + 8 I,
    h_AB = h_BA = 12 I.                                    (9.9)

There are no matrix elements between label blocks. The spin operators within each block are spin-1/2 operators. Thus the active AA block is (9.2) with J = lambda = 1.

[Proved, exact bond scale] The BB singlet has energy -4 and its triplets have energy 12. The AA eigenvalues are -1/4 twice and 1/4 plus or minus sqrt(5)/2, all strictly between -4 and 12. The mixed blocks have energy 12. Therefore h has spectrum contained in [-4,12] and raw half-spread 8. The BB singlet and zero-magnetization BB triplet have equal left and right marginals, so their expectations of every gauge divergence vanish. Their energy difference 16 proves M >= 8. Consequently

    M = 8.                                                 (9.10)

[Proved, uniqueness of the minimizing gauge] Remove an irrelevant scalar shift and write a competing representative as h+D_o, where D_o = o tensor I - I tensor o. Its expectations in the BB singlet and triplet are still -4 and 12. If it attains half-spread 8, its entire spectrum must lie in [-4,12]. A vector whose expectation equals a spectral endpoint is an eigenvector at that endpoint. In particular

    D_o |s_BB> = 0.                                         (9.11)

Decompose o into AA, AB, BA and BB blocks. The AB and BA output components of (9.11) occupy different label sectors, so each is zero separately. Full Schmidt rank of the BB singlet then implies o_AB = o_BA = 0. In the BB component, (o_B tensor I-I tensor o_B)|s_BB> = 0 implies o_B = c I: each traceless Pauli component maps the singlet to a triplet, and these three images are linearly independent. Thus o is block diagonal with scalar BB block.

[Proved, continued] On AB and BA the modified operators are respectively

    12 I + (o_A-cI) tensor I,
    12 I - I tensor (o_A-cI).

Both must have largest eigenvalue at most 12. Hence o_A-cI is both nonpositive and nonnegative, so o_A = cI. Thus D_o = 0. No nontrivial gauge divergence preserves the optimum M. This argument does not assume that the allowed gauge operator preserves the labels; its off-diagonal blocks were excluded by (9.11).

[Proved, uniform ground state] Every label configuration is an invariant subspace. For an AA domain with a internal bonds, the exchange terms are bounded below by -a/4, while the longitudinal terms telescope to its endpoints and are bounded below by -1. Thus its internal energy is at least -a/4-1.

[Proved] Two adjacent BB bonds are positive as a sum. Indeed, on three spins

    S_1 dot S_2 + S_2 dot S_3 >= -I,
    h_BB,12 + h_BB,23 >= 0.                                 (9.12)

To verify the first inequality, combine spins 1 and 3 into spin j_13 = 0 or 1. The operator is S_2 dot (S_1+S_3), with eigenvalues 0, -1 and 1/2. Pairing the internal bonds of any B domain then gives a lower bound -4 for the entire domain, with the possible unpaired bond supplying that -4.

[Proved] A periodic configuration containing both labels has an even number m >= 2 of mixed bonds, m/2 A domains, and m/2 B domains. If n_AA is its number of AA bonds, the total lower bound is

    E >= -n_AA/4 - m/2 - 4m/2 + 12m
       = -n_AA/4 + (19/2)m > -L/4.                          (9.13)

An all-B ring has nonnegative energy: sum (9.12) around the ring, counting every bond twice. An all-A ring has ground energy -L/4 because its longitudinal terms telescope and it is the ferromagnet. Therefore the full uniform chain has

    e0 = -1/4,

and exactly the all-A ferromagnetic ground space, with a nonconformal quadratic magnon branch. States containing B domains have a positive interface cost. No competing infrared branch was omitted.

[Proved, violation in the minimizing representative] Use the all-A coherent product state from (9.3) as a trial state of the full four-state weighted chain. Equations (9.4) and (9.5) now hold with J = lambda = 1 and M = 8. In particular

    liminf_(R -> infinity) D_R[u] >= epsilon pi/2,
    G_R[u] -> 0.                                           (9.14)

This disproves a universal quadratic-gradient lemma even in the unique representative minimizing the bond spread. It does not depend on a choice among degenerate minimizers. All constants and block dimensions are fixed before taking R -> infinity.

[Proved, gapped variant] Keep the BB and mixed blocks of (9.9), and replace AA by

    h_AA(g) = -(g/2)(Z tensor I+I tensor Z)
               + X tensor I-I tensor X,    0 < g <= 1.      (9.15)

Its eigenvalues are 0,0, plus or minus sqrt(g^2+4), all inside [-4,12]. The same proof gives M = 8 and a unique minimizing representative for every g. In an A domain the field contribution is bounded below by -g n_AA and the endpoint divergence by -2. Consequently a mixed label configuration satisfies

    E >= -g n_AA - 2m/2 - 4m/2 + 12m
       = -g n_AA + 9m > -g L.                              (9.16)

The all-B ring again has nonnegative energy. The uniform full chain therefore has the all-A all-up unique ground state, e0 = -g and gap 2g for rings of length at least three. Active spin flips attain this gap; mixed sectors cost at least 18, and an all-B sector lies at least gL above the ground energy.

[Proved] The exact weighted all-A ground state is a legitimate trial state of the full system. Equation (7.3), now with lambda = 1, gives

    liminf_(R -> infinity) D_R[u]/G_R[u] >= 1/(2g).           (9.17)

No exact evaluation of other weighted label sectors is needed for this lower bound. As g -> 0, it disproves any universal coefficient C M f(d) at fixed d = 4 and fixed M = 8, even for unique gapped uniform ground states and unique minimum-spread representatives.

[Proved, conformal variant at the level of the uniform lattice] Instead replace AA by

    h_AA(Delta) = S_1^x S_2^x + S_1^y S_2^y
                  + Delta S_1^z S_2^z + S_1^z-S_2^z,
    -1 < Delta < 0.                                        (9.18)

Its spectrum stays inside [-4,12] and its norm is less than 2, so M = 8 and uniqueness follow as above. Its periodic uniform restriction is the ordinary critical XXZ chain.

[Proved, exclusion of other label sectors for (9.18)] If an active open interval has s sites and ground energy E_open(s), product tiling of its ground state on the infinite active chain gives

    e_A <= [E_open(s)+||h_AA||]/s,
    E_open(s) >= s e_A-2.                                   (9.19)

The first inequality is variational, and uses that the expectation of each bond between adjacent tiles is at most ||h_AA||. The active energy density e_A is negative, as already follows from an all-up trial state when Delta < 0. Each B domain has energy at least -4, hence also at least s_B e_A-4. Summing over domains and mixed bonds gives E >= L e_A+9m. Thus the full ground state and its gapless infrared are those of the all-A XXZ chain, with the additional label excitations separated by a positive interface cost.

[Proved within the standard XXZ infrared matching] The all-A weighted trial ground energy gives at least the susceptibility response (9.7), with lambda = 1. It follows that

    liminf_(R -> infinity) D_R/G_R >= chi(Delta)/2,
    chi(Delta)/2 diverges as 1/[4(1+Delta)].                 (9.20)

This establishes the obstruction even within the class with a unitary CFT infrared, at the same physics standard as the XXZ dictionary. Its conformal coefficient A = v tends to zero. Equations (9.14) and (9.17) are independent exact lattice obstructions and do not rely on (9.20).

[Proved, remaining scope] A representative with a norm slightly larger than the minimum can remove the active-sector divergence by transferring it to the massive interface blocks. For (9.9), subtracting the global on-site divergence with o = P_A S^z changes the mixed energies to 12 plus or minus 1/2, while removing the AA improvement. Thus a theorem allowing an infrared-adapted, controlled-norm representative is not refuted by this construction. Such a prescription is different from minimizing the bond spread and needs its own matching proof.

## 10. Other nonconformal sectors and what z can control

[Proved, another exact z = 2 check] A hard-core-boson chain at its vacuum threshold has the local bond

    h_i = t[n_i+n_(i+1) - b_i^dagger b_(i+1)
                         - b_(i+1)^dagger b_i],    t > 0.    (10.1)

Its two-site eigenvalues are 0,0,2t,2t. The empty state minimizes all bonds, hence D(w) = 0 for all positive weights. The one-particle dispersion is 2t(1-cos k) = t k^2 + O(k^4). This includes the free quadratic-band vacuum edge after the usual spin or fermion identification. Neither this canonical z = 2 example nor the canonical ferromagnet violates the lemma.

[Heuristic, explicit scaling assumptions] In a one-dimensional scale-invariant infrared with dynamical exponent z, assume hyperscaling, a well-defined conserved energy density of dimension Delta_h = 1+z, and a regular static vacuum response. Fourier scaling of the integrated two-point function then gives

    S_IR(q) proportional to |q|^(2 Delta_h-(1+z))
                         = |q|^(1+z),
    D_IR per period proportional to R^(-z).                  (10.2)

For z > 1 this is smaller than the two-gradient allowance R^-1; for z = 1 it has the same power. This says nothing universal about the coefficient. A local analytic q^2 term from microscopic modes can dominate (10.2). At integer resonances logarithms require separate analysis.

[Heuristic, limits of the scaling statement] Energy-density improvements change the dimension of the operator being sourced. A gradient of an operator of dimension Delta_O can instead contribute a static power |q|^(2 Delta_O+1-z), if the associated response exists. Ground-state degeneracy, zero-energy sectors, dangerous irrelevant couplings, failure of hyperscaling, and singular equation-of-state response can invalidate (10.2). The ferromagnetic example is a concrete reason why z by itself cannot bound Def.

[Proved, scope of the search] The examples checked here are the ferromagnetic XXZ edge, its critical approach, the free quadratic-band vacuum edge, and independent gapped spins. No exhaustive classification of gapless nonconformal chains is claimed. A universal theorem for every such chain is neither implied by CFT matching nor established by these examples.

## 11. What a repaired research target must state

[Assumed, possible sufficient target] Select one representative h_*(H) for each gauge class, with a specified control of its norm by M. Require the spectral stress matching S for chains with conformal sectors, and prove the inverse-current-moment bound (6.5) for that representative. This would imply the desired symmetry-free bound on A with the explicit factor 24 pi. It need not assert equality between A and the complete microscopic response.

[Assumed, stronger target needed for equivalence] To call the full finite-amplitude smooth lemma an equivalent formulation, additionally establish P, or identify and subtract a proved source-response term B_h while controlling all gradient-induced relevant operators. Subtraction changes the lattice observable and must be stated. Choosing a representative that minimizes a bond norm does not, by definition, minimize its infrared susceptibility or remove every source contact term. Section 9.1 disproves those properties for the minimum-spread prescription, even when the minimizer is unique.

[Proved, final status] The pure-CFT calculation and both conditional directions are complete. The unconditional all-bond lemma is disproved, including the prescription of using its unique minimum-spread representative. The gapped estimate has explicit gap and correlation inputs, and the nonconformal examples distinguish genuine failures from harmless z = 2 points. The proposed symmetry-free bound on A itself is not disproved or proved by this work.

## 12. Checks and sources

[Numerical] checks.py verifies the exact gapped formula, the ferromagnetic coherent-state identity, two-site spectra and witnesses, and the cancellation of the circle term for a global conformal deformation. It uses only scalar arrays of at most 4096 entries and eigensystems of size at most 16. It also checks the four-state bond extrema and the three-spin positivity input (9.12). No many-body diagonalization, fitting, or heavy numerical job was used. Results are recorded in checks.json. These checks supplement the analytic proofs above.

[Numerical] With g = J = 1, lambda = 0.7 and epsilon = 0.4, D_gapped/G at R = 2048 is 0.2449998826, approaching lambda^2/(2g) = 0.245. The ferromagnetic lower bound at that R is 0.4374132295, while G = 0.00160930094, so the lower bound on D/G is 271.8032521. The limiting coherent-state lower bound is 0.4398229715. The continuum circle cancellation residual is -6.25e-17. All assertions in the script passed with python3; the initial attempt using the unavailable executable named python exited with code 127.

[Assumed, local background] The task context is specs/_bound_background.md, specs/_smooth_lemma.md, REPORT.md sections 9.68 to 9.77, contrib/fable_bound_modular/README.md, and codex/work_C105/{STATE,RESULT,PROOF}. Their uniform-chain results are used with their recorded qualifications. The derivative and gauge caveats in C105 are retained rather than promoted to an unconditional equivalence.

[Proved in primary literature, used as specified above] Fewster and Hollands, Quantum Energy Inequalities in two-dimensional conformal field theory, https://arxiv.org/abs/math-ph/0412028 , Theorem 4.1, supplies sharp stress bounds. Moosavi, Inhomogeneous Conformal Field Theory Out of Equilibrium, https://arxiv.org/abs/1912.04821 , Lemma 5.3, supplies the weighted-stress conjugacy. Sirker and Bortz, The open XXZ-chain: Bosonisation, Bethe ansatz and logarithmic corrections, https://arxiv.org/abs/cond-mat/0511272 , supplies the XXZ dictionary. The lattice counterexamples and reductions in this note are derived explicitly here; they are not attributed to those papers.
