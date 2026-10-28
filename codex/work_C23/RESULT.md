# C23. Route A from the microscopic bond to the infrared target

## 1. Verdict and conventions

[proved, scope of the assembled results] Route A is not yet a fixed-local-dimension counterexample. There is an exact family of finite matrix cells producing the required fermion hopping at second order, an explicit constant-alphabet history clock, a mobile wrapper theorem for initialized positive frustration-free cores, and a conditional stable infrared target. There is no supplied repeated two-site bond that satisfies all these specifications at once. Consequently neither a complete function b(n), a numerical fixed d0, nor the normalized accumulation bond is established.

[heuristic, physics conclusion conditional on the missing interfaces] A completed version can have n payload bosons and one length boson, hence c=n+1. Nothing in the checked finite-dimensional pieces or in C13/C20 forces an instability for every finite n. The strict weak-coupling inequalities have solutions if the required finite-cell resources and structural accuracy can first be obtained. This is a conditional existence statement about inequalities, not the missing existence proof for a microscopic Hamiltonian.

[proved, restricted failures] Three shortcuts do fail. Positive static Bell checks with port-idle terms that annihilate the code give scalar source resolvents and cannot synthesize C19. Equalizing first-order port maps does not equalize second-order exchange. Finally, finite-volume band perturbation theory does not compare its extensive remainder with conformal gaps of order 1/L. None is a no-go theorem for route A with a richer interface.

[assumed, conventions] n is the number of complex logical flavours, N=2^n their on-cell Fock dimension, k=n+1 the number of accessible positive spectral energies, b the short physical cell length, and b+1 the long length. M is the number of cells on a physical ring of L sites, r=L-bM, nu=r/M, and a=b+nu. A work-core length is denoted ell and a number of clock vertices C. Neither ell nor b is identified with n. Couplings are fixed before L tends to infinity. Energy dimensions are retained unless a numerical example explicitly sets the matrix-cell gap Delta=1.

[proved, status convention] Every substantive claim below is tagged. “Proved” means a derivation here or in the named source, not human refereeing. “Numerical” specifies what was computed. “Heuristic” identifies standard local elimination or RG reasoning. “Assumed” identifies an interface not supplied by the sources. All source paths are relative to this file. The source snapshots were read during the C23 run ending 2026-10-28; concurrent work is identified separately below.

## 2. Parts list and microscopic interface ledger

### 2.1 The components

| Component | Alphabet or Hilbert space | Repeated bond terms or matrix definition | Status and source |
| :-- | :-- | :-- | :-- |
| C7 selector | Four letters #, L, X, R | Forbidden-pair penalty KQ; half-onsite marker chemical potential; head rewrite LX <-> XR with amplitude -t | [proved] [C7, construction after Theorem 4](../work_C7/C7_result.md) |
| Mobile marker | Same four selector letters | #L <-> R# with amplitude -u, decorated by the vacancy conversion below | [proved] [mobile_selector §§1–2](../work_C8/mobile_selector.md), [selector_wrapper §5](../work_C8/selector_wrapper.md) |
| Positive initialized core, required input to wrapper | q_core local states on ell occupied sites | Positive uniform h_pay and endpoint a_pay,b_pay, with complete FF kernel of dimension 2^n | [assumed for C19] [arrow_wrapper §1](../work_C8/arrow_wrapper.md) specifies the hypothesis |
| Occupation filter | One qubit per occupied core site | Twisted singlet projector h_phi and two endpoint pins, phi=2pi/(ell+2) | [proved] [length_filter §§1–3](../work_C8/length_filter.md) |
| Vacancy and arrow wrapper | H_0,H_1,H_2, and L_alpha,B_alpha,R_alpha with alpha=1,...,2q_core; dimension 6q_core+3 | Arrow violations, vacancy domain-wall energy, vacancy/particle exchanges, lifted core terms, endpoint pins and rewards | [proved conditional on the core hypothesis] [arrow_wrapper §§2–4](../work_C8/arrow_wrapper.md) |
| Selector plus wrapper | 4(6q_core+3)=24q_core+12 | Selector bond plus P_int h_wrap; no wrapper term across the cut into # | [proved conditional] [selector_wrapper §§2–4](../work_C8/selector_wrapper.md) |
| Helix equalizer | Two additional qubits per physical site, total 96q_core+48 | Internal twisted-singlet bonds, pins at cuts, compensated endpoint observables; marker move gets a local helix rotation | [proved first-order identities only] [helix_equalizer §§1–3](../work_C8/helix_equalizer.md) |
| Full-sweep history controller | 15 clock letters times work dimension w; C18 adds separate left/right qubit letters, total 15w+4 | Decorated rewrite projectors or hopping rewrites; legality, initialization and energy events | [proved controller, assumed C19 program] [C8 run 7 §37](../work_C8/RUN7_RESULT.md), [C18 §§2–6](../work_C18/RESULT.md) |
| C19 matrix core | Two qubit ports, two partners, flag dimension n+2, logical dimension 2^n; D_ref=16(n+2)2^n | H_cell=Delta R H_diag R, defined in §3 below | [proved finite matrix, spatially nonlocal] [C19 §§3–7](../work_C19/RESULT.md) |
| C19 separated-end alternative | Separate flags at the two ends; D_sep=16(n+2)^2 2^n | Positive h_L+h_R, each touching one port and the shared logical register | [proved finite matrix] [C19 §6](../work_C19/RESULT.md); spatial locality still missing |
| Physical intercell link | Right qubit port of one cell and left port of next | lambda v, v=X_R X_L+Y_R Y_L+Z_R Z_L, norm 3 before wrapper compensation | [proved for the matrix cells] [C19 §6](../work_C19/RESULT.md) |
| Logical target | n fermion modes per ordered cell; this is an effective, growing alphabet | Complex hopping and chemical potentials, §3 | [proved second-order synthesis; heuristic stable IR] [C19 §7](../work_C19/RESULT.md), [C20 §§1–5](../work_C20/RESULT.md) |

[proved, explicit selector bond] The allowed pairs are #L,#X,LL,LX,XR,X#,RR,R#. With M_#=|#><#|,

    h_sel = KQ + (mu/2)(M_# tensor I + I tensor M_#)
                 -t(|LX><XR|+h.c.).                                  (1)

All coefficients are independent of the physical ring length. The arrow description #=(R,L), L=(L,L), X=(L,R), R=(R,R) proves conservation of each legality projector under both head and marker moves. Marker number is exactly conserved after adding motion. Payload parity is a separate obligation.

[proved, explicit wrapper bond] Let n_v project onto vacancies and let Pi_bad test mismatched arrow vertices. The charge bond is U(n_v I+I n_v-2n_v n_v) together with -v times the three rewrites H_0 L_alpha <-> L_alpha H_1, H_1 B_alpha <-> B_alpha H_1, H_1 R_alpha <-> R_alpha H_2. The full wrapper bond is

    h_wrap = 6U Pi_bad + h_charge + h0_tilde
             +(a_L-2U P_L) tensor I + I tensor (b_R-2U P_R),          (2)

where h0=h_pay+J_f h_phi includes the occupation-filter qubit. Choose 0<v<U/(ell+1). The selected occupied length is ell, and the two physical lengths are b=ell+1 and b+1=ell+2. Their complete codes are chi_p tensor V tensor g_ell, p=1,2, in ordered occupied coordinates. This is not just the two-placement ansatz from [padding_transport](../work_C8/padding_transport.md), which by itself has no local parent.

[proved, explicit assembled skeleton] Let P_int select allowed selector bonds whose right letter is not #, and P_cut select X# and R#. The skeleton bond is

    h_base = h_sel tensor I + P_int tensor h_wrap
                  + internal helix bonds + cut-conditioned helix pins,
    h_move = -u[|R#><#L| tensor (C_new tensor I)+h.c.],
    h_link = lambda P_cut tensor v_comp.                            (3)

Here C_wrap=|H_2><H_0| acts on the transferred site. With equalizer pitch phi, U_b=diag(1,exp(i b phi)) and C_new=exp(-i b phi) C_wrap tensor U_b tensor U_b. Every displayed operation is on one physical bond. Equation (3) is only a complete route-A bond after h_pay, its endpoint terms and v_comp with the required resolvents have been supplied.

[proved, controller specification] The 15 letters are L,0,1,C0,C1,R0,R1,I,Z,A0,A1,D0,D1,E0,E1. The forward rules are C1 x -> L Cx, C0 -> R1, L Rx -> R0 x, IZ -> LI, right I -> R0, left Rx -> Dx, Da x -> Aa Dx, right Dx -> Ex, Aa Ex -> Ea x, left Ex -> Cx. Each is accompanied by its adjoint and may carry a fixed work unitary. In positive form its term is (J_c/2)(|a><a|+|a'><a'|-|a'><a|U-h.c.). In hopping form retain just the last two terms with hopping coefficient -J_c. These are different parents and have different ground weights. C18 supplies the port boundary rules and full malformed-sector proof. The clock does not supply the work collision table that computes the C19 filters.

[proved, alphabet accounting limitation] A fixed w would make 15w+4 independent of n. Formally wrapping a core with q_core=15w+4 would cost 1440w+432 states. This multiplication is not a construction: C18's whole skeleton has not been shown to meet the wrapper's positive initialized-core hypothesis with the required spectral data. The alternative counts 9w+4, 27w+4 and 45w+4 in C19 refer to different bare clock/selector choices. Neither the historical estimates d0 approximately 45 or 130 in [minimal_d0](../../contrib/minimal_d0/README.md), nor D_ref, is an established microscopic d0.

### 2.2 Which interfaces actually match

| Interface | What has been checked | What remains assumed or fails |
| :-- | :-- | :-- |
| Finite code -> two ports -> second-order logical hopping | [proved] Exact kernel, parity, source measures and all second-order chain terms in C19 | Spatial NN realization with fixed alphabet |
| Work gates -> history clock | [proved] A specified fixed local unitary decorates the legal path, with illegal sectors excluded | The required spectral energy program, initialization and precision schedule |
| Port-idle history -> source resolvent | [proved] C18 weighted Schur-complement estimate under its exact-code hypotheses | Non-FF code splitting and any first-order port error if those hypotheses are relaxed |
| Static positive Bell events -> C19 | [proved obstruction] The source resolvent is scalar | This particular implementation cannot work |
| Initialized positive FF core -> wrapper | [proved conditional] Complete two-length code, all unwanted sectors excluded | A non-FF energy simulator is not automatically a valid input; richer positive events might be |
| Wrapper -> scalar marker transport | [proved conditional] Positive overlap times identity on both payloads | The C19 spatial code must meet the same hypothesis |
| Equalizer -> first-order port maps | [proved] Equal at both lengths | C19 has zero such maps; their equality says nothing about exchange |
| Equalizer -> second-order source measures | Not checked | Four length-pair resolvents, including excursions in helix and vacancy tracks |
| Matrix target -> C20 | [proved algebraically] Arbitrary phases/fields allow jointly nonresonant momenta | Physical translation/vertex gluing and dressed stiffness for the completed mobile cell |
| Projected mobile model -> finite-coupling IR | [heuristic] C13/C20 give a weakly coupled Gaussian phase under explicit conditions | Local massive-sector elimination for the actual parent and its coefficients |
| Mobile phase -> all-ring unique ground state | [proved in a positive-hopping subclass, heuristic for displaced parity-only target] | Complex amplitudes remove the supplied Perron-Frobenius proof; small rings and sector crossings remain |

[proved, length convention] C18's clock-selected cell is P_L [ell bulk letters] P_R, of length ell+2. C8's scalar-marker wrapper cell includes its first # and has b=ell_core+1. Their marker moves cannot be equated without a new rule table. Also adjacent raw clock lengths generally execute different computations. The two mobile lengths must contain the same n-flavour occupied core, not n and n+1 flavours. The vacancy wrapper accomplishes this for an eligible core.

[proved, attribution] G1 established exchange through small ports, and C9 established arbitrary finite matrix synthesis. Their symmetric SU(N) targets are not used. [graded_route_a](../../contrib/graded_route_a/README.md) repairs the parity compatibility of C9 and separates the ends, but does not compile a fixed alphabet. The new [exact_route](../../contrib/exact_route/README.md) preparation, [sweep interface](../../contrib/exact_route/SWEEP_INTERFACE.md), and [torus initialized-core construction](../../contrib/torus_bottleneck/EXISTENTIAL_CORE.md) supply useful acceptance and controller methods. Their sparse first-order family has growing local dimension O(q log q), as stated in [PROPOSED_COUNTEREXAMPLE](../../contrib/research_continuation/PROPOSED_COUNTEREXAMPLE.md). It does not complete fixed-d route A. Its [scale hierarchy](../../contrib/research_continuation/SCALE_HIERARCHY.md) makes every unwanted term higher order than its first-order target. That last property does not transfer to a second-order target.

## 3. The exact matrix cell and its logical couplings

[proved, C19] On n qubits set gamma_(2a-1)=Z_1...Z_(a-1)X_a, gamma_(2a)=Z_1...Z_(a-1)Y_a, P=product Z_a and eta_u=i gamma_u P. Global Jordan-Wigner strings give i Gamma_(j,u)Gamma_(j+1,v)=eta_(j,u)gamma_(j+1,v). The ordinary tensor-product logical bond is therefore

    K = sum_(a,u,v=1,2) J_a[u,v] eta_(2a-2+u) tensor gamma_(2a-2+v),
    J_a = -(t_a/2) [[sin vartheta_a, cos vartheta_a],
                    [-cos vartheta_a, sin vartheta_a]],
    F = -sum_a (mu_a/2) Z_a.                                         (4)

It represents

    H_* = sum_(j,a) [-t_a(exp(i vartheta_a)c†_(j,a)c_(j+1,a)+h.c.)
                         +mu_a(n_(j,a)-1/2)].                        (5)

To give band center theta_a and half-width k_a choose vartheta_a=-theta_a, mu_a=2t_a cos k_a. The dispersion is -2t_a cos(p-theta_a)+mu_a, with roots theta_a +/- k_a and speed magnitude 2t_a sin k_a. This sign conversion distinguishes C19's hopping phase from C20's band-center angle.

[proved, explicit filters] Put E_e=e, e=1,...,k=n+1, C_ef=1/(e+f), p_e=1/k. With j=0,...,n-1,

    Q_j(x)=sum_(l=0)^j binom(j+2,l)binom(j+1,j-l)x^l(x-1)^(j-l),
    h_j=(j+1)/[(2j+4)(j+3)],
    x(1-x)Q_j(x)/sqrt(h_j)=sum_(e=1)^k X_(e,j+1)x^e.

C19's Jacobi orthogonality proof gives sum_e X_ea=0 and X^T C X=I. Thus these are explicit finite functions of n, without an arbitrary exponentially large matrix whitening.

[proved, weights and visibility] Use Bell pairs on port/partner registers and a flag-zero logical code G. Its six Pauli sources are mutually orthogonal and G†F_sG=0. In X,Y,Z syndrome order set A_1=(|X><Z|+|Z><X|)/sqrt(2), A_2=i(|X><Z|-|Z><X|)/sqrt(2), and

    D_e^R=sum_(a,u) X_ea A_u tensor eta_(2a-2+u),
    D_e^L=-sum_(a,u,v) X_ea A_u^T tensor J_a[u,v] gamma_(2a-2+v),
    M_e^R=p_e I+rho D_e^R,
    M_e^L=p_e I+rho D_e^L+sigma d_e |Y><Y| tensor F,
    d=(1,-1,0,...), g_f=1/[2(k+2)], sigma=-rho^2/g_f,
    r_max^2=2 max_e sum_a X_ea^2,
    rho=min{1/(4k r_max), sqrt[g_f/(4k ||F||)]}.                    (6)

The second entry is omitted when F=0. For t_a<=2 all weights are at least 3I/(4k), and their sums are identity. For the t_a=1 target below, ||F||<=n, so the n-only conservative choice is

    rho_n^2=min{1/(16 k^2 r_max^2), 1/[8k(k+2)n]}.                  (7)

The logical visibility rho_n^2 is dimensionless. It is not a port occupation probability or a clock-readout fraction.

[proved, Hamiltonian and parity] Let S collect the flag-zero sources and S_e their copies in flag e. W=sum_e S_e sqrt(M_e^L direct_sum M_e^R) is an isometry orthogonal to S. The reflection R=I-SS†-WW†+SW†+WS† fixes G and exchanges S,W. H_diag is zero on G, e on flag e>=1 and k+1 on the flag-zero complement. Then

    H_cell=Delta R H_diag R,
    gap(H_cell)=Delta,
    ||H_cell||=(n+2)Delta,
    S†P_eS=M_e^L direct_sum M_e^R.                                 (8)

The physical product grading is Q=Z_pL Z_sL Z_pR Z_sR P. Both H_cell and the intercell Pauli bond preserve it. The flag is even. A spatial program must realize this as a repeated on-site involution and preserve it term by term.

[proved, complete second order] Insert the two source spectral measures in -PVH_0^(-1)VP. Because tr(A_u A_v)=delta_uv, the transpose on D^L is necessary. Jacobi orthogonality removes cross-flavour terms; the Y entry supplies the field. The result is

    H_eff,bond^(2)=(lambda^2/Delta)[a0 I+rho_n^2(K+I tensor F)],
    a0=-3 p^T C p,
    W_log=lambda^2 rho_n^2/Delta.                                  (9)

Each cell receives its field once on a ring. On an open chain this convention omits the first cell's field. Distinct bonds have no cross contribution at second order: an outer cell is acted on once and cannot return to its code. The result is the full chain target through that order, not just a two-cell spectrum match.

[proved, bandwidth distinction] A single flavour has full one-particle bandwidth 4t_a W_log. With t_a=1 the cell-coordinate Fermi velocity is 2W_log sin k_a and physical velocity a times that value. The many-flavour local energy norm grows with n; for example ||K+I tensor F||<=sum_a(t_a+|mu_a|/2)<=2n. The scalar energy a0 lambda^2/Delta is much larger than W_log when rho_n is small. In mobile cells a common scalar per cut shifts the marker chemical potential and must be retuned. It cannot be ignored as an L-independent constant when M varies.

## 4. The stable target, including the length fluid

[assumed, precise continuum input] Use C20's product compactification in cell coordinates. A local vertex is

    V_(m,A,B)=exp i[2m phi_0+sum_a(A_a phi_a+B_a theta_a^B)],
    m in Z, A_a=B_a mod 2, sum_a B_a=0 mod 2,
    q/(2pi)=-m nu+sum_a(B_a tau_a-A_a kappa_a) mod Z,                (10)

where tau_a=theta_a/(2pi), kappa_a=k_a/(2pi). Exact marker conservation excludes the length dual field. Exact payload parity does not exclude pairing. At the free point the dimension is m^2+sum(A_a^2+B_a^2)/4.

[proved conditional on this lattice, C20] In a neighborhood with scaling metric G >= (1-epsilon_G)I and epsilon_G<1/3, only charge vectors of squared norm 4 or 8 can be relevant or marginal. The rest have dimension at least 3(1-epsilon_G)>2. The number of conjugate vertex pairs is

    M_n=1+8n+32 binom(n,2)+96 binom(n,3)+128 binom(n,4),
    delta_n=1/(4M_n),                                               (11)

and a nonempty open set has every normalized mismatch greater than delta_n. Thus the angular margin is 2pi delta_n=O(n^-4). The maximum l1 coefficient norm is eight. Moving every normalized coordinate by less than delta_n/16 leaves a margin greater than delta_n/2.

[proved, an additional velocity bound] The C20 slab exclusion uses at most half the full parameter-domain volume. Restrict also nu to [1/8,7/8], removing one quarter, and kappa_a to [1/(32n),1/2-1/(32n)] for every flavour, removing at most one eighth by the union bound. A positive fraction remains. For t_a=1 this guarantees

    v_payload/W_log >= s_n=2 sin[pi/(16n)],
    v_length/J_len >= s_0=2 sin(pi/8).                              (12)

This proves simultaneous positive velocities and the stated all-n detunings. It does not rely on the geometrically crowded example in C19 §2. Infinitesimal irrational displacements can be made inside this open region. Excluding all arbitrarily high-order rational relations is unnecessary for stability.

[heuristic, RG requirements] If y is a dimensionless non-Gaussian coupling and q its reduced mismatch, its growth stops near length 1/|q|. A sufficient weak-coupling requirement is |y||q|^(-(2-d)) much smaller than one for a relevant primary of dimension d, and C|y|log(1/|q|) much smaller than one for a marginal one. At epsilon_G=1/4 the smallest allowed nonzero dimension is at least 3/4, so the worst power is 5/4. Derivative descendants have the same mismatch and are no more dangerous. Neutral currents shift the background; neutral current products change the Gaussian stiffness.

[assumed, conversion constants] A microscopic local-norm error is not automatically its bosonized y, its stiffness error, or its Fermi-point displacement. Denote a bound on these fixed-n conversion constants by C_RG(n,b)>=1. The sources give no universal numerical value. A deliberately small bookkeeping tolerance is

    epsilon_n=delta_n^2 s_n/[10^4(n+1) C_RG(n,b)].                    (13)

At C_RG=1, epsilon_n/delta_n^(5/4) and epsilon_n log(1/delta_n) are small in the examples below. The factor n+1 also allows a row-sum bound on stiffness corrections. Equation (13) is a conservative physics tolerance, not a proved microscopic threshold. Require separately G>=3I/4 and positive full stiffness H, including density-current blocks. The exact test is H>0, equivalently A>0 and B-C^T A^(-1)C>0 for H=[[A,C],[C^T,B]].

[proved within C13's reference dictionary] Physical translation can be checked using y(x)=x/a-phi_marker/pi and the marker comb sum_D exp(2pi i D y). A payload monomial of chiral charges z_i contributes S=sum_i p_i z_i/(2pi). Its physical momentum and free dimension are

    k_phys=2pi(D+S)/a mod 2pi,
    d=(D+S)^2/a^2+sum_i z_i^2/2.                                   (14)

Uniform resonance means D+S=h(b+nu), with h integer. Therefore sum_i p_i z_i-2pi h nu=0 mod 2pi and d_res=h^2+sum_i z_i^2/2. This agrees with the C20 cell resonance test after a consistent sign convention. The small physical marker Luttinger parameter 1/a^2 does not make infinitely many resonant operators relevant.

[assumed, remaining gluing interface] Equation (14) uses the scalar-transport reference quotient. A completed route-A parent must preserve its local operator lattice, physical translation character and ring twists, or supply their replacement. C20 explicitly leaves this input open. A boundary-state-dependent transport unitary, parity twist or change of compactification must be included before applying (10). The clock itself does not prove this interface.

## 5. Every available scale, and the missing resource functions

### 5.1 Selector and auxiliary gaps

[proved, bare C7 crossing] Put a_l=2t[1-cos(pi/l)]. At the b,b+1 crossing,

    mu_c=2t-[(b+1)a_b-b a_(b+1)],
    e_c=a_(b+1)-a_b,
    g_sel=min{a_(b-1)-2a_b+a_(b+1),
              a_b-2a_(b+1)+a_(b+2)} >0,
    omega_l=2t[cos(pi/l)-cos(2pi/l)].                               (15)

Discrete convexity makes the nearest rejected lengths the cheapest. Thus g_sel is the minimum wrong-length excess above e_c times length, rather than the energy-density difference. It is O(t b^-4), while omega_l=O(t b^-2). At a crossing the b and b+1 states are retained, not assigned a selector gap.

[proved, loaded wrapper selection] Use t=1 in the definitions a_l, mu_c,e_c for the following formula, and set

    eta_b=min{(a_(b-1)-2a_b+a_(b+1))/2,
              (a_b-2a_(b+1)+a_(b+2))/3},
    t_sel=2U[1/eta_b+b/mu_c+1/|e_c|],
    mu=t_sel mu_c+4U-(b+1)e_1+b e_2,
    e=t_sel e_c+e_2-e_1,
    K_sel=4t_sel+9U.                                               (16)

Here the two wrapper ground energies are -4U+e_p with 0<=e_p<=U. The wrong-length excess is at least (t_sel eta_b-U)(1+|l-b|), and illegal selector sectors cost at least U per violation. This is [selector_wrapper §§3–4](../work_C8/selector_wrapper.md), and includes arbitrarily long competing cells. Using the bare selector crossing without these energy shifts would not select the intended loaded cells.

[proved, auxiliary-resource formula] With ell=b-1 define kappa_*=1/[(ell+2)(2ell+3)], beta_*=sin^2(pi/(ell+2)) and

    epsilon_filter=kappa_* beta_*/[(1+sqrt(kappa_*))^2+beta_*],
    J_f>=2U/epsilon_filter,
    g_helix=J_h/[(b+1)(2b+1)].                                    (17)

The first controls the occupation filter. The last is the smaller gap bound for the two single-pinned equalizer tracks at strengths J_h. If gamma is the initialized occupied-core gap, d_p is the vacancy-charge gap, p_min its least required adjacent-particle probability, and B_core=(ell-1)||h0||+||a||+||b||, the source gives

    g_wrap,p >= min{
       d_p p_min gamma/[(sqrt(d_p)+sqrt(B_core))^2+p_min gamma],
       U-v(b+p-1)}.                                               (18)

All d_p,p_min and endpoint probabilities are finite eigenproblems once the parent is given, not known functions of n for C19. A lower bound for the unperturbed full mobile gap is

    g>=min{g_wrap(1),g_wrap(2),omega_b,omega_(b+1),
            t_sel eta_b-U,U,|e|b(b+1),g_helix}.                     (19)

Here g_wrap(p) denotes the two separate wrapper gaps. This gap is above the entire unperturbed tiling/code manifold, not above its eventual conformal ground state.

### 5.2 Clock scale and source accuracy

[proved, full-sweep clock] For ell clock bulk sites,

    C(ell)=(2ell+4)2^ell+ell-4.                                    (20)

The positive path has clock gap J_c[1-cos(pi/C)]. With the supplied initialized commuting pins one may use a conservative bound Gamma_init>=J_c/(8C^2), provided the legality and pin coefficients are chosen as in the source. Additional acceptance machinery needs its own bound. In the hopping version, with c=C+1,

    omega_clk=4J_c sin(3pi/(2c))sin(pi/(2c)),
    f_min=2 sin^2(pi/c)/c,
    Gamma_in>=omega_clk eta_pin f_min/(omega_clk+eta_pin).          (21)

With eta_pin proportional to J_c this conservative bound is O(J_c C^-5). An endpoint pin actually has a shift O(J_c C^-3), so substituting the O(J_c C^-2) propagation gap for every initialized gap is invalid. A clock-selected length uses a_ell=2J_c[1-cos(pi/(C(ell)+1))] and period ell+2. Its midpoint length margin is (a_(ell-1)-2a_ell+a_(ell+1))/2 and its competing density margin is that quantity divided by ell+3. This exponentially small alternative is not (15).

[proved under the port-idle hypotheses] If D is the energy-program perturbation, ||D||<=B_D, and the initialized clock gap is Gamma>B_D, the clock Schur correction is bounded by

    E_clk=B_D^2/(Gamma-B_D),
    ||R_sim(z)-R_target(z)||
       <= ||F_a||||F_b|| E_clk/[(Delta+z)(Delta+z-E_clk)].           (22)

This is C18 §6, including nonuniform hopping weights f_t. Energy events must average with those weights. Multiplying an endpoint event by C is correct only for uniform histories; hopping endpoint weights can be O(C^-3). With timewise code annihilation the code and its zero port compression are exact. Otherwise an O(B_D/Gamma) rotation and code splitting must be controlled as well. The phrase O(1/J_c) conceals C, initialization and B_D factors.

### 5.3 Marker visibility and equalizer norms

[proved conditional on the wrapper core] The marker amplitude is

    J_len=beta u,
    beta=A_b^2 beta_wrap cos(phi/2)^(2b),
    A_b=sum_(i=1)^(b-1) [2/sqrt(b(b+1))] sin(pi i/b)sin(pi i/(b+1)),
    beta_wrap=[sum_(s=0)^ell chi_1({s})chi_2({s,ell+1})]^2>0.       (23)

The cell-coordinate length bandwidth is 4J_len and velocity 2J_len sin(pi nu). Physical velocities gain a factor a. For phi=1/b the helix factor is exp[-1/(4b)+O(b^-3)]. It need not be exponentially small.

[proved, compensator cost] If r_1,r_2 are endpoint occupation probabilities, the positive helix compensator has norm at most (1/r_1+1/r_2)/sin^2(phi/2). A smaller Hermitian choice has norm (1/r_1+1/r_2)/2+|1/r_1-1/r_2|/[2sin(phi/2)]. Thus the physical port-bond norm B can be as large as 3||M_L||||M_R|| times the original lifted-port norms. This cost must be included in perturbation theory. It is not legitimate to set B=3 after compensation without checking it.

[assumed, second-order visibility] For a completed spatial core define j(n,b)>0 by its common intended non-scalar second-order coefficient W_log=j lambda^2. For the bare C19 cell j=rho_n^2/Delta. If each physical source is rescaled by a common factor alpha while its resolvent is otherwise unchanged, j=alpha^4 rho_n^2/Delta. In a real wrapper the resolvent need not change by a scalar at all. In that case j and the residual matrix must be computed, not inferred from r_p or history readout fractions.

## 6. Leading corrections and which ones can be dangerous

### 6.1 C19's cubic term

[proved, finite-bond bound] For the unwrapped matrix cells, PVP=0 and

    K3=PVQ H_0^(-1) QVQ H_0^(-1) QVP.

The first and last source states have both cells excited, with energy at least 2Delta. Consequently ||H_0^(-1)VP||<=3/(2Delta), while ||QVQ||<=3. Hence

    ||lambda^3 K3||<=27 |lambda|^3/(4Delta^2),
    ||lambda^3 K3||/W_log<=27|lambda|/(4Delta rho_n^2).               (24)

This holds for the C20 target as well as C19's original example. No fitted coefficient is needed. The bound controls the full cubic operator. A Hilbert-Schmidt projection onto Gaussian operators is not an operator-norm contraction, so it should not be used to claim the identical bound for the norm of its projected residual. Each normalized Majorana-monomial coefficient is bounded by (24); bounding a sum of such coefficients adds a finite n-dependent factor.

[proved, support check] On the infinite nearest-neighbor line a cubic product involving more than one distinct port bond has an outer cell hit exactly once. Its excited component cannot be removed by the other bonds or by the decoupled-cell resolvent, and its code component is zero. It vanishes. The only cubic bulk terms are repeated insertions on a single bond. A periodic three-cell ring is exceptional because its three bonds form a triangle. Mixed marker/port processes and a non-decoupled wrapper do not inherit this argument without verification.

[numerical, inherited evidence] C19 §10 finds non-Gaussian operator-norm residuals 2.16345e-7 and 1.55284e-9 in the cubic coefficient at n=1,2, Delta=1. These are two-cell full spectral computations for C19's original geometric target, with reflection-cell dimensions 96 and 256. They are not coefficients for the newly chosen C20 scaffolds and are not thermodynamic central-charge measurements. The stronger n=3 information here is the analytic bound (24), not new full ED.

[heuristic, infrared decision] Cubic terms do not merely renormalize hopping at the lattice scale. Their neutral quadratic pieces shift t,mu and the centers, their forward scattering changes the Gaussian matrices, and their non-neutral bilinear/quartic pieces can be relevant or marginal before oscillations cut them off. At the jointly detuned target, if their dimensionless strengths satisfy §4 and the dressed momenta stay within delta_n/16, these pieces remain perturbative and leave only Gaussian renormalizations and irrelevant operators in the infrared. If the centers are symmetric, an arbitrarily small allowed opposite-branch pairing can gap a flavour. The routeA_numerics [README](../../contrib/routeA_numerics/README.md) documents this failure with a D=32 matrix cell, two-/three-cell effective extraction and FES up to chi=128. That warning is applicable to the symmetry assumptions, not a numerical test of the present displaced target.

### 6.2 Clock errors

[proved conditional bounds] Let s be a code-splitting error, eta_port a first-order interface error defined so its bond contribution is bounded by |lambda| B eta_port, and dR a source-resolvent error bounded using (22). Sufficient requirements are

    s << epsilon_n W_log,
    |lambda| B eta_port << epsilon_n W_log,
    lambda^2 B^2 E_clk/g^2 << epsilon_n W_log.                     (25)

The last inequality is a conservative resolvent-identity estimate on a relevant gapped virtual block; fixed contraction constants can be included in B. If each individual port compression is O(eta), a factorized link may improve the middle estimate to O(lambda eta^2). The conservative form avoids assuming that factorization for a perturbed mobile code.

[heuristic, decision] Exact parity and marker conservation survive elimination if preserved by every term. Clock errors then lie in the same local operator classes as §4. They can be harmless after (25), retuning and the stiffness test. Without these inequalities, code splitting or a first-order anomalous term can dominate a target of order lambda^2 rho_n^2. Sending lambda to zero at a fixed splitting makes this problem worse.

### 6.3 Residual length dependence after equalization

[proved, quantity that must be computed] Let V_p encode the full length-b+p-1 code, and R_p,ab^s(t)=V_p†F_a^s exp(-tH_p)F_b^s V_p, with cell ground energies subtracted and excited projections understood. The four second-order matrices are

    K2_(p,q)=-sum_(a,b) integral_0^infinity
                     R_p,ab^R(t) tensor R_q,ab^L(t) dt.             (26)

After removing chosen scalar/chemical terms and defining the common target j H_bond, let d2=max_(p,q)||K2_(p,q)-K2_target||. Its contribution is lambda^2 d2, so its ratio to W_log is d2/j. There is no factor lambda left. First-order helix identities constrain V_p†F V_p, not (26).

[proved, a minimal obstruction to scalar compensation] Even suppose the only length dependence is a single accessible energy d_p>0. The resolvent factor is c_pq=1/(d_p+d_q). Arbitrary separate endpoint rescalings multiply rows and columns of this two-by-two matrix. To make all four entries the same requires c_00 c_11=c_01 c_10. Instead

    c_00 c_11-c_01 c_10
       =(d_0-d_1)^2/[4d_0 d_1(d_0+d_1)^2].                         (27)

It vanishes only when d_0=d_1. For d_0=1,d_1=2 it is 1/72. Thus equality of first-order zero maps, or compensating their occupation factors, cannot establish second-order equality even in this scalar example. Actual helix insertions also excite ancillary states, adding more poles.

[proved, available bound and its limitation] If the physical link norm is B and the excluded gap is g, ||K2_(p,q)||<=B^2/g. Comparing two such matrices gives at most 2B^2/g. Relative to j this is generically enormous, since j contains rho_n^2. It is not a useful small-error certificate. A possible sufficient target is d2<epsilon_n j after allowed retuning, but no source proves it for the C19 wrapper.

[heuristic, decision and possible repair] This mismatch need not kill the phase. Length-dependent scalar terms can be included in the length Hamiltonian; a common scalar per cell is a marker chemical potential. Non-scalar terms can average to changed payload couplings and produce Gaussian mixing. They can also generate resonant pairing or a negative stiffness. Their actual matrices are required to decide. Choosing a faster length fluid, rather than J_len=W_log, may suppress fluctuating mixed interactions and permit a retuned target without exact equalization. It does not by itself identify the averaged payload Hamiltonian. Equation (27) kills an inference, not the general use of mobile cells.

### 6.4 Marker excursions, mixed processes and leakage

[proved at finite perturbative order] For a gap g and link norm B, conservative local second-order scales are

    E_uu <= 8u^2/g,
    E_uλ <= 16B |u lambda|/g.                                      (28)

For a decoupled-cell chain the factor eight counts same-bond and ordered overlapping-bond products in the interaction norm sum_(X containing a fixed cell)||h_X||; expanding (u+B lambda)^2 gives the mixed factor. More general mobile-parent resolvents may require additional fixed-cell constants. These estimates are used as upper budgets, not claims that the mixed term is nonzero. Scalar transport after projection alone does not prove its vanishing outside the code.

[heuristic, decision] E_uu shifts length chemical potential, hopping and density interactions. E_uλ can give payload-dependent marker hopping and new mixed vertices. Their neutral current products renormalize stiffness, but their oscillatory parts must obey the mixed C20 exclusions. If u is held of order lambda while the payload begins at lambda^2 rho_n^2, these errors need not be small. Matching marker hopping to the logical scale is a sufficient, though not necessary, way to avoid that problem once structural exchange is controlled.

[proved on finite clusters] First-order leakage amplitude is bounded by O((B|lambda|+|u|)/g), with probability of its order squared. The C19 two-cell source estimate is the sharper 3|lambda|/(2Delta). This bounds local dressing on a fixed cluster, not the norm difference between many-cell ground states. A large tensor product can have order-one total overlap loss from arbitrarily small per-cell dressing.

[heuristic, physical input] A positive dressed cost for excluded local excitations and frequencies omega<<g permit local integration with retarded kernels decaying on time g^(-1). Their low-frequency expansion produces the same parity-even local operators and derivative corrections. It can preserve the Gaussian phase under the foregoing bounds. Neither a finite-cluster SW expansion nor the inherited data prove an exact volume-uniform block diagonalization; [ir_control §5](../../contrib/ir_control/README.md) explicitly discusses this limitation. At the requested physics standard local massive-sector integration is a legitimate stated input once an actual parent and local bounds exist. It is not a substitute for the missing parent.

## 7. A conditional all-n scale ordering, and why it is not yet the full ordering

[assumed, structural input] Fix n and a finite b for which an actual core has been supplied. Suppose its dressed local eliminated gap is bounded below by g>0, its physical link norm is B, its marker overlap is beta>0 and its intended second-order coefficient is j>0. Suppose code splitting, port errors and residual matrices (26) can meet (25) and d2<epsilon_n j, with the same physical resources. Fixed geometry constants can be absorbed in the error budgets below. These are the unsatisfied premises, not conclusions of this section.

[proved, elementary conditional feasibility] A sufficient explicit choice is

    lambda=(epsilon_n/1000) min{
        g/B, g^2 j/B^3, beta g/B, beta sqrt(g/j)},
    W_log=j lambda^2,
    u=W_log/beta,             J_len=W_log.                         (29)

Every quantity is positive for finite n,b. The units in each entry of the minimum are energy, since j has inverse-energy units. In particular B lambda/g<=epsilon_n/1000. Allow a conservative cubic local budget 100B^3|lambda|^3/g^2. Dividing the unwanted terms by W_log gives

    E3/W_log <= epsilon_n/10,
    E_uλ/W_log <= 16 epsilon_n/1000,
    E_uu/W_log <= 8 epsilon_n^2/10^6.                              (30)

For the last inequality use lambda<=epsilon_n beta sqrt(g/j)/1000. The finite-range cubic constant can be checked on a decoupled-cell line. The canonical third-order expression has one PV R V R VP term and an anticommutator with two half-weight folded terms, so each ordered triple has norm at most twice the product of its bond norms divided by g^2. There are 44 connected ordered triples whose union contains a fixed cell: two single-bond triples, 3 times 6 triples using two adjacent bonds, and 4 times 6 triples using three consecutive bonds. Thus 88 is an upper bound and 100 is sufficient in this setting. Extending this connected-support/resolvent argument to the moving parent is the stated local-elimination hypothesis. Mixed cubic terms carry additional u/(B lambda) or an additional small perturbation ratio. Equation (24) is the sharper rigorous C19 bound.

[proved, remaining clock arithmetic] If Gamma=J_c gamma_clk for a known initialized-clock bound and B_D is independent of J_c, define

    E_req=min{epsilon_n W_log/100,
              epsilon_n j g^2/(100B^2)},
    Gamma=B_D+max{B_D^2/E_req,
                  400B |lambda| B_D/(epsilon_n W_log)},
    J_c=Gamma/gamma_clk.                                         (31)

Then E_clk<=E_req. A code rotation bounded by 4B_D/(Gamma-B_D) also obeys the middle condition of (25) with a factor 100 of slack. Any intrinsic circuit-approximation error must be added separately. Exact-code interfaces omit the splitting and rotation requirements and can be much cheaper. All coefficients remain fixed as L grows.

[proved, what this establishes] At the level of independently supplied resources, there is no algebraic contradiction between n^-4 detunings, small rho_n, small g, marker motion, cubic errors and finite clock coefficients. Very small positive lambda suffices for every n. The required accuracy for a second-order structural error, however, is d2/j, and cannot be improved by that choice of lambda.

[assumed, unresolved resource compatibility] The resources are not actually independent in the proposed spatial construction. Increasing accuracy may enlarge b and C; increasing J_c changes ||h0|| and the wrapper bound (18); a non-FF code can change its eigenvectors and splitting. Equation (31) therefore cannot be called a proof of consistency for the full route. For example, a crude certificate g_wrap proportional to 1/J_c, combined with a generic splitting proportional to 1/J_c, cannot establish splitting much smaller than lambda^2 while also requiring lambda much smaller than g_wrap. This is a failure of that set of sufficient estimates, not a physical no-go theorem: the actual gap may be larger, or exact annihilation may eliminate splitting. A common-kernel positive clock term allows a monotonicity lower bound from a fixed J_c instead, but a general non-FF simulator does not inherit it.

[heuristic, alternative scale option] If a length-only second-order term is O(lambda^2/g), it need not be suppressed below W_log. One may arrange lambda^2/g << J_len << lambda j g/B, so the length correction is small relative to its own bandwidth and mixed corrections remain small relative to W_log. At each fixed j,g,beta such a window exists for sufficiently small lambda, with additional u^2/(gW_log) control. This illustrates why (29) is sufficient rather than necessary. It does not cure an unknown averaged payload matrix or provide the missing bounds.

## 8. Numbers for n=2 and n=3

### 8.1 Certified finite-matrix and resonance data

[proved and finite arithmetic checks] The script [check_assembly.py](check_assembly.py) checks sum-zero and all Jacobi Gram entries using rational arithmetic. For n=3 the new third column is sqrt(40/3)(3,-21,39,-21); the first two columns are the n=2 columns padded by zeros. The values are

| Quantity | n=2 | n=3 |
| :-- | --: | --: |
| N=2^n | 4 | 8 |
| Accessible energies, in units Delta | 1,2,3 | 1,2,3,4 |
| D_ref | 256 | 640 |
| D_sep | 1024 | 3200 |
| r_max^2 | 1200 | 41160 |
| rho_n^2 | 1/172800 | 1/10536960 |
| rho_n^2 decimal | 5.7870370e-6 | 9.4904033e-8 |
| C20 conjugate candidates M_n | 49 | 217 |
| Normalized margin delta_n | 1/196 | 1/868 |
| Angular margin 2pi delta_n | 0.032057068 | 0.0072386928 |
| Allowed normalized coordinate shift delta_n/16 | 0.00031887755 | 0.000072004608 |
| epsilon_n with C_RG=1 | 1.7009777e-10 | 4.3403933e-12 |

[proved by integer enumeration for these two cases] A rational scaffold with common denominator 1000003 is specified by numerators in the order (nu,tau_1,...,tau_n,kappa_1,...,kappa_n):

    n=2: (728872,330410,394497,127563,292196)/1000003,
    n=3: (638850,31803,808569,515884,280072,166260,327734)/1000003.

The minimum normalized mismatches over the complete shells are 16005/1000003 and 7486/1000003, respectively. The check uses integer residues, not floating-point near-zero tests. Both exceed (11). An optional change nu -> nu+sqrt(2)delta_n/1000 gives an irrational density without spoiling the margins; its normalized displacement is much smaller than delta_n/16. High-order commensurabilities in the unshifted rational scaffold are irrelevant near the stated free metric.

[numerical evaluations of the explicit scaffold, not ED] With t_a=1 and vartheta_a=-theta_a, the parameters of (5) are

| n | theta_a, radians | k_a, radians | mu_a | 2 sin k_a |
| :-- | :-- | :-- | :-- | :-- |
| 2 | (2.076021,2.478690) | (0.801500,1.835916) | (1.391260,-0.524050) | (1.436800,1.930122) |
| 3 | (0.199824,5.080374,3.241385) | (1.759739,1.044639,2.059207) | (-0.375641,1.004428,-0.938446) | (1.964407,1.729487,1.766159) |

[proved, conservative cubic condition] For a single bare matrix bond, demanding its cubic error be at most epsilon_n W_log is satisfied by |lambda|/Delta <=(4/27)epsilon_n rho_n^2. These upper values are approximately 1.46e-16 for n=2 and 6.10e-20 for n=3 when C_RG=1. They are sufficient norm bounds, not optimized physical thresholds. No measured CFT is inferred from them.

### 8.2 Length-dependent scales and an arithmetic feasibility example

[numerical evaluations of proved bare formulas] To make the b dependence concrete, take b=8 and b=9 as illustrative lengths, respectively. These are not compiled lengths for the n=2 and n=3 target. With bare-selector t=1 and equalizer J_h=1,

| Quantity | b=8 | b=9 |
| :-- | --: | --: |
| Wrong-length excess g_sel | 0.0088983855 | 0.0058548764 |
| Loaded-selector coefficient eta_b | 0.0029661285 | 0.0019516255 |
| Smaller head gap omega_(b+1) | 0.34729636 | 0.28407904 |
| Equalizer gap bound 1/[(b+1)(2b+1)] | 0.0065359477 | 0.0052631579 |
| A_b | 0.97545761 | 0.98035995 |
| C(b-1) | 2307 | 5124 |
| Conservative positive-clock gamma_clk=1/(8C^2) | 2.3486312e-8 | 4.7609297e-9 |

[assumed, deliberately explicit artificial resource assignment] The following table evaluates (29)–(31) under g equal to the smaller bare/helix bound above, j=rho_n^2 with Delta=1, B=3, beta=A_b^2 cos(1/(2b))^(2b)/4, B_D=C(n+2), gamma_clk=1/(8C^2) and C_RG=1. Thus it assumes the ideal padding overlap 1/4, no loss of spectral visibility, a wrapper gap at least g and no extra compensator norm. These assumptions have not been verified together and need not hold for C8's actual finite-v wrapper. The assignment tests algebraic consistency only. It is not a parameter proposal for a constructed bond.

| Scale or bound | n=2, b=8 | n=3, b=9 |
| :-- | --: | --: |
| Assumed beta | 0.23055591 | 0.23369057 |
| lambda | 1.5574290e-24 | 4.2261293e-28 |
| W_log=J_len | 1.4036951e-53 | 1.6950021e-62 |
| Marker coefficient u | 6.0883067e-53 | 7.2531899e-62 |
| Per-flavour and length bandwidth 4W_log | 5.6147804e-53 | 6.7800083e-62 |
| Conservative cubic budget | 2.3876541e-64 | 7.3569757e-75 |
| Mixed u lambda budget | 6.9636583e-73 | 2.7955462e-85 |
| Marker u^2 budget | 4.5370594e-102 | 7.9965322e-120 |
| Clock Schur budget E_req | 2.3876541e-65 | 7.3569757e-76 |
| Sufficient J_c for the assumed B_D | 1.5185494e80 | 1.8739895e92 |
| Initialized clock gap Gamma | 3.5665126e72 | 8.9219325e83 |
| Leading leakage-amplitude budget (3lambda+u)/g | 7.1485992e-22 | 2.4088937e-25 |

[proved arithmetic, conditional interpretation] Every comparison in (30) passes, and the relevant growth budgets epsilon_n/delta_n^(5/4) are 1.2474e-7 and 2.0449e-8. The marginal budgets epsilon_n log(1/delta_n) are 8.9780e-10 and 2.9368e-11. The huge clock coefficients illustrate why absolute and bond-normalized gaps must be distinguished. Making such a coefficient large does not prove that the assumed wrapper gap remains fixed.

[proved, limits of the numerical answer] Actual values of b(n), d0, g_wrap, j, B, beta and required B_D for the route-A spatial parent are absent. Therefore no honest table can give all physical n=2 or n=3 scales of that parent. The exact matrix/scaffold table and the explicitly conditional resource table are the available answer, rather than invented compiler data. The numerical work here took less than a second, used NumPy and Python rational/integer arithmetic, and simulated no thermodynamic chain.

## 9. The normalized accumulation point

[proved, bare selector] Remove additive bond scalars and normalize before invoking compactness. For the four-state C7 family with bounded t,K, mu_c tends to 2t and

    a_b=pi^2 t/b^2+O(t/b^4),
    2t-mu_c=3pi^2 t/b^2+O(t/b^3),
    e_c=-2pi^2 t/b^3+O(t/b^4).                                    (32)

At mu=2t the normalized limiting bare bond is proportional to

    KQ+t(M_# tensor I+I tensor M_#)
            -t(|LX><XR|+h.c.),                                    (33)

with the chosen scalar subtraction understood. Marker-free all-L and all-R states have zero energy, while one length-L cell costs 2t[1-cos(pi/L)] approximately pi^2 t/L^2. The head dispersion near its minimum is tk^2. This is a z=2, nonconformal finite-size scale. It is not an infinite-central-charge finite-velocity CFT. C7 proves this limiting selector statement; C13 §9 keeps its scope explicit.

[proved, length must diverge] A core on ell sites with fixed alphabet q_core has dimension q_core^ell. Hosting a 2^n-dimensional code requires ell>=n log(2)/log(q_core). Thus b must diverge in a fixed-alphabet realization; treating n as growing inside one fixed finite cell is impossible. This dimensional bound is necessary, not a compiler runtime or accuracy estimate.

[heuristic, natural full-route scenario] If a completed normalized h_n retains the selector kinetics while the spectral core and link scales vanish, cell density tends to zero, logical velocities a W_log tend to zero, and the limiting spectrum can contain the dilute z=2 selector branch plus collapsed logical sectors. Such a limit is consistent with unbounded c_n at every fixed n because the conformal onset and energy scales need not be uniform in n.

[assumed, unknown actual limit] No full h_n on one common fixed alphabet has been supplied, and no bound forces the largest coefficient to remain the selector hopping. The clock, occupation filter or penalties can dominate normalization. In the artificial hierarchy of §8, J_c is enormous, so assuming (33) survives normalization would be unjustified. The limit might instead retain positive clock constraints with an extensive zero-energy space, static marker sectors, or several different vanishing scales. The exact equal-hopping d=4 selector is another family, with threshold |mu|=4a_hop, not mu=2t. These limits cannot be interchanged.

[proved, compactness conclusion] If a fixed-d family is eventually constructed, normalized bonds have a convergent subsequence. That fact does not identify the limiting bond or its dynamical exponent. The continuum c-theorem alone also does not rule out a finite-c lattice accumulation point without a uniform crossover/basin assumption, as C13 §9 explains. Thus the full accumulation point requested in the task is currently undetermined. The only identified limit is the bare selector (33).

## 10. Rings, conformal interpretation and remaining checks

[proved, unperturbed geometry] C8's b/(b+1) tilings exist for every L>=b(b+1), including primes. The wrapper theorem excludes all other sectors at zero intercell coupling under its hypotheses. It does not settle the finitely many smaller rings or uniqueness after adding complex payload hopping.

[heuristic, expected infrared spectrum] With the verified local gluing, positive stiffness and the weak-coupling conditions, the mobile theory has n+1 right-moving and n+1 left-moving bosonic branches. It has entanglement central charge n+1, with generally different velocities. A velocity-weighted Casimir energy alone measures a sum of c_a v_a, not c_total. Entropy and independently resolved velocities or low-temperature free energy are appropriate complementary estimators.

[proved, finite-ring qualification from C13/C20] For irrational nu, finite rings use integer sectors with r/M tending to nu. Charge/current zero-mode offsets vary with L. The resulting towers need not approach a single untwisted vacuum spectrum along every integer L, even when every finite-ring ground state is unique. Generic uniqueness for all L requires avoiding a countable set of sector-crossing loci, not just belonging to an open thermodynamic phase. Complex displaced hopping lacks the supplied positive-matrix uniqueness proof. If the conjecture demands a fixed isolated CFT identity sector on every integer ring, this incommensurate construction does not meet that stronger interpretation. Under the usual Luttinger-liquid convention retaining shell offsets, it remains viable.

[numerical, inherited scope] C13 §8 calibrates soluble n=0,1,2 projected models at M=40,64,96,128,192 with cell entropy and shell-corrected velocity-weighted Casimir estimates. At M=192 the values for n=2 are 3.0006723 and 3.0000937. These are not measurements of the C19 spatial parent. Its generic physical-graph interacting tests at L=14,17,19 already have substantial free entropy-fit bias, so they do not establish the generic interacting c. No new CFT measurement was attempted here.

## 11. Remaining gaps, ranked and attached to concrete calculations

1. [assumed, decisive] A fixed-alphabet local spectral core. Supply the actual repeated-bond collision/energy rules, complete initialized code, product parity and malformed-sector exclusion. Compute the boundary spectral measures, not only eigenvalues. C19's target is (6)–(8), with n+1 accessible energies. The pole-count requirement 2n<=4(k-1) excludes a fixed-pole repetition trick. The smallest useful calculation is a spatial n=2 instance with the complete source resolvent compared to C19, followed by a proof that its rule alphabet stays fixed as n grows.

2. [assumed, equally structural] Compatibility with the positive mobile wrapper and its excited spectrum. Either implement C19 using richer positive port-touching events that evade the static-Bell obstruction while retaining an initialized FF core, or extend the wrapper proof to the proposed non-FF core. Calculate both complete length codes, all other occupied/tag sectors and a bound like (19). A finite matrix reflection alone supplies neither. For a signed energy program, explicitly bound the most negative energy density at every competing length before retuning the selector.

3. [assumed, potentially decisive for the present assembly] Four second-order mobile exchanges and marker/port cross terms. Compute K2_(1,1), K2_(1,2), K2_(2,1), K2_(2,2) from the actual full parent, along with the u lambda term and scalar shifts. Separate common retunable terms from residual length dependence. Test d2/j, rather than d2 or lambda alone. If scalar compensation fails as in (27), try common source measures, additional controlled virtual terms, or a faster length mode with a directly calculated averaged target. This is the next calculation that decides whether the existing equalizer can be used at second order.

4. [assumed, high importance] Close the resource inequalities for one family. Provide b(n,accuracy), w, B_D, initialized gap, wrapper gap, B, beta and j from the same rules. Substitute them into (25) and (29)–(31), including their dependence on J_c. An exact-code construction removes the most severe code-splitting requirement. Otherwise solve the coupled inequalities explicitly; an independent choice of g and J_c is insufficient. This calculation also gives the actual n=2,3 parameter tables and normalization.

5. [heuristic with unresolved microscopic coefficients] Finite-coupling infrared transfer and physical gluing. Derive the local translation action and winding/parity lattice of the completed mobile code. Expand its leading two-/three-cell corrections in Majorana monomials, retune momenta into the strict C20 region, and compute the full Gaussian stiffness including density-current mixing. Check the slowest mode and all norm-4/norm-8 mismatches. Use local massive-sector integration as the stated physics input; do not compare an extensive finite-size SW error with 1/L levels.

6. [assumed, secondary to existence but necessary for the requested full statement] All-ring uniqueness and the accumulation Hamiltonian. Check that the compiled rules add no projective or antiunitary degeneracy, address small rings, and specify shell-dependent conformal sectors. Then normalize the actual h_n and compare coefficient ratios. Only then decide whether its limit is (33), a constraint-dominated zero-mode space, or another gapless model. After a microscopic family exists, run enough cells for a calibrated physical entropy estimate and an independent velocity-resolved estimator.

[assumed pending the completed source derivation, concurrent reports] At the final read, [C21 PROGRESS](../work_C21/PROGRESS.md) reports a stronger boundary-moment obstruction: C19's left/right first moments do not commute for its original target, whereas disjoint boundary moments compressed into a product subspace do commute, even after a global logical change of basis. It also reports delayed Jacobi filters x(1-x)^s P_j^(2s,1)(2x-1), with k=n+s poles and scalar moments through order s-1, and a local pole ladder whose unitary-decorated tree returns remain scalar. Its exact rational filter checks cover n=1,2,3,5,8 and s=1,2,4,8. These are reported finite advances, not a spatial compiler. [C22 PROGRESS](../work_C22/PROGRESS.md) is testing a looped excitation bus after identifying the same tree-return obstruction. Neither has supplied a successful replacement core at this snapshot. These calculations were not duplicated here, and the numerical hierarchy above continues to use C19's original filters, not C21's delayed ones.

[proved, final accounting] The integration identifies no incompatibility that kills every route-A implementation. It does invalidate treating the current pieces as an established counterexample. The new finite results are the explicit n=2,3 target checks, a conservative cubic control at those margins, and the second-order scalar-equalization obstruction. The outstanding finite-dimensional tasks are more important than refining the RG estimates. All new files are inside work_C23. Project STATE.md and other agents' files were left unchanged because this task explicitly forbids writes outside this folder.
