# C54. Adversarial audit of input S

[proved, scope] This audit supplies an amplitude ledger, an exhaustive low-dimension vertex test conditional on S's operator lattice, and explicit failures of weaker substitutes for S. It does not prove or disprove the microscopic Gaussian phase of the scheduled C39/C44 chain. In particular, an upper bound on a possible coefficient does not prove that coefficient is generated. No nonzero coefficient is inferred from a symmetry allowance alone.

[proved, provenance] The inputs are C42 PROOF (7)-(11), C44 PROOFS (3)-(5), (9), (13)-(24), C45 PROOF Sections 3-8, C46's thermal comparison, and C47's demonstrator reduction. These finite results are imported, not independently re-proved here. The formal S statement in contrib/expert_note/note.tex was read. The newer matching, frozen-defect, physical-ensemble and counterexample-audit campaigns were also read. C48 and C53 had no progress file when inspected; C49-C52 had active progress files. This task does not duplicate their matching or constructive RG programs.

[proved, tags] Statements labeled conditional use the stated compact-boson/matching input. Numerical checks are small algebra controls or explicitly attributed published-in-repository diagnostics, not microscopic simulations. An assumed XXZ formula is marked separately. All writes are confined to this folder.

## 1. Units and the complete amplitude budget

[proved, imported normalization] Fix n before all limits. Write

    p = 2n+4,  V0 = q^(2n+2)/S_q^2,
    J0 = 2 V0 g^2 D^(-p),  c_n = J0/s,
    W = lambda^2 c_n,  J_len = W,
    u_mov = W/beta_n,  beta_n > 0.

The last normalization is C45's schedule. C44 also permits other positive length hoppings. A velocity in the following discussion is in cell coordinates unless stated otherwise. Its physical-site value is multiplied by b+nu, b=n+3. That coordinate factor does not alter scaling dimensions.

[proved, imported bounds] The normalized quadratic finite-D error and complete non-Gaussian conversion error obey

    E_D <= 2n B0(1+eta)/(D r0^(p+1)),
    E_ng <= 3(1+eta)(2+64n^2) g^2 D^p,
    r0=1/16384, B0=128.                                      (1)

The finite-g quadratic error is bounded by 48n(1+eta)g^2 D^p. These are errors relative to J0 before wrapping, hence relative to W afterwards. The full g^4 conversion remainder in the unscaled bond is at most 6(1+eta)V0 g^4. This retains the endpoint visibility.

[proved, imported equalization bound] With C44's x,y,L_c,L_b,a and t=s/delta_fast,

    epsilon_eq = x^2 zeta^2 L_b/a^2 + 2t L_c^2 x^2/a^2
                 +(2/a)(2t L_c x y+t^2 L_c^2 y^2)+t y^2,
    ||s K2_pq-K_bare|| <= epsilon_eq.                         (2)

The exact table decomposition is

    E_pq = E_0+l_i E_L+l_(i+1) E_R+l_i l_(i+1) E_LR,
    ||E_0|| <= epsilon_eq, ||E_L||,||E_R|| <=2epsilon_eq,
    ||E_LR|| <=4epsilon_eq.                                  (3)

Normalized trace extracts the scalar table without increasing its norm. Subtracting that trace increases the remaining bound by at most two. Consequently

    |V_len| <= 4W epsilon_eq/J0,
    ||mixed traceless LR term|| <= 8W epsilon_eq/J0.          (4)

The scalar V_len multiplies l_i l_(i+1). It is not a marker chemical potential. Common scalar and linear terms can be retuned; this one cannot.

[proved, normalization audit] Comparing epsilon_eq only with g^2 fails. Its correct divisor is J0, which includes D^(-p) and V0. Likewise reducing lambda cannot decrease E_ng or epsilon_eq/J0. C45 explicitly chooses g, the vacancy ratio, and s/delta_fast first. This candidate failure therefore attacks a weaker schedule, not C45's actual one.

[conditional on massive matching] Put r_SW=C_SW lambda^3. The low-energy correction from the remaining full-space residual is bounded by

    r_LE/W <= A_dec,n C_SW lambda/c_n.                       (5)

Before matching only r_SW is proved; replacing it by an exact retained-action bound without S or another matching theorem would be circular. A convenient conservative coefficient budget is

    E = C_b [E_D+E_ng+16 epsilon_eq/J0]
          + A_dec,n C_SW lambda/c_n,
    C_b=16(n+1)^2 exp(48(n+1)).                              (6)

One can increase C_b to include a chosen normal-monomial convention. It pays for spatial weights and coefficient conversion. Equation (6) deliberately overcounts the table entries and is not a claim that C45's numerical allocation sums to epsilon_n/4. If needed shrink its static budgets by the fixed extra factor 16. All choices still have nonempty positive intervals. Every additional non-Gaussian coefficient in that norm is at most W E. Numerical lower bounds or signs are unavailable for most individual vertices.

## 2. The low-dimension operator test

[assumed, operator lattice] Use the lattice explicitly required by S, rather than claiming it follows from U(1) alone. A local primary is labeled by integer h and z in Z^(2n), with

    sum_i z_i=0,
    d0(h,z)=h^2+(1/2)sum_i z_i^2,
    Q(h,z)=sum_i z_i p_i-2pi h nu modulo 2pi.                 (7)

Here h is a length density harmonic. It is not a length creation charge. Exact marker number excludes the latter dual-field vertices. Payload charge conservation imposes sum z=0. Identity transport and locality are needed to exclude additional fractional or twist operators; neither total charge alone nor a finite gap proves that exclusion.

[proved within (7), exhaustive classification] If d0<=2, then |h|<=1. For h=0, the nonzero z have squared norm two or four. Norm two gives one +1 and one -1. At norm four an entry +/-2 would already exhaust the squared norm and violate zero total charge, so there must be two +1 and two -1 on four distinct branches. For |h|=1 the squared norm is at most two, hence z=0 or a bilinear difference. These exhaust the list. A neutral vector with repeated incoming and outgoing labels first cancels those labels and reduces to one of these cases or zero.

[proved within (7)] Descendants increase dimension by a positive integer. The zero primary has currents of dimension one and quadratic current products of dimension two. They change fillings or the Gaussian stiffness. A derivative of a dimension-one nonzero primary retains that primary's nonzero momentum. Thus descendants add no untested marginal resonance. Higher current products and band curvature have dimension at least three.

[proved, dimension and momentum ledger within (7)]

| Mechanism | Free dimension | Momentum or exact rule | Size of additional coupling |
| :-- | :-- | :-- | :-- |
| Payload pair creation, e.g. psi_aR psi_aL | 1 | Nonzero total payload charge | Exactly zero |
| Interbranch bilinear psi_i^dagger psi_j, i != j | 1 | p_j-p_i, excluded by endpoint margin | At most W E after coefficient conversion |
| Pair transfer aR,aL to bR,bL | 2 | 2(theta_a-theta_b), covered by C42 pair sums | At most W E; conversion part bounded by W C_b E_ng |
| Other four-distinct-branch quartics | 2 | Difference of two distinct endpoint pair sums | Same bound |
| Repeated-branch spinless umklapp, 2R to 2L | 4 | 2(p_R-p_L); C42 also detunes it | Same coefficient-norm bound |
| Pure length harmonic h=1 | 1 | 2pi nu, bounded away from zero in the chosen interval | At most W E if generated |
| h=1 times a payload bilinear | 2 | p_i-p_j-2pi nu, covered by mixed margin | Static mixed part at most 8W epsilon_eq/J0 in operator norm |
| Pure length harmonic |h|>=2 | h^2>=4 | May be resonant at rational filling | Same norm bound, but irrelevant in the stated Gaussian neighborhood |
| Density-current and current-current interactions | 2 | Zero momentum, allowed | Length density coefficient obeys (4); complete coefficient sum obeys (6) |
| Other primaries and nonlinear descendants | >=3 | No further number-theoretic exclusion required | Collective norm, not a separate unweighted bound per vertex |

[proved within the free dictionary] The repeated-branch umklapp is not a dimension-two local four-fermion operator. For one chiral branch psi(x)psi(x+a)=a psi(x)partial_x psi(x)+O(a^2). Antisymmetrization requires one derivative in each repeated incoming/outgoing pair. Four fermion dimensions contribute two and the two derivatives contribute two more. Equivalently z=(2,-2) gives d0=4. For n=1 the tiny repeated-pair quartic margin therefore protects an operator that is already irrelevant near the free point. Treating its tiny margin as an unavoidable marginal instability would be a false alarm.

[conditional, interacting dimensions] In diagonal nonchiral notation, pair transfer has dimension K_a^(-1)+K_b^(-1), interflavor density backscattering has dimension K_a+K_b, one-flavor umklapp has dimension 4K_a, and a length harmonic times that flavor's backscattering has dimension K_len+K_a. General current mixing requires the full matrix, not these diagonal formulas. C45's relative stiffness bound gives

    sqrt((1-eta_G)/(1+eta_G)) d0 <= d <=
              sqrt((1+eta_G)/(1-eta_G)) d0,
    eta_G <= 1/4.                                           (8)

Thus a reference-dimension-two quartic may have dimension below two. The nonresonance condition is doing real work. Every remaining reference-dimension-three primary instead has

    d >= 3 sqrt(3/5) = 2.3237900077...,

and d0>=4 implies d>=3.0983866769.... The threshold for protecting all d0>=3 vertices is eta_G<5/13, obtained by solving 3 sqrt((1-eta_G)/(1+eta_G))>2. The reserved 1/4 is strictly inside it.

[proved, algebra control] checks.py exhaustively generates all vectors of squared norm at most six and zero sum at n=1,2,3,4, and h=0,+/-1. The numbers of nonzero primaries with d0=1 are 4,14,32,58; with d0=2 they are 4,30,150,532. In general, with N=2n, these counts are N(N-1)+2 and 6 binomial(N,4)+2N(N-1). The enumeration is exhaustive at its stated cutoff because sum |z_i|<=sum z_i^2<=6 and sum z_i=0 imply at most three positive and three negative units, which its multiset construction generates. The analytic classification is valid at every n; the enumeration only checks its implementation. The script does not enumerate arbitrary microscopic operators or establish (7).

## 3. Near resonances, slow velocities, and the order of choices

[proved, imported finite lower bounds] Let delta_pay be C42's certified margin. Its physical-symbol speed bound is

    v_cert = (t_*/2)^n delta_pay^(2n-1)/(B+1)^(n-1)>0.        (9)

C44 supplies the conservative length margin

    delta_mix >= pi/[4(5(4n+1)^4+1)].                        (10)

Use positive root patch data and v_*=min(v_cert,1), and reserve a positive delta below both dressed margins. C45's scale is

    e_det=min(e_lin,v_* delta/32),
    ell_*=log(e_top/e_det),  E_det=W e_det.                  (11)

In particular a momentum mismatch alone is not an energy threshold. Multiplication by the smallest speed is necessary. The empirical velocity behavior cited in the task is not substituted for (9); it is neither an all-n lower bound nor the physical unit W.

[conditional, relevant versus marginal growth] A nonzero momentum vertex of dimension d has dimensionless strength y(l)=y(0)exp((2-d)l) at linear order, until its oscillation is resolved. For a dimension-one bilinear, root stability requires its energy coefficient much smaller than W v_* delta. The equivalent free mass test is |m|/(W v_*) << delta. C45's root inequalities impose this directly. A marginal quartic instead has the conservative flow y(l)<=y(0)/(1-A_beta y(0)l) up to ell_*. Hence A_beta y(0)ell_*<1/2 prevents an exponentially small BKT scale from being reached before detuning. If d drifts slightly below two, its linear growth must also be included in the running estimate. In a weak Gaussian neighborhood 2-d=O(E/v_*) and its accumulated growth is bounded by exp[O(E ell_*/v_*)]; this is another finite factor covered by choosing E ell_*/v_* small. It is not justified by treating every quartic as exactly marginal at its dressed value.

[proved conditional on the RG estimates] C45 chooses epsilon_n below constants times v_*delta and v_*/(1+ell_*), with all flavor sums in its matching constants. The static bounds in (1)-(4) can be placed below any such strictly positive tolerance. Choose D, then g, then the vacancy ratio, then s/delta_fast. The actual charge gap can be as small as U rho^(2n+4); the explicit lower bound is

    delta_ch,p >= w_vac [rho/(k+p-1)]^(2kp)
                       /[binomial(k+p,p) kp], k=n+2.       (12)

It is positive after fixing rho, and is used before choosing s and lambda. The finite core gap and filter contribution enter C44's complete Delta. Using U in place of that Delta would invalidate this argument.

[proved, schedule logic] Once those parameters and beta_n are frozen, v0,f0,C_SW are finite. Taking

    lambda <= epsilon_n c_n/(16 A_dec,n C_SW)

makes (5) small. Extra restrictions ensure lambda v0<<Delta and W<<Delta. There is no lower restriction on lambda competing with these upper restrictions. Even arbitrarily fast growth of constants with n creates no contradiction for the quantifiers

    for every finite n, choose fixed nonzero parameters,
    then L -> infinity and T -> zero.

It does preclude a uniform finite-energy or finite-size test as n grows. Taking n->infinity first is not an argument for or against the proposed fixed-n phases. Uniformity of matching constants in the required auxiliary parameters remains a hypothesis; finite C42/C44 estimates alone do not supply it.

## 4. Higher-order marker and readout processes

[proved, imported selection rule] A single port insertion creates two bus excitations. Marker motion preserves bus versus vacuum token type. The resolvent also preserves it. Therefore every word with exactly one port insertion, including lambda u_mov, has zero code compression. This is stronger than power counting. Pure marker words preserve the occupied payload data and are scalar on it.

[proved, finite-support bounds] Let v_c bound the unit port on the support, let the unit marker operator have norm at most one, and put f_R=c_n mu S_q/q^2 for the physical readout coefficient divided by lambda^2. C44's resolvent estimates give the following upper bounds, divided by W:

| Process | Relative upper bound |
| :-- | :-- |
| Two marker moves | lambda^2 c_n/(beta_n^2 Delta) |
| Two ports and one marker move | 4 lambda^2 v_c^2/(beta_n Delta^2) |
| Marker and readout | 2 lambda^2 |f_R|/(beta_n Delta) |
| Two readouts | lambda^2 f_R^2/(c_n Delta) |
| Two ports and readout | 4 lambda^2 v_c^2 |f_R|/(c_n Delta^2) |
| Four ports | C4 lambda^2 v_c^4/(c_n Delta^3) on a fixed cluster |

The cluster constant C4 is not an all-volume convergence constant. These estimates show why none of these fixed orders has an unavoidable order-one ratio to W once lambda is chosen last. They do not establish convergence of the sum of all orders. Their relevant, marginal and irrelevant components have the dimensions in Section 2. An order in lambda is not a scaling dimension.

[proved at second order] C46 bounds the source norm by B_src=3(1+eta)V0(1+M_h^4). For |omega|<=Delta/2 a resolvent identity bounds the difference of static and frequency-dependent second-order kernels by 2lambda^2 B_src |omega|/Delta^2. At |omega|<=W e_top, its ratio to W is at most

    2lambda^2 B_src e_top/Delta^2.                          (13)

A frequency derivative raises the continuum dimension by one. A derivative bilinear can renormalize the kinetic term, while derivative quartics are irrelevant near the reference. Higher derivatives cost further powers of E/Delta. No claim about arbitrary-degree retarded kernels follows from this second-order check alone.

## 5. An explicit length-fluid failure outside the scheduled window

[proved] Consider just the scalar length model at fixed long-cell number r,

    H_len=-t sum_i (a_i^dagger a_(i+1)+h.c.)
                  + V sum_i l_i l_(i+1), t>0.              (14)

For V=-A with A>=2t, subtract V r. The two-site summand is zero on 00,11 and has matrix

    [ A/2  -t ]
    [ -t   A/2 ]

on 01,10. Its eigenvalues are A/2 +/- t, both nonnegative. Thus E_0(M,r)>=V r. A basis state containing a single droplet of r occupied sites, 0<r<M, has energy V(r-1) and zero hopping expectation. Consequently

    V r <= E_0(M,r) <= V r+|V|,
    lim_(M->infinity,r/M->nu) E_0/M = V nu.                 (15)

The energy is affine in density and has zero inverse compressibility, contrary to the strictly convex length-fluid energy required in C45. For A>2t the grand-canonical problem away from coexistence is empty or full. This is the phase-separated regime, not the selected interior c=1 length liquid. The strict case is sufficient here; the threshold A=2t has additional degeneracy.

[numerical, imported C44 diagnostics] C44's table at rho=.16 and s/min(delta_ch)=.003 reports scalar E_LR=-.00069091 at n=1 and -.00034838 at n=2. It reports bare traceless bond norms 1.97e-8 and 4.30e-10. If one sets the scaled length hopping to those respective norms, the scalar truncation has V/(2t)=-17535.8 and -405093.0. It satisfies (15) by many orders of magnitude. These are rounded diagnostic numbers, not certified microscopic phase boundaries. The full unequal-bond payload couplings have not been bounded relatively at these diagnostic points, so (15) is not a theorem about their full microscopic Hamiltonians.

[proved, scheduled exclusion of this failure] Under C45 equalization, epsilon_eq<=J0 epsilon_n/(16C_b), so with t=W,

    |V_len|/(2t) <= 2epsilon_eq/J0 <= epsilon_n/(8C_b).     (16)

This is far below one. The weak scalar interaction contains a smooth density product of dimension two and at commensurate half filling a cos(4phi_len) term of dimension four at the free reference. Equation (8) keeps the latter irrelevant. The diagnostic failure therefore establishes why equalization must precede marker-hopping normalization. It does not refute the final schedule.

## 6. Instantons, defects, and moving zero modes

[conditional, length instantons] A local phase slip is represented in the compact density field by h in (7). Exact marker number forbids length dual vertices but does not forbid these slips. At h=1 the momentum is 2pi nu, or 2pi nu plus a payload momentum if dressed. These are the excluded cases in Section 2. At rational nu=a/q0 in lowest terms, the first pure translation-invariant harmonic has |h|=q0 and dimension q0^2 K_len in diagonal notation. In the interior q0>=2, so (8) gives d>=4sqrt(3/5)>2. At irrational nu high harmonics can be arbitrarily close to commensurate, but their dimensions grow quadratically. They do not require an all-order Diophantine bound once the collective interaction norm and (8) hold. No independent microscopic instanton fugacity is supplied by C42/C44. Bounding it by (6) presupposes the local matched action. A defect event not represented in that action remains part of the missing matching problem.

[proved, reference defects] C44 gives wrong-length energy at least gamma_len(1+|ell-b|), gamma_len=t_sel eta_b-U>0, and positive illegal-sector costs. C45's block-diagonal reference has D>=E_P+delta N_Q with delta>=Delta/2 after its choices. The reference thermal contribution is exponentially small as established in C46. Its heat trace is finite at sufficiently low temperature because the wrong-length energy grows linearly while the number of length-ell states is at most exponential in ell. This verifies the available mechanism for the heat-trace hypothesis, not every all-word bounded-interaction hypothesis of the newer matching theorem.

[proved, weak consequence for the perturbed chain] Work in the unrestricted rooted ensemble where the imported threshold and its code trial state apply. From ||R||<=M r_SW,

    E_0(H)<=E_P+M r_SW,
    <D>_H>=E_P+delta <N_Q>_H,
    <R>_H>=-M r_SW.

Therefore

    <N_Q>_H/M <= 2r_SW/delta <=4C_SW lambda^3/Delta.        (17)

This excludes an order-one density of bare defects in that rooted regime. Transferring the statement to a fixed physical length also requires the tangent and rounding comparison; it is not automatic from this inequality. It neither proves zero dressed defect density nor a positive dressed excitation threshold. Virtual dressing itself can give a small nonzero bare density. A dilute phase competing at a nonperturbative scale is not excluded just by (17). C49 is pursuing sharper leakage/spectral estimates; this audit does not treat (17) as their replacement.

[proved, why a defect has no assigned Gaussian dimension here] An excluded-cell excitation has a positive reference energy, not a massless Gaussian primary. C42/C44 supply energy and activity bounds, not a scaling dimension for its creation field after dressing. Assigning it dimension infinity or declaring its exact fugacity zero would hide the unproved step. Only its local neutral descendants in a successfully matched action have the dimensions in Section 2. The new frozen-defect pressure theorem concerns a different frozen-sector comparison and cannot be silently substituted for the full mobile one.

[proved, exact marker selection rule] Marker number M commutes with the physical Hamiltonian. Different M sectors can cross but cannot hybridize. A term that creates an avoided crossing between them is identically zero at every order. Payload flavor charges are not all exact; their non-Gaussian transfer operators are instead governed by the momentum ledger.

[conditional, finite-size zero modes] A Gaussian charge sector has energy A(m-rho L)^2/L plus oscillator and other charge terms, with A>0. For the two adjacent integers around rho L, the splitting is A|1-2{rho L}|/L. For irrational rho its coefficient has subsequences tending to zero. Thus an all-ring lower bound on every positive gap by a fixed multiple of 1/L is false even for a Gaussian reference. It is not required to count its oscillators or its central charge. In the mobile chain the joint quotient

    (x,l,z) ~ (x+b+l_0,S l,S z)

imposes twists and correlated zero-mode quantum numbers. C44's formula includes these. The simple parabola illustrates the obstruction; it is not asserted to be the exact shell energy of every multiband wrapper ring. Crossing sectors do not provide a local relevant operator or remove a boson.

[conditional, winding tails] A quasi-local interaction with support-weight exponent a and norm r assigns total norm at most M r exp(-aM) to terms whose minimal connected support wraps the M-cell ring. Relative to the slow oscillator spacing W v_*/M this is at most

    (r/(W v_*)) M^2 exp(-aM) -> 0.                         (18)

This excludes an independent finite-ring winding obstruction once the assumed effective quasi-locality is available. A local bond crossing the chosen root seam is not such a term: it is a short cyclic interaction with the quotient identification retained. C44 correctly does not claim exact evenness in lambda on odd rings. Fixed-order winding words and (18) are consistent. Equation (18) is not a derivation of exact all-order spatial decay after massive integration.

## 7. Small demonstrators

[proved, imported C47 algebra] The ten-letter frozen demonstrator has

    h2=constant -(9/32)(sigma+ sigma-+h.c.)+(243/11648)ZZ.

The staggered rotation gives J=(9/16)lambda^2 and Delta_XXZ=27/182. The density interaction divided by hopping is 27/91, independent of lambda. These numbers are checked with exact fractions in checks.py. Its cell gap is 2/5 and its residual is bounded by C_a lambda^3. Relative to its hopping, the residual bound is (32 C_a/9)|lambda| before the matching constant. C47's imported matching bound includes an additional factor four.

[assumed benchmark, not a theorem application] Adopt the zero-magnetization XXZ Gaussian parameter

    K_XXZ=pi/[2(pi-arccos(Delta_XXZ))].                     (19)

This external exact-solution formula is used only as an explicit benchmark assumption; it is not derived here and no published infrared theorem for the retarded microscopic action is invoked. At 27/182 it gives K=0.91341427195 and half-filled umklapp dimension 4K=3.65365708782. One-cell translation excludes the staggered dimer field of dimension K; U(1) excludes pairing. Uniform umklapp is allowed but irrelevant at this reference. Under the relative Gaussian metric bound eta_G<=1/4 its dimension remains at least 4K sqrt(3/5)>2. A larger possible tolerance is eta_G<(4K^2-1)/(4K^2+1)=0.5388839..., obtained by solving that inequality.

[proved, theorem-scope limitation] Taking lambda smaller does not move Delta_XXZ toward zero. Thus an unspecified weak-XX theorem cannot be applied by citing small lambda. The newer demonstrator tuning changes the fixed microscopic port ratio first and makes the reference anisotropy arbitrarily small. That addresses this quantifier issue. It does not itself prove stability for arbitrary-degree retarded corrections. No instability of the original ten-letter demonstrator is established here.

[proved, inapplicable C42 margins] This frozen demonstrator and the six-letter mobile toy do not use the two-bus C42 pencil. Assigning their perturbations C42's momentum margins would be incorrect. Their own symmetries and actual couplings must be used. In the mobile toy there are two squeezed payload sites per cell; a uniform squeezed dimerization is at zero cell momentum. The local mobile_cells analysis excludes it using reflection, not C42 nonresonance. Its reflection-broken control D has reported c near one rather than two. In the free dimerized XX reference, alternating hoppings t1,t2 have bands +/-|t1+t2 exp(ik)|, a gap 2|t1-t2|, and a mass operator of dimension one. This is a concrete mechanism in that control, not an allowed perturbation of the reflection-preserving toy or of the scheduled C39 operator lattice. The finite-coupling toy numbers are imported numerical evidence, not new ED here.

[proved, scope of the n=2 selector issue] The separate N2_SELECTOR_AUDIT finds that sums of shared-port run lower bounds were treated as saturable selected-sector energies. Their difference does not establish a dressed per-defect margin. It also exhibits exact zero-coupling zero modes omitted in saved spectra. These are existing failures of a demonstrator's sector-selection evidence, not a computed relevant perturbation or a refutation of C42. This audit does not recalculate those settled findings or use the suspect data to infer a phase.

## 8. Limits and what weak matching cannot establish

[proved, order of limits] At fixed finite n and fixed nonzero couplings, E_det in (11) is positive but may be extremely small. The required infrared regime is L much larger than the fixed crossover length and T much smaller than E_det, with the thermodynamic limit taken first for pressure. Taking lambda->0 at fixed L or fixed T instead can show an effectively degenerate code with unresolved modes. These two experiments do not commute with the claimed infrared limit. Their noncommutation is not an inconsistency in a fixed-coupling statement.

[proved, weak pressure bound fails] C46's direct comparison yields, schematically,

    |f_H(T)-f_P(T)| <= r_SW+C T exp(-delta/T).               (20)

Even after subtracting the two ground energy densities the bound is at best 2r_SW plus the thermal term. Division by T^2 diverges at fixed lambda. A prethermal remainder r_pre>0 has the same defect. A global eigenvalue error M r_pre becomes larger than W v_*/M when M^2 r_pre exceeds W v_*. Neither estimate establishes the coefficient of T^2 or the asymptotic 1/L spectrum.

[proved, explicit countermodel to the weak inference] Add a staggered potential m(-1)^j n_j to a free nearest-neighbor fermion chain. The local perturbation norm is |m|. Its two-band spectrum is

    E_+/- (k)=+/-sqrt(4t^2 cos^2 k+m^2).

Any nonzero m opens a single-particle band gap 2|m|. The continuum perturbation has dimension one. A small fixed energy-density error therefore does not preserve a nonzero thermal central-charge coefficient. This countermodel breaks the one-cell translation/nonresonance conditions of C39 and is explicitly not a counterexample to full S. It proves that replacing S by only a small norm or a prethermal estimate is insufficient.

[proved, source-independent scalar information] Bounds on normalized correlators alone do not determine a vacuum free energy. Algebraically multiplication of a generating functional by exp(a M/beta) changes f by -a T^2 without changing normalized source derivatives. This is an information test, not a claim that every such factor comes from a fixed local Hamiltonian. To infer c=n+1 one needs the full Gaussian current/stress limit or scalar thermal estimates with o(T^2) error and identified velocities. C46 and C51 already isolate this obligation; no result in this audit removes it.

[proved, dimensional normalization of moment bounds] Time moments must be stated with their units. For an auxiliary two-level system with gap Delta and perturbation r O sigma_x, with a logical occupation projector O, the second-order external-field excursion has scalar kernel r^2 exp(-Delta t). Its second unscaled moment is 2r^2/Delta^3. Choosing W=Delta^2 and r=epsilon W gives W/Delta->0 and r/Delta->0 as Delta->0, but the second moment divided by r is 2epsilon/Delta, unbounded. Thus a literal bound on unscaled second time moments by A r with A independent of Delta is false for that elementary kernel. Gap-weighted moments, with Delta^m t^m, avoid this problem. The displayed frequency estimate in formal S uses omega/Delta and is consistent with that normalization. This is a formulation requirement, not a physical instability or a refutation of that frequency estimate.

## 9. What the newer matching work does and does not settle here

[proved, algebra using the imported abstract theorem] In INFINITE_CELL_PROGRESS the gap is denoted g_aux, tau=2L_weight/g_aux, kappa=g_aux/8, v_match=1/(4tau), and B=100 exp(2kappa tau)/tau. Its leakage Lipschitz factor is

    4B/v_match=1600 exp(L_weight/2),                        (21)

independent of g_aux and of W. This checks that its stated zeroth-moment bound has no automatic inverse-bandwidth obstruction. Unscaled time moments still carry the powers of 1/g_aux described above. Its smallness condition on full coefficient interactions is proportional to g_aux and must be imposed after fixing the wrapper gap.

[proved, application checklist] For the selected rooted cells, C44/C45 provide a finite retained code of rank 2^(n+1), positive reference gap, exact charges, locally block-diagonal reference, and a finite-order exponentially weighted residual. Linear wrong-length costs provide the candidate excited heat trace, as in C46. The finite-cell theorem applies most directly to the frozen ten-letter sector as checked in C47. For the full mobile chain this audit has not verified the bounded coefficient decomposition on every malformed invariant sector, all length-tilt/source hypotheses, or preservation of the intended local operator lattice through the exact matching. The physical-ensemble campaign supplies interacting pinning and fugacity theorems, so the older limitation to common bonds should not be repeated as though no extension existed. Their complete microscopic application and subsequent scalar/current RG remain conditional.

[proved, final boundary] No coefficient, instanton or level crossing found here defeats the complete C45 schedule under S. There is a precise destructive scalar attraction in a weaker unscheduled length model, an exhaustive conditional low-dimension exclusion, and several precise failures of weaker matching/finite-size/thermal claims. There is no unconditional proof that every microscopic nonperturbative process is harmless. Claiming one would amount to assuming the step this task was asked to attack.
