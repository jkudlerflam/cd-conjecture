# C66. The n=1 wrapper action and its conditional infrared comparison

[proved, scope] This note makes the two-channel specialization explicit. It does not prove the global routing left open in C61 or the comparison left open in C50. Every theorem below is either finite algebra or an implication with stated analytic hypotheses. In particular, setting n=1 does not establish input S.

[proved, tag convention] Each tag covers the following prose and equations until the next tag or heading. A proved statement importing a model or a theorem is an implication from that named input, not a new certification of the input.

[assumed, imported inputs] Use C39/C42's code and quadratic reference, C44's four projected bonds and identity transport, C45's fixed-coupling schedule, and C48's massive matching with C56's connected-network refinement. These inputs are not independently reproved. The separate Jordan-Wigner construction in ../../contrib/thermal_ensemble_20261028/exact/SEPARATE_JW.md supplies two local fermion species for interactions neutral under both charges. C58 supplies the conditional thermodynamic ensemble transfer. The present argument works in ordered-cell coordinates, in the infinite-volume bulk, with a fixed infrared cutoff removed last. It does not supply a finite-ring spectral theorem.

## 1. The actual two-species projected Hamiltonian

[proved from the imported code] At n=1 a cell has length 4 or 5 and one payload qubit. Write c_j for the payload fermion, d_j for the length fermion, n_j=c_j^dagger c_j, and m_j=d_j^dagger d_j. The short and long projectors are 1-m_j and m_j. Use separate Jordan-Wigner strings and a Klein factor, as in SEPARATE_JW, rather than a single interleaved string that would incorrectly make one hopping parity dependent. Both total charges are conserved. At fixed marker number M the retained physical length is L=4M+sum_j m_j. Boundary parity and cyclic gluing affect sectors and are not erased by this bulk notation.

[proved, definition of second-order coefficients] Let K2_ab, a,b in {0,1}, be C44's source-resolvent bond on the two payload qubits, with a=0 denoting its p=1 length and a=1 its p=2 length. Define the physical energy operators

    Q_ab = lambda^2 K2_ab,
    A=Q_00,  B=Q_10-Q_00,  C=Q_01-Q_00,
    D=Q_11-Q_10-Q_01+Q_00.                              (1)

The retained second-order Hamiltonian, including its readout and scalar fields, is

    H_proj = sum_j [-j_l(d_j^dagger d_(j+1)+h.c.)
                   + A_j + m_j B_j + m_(j+1) C_j
                   + m_j m_(j+1) D_j]
             + h_p sum_j n_j + h_l sum_j m_j + M c_0.   (2)

A_j,...,D_j act on payload sites j,j+1. This is an identity, obtained by expanding the four length projectors. It is not an average over length configurations. Additional orders of the microscopic elimination are represented by the kernels in Section 4, not declared absent from the microscopic model.

[proved, complete payload bond basis] Every Hermitian, payload-number-conserving two-qubit operator O has the unique form

    O = a_O I + b_O (n_j-1/2) + c_O (n_(j+1)-1/2)
        + U_O (n_j-1/2)(n_(j+1)-1/2)
        + tau_O c_j^dagger c_(j+1)
        + conjugate(tau_O) c_(j+1)^dagger c_j.          (3)

There are four real diagonal entries and one complex off-diagonal entry in the one-particle block, hence six real coefficients. The diagonal polynomials in (3) are a basis of the four diagonal entries. This proves completeness. In particular a_O=tr(O)/4. All coefficients in (1)-(3) can be extracted from the finite C44 matrices without a fit.

[proved, scalar attraction and mixed operators] Write a_ab=tr(Q_ab)/4. The scalar length interaction is

    V_l=a_11-a_10-a_01+a_00=a_D.                        (4)

Only the scalar constant and the linear length terms can be put into c_0,h_l. The term V_l m_j m_(j+1) remains. The traceless parts of B,C,D generate length-payload density interactions, length-dependent payload hopping and higher-field interactions. For example m_j U_B(n_j-1/2)(n_(j+1)-1/2) has six-field terms, and m_j m_(j+1)U_D(n_j-1/2)(n_(j+1)-1/2) has eight-field terms. Even (2) is not solely a quadratic Hamiltonian plus six contact quartics.

[proved, schedule bound] C44/C45 give

    s K2_ab=K_bare+E_ab,   ||E_ab||<=d_eq,
    W=lambda^2 J0/s.

Thus ||Q_ab-lambda^2 K_bare/s||<=W d_eq/J0. Since |tr E/4|<=||E||,

    |V_l| <= 4W d_eq/J0.                               (5)

The traceless part of D has norm at most 8W d_eq/J0. Bounds for B,C follow with factors 4. The complete coefficient norm has the additional fixed two-cell conversion constant C_b of C45. Reducing lambda leaves the ratios V_l/W and these second-order errors/W unchanged. Equalization, including its scalar part, must be done before lambda is chosen.

## 2. Scheduled dispersions and all four velocities

[proved from the supplied reference] For n=1, p=6, q=(w u)^2, w=1/4, and

    S_q=1+q^2+q^4,  J0=2 g^2 q^4 D_aux^(-6)/S_q^2,
    W=lambda^2 J0/s.                                  (6)

D_aux is C45's high-rail energy D, distinguished from the operator D in (1). Its one-vertex stub has A_stub=4, so the supplied Gram formula gives f_0(t)=t exp(-4t) and

    Gamma(w)=Gamma(wv)=integral_0^infinity t^2 exp(-8t)dt
                        =2/8^3=1/256=:gamma_0.         (7)

No large-n velocity estimate is needed for this scalar specialization. The normalized limiting payload dispersion is

    e_p(k)=-2 gamma_0 cos k+2 eta gamma_0 sin k+mu_p
          =-2 t_p cos(k+phi)+mu_p,
    t_p=gamma_0 sqrt(1+eta^2),  phi=atan(eta).           (8)

Here t_p is in W units, and mu_p is the scheduled readout, not an RG counterterm. If |mu_p|<2t_p, put k_p=arccos(mu_p/(2t_p)). Then

    p_pR=-phi+k_p,  p_pL=-phi-k_p,
    v_pR=+v_p,  v_pL=-v_p,
    v_p=2t_p sin k_p=sqrt(4t_p^2-mu_p^2).              (9)

The uniform phase can be gauged away in the bulk. On a ring this changes the boundary twist. Formula (9), rather than a congruent symbol's slope, is the actual scheduled reference velocity.

[proved from C44/C45 transport] C44's positive overlap beta_1 is the scalar multiplying identity transport. C45 chooses the physical marker coefficient t_m=W/beta_1. Consequently

    j_l=beta_1 t_m=W,
    e_l(k)=-2 cos k+2 cos(pi nu),
    p_lR=+pi nu,  p_lL=-pi nu,
    v_lR=+v_l,  v_lL=-v_l,  v_l=2 sin(pi nu).           (10)

The four physical speeds per cell are W times (v_p,-v_p,v_l,-v_l). A prescribed interior nu is chosen by the marker scalar and physical-length stationarity, C45 (34), not by adding a new independent microscopic length field. The grand-canonical h_l in (2) is useful bookkeeping in a fixed-M description.

[proved, distinction of references] The C45 choice has j_l/(W t_p)=1/t_p. For small eta and mu_p, v_p is near 1/128 whereas v_l is of order one. Finite D_aux, finite g, equalization and matching perturb these reference values. C45 chooses their full norm small relative to the actual v_min=min(v_p,v_l), root patches and momentum margins. Equations (8)-(10) are not formulas for the fully dressed speeds.

[numerical, imported comparison] verify_C44 instead uses units of the projected payload hopping |t| and sets j_l/|t|=1. It uses half-filled payloads and nu between 0.29 and 0.38. Its M1 and M2 tests are not the velocity hierarchy j_l/(W t_p)=1/t_p of (10). Its calibrated c=2 is evidence for those projected models, with the omissions listed in its README. It does not test the all-order retarded action at the scheduled four velocities.

## 3. Six density couplings and an explicit first-order dictionary

[proved, local classification] Write c_x approximately as sum_(sigma=R,L) exp(i p_psigma x) psi_psigma and d_x similarly with p_lomega. There are four chiral fields. Exact payload and length U(1) forbid cross-species quadratic terms and pair conversion between the species. Translation excludes same-species R-L bilinears when their Fermi points differ modulo 2pi. After excluding mixed exchange momenta, a nonzero marginal contact quartic is a product of two distinct branch densities rho_a=psi_a^+ psi_a^-.

[proved, six entries] In the order (pR,pL,lR,lL), define the real symmetric matrix

    G = [ 0       g_p       g_RR     g_RL
          g_p     0         g_LR     g_LL
          g_RR    g_LR      0        g_l
          g_RL    g_LL      g_l      0     ].            (11)

Then the local quartic action is

    integral dz [g_p rho_pR rho_pL + g_l rho_lR rho_lL
          + sum_(sigma,omega) g_sigma,omega
                                  rho_psigma rho_lomega]
      = (1/2) integral dz rho^T G rho.                 (12)

The two within-species and four mixed couplings exhaust the contact terms. Same-branch Grassmann quartics vanish. Their point-split current analogs can renormalize quadratic velocities, which are retained separately. The six entries are the allowed coordinates, not a claim that all are independently nonzero in the bare bond.

[proved, needed margins] In this two-species case the mixed exchange mismatches are 2(k_p-k_l) and 2(k_p+k_l), modulo 2pi, where k_l=pi nu. Keep their distances positive on the root box. Pure same-species umklapp has transfer 4k_a and two Pauli derivatives if its momentum is resonant, hence free dimension four. A strict unrestricted pair-sum certificate may instead detune it, as in C42. Exact half filling in the numerical verifier must not be described as satisfying that unrestricted certificate. Under C45's small Gaussian metric perturbation these derivative/vertex operators remain irrelevant; their all-scale bound still needs Lemma N. One may select nu=1/3 for the near-half-filled n=1 reference to keep the mixed margins open, or use C45's certified density choice.

[proved, stronger reference hypothesis qualification] If using a theorem that tests all pair sums without using the two separate charges, also impose C46's zero-pair test and its full finite exclusion list for nu. Here p_pR+p_pL=-2phi is nonzero modulo 2pi for the scheduled eta>0, whereas p_lR+p_lL=0. Separately gauging away phi is harmless for the actual charge-preserving action, but makes those two sums coincide. Thus one cannot use that gauge and then claim an unrestricted pair-sum certificate without modifying the reference hypothesis to respect the two charges. The explicit nu=1/3 example above certifies the stated symmetry-allowed mixed margins, not every unrestricted pair-sum test.

[proved, first-order matching prescription] Here is an explicit coefficient dictionary for (2), useful even when the finite bond entries are left symbolic. Use Gaussian normal ordering at the free Fermi seas and then localize the normal-ordered quartic. Put

    nbar=k_p/pi,  kappa_p=exp(-i phi) sin(k_p)/pi,
    kappa_l=sin(k_l)/pi,  M_l=nu^2-|kappa_l|^2.

For a payload bond O of (3) define

    E_p(O)=a_O+(b_O+c_O)(nbar-1/2)
          +U_O[(nbar-1/2)^2-|kappa_p|^2]
          +2 Re(tau_O kappa_p),
    H_psigma(O)=b_O+c_O+2U_O(nbar-1/2)
          -2U_O Re(conjugate(kappa_p) exp(i p_psigma))
          +2 Re(tau_O exp(i p_psigma)),
    H_lomega=2nu-2 Re(conjugate(kappa_l) exp(i p_lomega)). (13)

E_p(O) is its free expectation. H_psigma(O) is its coefficient of the slowly varying normal-ordered branch density. H_lomega is the corresponding coefficient of m_j m_(j+1). In the phase gauge of (10), H_lR=H_lL; no symmetry setting this equality for a general renormalized four-branch action is assumed.

[proved, projected first-order values] In W units the six coefficients of the Gaussian-normal-ordered static action are

    g_p^(1) = [4 sin^2 k_p/W]
                    [U_A+nu(U_B+U_C)+M_l U_D],
    g_l^(1) = [4 sin^2 k_l/W] E_p(D),
    g_sigma,omega^(1) = [H_psigma(B+C)
                              +H_lomega H_psigma(D)]/W. (14)

The scalar contribution to g_l^(1) is 4(V_l/W)sin^2 k_l. Traceless mixed payload terms also contribute to E_p(D); scalar trace subtraction alone is not the full Hartree-Fock length interaction. The four mixed entries in (14) happen to be independent of omega at this static order. Keeping six coordinates in the RG remains the safe general formulation.

[proved, derivation of (13)-(14)] For one species, Wick's identity gives <n_j n_(j+1)>=nbar^2-|kappa|^2. Vary the occupation of a narrow momentum packet at p by dN per unit length. Then dnbar=dN and dkappa=exp(ip)dN. Differentiating the expectation in (13) gives H. The coefficient of two distinct packet densities at p_R,p_L in nbar^2-|kappa|^2 is 2-2cos(p_R-p_L)=4sin^2 k. This is also the antisymmetrized direct-minus-exchange quartic coefficient, so it fixes the normalization of (12). Expectations factor between species. Applying the two packet derivatives to E_p(A)+nu E_p(B+C)+M_l E_p(D) gives (14), including contractions from the six- and eight-field endpoints. No extra combinatorial factor is left because (12) sums each unordered pair once.

[assumed, full matching data] The exact matched g_ab are the localized four-point kernels after massive matching, ultraviolet initialization and fixed crossover, in a specified normal-ordering and field-normalization convention. Equation (14) is the linear Gaussian-normal-ordering dictionary for the projected interaction; it is not their all-order value. If the complete perturbation has an analytic small-norm finite-step expansion, nonlinear finite-step corrections are O(norm^2), while additional supplied retarded endpoints can contribute at first order in their own norm. In particular a six-field kernel can contract to a quartic at first order, as C50 Section 6 proves. No numerical values of the six fully matched g_ab are supplied by the cited numerical c calibration.

## 4. The matched retarded action

[assumed, matching interface] Normalize time by t'=Wt. With lattice Fourier measure dk0 dk/(2pi)^2, retain both complete lattice bands during ultraviolet initialization. The free covariance is diagonal in species,

    C_a(k0,k)=1/[-i k0+e_a(k)],  a=p,l.                (15)

Use an equal-time prescription consistent with the Hamiltonian normal ordering, retaining induced quadratic counterterms. After one fixed crossover below e_det of C45 (24), resolve the four real smooth patches about trial roots p'_a. Write e_a(p'_a+q)=e_a(p'_a)+v_a q+r_a(q), r_a(q)=O(q^2). The notation a now denotes a branch. A concrete matched action is

    S_match = S_scalar + sum_a integral_k psi_a^+(k)
          [Z_a(-i k0)+Z_a u_a q+r_a(q)+E_0 nu_a] psi_a^-(k)
          + (1/2) integral_z rho(z)^T G rho(z)
          + R_2 + R_4 + sum_(m>=3) V_(2m) + S_sources, (16)

with the patch cutoffs implemented in the Gaussian measure. At zero interaction Z_a=1 and u_a=v_a from (9)-(10). Local frequency and spatial terms can equivalently be absorbed into a running covariance. Scalars are kept as pressure densities. E_0>0 is the fixed crossover energy in W units. Relevant nu_a are four bookkeeping counterterms, not four new microscopic controls.

[proved, what 'retarded' includes] For each even degree p,

    V_p=(1/p!) sum_(branch,charge labels) integral K_p(z_1,...,z_p)
                                       psi_1 ... psi_p. (17)

Integrals include lattice sums and relative-time finite measures. Kernels can have instantaneous delta components, nonzero relative times and all even degrees. Each kernel has equal creation and annihilation count within each species. R_2 is the matched quadratic kernel after its constant and first derivatives are localized; R_4 is the quartic kernel after its allowed zero-momentum density values are localized. Nonresonant terms are retained with their Fermi phases until the detuning scale permits their irrelevant bound. They are not set to zero as functions on the lattice.

[proved, normalization] For a degree-p coefficient density the t'=Wt substitution multiplies the kernel by W^(-p); integration over p-1 anchored relative times leaves a factor W^(-1). For measures use this pushforward identity. A physical time-tree rate k_T becomes k_T/W. Thus W<=Delta preserves a positive normalized temporal rate supplied by C56. The free curvature r_a is present at zero interaction and must not be charged an interaction prefactor.

[assumed, thermal version of the matching interface] Equations (16)-(18) describe the zero-temperature endpoints. At finite inverse temperature beta_phys, use beta'=W beta_phys, periodize the zero-temperature kernels and retain the thermal auxiliary corrections, including their scalar and source components. C48 bounds a nonzero thermal baseline by a fixed prefactor times exp[-beta_phys Delta/(4T_1)] in its original norm. The corresponding anchored normalization divides the energy norm by W. At fixed couplings this remains activated as beta_phys tends to infinity; it is not identically zero when leakage is zero. No source-independent thermal term is removed by writing a normalized fermion correlation. The later scalar infrared comparison is still a separate problem.

[assumed, sufficient endpoint norm] Retain C56's connected network, not only an exponential diameter. With fixed field radius H, one sufficient prepatch norm uses exp(T), where T is the terminal-tree length in a fixed positive space-time metric. A strong C61 endpoint hierarchy is

    H^p sup_anchor sum_labels integral exp(T) d|K_p^irr|
                   <= A epsilon^max(2,p/2-1), p>=2.    (18)

This is a sufficient input with a perturbative-degree reparameterization, not an inference from 'n=1'. Alternatively carry the actual forcing size separately. C56's refinement provides the terminal-tree geometry conditional on C48's polymer construction, with a_T=alpha/8 and k_T=L_0/(128 tau). Increase the microscopic support/radius budget before choosing lambda to pay fixed coefficient and moment losses. Smooth compact patch conversion uses polynomial tree norms or explicitly retained conversion legs. It is not asserted to preserve exponential spatial decay.

[assumed, interface audit] The code rank here is four, the two exact charges are those above, identity transport is C44's, and the positive auxiliary threshold is C45's minimum over the actual wrapper sectors, including the charge doublet. These are the n=1 substitutions in the imported matching theorem. Auxiliary heat-trace and full-space completion/ensemble hypotheses are inherited from the C48/C58 and thermal-ensemble application, not newly checked by diagonalizing a truncated cell. UV initialization uses the matrix-band extension in SINGLE_CHANNEL_RG Section 5.2 with both bands retained. Its finite-step result does not supply an infinite RG iteration.

## 5. Velocity dependence of localization and mixed trees

[proved, free line estimates] For a smooth annular cutoff and E_h=E_0 2^h, a branch with signed speed u has the bound

    |partial_t^a partial_x^b C_(u,h)(t,x)|
      <= C E_h^(1+a+b)/|u|^(1+b)
                   (1+E_h|t|+E_h|x|/|u|)^(-M),
    ||C_(u,h)||_1 <= C/E_h,
    Gram determinant of r lines <= (C E_h/v_min)^r.   (19)

Rescale k0=E_h q0 and k=E_h q1/u. The measure contributes E_h^2/|u|, the denominator E_h^(-1), and each spatial derivative E_h/u. The smooth annular multiplier is Schwartz. Integrating its absolute value gives the first two claims. Factoring its Fourier multiplier into two square roots gives Hilbert norms squared at most C E_h/|u|. A direct sum over the four branch Hilbert spaces and the Gram inequality give the determinant bound, also for positive BBF interpolation matrices. This proves (19).

[proved, mixed-channel interpretation] In a mixed graph each covariance still carries one fixed branch label. Different branches interact at endpoints, not through a free interspecies propagator. Choose one common localization metric and bound (19) using v_min and v_max. Its constants may be large when v_p is small, but no inverse |v_p-v_l| is introduced. Coalescing positive speeds is regular; a speed tending to zero is not. Do not partial-fraction a product of different branch propagators and introduce a spurious inverse velocity difference.

[proved, finite-channel combinatorics] If kernels use the collective label norm in (18), summing branch labels is already paid. Otherwise four choices per field produce at most 4^p, absorbed by one fixed field-radius reserve. The number of branch labels changes constants, not remainder orders or scale dimensions. It does not bound external degree, since arbitrarily many fields at distinct times can belong even to one species.

[proved, conditional lifting of a routing theorem] More explicitly, suppose a one-species routing proof uses only kernel moments, the uniform line/Gram bounds (19), exact Taylor identities and charge-compatible local remainders, rather than an identity special to a scalar interaction. Fix every field label in a four-branch tree. The same spatial proof then applies with the worst fixed speed constants. If the endpoint degrees are p_i, every internal propagator identifies two endpoint-derived field labels. There are at most sum_i p_i field labels to choose, so even ignoring those identifications the count is at most 4^(sum_i p_i). This is a fixed exponential endpoint loss, with no factor depending on the number of unary scale steps. A sum over the retained Gram determinant is handled by the direct-sum representation in (19). Thus such a routing theorem would lift to this action. Reference signed cancellations still require their multichannel theorem and are not consequences of this absolute-value argument. This implication identifies what a future one-channel routing proof would need to contain.

[assumed, normalization] On the running history domain keep |u_a|>=v_min/2 and the Z_a single-step ratios close to one. The anomalous products of Z ratios must remain in every marked tree. Dropping them would not prove the needed history derivative bound. Constants below depend on this fixed domain, actual speeds, patches and detuning margins, not on the infrared depth. No uniform v_min->0 statement is made.

## 6. Corrected Lemma N for these four branches

[assumed, N-action_2] For the initial action (16)-(18), the normalized source-free renormalized tree expansion, with one and two history marks, converges uniformly in infrared depth and on a compact box of trial Fermi points. A single field-radius reserve pays all degree sums. In every term the spatial routing, derivatives, endpoint moments and BBF anchor integrations use one and the same representation. Relevant/marginal localization has

    z(2)=2, z(4)=1, z(p)=0 for p>=6,
    d_R(p)=p/2-2+z(p)>=p/6.                            (20)

All engineering derivatives are counted once. The four-branch version of C61 requires the same simultaneous routing as its one-channel statement, now with a branch label attached to each line.

[proved, unchanged cut criterion] In a reduced scale tree let ell_v=h_v-h_parent(v)>0. Pair each coordinate moment at scale t to a localization derivative at an ancestor scale s<=t. Let F_v count paired paths crossing edge v. Then

    sum_derivatives a_f s_f - sum_moments b_l t_l
                         =-sum_v F_v ell_v.            (21)

For one pair its contribution is s-t, the negative sum of the traversed edge lengths; summing proves (21). If F_v>=z_v and the total routing coefficient mass is at most C^(sum_i p_i+N_end), it follows that the scale factor is bounded by C^(sum_i p_i+N_end) product_v 2^(-z_v ell_v). Four branch labels do not enter this telescoping identity. Consumed paths must be concatenated, since (u-t)+(s-u)=s-t, not erased.

[proved, local algebra retained] With R_m=1-T_(m-1), P_a a homogeneous external derivative polynomial of degree a, and |beta|=b, the exact normal form is

    partial^beta R_m(P_a F)
      =sum_(delta<=beta) binom(beta,delta) (partial^delta P_a)
         R_max(0,m-b-a+|delta|)(partial^(beta-delta)F).   (22)

Differentiating a Taylor polynomial lowers its order by b. Multiplying by a homogeneous degree-c polynomial lowers the residual Taylor order by c. Leibniz gives (22). For d=|delta| in a nonzero term, r=max(0,m-b-a+d),

    J=a+b-d+r=max(m,a+b-d)>=m,
    B=b-d+r<=2 when b,m<=2.                            (23)

If r=0, B=b-d<=2. If r>0, B=m-a<=2. The nonzero term has a>=d. This proves the counts for either species and for mixed quartics. It does not construct their paths in the BBF spatial tree.

[proved, weighted-line rule] Once exact identities justify a factor X^b C_(u,h)(X), keep it intact under later derivatives and evaluation at shifted positions. Every Leibniz term X^(b-d) partial^(a-d) C_(u,h) has integrated scale E_h^(a-b-1), with fixed speed-dependent constants from (19). Shifting only the propagator while retaining an original displacement gives x-y=X+t(x-p), whose internal offset is not paid by a moment of X. This is C61's remaining geometric issue in both species, not a new two-channel cancellation.

[assumed, N-comparison_2] In the same absorption and cutoff scheme, trees with an initial irrelevant endpoint or a lattice/reference covariance difference have short memory. Constant-history reference cancellation is a signed-sum input. History marks preserve the additional interaction factor after Z normalization. Include relevant-history control and continuity as specified in Section 7. These assumptions are separate from (20)-(23). Neither positivity of the Gaussian stiffness nor the finite number of six couplings proves them.

[assumed, N-source_2] For each specified physical or projected source O, let C_O(u;z) and C_O^ref(u;z) include its matched amplitudes and oscillations. Define

    D_O(u;z)=C_O(u;z)-C_O^ref(u;z),
    F_O^free(z)=D_O(0;z),
    D_O^int(u;z)=D_O(u;z)-D_O(0;z).                    (24)

The Gaussian subtraction uses the same roots, speeds and bare source as the action at that point of the parameter box. Require a nonzero matched source overlap and summable Gaussian-subtracted normalization differences. Separately require, for an identified scaling component,

    D_O^int=F_O^cross,int+sum_(h<=0) R_O,h^int,
    |R_O,h^int(z)|<=C ||u|| E_h^(2X_O(u)+alpha)
                       /[1+(E_h ||z||_0)^M],          (25)

with alpha>0 and M>2X_O+alpha uniformly in the domain. The fixed norm ||z||_0 may use any one reference speed since v_min>0 is fixed. Different propagation rays enter constants and scaling functions. Multiple leading source components must be resolved separately. Bare source moments do not have an interaction smallness factor.

[proved, free subtraction check] A source containing only payload fields in a decoupled two-species theory has the same normalized correlation as in the one-species theory: the length partition function cancels. C60's free pair example therefore survives unchanged when a spectator length channel is added. At k_p=pi/3, x=3m it has

    F_pair^free(x)=3/[4pi^2 x^2(x^2-1)].                (26)

For completeness, the equal-time free kernel is s(x)=sin(k_p x)/(pi x). Wick contraction of the adjacent-site pair gives s(x)^2-s(x-1)s(x+1). At x=3m this is 3/[4pi^2(x^2-1)]. Subtract its continuum leading term 3/(4pi^2 x^2) to obtain (26). This is a consistency test of the general two-channel lemma, not the prescribed near-half-filled parameter of (8). A mixed source can similarly contain free lattice descendants from either factor. Only D_O^int vanishes at u=0. Splitting the shell sum in (25) at E_h||z||_0=1 yields O(||u|| ||z||_0^(-2X_O-alpha)) by two geometric sums. The separate free term must still be compared with the proposed leading decay. Sources for length response also require the source-uniform estimates used in C58.

[assumed, unresolved part] No simultaneous all-tree routing proving the cut coverage, bounded derivative/moment loads and coefficient mass has been established here. C61 does not provide it even for one channel. The identities (21)-(23) neither prove nor refute it for the two-channel action.

## 7. C50 (7.1), with all six couplings and relevant terms

[proved, explicit norm choice] Let g_h contain the four normalized velocity corrections and the six real entries in (11), but not Z_h itself. One fixed collective norm is

    ||g||_v = sum_a |delta u_a|/|v_a|
                 +sum_(a<b) |g_ab|/[2pi sqrt(|v_a v_b|)]. (27)

This pays the slow payload velocity and the mixed geometric mean. Define B_h=||g_h-g_(h+1)||_v and suppose ||g_h||_v<=eta on the proposed domain. Denote ultraviolet/comparison forcing by d_for. It need not be epsilon^2 from (18), because the reference's own ultraviolet kernels can be O(eta). Curvature alone cannot generate an interaction beta function at zero interaction, but need not carry an interaction factor as a covariance remainder.

[assumed, actual comparison estimates] With 0<theta<kappa, assume the normalized lattice and reference beta vectors satisfy, on the same history domain,

    ||beta_h^lat(G,a;K)-beta_h^ref(G)||_v
                <= A_1 d_for 2^(theta h)
                    +b eta sum_(k>=h) 2^(-kappa(k-h)) |a_k|,
    ||beta_h^ref(g,...,g)||_v <= A_2 d_for 2^(theta h),
    ||partial_(g_k) beta_h^ref|| <= b eta 2^(-kappa(k-h)). (28)

Here a_h is the four-component dimensionless relevant history. A history derivative includes all induced Z-ratio derivatives. Equations (28), not a Gaussian diagonalization, are the missing replacement-action estimates. The multichannel constant-history theorem exists already: Mastropietro-Porta Theorem 8.5 gives decaying velocity and quartic beta functions on constant histories. Its scheme must be matched to (16). Porta-Scola-Singh Assumption 2.1 allows unequal nonzero speeds, and Remarks 4.12-4.13 describe short memory and covariance continuity for its tree expansion. Neither cited statement by itself establishes (28) for the supplied arbitrary retarded endpoints. Primary sources were checked at https://arxiv.org/html/2108.13708v3 and https://arxiv.org/html/2509.08665v1.

[proved, relevant-history closure conditional on its bound] Use the convention a_(h-1)=2a_h+b_h^a. A bounded critical trajectory has

    a_h=-sum_(k<=h) 2^(k-h-1) b_k^a.                  (29)

This follows by solving the recurrence from a lower cutoff and taking its bounded limit. Suppose

    |b_h^a|<=A_a d_for 2^(theta h)
                   +b_a eta sum_(k>=h) 2^(-kappa(k-h))|a_k|.

In the norm ||a||_theta=sup_(h<=0)2^(-theta h)|a_h| set

    C_rel=1/[2(1-2^(-(1+theta)))],
    C_mem=1/[1-2^(-(kappa-theta))].

Summing (29) and the preceding inequality gives

    ||a||_theta <= C_rel A_a d_for
                        +C_rel b_a eta C_mem ||a||_theta.

Thus if C_rel b_a eta C_mem<=1/2, every such critical trajectory obeys ||a||_theta<=2C_rel A_a d_for=:K_a d_for. Existence, uniqueness and history-difference continuity are still analytic hypotheses; C56's conditional sequence construction can supply them once its G1-G2 estimates hold.

[proved, reduction to C50] Freeze the reference history at g_(h+1), interpolate linearly to G, and integrate its derivative in (28). Telescoping gives

    ||g_k-g_(h+1)||_v <= sum_(j=h+1)^k B_j.

Insert the relevant bound from (29). With A_eff=A_1+A_2+b eta C_mem K_a one obtains

    B_h <= A_eff d_for 2^(theta h)
       +b eta sum_(k=h+1)^0 2^(-kappa(k-h))
                                    sum_(j=h+1)^k B_j. (30)

This is C50 (7.1) for the ten marginal coordinates of the two-channel action, allowing different forcing and history sizes. Setting d_for=eta=epsilon reproduces C50's convention. The four mixed density coordinates appear in both norms and all marks; they cannot be dropped because they are small at one projected point.

[proved, summable closure] The double geometric constant is

    C_(kappa,theta)=2^(-(kappa-theta)) /
              [(1-2^(-kappa))(1-2^(-(kappa-theta)))].    (31)

If b eta C_(kappa,theta)<=1/2, downward induction in (30) gives

    B_h<=2A_eff d_for 2^(theta h),
    ||g_h-g_*||_v <= 2A_eff d_for 2^(theta h)
                                      /(1-2^(-theta)). (32)

Indeed substitute the induction hypothesis with k-h=d and j-h=r. The convolution is 2A_eff b eta d_for 2^(theta h) times sum_(d>=1)2^(-kappa d)sum_(r=1)^d2^(theta r). Reversing the sums gives (31), so the convolution is at most A_eff d_for 2^(theta h). Summing the increments proves (32). Require the total drift 2A_eff d_for/(1-2^(-theta)) to fit inside the reserved history, speed and stiffness domains. This closes the conditional bootstrap; it does not construct the missing beta functions.

[proved conditional on C56, physical roots] For trial roots in a compact product of four patches, use the exact initialization a_0^phys(p')=A_init(p')[e(p')-mu]+u_init(p'). C56's continuous critical map gives the added-and-subtracted counterterm

    nu_phys(p')=A_init(p')^(-1)[u_init(p')-a_0(p')].

If its sup norm is <=C_nu d_for and patch slopes are at least v_min/2, the inverse dispersions map the root box to itself when 2C_nu d_for/v_min is below its radius. Brouwer gives roots with displacement <=2C_nu d_for/v_min. A quartic mismatch changes by at most four times this displacement; keep that below half the original margin. The four-dimensional map solves bookkeeping at fixed physical Hamiltonian. It does not require four chemical potentials. Continuity and the small bound still depend on the unproved normalized RG hypotheses.

## 8. Attraction and the two-channel stability reserve

[proved, exact length-only obstruction] Consider the isolated length chain with physical hopping j_l>0 and attractive nearest-neighbor coefficient V_l, after linear terms are absorbed. For V_l<=-2j_l the two-site operator

    h_(j,j+1)-(V_l/2)(m_j+m_(j+1))

has eigenvalues 0,0,-V_l/2-j_l,-V_l/2+j_l in the empty, full and one-particle sectors. It is positive. Summing gives H_l>=V_l N_l. A state with N_l consecutive occupied sites on a large ring has zero hopping expectation and interaction V_l(N_l-1), provided 0<N_l<M. At fixed interior density these bounds imply e_l(nu)=V_l nu in the thermodynamic limit. The energy is linear, not strictly convex. At V_l=-2j_l the local symmetric one-particle vector also has zero energy in the shifted operator. This is the ferromagnetic boundary, not a strictly stable two-mode Gaussian point.

[proved, necessary reserve for the attractive reference] To exclude this mechanism in the pure length t-V submodel one needs

    -2j_l < V_l <=0, equivalently |V_l|<2j_l.           (33)

The positive-V statement is not the same claim; (33) is specifically about attraction. The argument proves failure at and beyond the negative threshold. It does not prove that all mixed retarded perturbations are harmless throughout the entire interval (33). No sharp interacting phase boundary for the full two-species model is claimed.

[proved, scheduled safety margin] With j_l=W, (5) implies |V_l|/(2j_l)<=2d_eq/J0. Choosing d_eq/J0 well below 1/2 excludes the isolated collapse. C45 actually chooses it much smaller to control the complete weighted coupling norm, including payload and mixed terms, relative to v_min. Higher-order matched interactions and the contribution E_p(D)-V_l also require a reserve. This is why dropping the scalar four-length table or taking lambda small alone is invalid.

[proved, mixed Gaussian stability test] In canonical branch-current normalization the candidate forward Hamiltonian has energy density

    H_G=(1/2) rho^T A_G rho,
    A_G,0=2pi diag(v_p,v_p,v_l,v_l),
    (A_G)_(a,b)=g_ab for a!=b,                         (34)

with diagonal velocity shifts included. This convention is for the Hamiltonian current coefficients; an RG reference using other density normalizations needs its explicit finite conversion. If

    ||A_G,0^(-1/2)(A_G-A_G,0)A_G,0^(-1/2)||<=r<1,

then A_G is positive by congruence. For its payload/length blocks, positivity is equivalently A_pp>0 and A_ll-A_lp A_pp^(-1)A_pl>0, by completing the square. Thus mixed terms can reduce length stiffness even when the isolated inequality (33) holds. Their natural cost is g_pl/(2pi sqrt(v_p v_l)); the slow payload mode matters. The stronger reserve r<=1/4 used by C45 also controls the scaling metric of the allowed operator lattice, whereas positivity alone does not.

[proved, no separation singularity in the Gaussian reference] With S=diag(1,-1,1,-1), the mode generator is S A_G/(2pi). It is similar to the symmetric matrix A_G^(1/2) S A_G^(1/2)/(2pi). If A_G>0, this matrix has two positive and two negative eigenvalues by inertia, none zero. Coinciding speeds do not invalidate this argument. Mixed density-current terms can make left and right magnitudes unequal. The collaborator's MULTIVELOCITY_REVIEW gives the corresponding local mode reconstruction in its density-covariance convention, so that reconstruction is not an extra missing two-channel theorem. Identifying the microscopic scaling fields with it still needs the source comparison.

[numerical, imported diagnostics] In verify_C44, M1 has V_l/|t| about -0.32 and j_l/|t|=1, inside (33); its calibrated finite-entanglement estimates are 1.97 and 2.00 and interval estimates 2.03 and 2.00. C1 has V_l/|t| about -6.3 with the same length hopping, outside (33); its selected length sector freezes and the calibrated result is c=1. M2 multiplies traceless mixed terms by 300 and remains near c=2 in those finite tests. These are observations of the projected chain. At finite diagnostic s/delta the much larger attraction is another reason to obey the schedule's equalization order. None of these observations proves (28) or bounds its constants.

## 9. What is new at two channels and what remains

[proved, comparison ledger] Relative to a one-channel action, this case adds a second free dispersion, four mixed marginal coordinates, a second exact charge and the associated nonresonance conditions. The four relevant coordinates and separate field normalizations must be included in the fixed-parameter root problem. The slow-speed and geometric-mean weights in (27), and the coupled stiffness Schur complement, are essential. The scalar attraction (4) is an additional physical instability test. The C45 hierarchy pays for all of these at fixed nonzero couplings before the infrared limit.

[proved, unchanged mathematical obligation] The scale dimension, Taylor normal form, paired cut identity, arbitrary external degree, normalization marks and need for a Gaussian source subtraction are already present in one channel. Mixed branch labels do not solve or refute the global routing. They also do not create an inverse difference of velocities. C50's summation closes in the same way once the four-branch version of (28) is proved.

[assumed, outstanding theorem] Establish N-action_2 with simultaneous routing, the marked normalized comparison (28) and the corresponding critical/source continuity estimates for the actual retarded kernels. This would give the conditional flow and source implications written above. The scalar thermal estimate P/TP remains a separate input for an interacting T^2 coefficient; no claim about it is proved here. With a complete compact two-boson Gaussian limit one obtains c=2 in the project's sum-of-modes convention, using C51, but neither that microscopic limit nor a common Lorentz speed follows from the present partial result.

## 10. Checks and provenance

[numerical] checks.py checks the four-bond expansion, the first-order six-coupling dictionary by independent polynomial differentiation of Gaussian energy, scheduled slopes and speed weights, the exact attraction threshold matrix, the C50 convolution closure and a four-branch source free-subtraction control. These are small algebra/normalization checks, not a numerical RG or a simulation of the microscopic wrapper. The results are in CHECKS.md.

[proved, scope record] All new files and Python execution stay in work_C66. No network commands, subagents or large computations were used. The two cited primary RG papers were consulted through the browsing tool. The parent project STATE.md was not modified because the task prohibits writes outside this folder. Local STATE.md records the handoff. There is no independent review of this note.
