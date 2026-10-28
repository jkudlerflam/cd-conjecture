# C45. Local reduction and the infrared budget

[proved, scope] The new results here are a finite-order local change of basis with an explicit remainder bound, a ground-relative auxiliary threshold for its reference Hamiltonian, and a consistent all-n hierarchy using C42's finite bounds and C32's wrapper estimates. The inference of a Gaussian infrared phase uses the single standard continuum input S below. It is not a proof of C40's exact connected-kernel hypothesis from a bare cell gap. The finite wrapper/exclusion and identity-transport statements are supplied by C44's completed audit of C8/C32 for this core. C44 finished while C45 was in progress; its RESULT.md and PROOFS.md were read before finalizing this note.

[proved, provenance] All requested inputs were read. The local SW source reviews in `../../contrib/counterexample_campaign/route_a/M2_REVIEW.md` and `../../contrib/counterexample_campaign/infrared/SOURCE_REVIEW.md` suggested using a fixed truncated generator instead of an all-order series. The proof in Section 3 is given directly, so no unread primary theorem or network lookup is needed. No files outside this folder were changed.

## 1. Fix the normalization before choosing a coupling

[proved, inherited from C39/C42] For a fixed positive integer n, put

    p = 2n+4,
    q = (wu)^2,        w = 1/4,
    S_q = sum_(j=0)^(n+1) q^(2j),
    V0 = q^(2n+2)/S_q^2,
    J0 = 2 V0 g^2 D^(-p).

Here u is C39's geometric readout parameter, not marker hopping. Below the marker hopping is called t_m, and the wrapper vacancy hopping is w_c. The normalized physical quadratic reference is

    h0(k) = -2 Gamma(w) cos k + 2 eta Gamma(wv) sin k
                         + mu D_w E_u D_w,
    D_w = diag(w^(2j)),    E_u = diag(u^(4j)).                 (1)

[proved, inherited] C42's dovetailed rational search returns eta,u,v,mu, root brackets, a positive pair margin delta_pay, an endpoint tolerance tau, and t_*>0 for every n. Its finite-D, finite-g bounds preserve all the brackets. With B=1+2 tr Gamma(w)+2 eta tr Gamma(wv)+mu, a common lower speed bound is

    v_cert = (t_*/2)^n delta_pay^(2n-1)/(B+1)^(n-1)>0.         (2)

This is the velocity of the uncongruent physical symbol in J0 units. The singular congruence used to locate zeros must not be used to read off velocities.

[proved, inherited wrapper normalization] If s multiplies the positive core in the wrapper, the common second-order port Hamiltonian is lambda^2 K2_bare/s. Thus the mobile payload energy unit is

    W = lambda^2 c_n,             c_n = J0/s > 0.              (3)

All symbols on the right are fixed before choosing lambda. The slowest payload speed is at least W v_cert before the remaining errors. J0 contains both the endpoint visibility and the high-rail propagation factor. Keeping only g^2 would miss factors q^(2n+2) D^(-p)/S_q^2. C40 permits putting them into its normalized velocity, but a wrapper tolerance expressed in J0 units must retain them explicitly.

[proved, scalar accounting] The length-independent scalar in K2_bare is a constant per ordered cell. It can be absorbed into the marker chemical potential. For four length pairs, a scalar table a_(p,q) also contains

    a_22-a_21-a_12+a_11,

which multiplies a length-density interaction. A marker chemical potential cannot remove it. The complete four-bond difference, scalar included, will be required to be small compared with J0, not just with the larger scalar source energy.

[proved, distinction from C26] C26's all-finite-order rho^2 factor uses a reflection family with a spectator logical register at rho=0 and a separate sign-changing symmetry. Neither its numerical fourth-order constants nor that factor is an estimate for this wrapper. No such improvement is needed here. We allow an unsuppressed higher-port remainder and pay the full factor s/J0 when choosing lambda.

## 2. The mobile Hilbert space and the finite wrapper interface

[proved, inherited coordinate identity] Root a configuration at a marker and order its cells. A cell space is the direct sum K=direct_sum_ell K_ell over lengths and internal states. The unperturbed Hamiltonian, after the coexistence energy per physical length is subtracted, has the form

    A = sum_j A_j,     A_j >= Delta Q_j,
    P_j = 1-Q_j,       dim P_j = 2^(n+1).                     (4)

P_j contains both selected lengths b=n+3 and b+1 and their complete 2^n payload codes. This representation can include all strings between consecutive markers, not only the legal low-energy strings. C44's cut-at-violations proof applies inside such a string: a selector violation costs at least U after its at-most-two fragment losses, wrong lengths obey its linear secant bound, and the other selected-length cell sectors have the stated positive gaps. A length-one ## interval is itself penalized. Thus these excluded states belong to Q_j in (4); they are not silently discarded. The sector with no markers is treated separately below. A marker move transfers one boundary site between neighboring factors. Ports act on the same two factors. They preserve total physical length, total payload U(1), marker number and the cyclic root identification. Moving boundaries therefore do not make the retained projectors time dependent in these coordinates.

[proved, inherited from C32 and the completed C44 audit] The lower gap uses the minimum of the selector, filter, helix, charge and core bounds. Its charge splitting for the two-hole length is of order U(w_c/U)^(2n+4), with a positive finite coefficient. A lower bound for the charge/core combination is

    delta_ch d_min gamma /
       ((sqrt(delta_ch)+sqrt(M_core))^2+d_min gamma).          (5)

Use the actual scaled core/filter gamma and M_core, and minimize this with the other exclusions. This defines a positive Delta at fixed nonzero wrapper parameters. Neither U nor the unwrapped core gap is a permissible replacement for Delta.

[proved, bounded perturbations in this representation] In each well-formed initialized cell the arrow grammar leaves one token, so the conditional vacuum readout has norm at most one even if its length is wrong. Boundary port and transport terms have bounded incoming/outgoing degree and fixed compensator norms. The weak marker chemical potential is bounded per cell. Hence, within this invariant well-formed representation, write

    H = A + lambda V + lambda^2 F,                            (6)

with bounded one-cell and nearest-cell interactions V,F. F includes the retained readout, marker motion t_m=lambda^2 c_n/beta_n, and the chemical adjustments. V is the compensated port. Its norm is the full physical norm, not its small projection. Beta_n is the positive projected marker overlap. Pure marker motion is the identity on the transported payload.

[proved, inherited finite wrapper contract F_wrap] The application to the complete physical bond uses C44 PROOFS Sections 1-4 and 8: comparison with wrong occupation, illegal arrow and marker-free sectors, the positive bounds just specified, and the exact identity-transport quotient. C44 proves these for C39 for every n. The local theorem below includes the mark-rooted excluded cell states once the just-described one-cell bounds are used. It is not an independent enumeration of the 4464-letter physical bond. The marker-free comparison below is separate.

[proved, physical constraint caveat] The tensor-product calculation can be done before restricting total length. Every transformation constructed below commutes with that restriction. The ground threshold in Section 4 is first an inequality relative to the unrestricted rooted tensor-product ground energy at a chosen chemical potential. It does not by itself identify the minimizing marker number on a fixed physical ring. That last comparison uses F_wrap and the compressibility/rounding argument in Section 8, not an extensive bare-Q inverse.

## 3. An explicit volume-independent local SW bound

[proved, interaction norm] For interval interactions define

    ||Phi||_a = sup_j sum_(X contains j) exp(a |X|) ||Phi_X||,

where |X| is the number of ordered cells. A one-site term can be kept separately. If the largest supports of B,C have sizes b,c, respectively, direct summation of overlapping pairs gives

    ||[B,C]||_a <= 2(b+c) ||B||_a ||C||_a.                    (7)

Indeed a commutator is zero for disjoint supports, its norm is at most twice the product, and exp(a|X union Y|)<=exp(a|X|)exp(a|Y|). For a fixed marked cell, count first the terms with the mark in X and then those with it in Y. The sums over the overlap cost at most |X| and |Y|. This proves (7) without a volume factor.

[proved, local inverse] For an interaction B_X, let P_X=product_(j in X) P_j, Q_X=1-P_X and R_X=Q_X A_X^(-1)Q_X. Define

    B_X^d = P_X B_X P_X + Q_X B_X Q_X,
    B_X^o = P_X B_X Q_X + Q_X B_X P_X,
    L_X(B_X) = R_X B_X P_X - P_X B_X R_X.

Then L_X is anti-Hermitian for Hermitian B_X,

    [L_X(B_X), A_X] = -B_X^o,
    ||L_X(B_X)|| <= 2||B_X||/Delta.                           (8)

Both diagonal and off-diagonal maps have norm at most 2, a loose bound used below. These statements require only (4). They also hold on the direct-sum cell space: the reduced inverse is bounded, although A_X need not be bounded above.

[proved, construction] Apply these operations separately to each declared support, and set

    S1 = L(V),
    G2 = F + [S1,V^d] + (1/2)[S1,V^o],
    S2 = L(G2),
    S(lambda) = lambda S1 + lambda^2 S2,
    D(lambda) = A + lambda V^d + lambda^2 G2^d.               (9)

S1 has range two and S2 has range three. Expansion of the conjugation through order two, using [S1,A]=-V^o, shows

    exp(S) H exp(-S) = D + R,
    P D P = lambda^2 (P F P - P V Q A^(-1) Q V P),            (10)

since PVP=0. P denotes the all-code projection. The second expression is a local second-order coefficient, not the use of a global perturbed-band resolvent. On the present core cross-boundary port contractions vanish, so it is the sum of the complete four length-pair source bonds and PFP. The entire finite-g K2 is included.

[proved, quantitative constants] Put

    v0 = max(1, ||V||_a),      f0 = ||F||_a,
    G = f0 + 48 v0^2/Delta,
    C_SW = 2^16 (v0 G/Delta + v0^3/Delta^2).                  (11)

Equation (7) gives ||G2||_a<=G, ||S1||_a<=2v0/Delta and ||S2||_a<=2G/Delta. If

    0<lambda<=min(v0/G, Delta/(192 v0)),                     (12)

then the exact remainder in (10) satisfies

    ||R||_a <= C_SW lambda^3.                                (13)

The constants contain the compensators, both buses, all physical marker/readout coefficients and the actual Delta. They have no dependence on cell count. There is no claim that they are moderate or uniform in n.

[proved, convergence argument for (13)] This is a convergent BCH estimate for the fixed polynomial S, not a convergent all-order SW expansion. Write

    z=(lambda v0+lambda^2 G)/Delta,
    s0=||S||_a<=2z,       z<=1/96,
    C=[S,A]=-lambda V^o-lambda^2 G2^o,
    B=lambda V+lambda^2 F.

Then ||B||_a<=Delta z and ||C||_a<=2Delta z. For a range-three B or C, the kth nested commutator with S has support size at most 3+2k. Repeated (7) implies

    ||ad_S^k B||_a <= 12^k k! s0^k ||B||_a,

and the same for C. Use the identities

    exp(ad_S) A = A+C+sum_(k>=1) ad_S^k C/(k+1)!,
    exp(ad_S) B = B+sum_(k>=1) ad_S^k B/k!.

The terms with k>=2 are bounded by

    3 Delta z (24z)^2/(1-24z)
                  <= 3456 Delta z^3
                  <= 27648 lambda^3 v0^3/Delta^2.             (14)

The terms of degree at least three in [S,B]+[S,C]/2 have norm at most

    96 lambda^3 v0 G/Delta + 48 lambda^4 G^2/Delta
                    <=144 lambda^3 v0 G/Delta.

Combining these proves (13). The second identity uses only bounded local interactions. The first uses the bounded commutator [S,A], so it remains valid on the form domain for the possibly unbounded direct-sum A. It can also be proved first with a length cutoff and then with bounds independent of that cutoff.

[proved, locality and symmetry] All generators are connected interval interactions and commute with the exact charges. The unitary is a finite-time local evolution and has decaying spatial tails. Restriction to a fixed physical length and the cyclic quotient commute with it. A process on r neighboring selected cells acts on O(rb) physical sites, so finite n preserves spatial summability in physical coordinates for the retained kernels. The transformation does not create an extra gapless degree of freedom or a new local fractional charge. Preservation of the physical operator lattice uses the identity-transport part of F_wrap.

[proved, treatment of non-well-formed internal cells] Within the well-formed core-arrow sectors the readout bound used in (6) is uniform over lengths. In an illegal core sector the sum of physical readout projectors can instead grow with the number of tokens. This is not a bounded operator on the all-length direct sum. Its coefficient in (33) is positive, and core-arrow violations and token number are invariant under the ports and transport. Retain this positive onsite contribution exactly in the illegal blocks, rather than including it in f0. The local generators connecting the code have zero matrix elements into those invariant blocks. The untouched positive contribution cannot weaken (15) or the sector comparisons. This qualification avoids incorrectly applying a bounded-interaction theorem to an extensive per-cell readout on malformed words.

[proved, coefficient-norm conversion] For the retained local space of dimension d_P=2^(n+1), expanding an operator on r cells in a normalized matrix basis costs at most d_P^(2r) in coefficient l1 norm. Conversion to fermion monomials and an external-leg weight also costs only exp(O(nr)). One can, for example, take a=16(n+1) and reserve part of this weight for two spatial moments and leg weight 2. Thus (13) bounds a suitable complete instantaneous matching norm, not just a few quadratic coefficients. It does not yet bound exact retarded kernels after the remaining massive fields are integrated out.

[proved, selection-rule refinement not used] The alternating-cell excited-token sign removes odd bulk port coefficients after compression. A single port insertion with any marker insertions cannot return to the code. Thus the first projected corrections can be organized as lambda^4, lambda^2 t_m and t_m^2. Equation (13) instead controls the full transformed Hilbert space, including remaining off-diagonal terms, and safely uses the weaker cubic bound. Winding processes on an odd ring are not local cubic bulk coefficients.

## 4. A massive threshold for the reference without an extensive band gap

[proved] Suppose D=A+sum_X W_X with each W_X locally block diagonal relative to P_X, as in (9). Let E_P be the ground energy of PDP on the unconstrained rooted tensor product and let J=||W||_0. Then

    D >= E_P I + (Delta-2J) sum_j Q_j.                       (15)

This does not require [D,Q_j]=0 for each j.

[proved, proof] Choose any density matrix rho_j supported in P_j and define the unital completely positive extension of a code operator O by

    E_j(O)=O on P_j, and tr(rho_j O) I on Q_j.

The product map E extends the full code algebra. Positivity gives

    E(PDP) >= E_P I.

On a support X, W_X and E_X(P_X W_X P_X) agree in P_X, are block diagonal between P_X and Q_X, and have norms bounded by ||W_X||. Their difference is therefore bounded below by -2||W_X||Q_X. Since Q_X<=sum_(j in X)Q_j, summation gives (15). This is also a proof for arbitrary finite code degeneracy; no many-body gap in PDP was used.

[proved, application] From (11),(12), J<=2lambda v0+2lambda^2 G<=4lambda v0<=Delta/48. Consequently

    D >= E_P I + (Delta/2) sum_j Q_j.                         (16)

The auxiliary sector of the reference has a positive threshold above its ground. High-energy code states can still overlap auxiliary energies. Equation (16) is deliberately a ground-relative inequality, not isolation of the whole extensive code band. Applying it directly to the original bare projection after adding R would be invalid because R dresses the vacuum locally.

[assumed, passage from (16) to the exact phase] Stability of the massive threshold and of the critical subsystem under the small local R is part of the standard continuum input S. Equation (16) supplies its massive reference. It does not prove by itself an exponential bound on a bare off-code source correlator, which can contain retained low-energy spectral weight after dressing.

[proved, marker-free sector] Port and marker-motion perturbations require a marker and vanish in a marker-free configuration. The readout added here is positive. C44's marker-free lower bound therefore stays nonnegative. The selected trial states have energy density e_0+O(q_phys), where e_0<0 is the unperturbed coexistence energy and q_phys=lambda v_phys+lambda^2 f_phys bounds the physical per-site perturbation. Choosing q_phys<|e_0|/4 retains a negative-energy selected phase. No gapless marker-free competitor is being integrated into the retained fluid.

[proved conditional on S, length constraint and remaining sectors] To compare different M at fixed L, subtract the tangent physical energy density of the selected phase. The required change from e_0 is O(W), after the common port scalar is canceled by the marker chemical potential in (34). The wrong-length cost in C44 is at least gamma_len(1+|ell-b|). It absorbs an O(W)|ell-b| tangent change when W is sufficiently small. A bounded O(W b) constant per defective cell is absorbed in Delta as well. Repeating (15) with these reduced positive bounds leaves, for example, Delta/4 per excluded cell. In the retained block S gives the positive compressibility and finite-size zero-mode behavior described in (34). The constrained good ground energy differs from the grand-canonical tangent value by O_n(W/M); for the finitely many small M the boundary cost is bounded by O_n(W) instead. Choose W times those fixed constants below Delta/16. This rules out an excluded-cell ground state for sufficiently large physical rings, even though the full code spectrum is extensive. It is the constrained-ground application of (15), not an assumption of a gap between the complete interacting bands. Stability after adding R uses S, as explicitly stated above.

## 5. One explicit continuum input and the connected low-frequency statement

[assumed, standard input S] Use weak-coupling Wilsonian matching and Gaussian phase stability for finitely many one-dimensional modes coupled locally to a massive complement. In the form needed here, assume the following single statement.

[assumed, S in quantitative form] A reference with (i) a positive local auxiliary threshold Delta/2, (ii) a retained number-conserving quadratic dispersion with finitely many simple Fermi roots and positive low-order momentum margins, and (iii) the identity-transport compact-boson operator lattice admits continuous local Wilsonian matching. If its retained interaction is epsilon-close to the quadratic reference in the weighted coefficient norm, and an additional symmetry-preserving full-space local perturbation R obeys ||R||_a/W<=epsilon, the matching constants for vertices, two frequency/spatial moments, background filling derivatives and the collective irrelevant norm are finite. They can be chosen uniformly for sufficiently small W/Delta and ||D-A||_a/Delta. Subsequent RG obeys the usual weak-coupling marginal majorant and irrelevant decay stated in Section 7. Scalar ground energies and allowed quadratic/current terms are matched, not discarded. The local operator lattice contains no extra fractional vertex.

[assumed, interpretation] This is the usual decoupling of massive local degrees of freedom plus Luttinger-liquid RG stability, stated with the uniformity that the hierarchy needs. It is the only continuum input. It includes smooth response of dressed fillings, not just an energy-density estimate. It is not a theorem proved in C45 or by the local SW source reviews. In particular, calling the reference gap positive without supplying the uniform matching and symmetry hypotheses would not establish S.

[proved conditional on S, concrete replacement for C40 M2] Choose the matching norm to sum connected irreducible retained vertices K_(X,m), with one time anchored at zero,

    sup_j sum_(X contains j) exp(a' |X|)
       sum_m 2^m/m! integral product dt_i
       (1+Delta sum_i |t_i|)^2 |K_(X,m)(t;0)|_1.              (17)

Take a'<a after the coefficient conversion above. Include scalar length terms, instantaneous delta measures and the first two filling derivatives. All propagation reducible through retained lines is subtracted. Within the low-external-frequency region of S, the change due to R in (13) satisfies

    ||delta S_R||_LE <= A_dec,n C_SW lambda^3,
    ||delta S_R(omega)-delta S_R(0)||_LE
           <= A_dec,n C_SW lambda^3 |omega|/Delta.           (18)

Here (17) refers only to massive matching with the full retained lattice fields left unintegrated. It does not include the later gapless band-matching loops, whose time scale can be 1/W. A_dec,n is one of S's finite constants. Equation (18), rather than a bound on an unrestricted imaginary-time interaction-picture operator, is the concrete connected low-frequency remainder criterion for this chain.

[proved conditional, comparison with the actual slow scale] Including the static core and equalization errors epsilon_static, a sufficient normalized bound is

    epsilon_total <= epsilon_static
              + A_dec,n C_SW lambda/c_n.                    (19)

Since c_n=J0/s>0 is fixed first, the last term can be made smaller than any prescribed tolerance by one fixed nonzero lambda. It is allowed to lose all inverse visibility, conversion and propagation factors. At retained frequencies |omega|<=W times a fixed normalized bandwidth, the frequency-dependent part in (18) has an additional small factor W/Delta. Taking lambda smaller improves both ratios.

[proved, optional untransformed second-order kernel check] At second order every source excursion has a fast denominator at least Delta. For |omega|<Delta/2 the resolvent identity gives a norm bound of the form 2 lambda^2 ||Z||^2 |omega|/Delta^2, with connected-support factors included. In W units at |omega|<=W E_top, it is at most 2 lambda^2 ||Z||^2 E_top/Delta^2. The static local transform has already absorbed its normalization effects, so this term need not be added twice. It gives an independent scale check that frequency dependence cannot leave an order-one error at the payload bandwidth as lambda decreases.

[proved, what has not been established] Equation (13) does not imply (18) without S. Nor do finite-order connected coefficients imply convergence of C40's all-order series. The active-defect count misses histories that create a separate defect pair and join it later and histories carrying a changed retained state after the local defect disappears. Bare-Q spectral bounds and subtraction of a scalar expectation do not fix this. C26's fourth-order normalization term illustrates the necessary retained-propagation subtraction. The present reduction avoids claiming that those exact cancellations were proved from a bare defect gap.

[heuristic, why S is the physics input here] After (9)-(16), the problem is an ordinary strict Gaussian basin with a separated massive local complement and an arbitrarily small quasi-local perturbation. Auxiliary interactions of order lambda may exceed W, but are small relative to Delta and live outside the retained block. They are not additional marginal couplings of order lambda/W. This is the familiar setting for integrating out massive fields in effective field theory. There is no additional compiler-specific assumption of an isolated extensive band. A rigorous proof of S for the full class would still be substantial work.

## 6. The allowed operator lattice and a mixed margin

[assumed within S, inherited dictionary] Let z_i be integer chiral payload charges, i=1,...,2n, and let h be the integer length-density harmonic. Exact payload U(1) imposes sum_i z_i=0, and exact marker conservation forbids a length dual-field vertex. In normalized cell coordinates the free primary dimension and momentum are

    D_free(h,z)=h^2+(1/2)sum_i z_i^2,
    Q(h,z)=sum_i z_i p_i-2pi h nu modulo 2pi.                 (20)

Cyclic gluing changes allowed zero-mode sectors, not this local lattice. This last assertion relies on the exact transport audit, not only on conservation of marker number.

[proved within (20)] The nonzero primaries of dimension at most two are payload bilinear differences, distinct-branch two-in/two-out quartics, h=+/-1 alone, and h=+/-1 times a payload bilinear. Repeated same-chirality fields have a Pauli derivative cost and do not supply an omitted low-dimension vertex. Equal incoming/outgoing branch multisets yield densities or current products. Total U(1) forbids anomalous pairing but not pair transfer; C42's quartic margin is essential.

[proved, density choice] In [1/4,3/4] exclude solutions of

    2pi nu = p_i-p_j modulo 2pi, i != j.

There are at most K=2n(2n-1) such points. Include the two interval endpoints when locating the largest complementary interval. Its length is at least 1/[2(K+1)]. Choose nu at its midpoint, or a rational certified sufficiently close to it. One can reserve

    delta_mix >= pi/[4(K+1)]>0,                              (21)

and pure h=+/-1 has margin at least pi/2. No all-order Diophantine condition is necessary. Higher harmonics and higher payload vertices have D_free>=3.

[proved within Gaussian theory] Let H0 be the free stiffness with all actual speed ratios and H its dressed value. If

    eta_G = ||H0^(-1/2)(H-H0)H0^(-1/2)|| <= 1/4,

the Gaussian scaling metric is bounded below by sqrt(3/5) times the free metric. Thus every remaining nonzero primary has

    D >= 3 sqrt(3/5) = 2+sigma,
    sigma = 3 sqrt(3/5)-2 > 0.                               (22)

This follows from matrix monotonicity of H^(-1) # (J H J), as derived in C40. It includes density-current and length-payload entries and does not assume comparable velocities. Positivity alone would not give (22).

## 7. RG down through all bandwidths and to finite-size energies

[proved, finite reference data] On the normalized physical dispersion choose disjoint root patches with branch derivative at least v_*/2, gap d_*>0 to other branches, radius r_*>0, and nonzero energy e_off>0 outside. Include the length band with hopping one and speed 2sin(pi nu). Here v_* can be taken no larger than min(v_cert,1). Let M1 bound the first derivative on the patches. These data are computable by finite interval bounds; all are positive at fixed n. Use half the margins in (21) and C42 to reserve a dressed mismatch delta>0.

[proved, root budget] For a C1 symbol error e, it suffices that

    e < d_*/4,  e < e_off/2,  e < v_* r_*/4,
    e(1+4 M1/d_*) < v_*/4,
    4e/v_* < delta/8.                                      (23)

The projector estimate and Hellmann-Feynman bound in C40 give the derivative condition. Monotonicity then gives one root per patch and displacement at most 4e/v_*. The off-patch gap excludes extra roots. Interacting self-energy and filling shifts satisfy the corresponding bounds through S, not through finite matrix perturbation alone.

[proved, fixed matching interval] In W units take

    e_lin = min(e_off/2, v_* r_*/8),
    e_det = min(e_lin, v_* delta/32),
    e_top = max(1, sup_k ||h0(k)||, 4),
    ell_* = max(0, log(e_top/e_det)).                         (24)

These are independent of lambda. All retained dispersions must be kept until their own linear regimes are reached. A narrow bandwidth is not a mass gap. Matching constants may contain powers of speed ratios; put them into S's A_n. Using only log(1/delta) would omit the hierarchy cost in ell_*.

[assumed within S, RG estimates] In a collective coefficient norm, let x be the marginal deviation, including forward terms, and Y the irrelevant norm. With finite A_beta,C_irr and positive sigma from (22), weak-coupling RG gives

    X0 = x0+C_irr(Y0+Y0^2)/sigma,
    x(ell) <= X0 + A_beta integral_0^ell x(t)^2 dt.

Thus x(ell)<=X0/(1-A_beta X0 ell). Matching and filling/stiffness shifts are bounded by A_n epsilon_total/v_*, with additional finite constants absorbed into A_n. Flavor sums are included in this collective norm.

[proved conditional on these estimates] Choose a tolerance epsilon_n>0 so that (23), the filling shifts, eta_G<=1/4 and

    A_beta X0 ell_* < 1/2                                   (25)

hold with strict reserves. This is an explicit finite list of inequalities. For example, first choose a positive q_patch satisfying the first four bounds in (23). Enlarge A_n>=1 to include A_beta, C_irr, 1/sigma, the filling constants and their products. Then the conservative choice

    epsilon_n <= min(1, q_patch/(1024 A_n),
                     v_* delta/(1024 A_n),
                     sigma v_*/[1024 A_n(1+ell_*)])          (26)

suffices after using the same enlarged constant for the complete drift bounds. A_n is a standard matching constant, not a number supplied by C42. Accordingly the final schedule is quantitative conditional on S, rather than a numerically certified value of lambda for arbitrary n.

[heuristic under S, all scales] Before W e_det, the finite majorant prevents any nearly resonant non-Gaussian coupling reaching strong coupling. At energies below W e_det the dangerous oscillatory vertices have no support in the low-energy patches. The only surviving marginal terms are Gaussian; all remaining nonlinear vertices obey (22). The bootstrap gives, for E<=W e_det,

    Y(E) <= C Y(W e_det) (E/(W e_det))^(sigma/2),

with finite integrated Gaussian and filling drift. The complete stiffness remains inside its reserved neighborhood, and each velocity remains positive. There are n payload bosons and one length boson.

[proved, meaning of control at 1/L] No fixed nonzero extensive remainder norm is smaller than 1/L for all L. Such a demand would be false even for a harmless velocity renormalization. The RG statement is that corrections to the renormalized Gaussian finite-size spectrum are irrelevant: their ratio to an O(v/L) oscillator spacing tends to zero with a positive power once L exceeds the fixed crossover length. Marginal forward terms change that spacing by a small fixed relative amount and are included in the dressed velocities. The local error never becomes a relevant mass at a later scale under (25),(22). This is the appropriate answer to the requested control down to 1/L.

## 8. A complete order of choices in n

[proved conditional on S and F_wrap] The following choices close the scale inequalities. Each step leaves a strictly positive interval for the next. Nothing depends on the eventual ring size.

[proved, step 1] Run C42's dovetailed search at the desired n. Obtain h0, tau, t_*, delta_pay and its rational root brackets. Choose nu and delta_mix by (21), compute the patch data, and select epsilon_n using (26). Set a=16(n+1) and, for translating a two-cell bond estimate into the weighted static interaction norm, use the loose constant

    C_b = 16(n+1)^2 exp(3a),
    theta = min(tau, epsilon_n/(16 C_b)).                     (27)

C_b pays for two incident bonds, hopping to many-body norm conversion and spatial weights. The subsequent coefficient/continuum conversion is already in A_n.

[proved, step 2, finite D] With r=1/16384 and B0=128 choose an integer

    D >= max(2/r, 8n B0(1+eta)/(theta r^(p+1))).              (28)

C42 bounds the normalized quadratic symbol and its first derivative by E_D<=theta/4. The derivative obeys the same estimate because only nearest-neighbor hopping is changed. This error is included in the static budget, not only in the bracket certificate.

[proved, step 3, nonzero g] Choose a positive dyadic g satisfying

    g<=1/4,
    g^2 <= theta/[192 n(1+eta)D^p],
    g^2 <= epsilon_n/
       [16 C_b (2+64n^2) 3(1+eta)D^p].                       (29)

The first two are C42's finite-g quadratic certificate with a reserve. The last bounds the complete non-Gaussian second-order remainder by epsilon_n/(16 C_b) in J0 units. The full conversion remainder, including its quadratic part, is also small. A common scalar is absorbed into the marker chemical potential. Compute J0 now and freeze g. Reducing lambda alone cannot suppress this second-order conversion error relative to W.

[proved, step 4, complete wrapper equalization] Use C32's quantities for the two-bus source, taking maxima over all four length pairs. In its notation x=||PZ||, y=||QZ||, L_c bounds ||W_core,pair||, a_c is a uniform positive lower bound for P W_core,pair P, and L_bond bounds the sum of the positive bare core-bond norms. Avoid confusing a_c with the spatial norm exponent a. C44 now makes these constants explicit. Set

    h_c,max=K_c+D+12,  a_c,max=K_c+D+5,
    M_c=(k-1)h_c,max+2a_c,max,
    L_c=2M_c,  L_bond=2(k-1)h_c,max,
    a_c=1-|g|,   x=sqrt(3(1+eta)V0).

For zeta<=1/2, a uniform compensator bound is

    M_h0=16/[3 sin^2(1/(2b))],     y<=x M_h0^2.

Thus the source norms needed here are bounded before any equalization limit. Let

    d_eq = J0 epsilon_n/(16 C_b),
    zeta = w_c(k+1)/(U-w_c(k+1)),    k=n+2.

Choose w_c/U>0, with zeta<=1/2, such that

    x^2 zeta^2 L_bond/a_c^2 <= d_eq/2.                       (30)

The compensated ground-source x is the bare source and does not diverge in this step. For an explicit choice let zeta0 be the minimum of 1/2 and sqrt(d_eq a_c^2/(2x^2 L_bond)), and take

    0<w_c/U<=zeta0/[2(k+1)(1+zeta0)].

If the numerator x or L_bond vanishes, omit its vacuous constraint. After fixing this ratio, calculate the fast gap delta_fast>0 and the finite compensator/source norms, including y. No asymptotic gap coefficient is needed: C44 proves for N_p=k+p, m_p=binomial(k+p,p), a_p=(w_c/U)/(N_p-1),

    delta_ch,p >= w_c a_p^(2kp)/(m_p kp)>0.

Minimize these two explicit bounds with J_h/[(b+1)(2b+1)] to obtain delta_fast. The occupation-filter penalty does not enter the source equalization bound because its ground factor is invariant, but it does enter the final Delta in (5).

[proved, step 4 continued] Put

    A1=2 L_c^2 x^2/a_c^2+4 L_c x y/a_c+y^2,
    A2=2 L_c^2 y^2/a_c,
    t=s/delta_fast.

Choose

    0<t<=min(1, a_c/(2 L_c^2), d_eq/[2(A1+A2)]).             (31)

C32's Schur-complement estimate then gives

    max_(p,q) ||s K2_(p,q)-K2_bare|| <= d_eq.                 (32)

This includes the scalar table and all noncommuting core terms. It does not replace the fluctuating adjacency operator by its mean without an error. In W units, (32) contributes at most C_b d_eq/J0<=epsilon_n/16. The large D and small g may force an extremely small w_c/U, which in turn forces a still smaller s. Their gaps remain positive. This has no circular dependence on lambda.

[proved, step 5, choose length and readout scales] Recompute Delta from the full wrapper after (30),(31), including the small charge tunneling splitting. Compute beta_n>0 and the compensated port norm. Set

    t_m = lambda^2 J0/(s beta_n) = W/beta_n,
    f_read = J0 mu S_q/(s q^2),
    H_read = lambda^2 f_read R_read.                         (33)

The readout's compression is the field in (1), equally on both lengths. Neither its off-code part nor the physical norm 1/beta_n may be omitted. All such terms belong to F in (6). Include the scalar/marker chemical adjustments in f0 as well.

[proved, marker-density choice at leading order] The physical constraint is L=M(b+nu). With the length hopping set to W, the free energy per ordered cell, in W units, has the form

    e(nu)=e_pay-(2/pi)sin(pi nu).

An adjustable scalar c per cell gives E=M[c+W e(L/M-b)]. Stationarity with respect to M at the chosen nu requires

    c=W[(b+nu)e'(nu)-e(nu)].                                 (34)

Its curvature is W(b+nu)^2 e''(nu)/M>0 in the chosen interior interval. Thus one marker chemical potential selects that density; no independent flavor chemical potentials are required. S controls the small change of this convex energy and of the equilibrium density. Nearest-integer rounding of M changes nu by O(1/M) and the energy by O(W/M), which is ordinary zero-mode behavior. F_wrap excludes states outside the selected-cell phase. This reasoning does not assert unique finite-ring ground states.

[proved, step 6, choose lambda last] Calculate v0,f0,G,C_SW in (11) using the physical operators fixed above. Let A_dec,n be the constant in (18). One sufficient fixed nonzero choice is any lambda below

    min( v0/G, Delta/(192 v0),
         epsilon_n c_n/[16 A_dec,n C_SW],
         sqrt(Delta/[16 c_n (b+1) e_top]) ).                       (35)

The last condition puts the retained bandwidths and the tangent energy per cell below the auxiliary threshold. Enlarge the fixed factor (b+1)e_top if necessary to include the finite constrained-ground/boundary constant in Section 4. Also impose lambda<=1 and

    lambda <= |e_0|/[8(v_phys+f_phys)],

using the actual physical per-site coefficient norms, to retain the marker-free comparison. This is another finite positive restriction. If the original frequency-kernel estimate is used separately, also require

    lambda^2 <= epsilon_n Delta^2/(32 B_time,n e_top),        (36)

where B_time,n bounds the finite second-order source/moment norm. This extra restriction is sufficient and unnecessary in the static-transformation formulation. Conditions (35),(36) are finite positive numbers for each n. Taking, for example, a positive dyadic strictly below their minimum fixes every coupling before L is increased.

[proved conditional, total budget] Equations (28),(29),(32),(35) put the quadratic approximation, conversion interactions, all scalar and payload equalization errors, and the complete local SW remainder inside the tolerance of Section 7, with unused reserve for matching and running. Mixed length-payload terms from the unequal bonds are already in (32). Higher marker excursions and readout leakage are included in R. If one instead resolves their orders, their sizes are t_m^2/Delta, lambda^2 t_m/Delta^2 and lambda^4/Delta^3 times fixed norm factors. Dividing by W makes each tend to zero with lambda. There is no leftover order-lambda^2 length interaction of uncontrolled relative size once (32) is enforced.

[proved, explicitness and limitation of this schedule] The reference search and D,g bounds are terminating finite algorithms from C42. Wrapper parameters are given by explicit inequalities and finite-cell spectral data. The final tolerance and lambda use the constants of the single standard continuum input S. They are not a new numerical certificate of the complete 4464-letter chain at these parameters. No practical lower bound on the crossover energy is claimed. The scheme establishes existence of fixed nonzero coefficients under the stated physics input, not feasible experimental energy scales.

## 9. Failure tests, checks and conclusion

[proved, genuine failures of weaker arguments] The following substitutions fail: comparing a global norm O(L lambda^3) to 1/L; dividing a wrapper error only by g^2; dropping scalar length interactions; using U as the smallest wrapper gap; reducing lambda to repair an order-lambda^2 g^4 interaction at fixed g; inferring exact low-frequency kernels by summing C40's finite-order C_r; importing C26's reflection suppression; or using only stiffness positivity with an extremely slow mode. None is used above.

[proved, conditional closure] No scale obstruction remains after C42's nonresonant reference and the inherited wrapper inequalities are supplied. Each subsequent error has an independently adjustable parameter, and the final local remainder is o(W). The hierarchy in n can be arbitrarily severe without contradicting fixed-n infrared stability.

[numerical, independent controls] `checks.py` tests (9),(10) on two and three sites with local dimension four, code dimension two and a conserved U(1), in full dimensions 16 and 64. The second-order identity errors are below 6.6e-17. For lambda=.02,.01,.005, ||R||/lambda^3 approaches .341 for two sites and .621 for three. The defect-threshold inequality has minimum eigenvalue no smaller than -3.3e-15, consistent with zero roundoff in its code block. Exact rational checks verify the BCH constants and a synthetic schedule. These are checks of the new algebra, not tests of the C39 wrapper or of central charge.

[heuristic under S and F_wrap, infrared conclusion] The wrapped chain has n+1 compact Gaussian modes at every fixed n for the schedule above, in the project's multi-velocity convention. Each nonchiral boson contributes one to c, so c=n+1. A common Lorentz velocity is not implied. The Casimir coefficient is velocity weighted and cannot by itself measure this unweighted c. Taking the thermodynamic/infrared limit for each n before comparing n therefore gives unbounded c at d0=4464, conditional on the finite wrapper contract and the stated conventional continuum input.

[heuristic, final physics-standard verdict] Yes. Combining C39/C42, C44's now-completed finite audit and the reduction here establishes c=n+1 at every fixed n at the project's theoretical-physics standard, using the explicit standard input S. There is no remaining independent all-n band or wrapper hypothesis awaiting construction in this argument. The claim remains conditional on ordinary massive-mode matching and Gaussian RG, not a rigorous microscopic CFT theorem.

[proved, status boundary] C45 does not prove exact C40 kernels without S, independently rerun C44's full wrapper calculations, measure microscopic c, or establish uniqueness on every ring. The complete chain has multiple velocities; the conclusion uses the same sum-of-mode-central-charges convention as the project. If the intended conjecture demanded a single common Lorentz velocity, that additional claim is not established here.
