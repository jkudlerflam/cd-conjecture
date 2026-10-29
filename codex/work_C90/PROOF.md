# C90. Infrared transfer for the vacancy mobile chains

## 1. Status, conventions, and the claim being transferred

[Definition] A tag governs the following paragraphs and displays until the next tag or section. [Proved] means a deduction given here. [Imported, proved] records a supplied result with its original scope, not a new independent verification. [Assumed] names an unproved input. [Numerical] denotes a finite arithmetic check. No numerical check below establishes an infrared limit.

[Proved, scope] The reusable conclusions are (i) C48's abstract massive matching, conditional on its representation, gap, heat trace, and smallness hypotheses, (ii) C58's thermodynamic operations, conditional on the appropriate uniform pressure comparison, and (iii) C51's identification of a complete compact Gaussian theory. The microscopic verification in C58 for d0=4464 is not a verification for a vacancy chain. In particular, the requested three remaining labels, N, TP, and dressed-defect comparison, are sufficient only if N is explicitly enlarged to include vacancy initialization and Gaussian completeness, and the defect comparison includes thermal traces. They are not sufficient with their narrower existing meanings.

[Definition] Write K_pen for the grammar penalty and K_gap for the Luttinger parameter. Let L be physical length, N_occ the number of material letters, M the number of material cells, and b their selected material length. In the good sector N_occ=bM. Set rho=N_occ/L and ell=L/M=b/rho. All couplings are fixed and positive before taking L to infinity and then T to zero. Density sources vary in a fixed neighborhood of the chosen zero-temperature source. They are not temperature-dependent counterterms in the physical Hamiltonian.

[Imported, proved] The relevant reference data are as follows. The d=24 entries refer to the original 23-letter material construction. The two d=14 constructions are different chains.

| Chain | Material letters | Vacancy colors | n | b | Proposed Gaussian count |
| :-- | --: | --: | --: | --: | --: |
| C88, d=13 | 11 | 2 | 13 | 15 | 14 |
| C83, d=14 | 12 | 2 | 13 | 14 | 14 |
| C85, d=14 | 11 | 3 | 13 | 15 | 14 |
| Original d=24 lift | 23 | 1 | 23 | 24 | 24 |

C88's material certificate has 26 simple Fermi points, root separation greater than 0.0156742281 and pair-sum separation greater than 3.2513114e-6. Its selected length is 15, not 14. Its finite material schedule is C88 PROOF (16), with g=2^-272, lambda=2^-608 and K_pen at least 2^2560. Its mobile restrictions on K_pen and J/t are additional. None of these finite certificates proves strong-norm vacancy initialization.

## 2. Coloring and the exact retained charge variables

[Imported, proved] C88 PROOF sections 1 to 4 establishes conservation of N_occ, M, Q, and the squeezed pin count F on the full finite-K_pen space. Ordered material cell lengths are conserved up to cyclic rotation. The zero-mismatch projection P is not invariant. The F=0 finite-K_pen material restriction equals C85's material restriction, but their finite-K_pen mobile Hamiltonians differ. C83 and C85 establish the corresponding coloring statements for their automata.

[Proved, coloring transfer] Suppose a material letter a has one outgoing automaton state o(a), an incoming set I(a), and a vacancy V_q copies q. Between consecutive material letters a and c with g>0 vacancies, the first vacancy must have q=o(a). Each subsequent vacancy has the same q. The final bond is legal if and only if o(a) belongs to I(c). Thus, for N_occ>0, a legal squeezed word and its occupancy pattern determine one legal coloring. This is a basis bijection and hence an isometry U_col. The rule is independent of the choice of root, so U_col intertwines physical rotation.

In C88, o=0 on L and markers, o=1 on R,D,X, I(R_0)=I(E)={0,1}, I(R_1)={1}, and I={0} otherwise. Each illegal squeezed junction requires a physical mismatch in its intervening interval. Distinct squeezed junctions have disjoint intervals, proving N_matbad<=N_bad. The two constant all-vacancy words are an exception to uniqueness. They belong to the competing sectors, not to a second local color fluid.

[Proved, projected operator] Let g_j be the number of vacancies after material letter j, A_j=1_(g_j=0), and E_j^-|g>=|g-1> for g>0, E_j^-|0>=0, E_j^+=(E_j^-)^*. A swap of a material letter and a neighboring vacancy transfers one gap unit with amplitude -t. Unique coloring makes the allowed forward and reverse amplitudes identical. A material bond acts only at contact. Consequently the legal compression, with the actual root gluing retained, is

    U_col^* P H P U_col
      = -t sum_j(E_j^+ E_(j+1)^- + h.c.)
        + J sum_j A_j h_mat,j - mu N_occ.                 (1)

Here h_mat is itself compressed to the legal material words. Material half-site energies in C88's literal bond are included in h_mat and multiplied by A_j. They cannot be replaced by independent physical onsite energies. Equation (1) does not say that the finite-K_pen chain equals its legal compression.

[Proved] The gap shifts have norm one. The gap Hilbert space is infinite dimensional because g_j is unbounded, but the shifts and contacts are bounded. The gap occupation is not a binary length bit. There are n payload modes per cell and one collective gap fluid, not n+b independent modes. Fixing sum_j g_j removes a zero mode, not the nonzero-momentum oscillator of that fluid. Color numbers are not independently conserved in C88.

## 3. The charge reference and K_gap=4

[Imported, proved reference identification] The gap Hamiltonian in (1), without the material interaction, is the squeezed representation of hard-core particle hopping on the physical chain. On an open chain, ordering the occupied positions gives a bijection with gaps plus boundary data. The hopping matrix elements are all -t. Jordan-Wigner fermionization of the physical occupancy variables gives dispersion -2t cos k. On rings, the material rotation eigenvalue and occupancy parity give boundary twists. Those twists affect the finite-size sewing, not the bulk energy density. This is the charge reference used in C81 PROOF section 5 and the d=24 vacancy notes.

[Proved, reference thermodynamics] Filling the free sea to k_F=pi rho gives

    e_p(rho) = integral_(-pi rho)^(pi rho) (-2t cos k) dk/(2pi)
             = -(2t/pi) sin(pi rho),
    v_p = 2t sin(pi rho),
    e_p''(rho) = 2pi t sin(pi rho)>0, 0<rho<1.             (2)

Put r=(L-N_occ)/N_occ=rho^-1-1. Energy per occupied letter and velocity in occupied-index coordinates are

    e_gap(r)=e_p(rho)/rho,
    v_gap=rho v_p,
    e_gap''(r)=rho^3 e_p''(rho).                          (3)

For the derivative identity, d rho/dr=-rho^2, so e_gap'=e_p-rho e_p' and differentiating once more gives rho^3 e_p''. The velocity relation follows from j=rho x at the homogeneous reference.

[Assumed, reference field dictionary] Use the Gaussian convention delta r=-partial_j phi_gap/pi and dimension(exp(2i phi_gap))=K_gap. In this convention the static energy density is (v_gap/(2pi K_gap))(partial_j phi_gap)^2, hence e_gap''=pi v_gap/K_gap. Identifying the full microscopic contact and vertex fields with this compact Gaussian dictionary is part of the vacancy initialization input below. The thermodynamic computation alone does not construct those fields.

[Proved conditional on that dictionary] Equations (2)-(3) imply

    K_gap=rho^-2,             K_particle=1,
    rho=1/2: K_gap=4, v_gap=t, v_particle=2t.             (4)

This is a change of coordinates and compactification. It is not a perturbation of size three away from a K=1 fermion reference in the same gap variables. A weak-fermion theorem for C58's binary length variable cannot be applied to E^+/- by identifying these operators with canonical boson or fermion creation operators.

[Proved, contact normalization] At separation one, the free sea has <c_j^*c_(j+1)>=sin(pi rho)/pi. The fermion contraction gives <n_j n_(j+1)>=rho^2-[sin(pi rho)/pi]^2. There are rho occupied letters per physical site. The contact fraction per occupied letter is therefore

    alpha(rho)=rho-sin^2(pi rho)/(pi^2 rho),
    alpha(1/2)=1/2-2/pi^2>0.                             (5)

Replacing A_j by alpha is the first-order reference projection. It is not an exact identity at fixed nonzero J. The connected fluctuation A_j-alpha couples to the whole material operator, including the clock and auxiliary states.

[Proved, coordinate check] In cell coordinates the charge speed is v_charge,cell=v_gap/b. Multiplying by ell=b/rho gives u_charge=v_particle, so u_charge=2t at half occupancy. Payload speeds in the reference are v_pay,cell=J alpha v_mat,cell and u_pay=(b/rho)J alpha v_mat,cell. Every speed in the final theorem is dressed; these equalities specify the reference only. For C88 the cell-to-physical length factor is 30. For C83 it is 28, for C85 30, and for d=24 it is 48.

## 4. Finite penalty and what its elimination proves

[Imported, proved] C88's local integer-penalty construction yields, with kappa=J K_pen and unpenalized local strength nu,

    U H U^* = kappa N_bad + D_2 + R_2,
    [D_2,N_bad]=0,
    ||R_2||_a <= C_3 nu^3/kappa^2.                        (6)

Here a can be a weaker fixed exponential support weight. The retained legal block contains the complete correction

    -[t^2 C + J^2 A + 2J^2 lambda^2 Z]/(J K_pen),         (7)

with the counters defined in C88 PROOF section 4. In particular C counts R_0 or E followed by a vacancy. Its payload-dependent field has scale t^2/(J K_pen), without a lambda^2 factor. C85 has a different color counter. C83 also has a color-induced bit field. The original d=24 lift has no two-color ambiguity term, but still has finite-coupling vacancy matching to perform.

[Proved, sufficient budget] A bound on the local norm of (7) is A2/(J K_pen), where A2=2t^2+3J^2+6J^2 lambda^2. To place (7) and R_2 below J tau/2 each, it is enough to require

    K_pen > 2 A2/(J^2 tau),
    K_pen > [2 C_3 nu^3/(J^3 tau)]^(1/2),                (8)

as well as the convergence condition for (6) and the material lower bound on K_pen. These are direct rearrangements of the two norm inequalities. At fixed positive J,t,tau they admit finite K_pen. For a band-sensitive application tau must be measured against the actual small payload stiffness and momentum tolerance, not against one material unit.

[Proved, parameter-order limitation] With t fixed, the color correction relative to J is of order t^2/(J^2 K_pen). Thus sending J/t to zero at a fixed K_pen eventually worsens this error. The phrase "take J small last" is not an adequate colored-chain schedule. One may choose J first and then K_pen for the algebraic bound (8). If a separate matching theorem has constants growing with K_pen, the joint feasibility must also be proved. For example, upper bounds requiring J/t<const/K_pen^2 would be incompatible at large K_pen with J/t>const'/sqrt(K_pen). We do not assert that the actual optimal constants have this behavior. We identify why their dependence cannot be ignored.

[Proved, remainder limitation] The norm bound (6) is a useful exact conjugation with a retained remainder. Discarding R_2 gives at best a free-energy error of order ||R_2||_0, independent of T. At fixed couplings it is not o(T^2). Nor does (6) give a volume-uniform global projection onto N_bad=0. A C48-type matching theorem must integrate this remainder with its scalar and sources, or an equivalent exact matching must retain it. The statement holds for every fixed finite order of penalty elimination.

## 5. C48 carries over as a theorem, not as an already verified application

[Imported, proved abstract theorem] C48 assumes a product cell space with a finite retained Fock code P, an onsite A satisfying

    A P=0, A>=Delta Q,
    Tr_Q exp(-s A)<=C0 exp(-s Delta), s>=s0,

bounded locally code-reducing B, a remainder R, and a local graded representation preserving the retained reference. It integrates the complete auxiliary operator rather than expanding its internal interactions as Grassmann vertices. The theorem is independent of d0=4464.

[Proved, finite-cell part] If a selected material cell of dimension D_cell has a separated 2^n-dimensional low subspace with gap Delta, its excited heat trace is bounded by (D_cell-2^n)exp(-s Delta). Small finite-K_pen splittings within the low subspace can be retained in B rather than falsely declaring an exactly degenerate code. Given the separated subspace, diagonalize its finite Hamiltonian, subtract its common reference energy, and define A to vanish on that entire subspace and equal the positive excited block on its complement. The difference on the code is a bounded retained term. Thus no countable heat-trace difficulty occurs for the selected finite material cell itself. This does not prove that it stays gapped in a mobile bath.

[Proved conditional transfer] Suppose a mobile preparation actually produces the preceding C48 hypotheses, uniformly on a source strip, with m finite logical modes per prepared block and local remainder r=||R||_a. The matrix-unit conversion on r_X cells costs at most (8 H^2)^(m r_X), independent of auxiliary dimension. Choose the input support exponent to pay this factor before fixing couplings. Define

    K_D,X = J_X B_X J_X^* - E_ref(P_X B_X P_X),

where E_ref is the isometric reference extension in that representation. Since B_X reduces P_X and E_ref agrees with it on auxiliary vacuum, each auxiliary coefficient of K_D,X annihilates that vacuum from both sides. Every vacuum word made only of K_D is zero. Analyticity of C48's convergent polymer map then gives

    ||W_infinity|| <= C_match r,
    ||W_beta-Per_beta W_infinity|| <= C_th exp(-c beta Delta). (9)

The constants include the chosen field and moment weights. The proof of the zero-temperature cancellation does not depend on the color number or K_gap. Fixed derivatives of bounded sources follow by Cauchy's formula on a smaller analytic ball. Scalar coefficients and the actual thermal endpoint are retained. The invalid embedding states have an activated pressure bound if their gap exceeds the incident interaction strength, as in C48 section 7. These invalid states are distinct from physical grammar, length, or pin defects.

For an explicit unchanged C48 prescription at field weight four, choose an output spatial rate alpha and put L_w=8(alpha+log(32)+10), Theta>=max(2L_w,Delta s0,log(C0)+2), and a=L_w+m log(128)+1. Its auxiliary coefficient strengths k,r must satisfy k+r<=Delta/(64 Theta). Then beta>=2Theta/Delta, temporal decay rate Delta L_w/(16Theta), leakage constant 2^22 exp(L_w), and a thermal bound of order [2^16 exp(L_w)Delta/Theta]exp[-beta Delta/(4Theta)] are available, with the moment losses of C48 section 4. These constants use the newly verified cell data, not the wrapper's old d0-dependent substitution. If the charge field remains a bosonic gap system rather than a finite Fock code, this particular prescription requires an extension and is not an application already covered by C48.

[Proved conditional source extension] If an auxiliary source q commutes strongly with A and |q|<=c_A A, then on |Re s|<1/(2c_A), Re(A+s q)>=A/2. In their common spectral representation,

    ||Q exp[-u(A+s q)]||<=exp(-u Delta/2),
    Tr_Q |exp[-u(A+s q)]|<=Tr_Q exp(-u A/2).

These are the contraction, passive-wire and trace estimates used by C48. They prove the same polymer bounds and holomorphy uniformly on an inner strip. Finite spectral cutoffs converge coefficientwise by strong convergence between bounded insertions and trace-class semigroups. Their absolute polymer sum is dominated by the same convergent majorant, so dominated convergence removes the cutoff. Cauchy's formula on a circle strictly inside the strip gives ||partial_s^k W||<=k! r_s^-k sup||W||. This supplies the routine tilted-matching argument when its relative bound is available.

[Proved, three failures of a literal substitution] First, the on-cell physical length is b+sum g_j, unbounded on the retained gap fluid. It cannot obey |q|<=c_A A with A zero on the retained space. The gap fluid must remain a gapless retained system or be replaced by an exact representation that keeps its full trace. A positive length tilt makes its trace finite, but does not turn it into a C48 massive auxiliary.

Second, the microscopic material excitation cost in (1) is contact-dependent. A gap configuration with every g_j>=1 has A_j=0 for all j, so the entire material term vanishes there. It admits arbitrary material excitations with zero contact energy. Therefore there is no pointwise operator inequality J sum_j A_j h_exc,j>=Delta sum_j Q_exc,j with Delta>0 on the retained gap space. A dynamically dressed threshold may exist, but replacing A_j by its positive mean in such an inequality is invalid. Likewise t is not small compared with the intended material gap of order J. C48's hypothesis verification needs a finite-J preparation, not just a finite-cell heat trace.

Third, for C88/C83/C85 the code charge is Q=2N_payload. The operator exp(i pi Q/2) gives logical parity on the code but has eigenvalues +i and -i on odd-Q auxiliaries. It is not a Z2 grading. The actual parity exp(i pi Q) is trivial on the code. The audited SEPARATE_JW theorem requires suitable integer parity charges on the entire cell, including auxiliary exchange. It cannot be invoked with Q/2 as an integer charge on this space. A different representation or a preceding elimination of half-charge auxiliaries is needed. This is an unverified interface, not a no-go theorem for such a representation. The original d=24 frozen N23_MATCHING application has Q=N_payload on its code and verifies this particular grading interface; its mobile preparation remains separate.

[Proved, physical locality qualification] A fixed interval of squeezed cells can contain arbitrarily many physical vacancies. Hence exponential locality in the squeezed cell index does not by itself imply exponential locality in physical distance. C49/C58's support bound used code lengths b or b+1, which are bounded before the physical lift; that argument does not apply to a code retaining arbitrary gaps. Pressure arguments below need bounded seams, which the literal Hamiltonian supplies. A physical current or charged-field identification additionally needs the correctly transformed sources and control of this fluctuating coordinate map. This requirement is included in I and G_mob, not inferred from a cell-index norm.

[Assumed, I, mobile initialization] For the chosen fixed finite parameters, construct an exact, source-uniform matching of the selected mobile trace and its physical observables to n payload channels plus one compact charge field. Require all of the following.

1. A justified elimination of the massive material and grammar excitations, or a hybrid preparation to which C48 applies, with the correct full-space grading or a replacement proof. The positive massive threshold is dressed and is not inferred from the frozen gap alone.
2. The K_gap=4 reference dictionary at rho=1/2, its compact charge lattice and physical harmonic sewing, with the whole gapless charge field retained. Low-frequency charge-payload matching at fixed J/t must control the full connected action, not just its two-point static kernel.
3. All-degree endpoint bounds in the strong spatial/time moment norm needed by N, including its degree hierarchy or an equivalent proved summability condition. The normalized error must be inside a stated positive RG tolerance relative to the smallest payload stiffness and the nonresonance margins. Finite-g, port, color-field and penalty-remainder terms are included. A finite matrix entry certificate does not establish this coefficient-norm bound.
4. The exact thermal scalar normalization, actual source representatives, and enough uniform source dependence to compare ensembles. If an enlarged trace is used, its comparison with the selected physical trace is bounded by C T exp(-delta/T), uniformly on the source interval, or by another proved o(T^2) bound with the same zero-temperature limit.
5. A consistent choice of fixed J,K_pen and other couplings satisfying both (8) and the initialization estimates. Alternatively, formulate and prove these estimates uniformly in the strong-penalty limit before choosing those finite parameters.

These are concrete outstanding initialization requirements in C81, C83, C85, C88 and the supplied d=24 vacancy analysis. C48 proves the implication (9) once its hypotheses are met; it does not prove I for these mobile chains. The supplied gap-free ultraviolet normal form controls a high-frequency remainder in a first spatial moment norm. Its retained low-frequency operator is still operator-valued, so it does not fill item 2.

## 6. C58's ensemble arguments that do transfer

[Proved, fixed-density pinning] Let p_i project onto a material letter, or onto any of the marker letters, on physical site i. On a sector with conserved R=sum_i p_i>=1, define

    V psi=R^-1/2 direct_sum_i p_i psi,
    K=direct_sum_i (p_i H p_i)|ran(p_i).

Then V^*V=1 and

    V^* K V=H-(1/(2R))sum_i[p_i,[p_i,H]].                (10)

For ||h||<=B only the two bonds incident on i contribute. The estimate ||[p,[p,h]]||<=4||h|| gives ||V^*KV-H||<=4BL/R. Scalar convexity in a basis diagonalizing the compression gives Tr exp(-beta V^*KV)<=Tr exp(-beta K). Translation invariance makes the latter L Z_pin. Applying the same scalar convexity to p_i H p_i and summing over i gives R Z>=L Z_pin. Thus

    (L/R) Z_pin <= Z <= L exp(4 beta B L/R) Z_pin.       (11)

This proof uses neither the wrapper's alphabet nor positivity of off-diagonal matrix entries. It applies to C88 in the fixed N_occ and M sectors, including its invariant good-sector restriction. At fixed K_pen, B is finite, even if very large. At R/L bounded below, the logarithm of the comparison divided by L tends to zero at every fixed T.

[Proved, seam] Cutting the physical bonds meeting a pinned root removes a bounded number of terms of norm at most B. Rooted rewriting is a basis isometry on its pinned space. Restoring the root move is a finite sum of partial isometries with the same literal amplitudes. Hence its norm is bounded by a finite B_seam depending on the fixed bond and alphabet, but not on L or the intervening gap lengths. For (1) the root gap transfer is a product of norm-one shifts. An exponentially summable retained interaction has the same conclusion when sum_(X crossing seam)||Phi_X|| is finite. Min-max gives a factor exp(+-beta B_seam) in the trace. Together with (11) this removes the root and cyclic seam in pressure. It does not identify independently periodic payload and charge fermion traces, and it gives no o(1/L) spectral error.

[Proved, trace-class gap tilt] On a selected block of b material sites there are b nonnegative gaps. With positive length tilt s,

    Tr_gap exp[-beta s sum_(j=1)^b g_j]
       = (1-exp(-beta s))^-b.                            (12)

Multiplying by exp(-beta s b) and the finite material dimension gives the free rooted trace. Bounded contact and shift interactions preserve trace-classness by the min-max trace bound. At fixed s in a compact subset of (0,infinity) the estimate is uniform. This is simpler than C58's all-string length heat bound, but (12) retains the gap entropy instead of eliminating it. For a fully enlarged model, the corresponding tilted auxiliary estimate must also hold, as required in I.

[Proved, pressure operations] Let p(T,s) be pressure per selected cell with tilt +s L_M and retained occupied source -mu b M. Positive lengths imply

    p(T,s+a)<=p(T,s)-b a, a>0.                           (13)

The lower bound b can be replaced by one when other positive cell lengths are included. With finite cut strength and the trace-class interval, the C58 blocking proof applies. To make the canonical step explicit, cut into q blocks of k cells. The removed interactions have norm at most C_cut(q+1). The fixed-total-length trace of the cut system is a convolution of positive block traces w_l. The upper coefficient bound follows from any tilt. For a lower bound keep finitely many lengths and prescribed rational frequencies of them; their multinomial weight gives the entropy term -sum r_l log r_l. Finite distributions approximate the tilted block trace and its first length moment because the tilt interval gives exponential tails. Each length bk+g, g>=0, is available here, so integer adjustments produce every compatible total without an extensive cost. Strict finite-difference slope brackets confine the tilt. Taking q to infinity and then k to infinity removes C_cut/k and proves

    f_cell(T,ell)=sup_s[-p(T,s)-s ell].                  (14)

This is a fixed-density statement for ell in the interior bracket. It does not extract coefficients from a bound on a generating function.

Summing over cell counts at fixed physical length gives the pressure-zero version. Positive cell lengths make subspaces with different cut prefixes orthogonal. Compress to those subspaces, delete the bounded crossing terms and apply scalar trace convexity. The resulting fixed-length sums are supermultiplicative up to exp(-beta C_cut). Blocking proves their exponential growth rate on compatible lengths. The positive identity sum_L exp(-beta s L)sum_M Z_(M,L)=sum_M Z_M(s), and (13), identify that rate with the unique interior zero of p(T,s). For the selected good cells, every L>=b is available with some M. Sectors with M=0 are handled by the defect hypothesis, not by this renewal identity.

[Proved, stability under pressure errors] If |p_1-p_2|<=d_T uniformly where both zeros lie, (13) gives |s_1-s_2|<=d_T/b by evaluating one pressure at the other's zero. If both canonical maximizers lie in that interval, (14) gives |f_cell,1-f_cell,2|<=d_T. Multiplication by M/L=1/ell gives the physical canonical bound. In particular d_T=C T exp(-delta/T) is o(T^2). These are the routine C58 transfers for every row of the table in section 1.

[Proved, source interval qualification] For the free charge reference at mu=0, the physical pressure zero at T=0 and rho=1/2 is s0=2t/pi>0. Indeed e_gap(1)=-4t/pi, e_gap'(1)=-2t/pi, and minimizing e_gap(r)+s(1+r) gives r=1 at s=s0 and zero minimum. For b gaps per cell the same zero is multiplied by b, not moved. The positive curvature (2) gives strict density brackets around half occupancy. A small uniform norm/pressure error preserves a finite-difference bracket if 2 error/(source spacing) is smaller than the reference margin. This is sufficient to place canonical saddles, but a differentiable dressed response and precise density tuning still belong to the source part of N and I. No derivative bound on arbitrary observables follows just from a small pressure error.

## 7. The missing physical defect comparison

[Proved, nontransfer] C58's completion removal uses invariant per-cell mismatch flags and a product onsite bound on bad cells and their successors. C88's reverse colored hopping and some material moves change N_bad. Thus those flags are absent, and the completed-wrapper proof cannot be copied. The exact pin count F is a different conserved quantity. Its conservation alone gives no energy or trace advantage to F=0. C88 proves frozen pin selection and assumes mobile pin selection S_pin.

[Assumed, D_E, ground comparison] At fixed compatible N_occ=bM and interior vacancy density, require a common-bulk-subtracted dressed inequality of the form

    E_sigma-E_good >= delta_def D(sigma)-C_boundary,
    delta_def>0,                                       (15)

with positive finite-defect costs wherever the boundary term leaves them undecided. D weights long wrong cells at least by their excess material length, and includes pin defects for C88, wrong selected lengths, and the relevant interfaces. Marker-free, all-vacancy, and density competitors are separately excluded at the chosen source. This is the extension of C88's S_pin needed for ground-state selection. A global O(J^2 L) error does not imply (15).

[Assumed, D_T, thermal comparison] For the thermal conclusion require, uniformly in the source window and at fixed small T after the thermodynamic limit,

    0 <= p_full(T,s)-p_good(T,s)
         <= C T exp(-delta/T).                          (16)

Here p_good uses the exact selected sector with its virtual grammar excitations and all scalar terms, not P H P. Alternatively require a finite-volume trace comparison with subextensive boundary factors that implies (16). Virtual defects that are not invariant sectors must already be included by exact matching or a dressed trace estimate. A ground-energy inequality such as (15) by itself does not bound the excited spectra or their entropy and is not (16).

[Proved, sufficient trace version of D_T] One usable stronger hypothesis is a decomposition or trace domination in which defects occupy disjoint marked intervals, the cost of an interval of weight l is at most A^l exp(-beta delta_0 l) relative to the common good bulk, and every defect has l>=1. Include internal multiplicities and allowed wrong lengths in A^l. Then the activity summed over lengths is

    z_beta=sum_(l>=1) A^l exp(-beta delta_0 l)
          =A exp(-beta delta_0)/(1-A exp(-beta delta_0)).  (17)

Ignoring disjointness only increases the sum. At most V possible anchor sites give Z_full/Z_good<=exp(o(V)/T)(1+z_beta)^V. Taking T/V times the logarithm yields at most T log(1+z_beta)=O(T exp(-delta_0/T)). A fixed conversion between V and cell or physical length only changes C. This proves the activated transfer from an all-spectrum defect bound with a length-weighted entropy cost. It also explains what must be added if "dressed-defect comparison" is intended to cover C58's thermodynamic use. No such comparison is established for C88 here.

## 8. The multichannel RG hypotheses change

[Assumed, harmonic dictionary from I] In physical particle normalization let N(x)=rho x-phi_p(x)/pi and dimension(exp(2ia phi_p))=a^2 at the reference. The selected cell coordinate is y=N(x)/b. For a payload operator of cell momentum q and dimension Delta_pay, the cell comb includes every integer harmonic Q_h=q+2pi l.

[Proved within the dictionary] Substitution of y(x) gives the branch momentum and dimension

    k_phys=rho Q_h/b mod 2pi,
    Delta_ref=Delta_pay+[Q_h/(2pi b)]^2.                 (18)

At rho=1/2 a uniform term satisfies Q_h=4pi b h, h integer. For h!=0 its charge dimension is 4h^2>=4. For h=0, the certified internal pair-sum margins exclude nonforward neutral quadratic and quartic resonances; the remaining forward quartics are retained in the Gaussian current form. Higher neutral fermion vertices have reference dimension at least three. Exact charges exclude operators changing net payload or occupied number in this retained operator lattice. This last phrase matters: the microscopic auxiliary grading problem is not solved by a low-energy charge selection rule.

For dimensions below three, |Q_h|/(2pi b)<sqrt(3). The distance to a nonzero physical reciprocal vector is at least pi(2-sqrt(3)); distance from a nonforward internal sum at the zero reciprocal vector is at least delta_internal/(2b). Thus

    delta_phys>=min(delta_internal/(2b),pi(2-sqrt(3))).  (19)

For C88 this is greater than 1.0837704e-7. The first marker harmonic has q=0, l=1, wavevector pi/15 and dimension 1/225. Its small dimension is harmless to a uniform Hamiltonian because it oscillates. Discarding it from the observable algebra would nevertheless be incorrect.

[Imported, proved Gaussian bound] In a relative neighborhood ||G_ref^-1/2(G-G_ref)G_ref^-1/2||<=eta<1/4, the supplied Gaussian metric comparison gives Delta_G>=sqrt((1-eta)/(1+eta))Delta_ref. In particular reference dimension three remains above 3sqrt(3/5)>2. This is a conditional stability margin, not a proof that the microscopic flow stays in that neighborhood. Constants may depend on all fixed velocity ratios.

[Assumed, N_mob] The needed multichannel version of corrected Lemma N is an all-scale theorem for the initialized mixed charge/payload action of I. It must provide normalized tree convergence, comparison with the full interacting Gaussian/current reference, continuous counterterms moving all relevant Fermi data, and source estimates on a common source window. It must preserve positive stiffness, the harmonic lattice and the margins (19), and give a differentiable ground pressure with an exposed interior density. The charge reference is K_gap=4, or its exactly equivalent physical K_particle=1 representation with the nonlocal squeezed operators and gluing accounted for. It is not the weak off-half-filled single-channel action in C55 or C59 without further work.

The source version must subtract the free lattice-to-reference descendant separately, as in C61 (N3a): only D_O(u)-D_O(0) may receive an interaction-small prefactor. C60's nonzero free pair tail rules out the unsubtracted version. C79 proves the source-free general-tree coefficient-moment estimate in its scheme; it leaves N-comparison and N-source open. Neither its proof nor the supplied one-step winding estimates prove N_mob for the phase-model charge sector.

[Assumed, G_mob, Gaussian completeness] For a Virasoro central-charge conclusion also require C51's G: all n+1 compact current sectors with their positive levels, their quadratic Hamiltonian and charged vertices are identified with the physical infrared theory, and exhaust it. A complete source theorem could include this as an explicit conclusion of N_mob. If N_mob means beta functions and selected current correlations only, G_mob is a separate input. Correct physical zero-mode gluing is included; finite-volume tower multiplicities would require the stronger spectral contract F of C51.

## 9. TP and the scalar summation transfer

[Assumed, TP_mob] After subtraction of the complete interacting reference at the dressed velocities and both ground energies, the scalar difference of the exact matched thermal action has an exact shell sum, including normalization determinants and the actual thermal histories. Uniformly on the source interval, with q_h=E_* 2^-h, require

    |delta_h(T,s)|<=C q_h^2 eta_h min(1,(T/q_h)^p),
    p>2, sup_h eta_h<infinity, eta_h ->0,
    |r_UV(T,s)|<=C E_*^2(T/E_*)^p.                      (20)

Constants carry the fixed length and energy units needed for pressure per cell. Convergence of eta_h is uniform in s. Massive thermal remainders from (9) are retained and separately bounded. Comparing only with the quadratic fermion determinant is insufficient when forward interactions survive. A trace-normalized multichannel compact reference, including its K_gap=4 charge component, is needed. This is the mobile version of C59/C62's TP, not a consequence of a two-time-moment endpoint bound.

[Proved, summation] Choose H with q_H<=T<2q_H. For h=H+j with j>=0, divide (20) by T^2 to bound the summand by C 2^-2j eta_(H+j). For j<0 the bound is C 2^((p-2)j) eta_(H+j). Both weight sequences are summable since p>2. On any fixed finite j-window eta_(H+j) tends uniformly to zero; outside it the uniformly bounded eta_h times the geometric tails is arbitrarily small. The UV contribution divided by T^2 tends to zero as well. Therefore the whole difference is o(T^2), uniformly in s. If eta_h=O(q_h^omega), geometric summation gives power min(omega,p-2), with one logarithm at equality. This is C51's scalar argument with no change from length to vacancy variables.

[Proved, Gaussian coefficient] One positive chiral oscillator of speed v has free-energy density

    T integral_0^infinity dk/(2pi) log(1-exp(-v k/T))
      = -pi T^2/(12v).

To evaluate the integral, expand -log(1-exp(-x))=sum_(m>=1) exp(-m x)/m and integrate to obtain sum m^-2=pi^2/6. Sum over the 2(n+1) chiral oscillators. Zero-mode sums on a fixed full-rank compact lattice contribute O(log L) to log Z at fixed T: positive definiteness bounds the lattice Gaussian above by C(1+LT)^(n+1), while a closest lattice point gives a lower bound exp(-C/(LT)). Shifts in a compact fundamental domain have uniform constants. Restricting finitely many conserved zero modes or finite-index cosets changes this polynomial bound, not the oscillator pressure. This is the Gaussian statement; the microscopic sewing identification remains in I and G_mob.

Thus N_mob's velocity and source output, TP_mob, and I give a matched cell pressure

    p_matched(T,s)=p0(s)+A(s)T^2+o(T^2),
    A(s)=(pi/12)sum_(a=1)^(n+1)
                     [1/v_(a,R)(s)+1/v_(a,L)(s)].       (21)

## 10. Conditional theorem for C88 and its proof

[Assumed, hypotheses of the theorem] Fix one literal C88 bond with n=13, b=15, the stated certified material data and additional finite J,t,K_pen,mu. Require I of section 5, D_E and D_T of section 7, N_mob of section 8, and TP_mob of section 9 on one source interval with interior saddles. Require G_mob for the central charge assertion. All estimates are at these same fixed couplings. This statement neither proves the existence of a jointly admissible parameter point nor assumes that the displayed C88 material numbers already meet every uncomputed strong-norm threshold.

[Proved conditional theorem, thermal part] Let s_* be the unique zero of p0 and ell_*=-p0'(s_*)>0. I and D_T transfer (21) to the full physical pressure with the same o(T^2) remainder. By (13), the finite-temperature zero s_T obeys s_T-s_*=O(T^2). Differentiability of p0 and continuity of A give

    0=-ell_*(s_T-s_*)+A(s_*)T^2+o(T^2),
    s_T=s_*+[A(s_*)/ell_*]T^2+o(T^2).                    (22)

Our Hamiltonian includes the chosen -mu N_occ, and no bulk centering has been subtracted in this notation. Its physical free energy is -s_T. If a centering e L is used, restore e throughout. Hence, with dressed physical speeds u_(a,R/L)=ell_* v_(a,R/L),

    f_phys(T)=f_phys(0)
      -(pi/12)sum_(a=1)^14[1/u_(a,R)+1/u_(a,L)]T^2
      +o(T^2).                                         (23)

At the tuned half-occupancy point ell_*=30. At fixed canonical rho=1/2 and N_occ=15M, (14) instead gives

    f_phys(T,rho)=sup_s[-p0(s)/ell-s]
                     -[A(s_rho)/ell]T^2+o(T^2),
    ell=15/rho, -p0'(s_rho)=ell.                        (24)

For completeness, this optimization needs no derivative of the thermal remainder. Strict exposure implies any optimizing s_T tends to s_rho. Evaluate the perturbed supremum at s_rho for a lower bound. For the upper bound, its zero-temperature objective at s_T is no greater than the maximum, while A(s_T)=A(s_rho)+o(1). Uniformity of the remainder proves (24). Fixing N_occ imposes one global charge constraint but does not delete the charge contribution to (23).

[Proved conditional theorem, Gaussian central charge] Under G_mob the positive quadratic form V of the 28 currents can be diagonalized without any equal-speed assumption. Put S=diag(1_14,-1_14), W=V^1/2 S V^1/2, choose O orthogonal with O^T W O=D, and set R=V^-1/2 O |D|^1/2. Inertia gives 14 positive and 14 negative eigenvalues of D. Direct multiplication gives R^T V R=|D| and R^T S R=S. The resulting right and left current algebras have positive levels normalized to one and positive propagation speeds.

For each normalized chiral current, define L_m=(1/2)sum_r :J_(m-r)J_r:. The current commutator gives [L_m,J_l]=-l J_(m+l). On a charge-x cyclic vector,

    ||L_(-m)|x>||^2=m x^2+(1/2)sum_(r=1)^(m-1)r(m-r)
                    =m x^2+(m^3-m)/12.

Subtracting 2m<L_0>=m x^2 leaves the Virasoro central term (m^3-m)/12. Different normalized currents commute, so their central terms add. Completeness identifies their sum with the full infrared theory. Therefore

    c_R=c_L=14, c_cd=(c_R+c_L)/2=14>12=d-1.             (25)

K_gap=4 changes vertex dimensions and the compactification lattice, not this central term. Neither vacancy color adds a current. Arbitrarily unequal nonzero speeds change (23), not (25). The conclusion uses the collection-of-conformal-sectors convention of C51; it does not assert one common Lorentz velocity.

[Proved, ground and spectral qualifications] D_E is what relates the selected phase to the full ground-state problem, including C88 pin sectors. For compatible particle counts, N_occ=15M. At a smoothly compressible half-density reference the nearest such count differs from L/2 by at most 15/2; Taylor expansion of L e_p(N_occ/L) about the stationary density gives O(t b^2/L) rounding energy, before the bounded-twist O(t/L) contribution. It vanishes at fixed positive defect gap, but it is not an all-ring spectral theorem. The bounded seams, scalar pressure argument and C48 local norm bounds do not prove o(1/L) matching, tower multiplicities, interval entropy asymptotics, or a unique finite-ring ground state. C51's extra spectral contract would be required for those claims.

## 11. Transfer ledger and checks

[Proved, conclusion of the audit] The ledger records theorem transfer separately from microscopic verification.

| Step | What transfers | What still needs verification for C88 |
| :-- | :-- | :-- |
| C48 resummation | Spectral auxiliary treatment, leakage-linear cancellation, all-degree kernels, massive winding and source analyticity | Mobile dressed gap, mixed retained charge representation, strong-norm initialization and full-space grading, I |
| Finite K_pen | Local integer-penalty conjugation and explicit correction, imported from C88 | Exact treatment of its remainder and a joint finite parameter schedule inside I |
| C58 rooting | Fixed-density pinning, bounded seam and cyclic-pressure equality | No new hypothesis for the finite literal bond at positive root density |
| C58 ensembles | Tilted trace, pressure zero, canonical Legendre transform and preservation of uniform o(T^2) errors | Interior source brackets and matched physical trace supplied by I, N_mob and D_T |
| C58 physical defects | Abstract activated comparison once defect traces are controlled | Wrapper flags do not transfer; D_E and D_T, including C88 pins and long cells |
| N | Moment estimates and conditional recurrence machinery after their norm hypotheses | Charge/payload comparison, counterterms, source theorem for the actual mobile action |
| TP | C51 shell summation and Gaussian oscillator coefficient | Uniform interacting-reference scalar comparison on actual thermal histories |
| C51 c | Positive-form diagonalization, Sugawara central terms and zero-mode qualifications | G_mob, completeness of the physical compact Gaussian limit |

[Proved, distinctions between dimensions] C83 and C85 use the same proof structure with their own b and color counters; C85's physical pin penalty differs from C88's removal and conserved F. The original d=24 construction has a valid full-space fermion grading in its frozen matching and no color-ambiguity self-energy. All three dimensions still need finite-J vacancy initialization, a mobile defect comparison, and the appropriate N/TP outputs. The d0=4464 wrapper has instead a binary length fluid with verified separate charges and an all-string onsite construction. Its completed application cannot supply those missing vacancy hypotheses by a change of numerical constants.

[Proved, three-label formulation] It is possible to summarize the conditional result using three labels by defining N_extended := I + N_mob + G_mob and D_extended := D_E + D_T, and using TP_mob as the third. This is an explicit enlargement of the open inputs. It is not a proof that the original narrower Lemmas N and TP plus C88's ground-energy S_pin imply (23)-(25). The expanded formulation states the intended research target without hiding the initialization, thermal or completeness obligations.

[Record] Sources read include the requested REPORT sections, C48/C58/C51/C88, referee_S and referee_imported, and the relevant C81/C83/C85, C61/C79, C59/C62 and d=24 matching/vacancy notes. Imported finite certificates were read and used, not rerun. The analysis uses no external literature search and makes no claim that the remaining inputs have been proved elsewhere. Local checks and their results are recorded in RESULT.md. Only work_C90 is modified.
