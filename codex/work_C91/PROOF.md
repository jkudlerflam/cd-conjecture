# C91. Finite-J vacancy initialization for the thirteen-state chain

[Proved, scope of this report] The requested finite-J initialization is not established. This report gives the exact legal coordinate representation, explicit volume-independent constants for its ultraviolet transformation, and three obstructions to using that transformation as the input to C48. It also gives a restricted grading obstruction and a quantitative constraint on the order of parameter choices. These are not a proof that the microscopic candidate fails. The missing statement concerns low-energy dynamics with the vacancy mode retained.

[Definition] A tag governs the following paragraphs and formulas until the next tag or section. Proved means derived here, or imported with explicit attribution and its original scope. Numerical means a finite floating-point control. Heuristic means a proposed field-theory interpretation. Assumed identifies an input that has not been proved. Definitions do not assert that the objects whose existence is requested have been constructed.

[Proved, provenance] C90 was read in full, including its check source and small check output. C88 and C48's proofs, handoffs and Python sources were read. All ten mathematical notes and the check source in the collaborator's vacancy folder, and the research notes and Python sources in counterexample_campaign, were read. Large parent data and build artifacts were not read. The campaign's initial endorsements are superseded by the withdrawn-counterexample status in its README and synthesis. Its finite matrices and continuum arguments do not supply the missing operator-valued matching estimate. No network access or parent-script execution was used.

## 1. The exact coordinate representation and its boundary

[Definition] Write N for occupied number, M for marker number, G for total vacancies, and K for the grammar penalty. In the selected material sector N=15M, F=0 and G=L-N. Write s=J/t. The notation K_gap is reserved for the charge stiffness. The material data are those in C88, with thirteen active bits, fourteen data sites, and material cell length fifteen. The literal bond is C88 equation (1).

[Proved, imported from C88 with the elementary identification repeated] At N>0, between two consecutive material letters a and b, every vacancy has the unique color o(a). The last pair is legal if and only if o(a) belongs to I(b). Thus a legal material word and its vacancy placements determine one legal physical word. This gives a basis isometry U_col onto the zero-mismatch space P. There are two exceptional constant-color all-vacancy words, outside this identification.

Choose a material root and set g_j equal to the vacancies following material letter j. The gap Hilbert factors are l2(N_0), with

    E^-|g> = |g-1> for g>0,    E^-|0> = 0,
    E^+ = (E^-)^*,            A = |0><0|.

Consequently

    E^- E^+ = I,    E^+ E^- = I-A,
    ||E^+|| = ||E^-|| = ||A|| = 1.                         (1)

The fixed-length space has sum_j g_j=G. The rooted tensor space before that constraint is useful for local estimates. It is not the whole physical ring by itself. With root position x, one must retain the identification

    (x; g_1,...,g_N; a_1,...,a_N)
      ~ (x+1+g_1; rotate(g); rotate(a)).                    (2)

The charge-only periodic matrix used in the numerical control below is the zero physical momentum representation. It is not an assertion of independently periodic charge and payload traces.

[Proved] A physical hop transfers one vacancy between adjacent gaps with amplitude -t. A material bond acts only at g_j=0. Hence the legal compression is

    U_col^* P H P U_col
      = t H0 + J V - mu N,
    H0 = -sum_j (E_j^+ E_(j+1)^- + h.c.),
    V = sum_j A_j h_j,                                    (3)

where h_j is the legally compressed unpenalized material bond. In particular the distributed material onsite energies are contact terms in (3). They are not independent of g_j.

This is locality in the material index. The map to g_j is a coordinate isometry, not a finite-depth local unitary on the original physical spins. Arbitrarily many physical vacancies can lie between two neighboring material letters. A locality statement in the material index must not be called a uniform physical-distance locality statement for sources.

[Proved, useful bounded extension] For local estimates on (3), implement each material move with the legality guards on its two exterior junctions and its internal junction. Its support is contained in [j-1,j+2], and it agrees with P h_j P on legal words. Use the product of the eleven-dimensional material-letter space and l2(N_0) as the enlarged local space. This extension leaves the legal sector invariant. A still smaller selected-length sector can be restricted afterwards. The guard factors have norm at most one. Summing the absolute norms of the literal move table and its diagonal terms gives the conservative bound

    ||v_j|| <= M0,   M0=2^13,   supp(v_j) subset [j-2,j+2], (4)

for D=1000, |g|,|t_d|,|lambda|<=1, 0<=mu_sel<=2 and 0<=f<=1. There are fewer than 64 unit-strength directed off-diagonal entries; the largest unpenalized diagonal is below 1010. The larger constant in (4) allows the guards and collection of terms. The chemical potential is scalar at fixed N and need not be included. The penalty JK is not included in M0.

[Proved, finite-K limitation] P is not invariant at finite K. For example, R_0 V_1 can hop backwards to either V_0 R_0 or V_1 R_0; the predecessor selects only one legal color. Therefore (3) is not an equality for the full finite-K Hamiltonian. Applying a theorem to (3) does not silently eliminate the other color or grammar states.

## 2. An explicit version of the cutoff-power-ten ultraviolet lemma

[Proved] The single-exponential method in the collaborator's SINGLE_EXPONENTIAL.md applies to (3). The following supplies one deliberately loose choice of its constants. All norms in this section are interaction norms for a specified decomposition, not extensive operator norms.

Define

    ||Phi||_(1,a) = sup_y sum_(X contains y)
                    (1+diam X) exp(a diam X) ||Phi_X||.

Take 0<delta<=1, M0>=1 as in (4), and put

    v=64 e,              T=delta^-2,
    R=ceil((2v+8)T),     C_R=2v+9,
    q=20 M0 R |s|/delta.                                 (5)

Then R<=C_R delta^-2. Work on intervals or rings with at least three gap sites; on a short ring replace a cutoff interval covering the ring by the entire ring. This can only reduce the counts below. Root seams can be kept as bounded boundary terms; the translation-covariant bulk statement is in material coordinates.

Let H0,I contain precisely the hopping bonds contained in a radius-R interval I_j. Define

    g_delta(u)=delta exp(-delta^2 u^2/4)/(2 sqrt(pi)),
    F_delta(u)=-(i/2) sign(u) erfc(delta |u|/2),
    S_j=integral_(-T)^T F_delta(u)
                         exp(iu H0,I_j) v_j exp(-iu H0,I_j) du,
    S=sum_j S_j,             U_delta(s)=exp(s S),
    w_j=integral_(-T)^T g_delta(u)
                         exp(iu H0,I_j) v_j exp(-iu H0,I_j) du.
                                                                  (6)

Each S_j is anti-Hermitian. Since integral_R |F_delta|=2/(sqrt(pi)delta),

    ||S_j|| <= 2 M0/delta,       ||w_j|| <= M0.            (7)

In every finite volume U_delta is an actual unitary, including for unbounded gap occupations. There are finitely many bounded summands. No occupation cutoff or gap-fluid spectral gap is needed.

[Proved, locality estimate used in the constants] For an observable O of norm at most M0 supported on at most five consecutive sites, and the sum H_boundary of the two hopping bonds crossing I_j,

    ||[exp(iu H0,I_j) O exp(-iu H0,I_j), H_boundary]||
       <= 128 M0 exp[-(R-4)+v |u|].                       (8)

Here every hopping bond has norm at most two. One elementary way to obtain (8) is the commutator integral inequality iterated over overlapping bonds. After the initial set, an edge has at most three overlapping edges; the initial set meets at most six. Each step contributes at most 2||h_edge||<=4 and a time-simplex factor |u|^n/n!. A path reaching a boundary must have at least R-4 steps. The count is bounded by 2*3^n, with the initial commutator bounded by 2||O||||h_boundary||. The tail of the exponential series is bounded using

    sum_(n>=d) (12|u|)^n/n! <= exp(-d+12 e |u|).

Summing the two boundaries and increasing the prefactor and velocity gives (8). The integral inequality can first be proved for finite local spaces by differentiating the commutator and conjugating away the terms internal to its support. All operators involved here are bounded, so the same integral equation and iteration apply directly on l2(N_0). This proof uses the norm bounds, not the local dimension.

[Proved, first-order identity] Distributionally F_delta'=-i delta_Dirac+i g_delta. Integration by parts gives, with H_boundary=H0-H0,I_j-H_out,

    v_j+[S_j,H0]=w_j+e_j,
    e_j=i F_delta(T) tau_T^I(v_j)
          -i F_delta(-T) tau_(-T)^I(v_j)
          +integral_(-T)^T F_delta(u)
                              [tau_u^I(v_j),H_boundary] du.       (9)

The sign follows from [tau_u^I(v_j),H0,I_j]=i d tau_u^I(v_j)/du. The error is supported in I_j enlarged by one site. The endpoints cost at most M0 exp[-1/(4delta^2)]. Equations (5), (7), (8) imply a boundary cost at most 256 M0 delta^-1 exp(-4/delta^2). In particular

    ||e_j|| <= 512 M0 delta^-1 exp[-1/(4delta^2)],
    ||sum_j e_j||_(1,0)
       <= 2^15 M0 R^2 delta^-1 exp[-1/(4delta^2)].        (10)

For the last estimate a site is in at most 2R+3 supporting intervals, each with 1+diam at most 2R+3. Their product is at most 49R^2. This accounts for incidence as well as range.

[Proved, nonlinear remainder] For a local original term B_j of norm B supported at radius at most two, an n-fold commutator with S has at most

    product_(k=1)^n (5 k R) = (5R)^n n!

potentially nonzero anchor sequences. Indeed the next support must intersect the preceding union, and its center is at distance at most 2kR from j. Each commutator costs at most 4M0/delta. Therefore

    ||ad_S^n(B_j)|| <= B n! (20 M0 R/delta)^n.            (11)

Assign this coefficient to the interval of radius 2+2nR about j. Its incidence times first-moment weight is at most 49n^2 R^2 for n>=1. The BCH factorial cancels n! in (11). If q<=1/2,

    sum_(n>=1) n^2 q^n <= 12q,
    sum_(n>=2) n^2 q^n <= 22q^2.

It follows that the nonlinear remainder after the terms H0+s(V+[S,H0]) is bounded by

    49R^2 [44q^2+12 |s| M0 q]
       <= 2^20 M0^2 C_R^4 s^2 delta^-10.                (12)

The first bracket uses ||h0,j||<=2. The second contribution is the transformed sV and retains its extra factor |s|. This verifies the power ten without an extensive-volume estimate.

[Proved, quantitative output] Define C_UV=2^20 M0^2 C_R^4. Combining (9)-(12) gives the exact finite-volume identity

    U_delta(s) (tH0+JV) U_delta(s)^*
        =tH0+J W_delta+R_UV,     W_delta=sum_j w_j,
    ||R_UV||_(1,0)
        <=C_UV t s^2 delta^-10
          +2^15 M0 C_R^2 J delta^-5 exp[-1/(4delta^2)].  (13)

All constants are independent of N,G and volume. For M0=8192, log_2 C_UV is approximately 79.9182. The hypothesis q<=1/2 is substantive; the numerical script reports its domain rather than treating moderate s as admissible.

An exponential weight is also available. Set

    kappa_x=(log 2)/(12R),    C_exp=2^28 M0^2 C_R^4.

The order-n support diameter is at most 6nR, so its extra weight replaces q by at most sqrt(2)q. The geometric sums at q<=1/2 are bounded by 256q and 256q^2, respectively. The error support has diameter at most 2R+2, whose additional weight is below two. Thus a sufficient bound is

    ||R_UV||_(1,kappa_x)
      <= C_exp [t s^2 delta^-10
                  +J delta^-5 exp[-1/(4delta^2)]].       (14)

With delta=s^z, 0<z<1/10, the first-moment remainder divided by J tends to zero, and q tends to zero. The rate kappa_x also tends to zero. This is a proved ultraviolet Hamiltonian statement. It is not a bound on a massive influence action.

## 3. Why the filtered contact is not the averaged material Hamiltonian

[Proved] For one contact term v_j=A_j h_j, H0 acts only on gaps, so the filtered operator is a gap operator times the same h_j. On the gap vacuum inside I_j, H0,I_j is zero and A_j is one. Thus

    w_j (|0,...,0>_I tensor u)
      =erf(1/(2delta)) |0,...,0>_I tensor h_j u.           (15)

Choose a sufficiently long canonical half-filled ring G=N. Put all its vacancies outside I_j, leaving the gap vacuum inside I_j. This vector belongs to the same fixed total-gap sector as the homogeneous sea. It proves, for any nonzero h_j on the chosen material sector,

    ||w_j-alpha_c h_j||
      >= |erf(1/(2delta))-alpha_c| ||h_j||,
    alpha_c=1/2-2/pi^2.                                  (16)

For 0<delta<=1/2 the lower coefficient exceeds 0.545. As delta tends to zero it tends to 1-alpha_c, approximately 0.702642. This is a canonical-sector operator-norm obstruction, not an artifact of admitting the zero-vacancy sector. Those test vectors have large energy above the sea. Accordingly (16) refutes a full-retained-space replacement, not a low-energy form or connected-kernel estimate.

The same distinction is visible in the ideal filter's energy basis. It multiplies (v_j)_(mn) by exp[-((E_m-E_n)/delta)^2]. It preserves every diagonal element, regardless of E_m. Restricting energy transfer cannot project onto small absolute energy. The strictly local filter differs from the ideal one only by the controlled truncation error on an anchored observable, but neither turns (15) into alpha_c h_j.

[Proved, no bare material gap] At half filling the basis vector with g_j=1 for every j has A_j=0 for every j. For any nonzero material excited vector u,

    <g,u|J sum_j A_j h_exc,j|g,u>=0.

There can therefore be no inequality

    J sum_j A_j h_exc,j >= Delta_* sum_c Q_c

on the entire retained gap space with Delta_*>0 and nonzero material excitation projectors Q_c. Adding an arbitrary occupation cutoff that keeps g=1 does not repair it. This does not refute a threshold measured above the interacting ground state of tH0+JV.

[Proved, leakage already at lambda=0] The frozen vacuum clock in a length-fifteen cell has components

    v_l=sqrt(2/16) sin(pi(l+1)/16),   l=0,...,14,
    h_clock=-sum_(l=0)^13 (|l><l+1|+h.c.).

It is an eigenvector of the sum, not of each edge. Let P_v=|v><v| and k_0=-(|0><1|+|1><0|). Then

    ||(1-P_v) k_0 v||^2
       =v_0^2+v_1^2-4v_0^2 v_1^2 >0.                   (17)

A gap configuration with only the corresponding internal material junction in contact realizes this single-edge term on the vacuum clock. At G=N take that gap to zero, a different gap to two, and every other gap to one. No marker contact is needed. The field is zero at lambda=0, and the g conversion has no vacuum-to-carrier matrix element. Hence the bare frozen-code leakage is present at order J, even at zero port. Numerically the coefficient in (17) is 0.1507149408.

This explains why the C48 strategy must first dress the zero-port mobile clock. Treating the contact fluctuation as a remainder beginning at lambda^3 is incorrect. Decreasing lambda does not reduce (17). A new mobile-code definition could remove this leakage, but no such construction has been verified here.

[Assumed, missing low-frequency theorem] A sufficient next theorem would retain the gapless fluid and construct dressed material degrees of freedom at lambda=0, with a positive threshold of order J, uniform in volume and in a density-source interval. It would also bound their operator-valued contact couplings in the coefficient norm after a justified representation change. A scalar static susceptibility, a finite-size charge projection, and (13) do not prove that statement. The collaborator's tagged-defect tests and continuum susceptibility calculation are explicitly weaker inputs.

## 4. The strong norm required by C48

[Proved, imported C48 statement] For thirteen retained fermion modes on r cells, C48's matrix-unit conversion at field weight four costs at most 128^(13r). If a graded local representation is available, its prescription for output spatial exponent alpha is

    L_w=8(alpha+log 32+10),
    a_in=L_w+13 log128+1,
    Theta>=max(2L_w,Delta_* t0,log C0+2),
    k+r <= Delta_*/(64 Theta).                            (18)

Here k and r are the auxiliary coefficient strengths in that representation, and the auxiliary heat trace satisfies Tr_Q exp(-u A)<=C0 exp(-u Delta_*) for u>=t0. The retained gap fluid is not part of this heat trace. Its onsite space is infinite and it is gapless.

[Proved, norm non-implication] A first spatial moment and an exponential tail at rate kappa_x do not imply the input norm in (18). For example, take one anchored interval term of length l with norm exp(-2 kappa_x l)/l^3, translated at every site. Its first-moment sum is finite. Its norm with weight exp(a l) is infinite for a>2 kappa_x. This is a counterexample to a general norm implication, not a claim that this particular tail equals the actual R_UV.

In cell coordinates one unit is fifteen material sites. Consequently (14) provides at most a spatial rate of order 15/R in those coordinates, while (18) asks for a fixed rate, already greater than 100 at alpha=0. Blocking at a length proportional to R does not establish the conversion either. It increases the number of retained payload modes and therefore the coefficient cost.

[Proved, explicit stronger Hamiltonian estimate and its cost] The same connected-commutator proof supplies an estimate at any prescribed material-site spatial rate a>=0 if

    q_a=q exp(6aR)<=1/2.

Indeed the declared diameter at order n is at most 6nR. Thus

    ||R_nonlinear||_(1,a)
       <=49R^2[44 q_a^2+12 |s| M0 q_a],
    ||J sum_j e_j||_(1,a)
       <=J 2^15 M0 R^2 delta^-1
                    exp[a(2R+2)-1/(4delta^2)].            (19)

This is an explicit bound, including when its right side is large. In particular the sufficient smallness condition is

    |s| <= delta exp(-6aR)/(40 M0 R).                    (20)

With delta=s^z and fixed positive a, the majorant q_a diverges as s tends to zero: log q_a contains a positive term proportional to s^(-2z), which dominates log s. Thus the polynomial cutoff schedule in the supplied power-ten lemma does not certify a fixed strong norm by this proof. It would be incorrect to report (14) as (18).

This conclusion is about these estimates. A different filter radius, a stronger initial Lieb-Robinson weight, or a different coefficient organization might give better strong Hamiltonian bounds. Moreover at fixed delta, (20) can be met by taking s sufficiently small. Neither observation controls the retained W_delta or supplies a dressed gap, so neither closes initialization.

[Proved, additional interface] C48 leaves a finite logical Grassmann algebra external to a product massive auxiliary trace. In (3), keeping l2(N_0) shifts external instead gives noncommuting operator coefficients with gap dynamics linking neighboring cells. A mixed operator-valued version, or a separate controlled charge representation, must justify the connected expansion and source matching. Calling the gaps additional fermion modes, replacing the shifts by canonical bosons, or putting the gaps in C48's massive heat trace changes the problem. Equations (1) rule out the first two operator identifications directly.

[Assumed, requested coefficient output] What is needed is an actual matched action with fixed field weight H>1 and fixed positive spatial and temporal weights, all external degrees summed absolutely, and at least two controlled time moments. The constants must be specified at the same finite J,K,lambda, normalized by the smallest payload scale. No such constants, and no numerical Delta_*, have been derived for C88. The temporal decay of the massive influence cannot be obtained from an ultraviolet Hamiltonian range estimate alone.

## 5. A restricted full-space parity obstruction

[Proved] The failure of exp(i pi Q/2) is not repaired by assigning arbitrary diagonal signs to the thirteen microscopic letters. Let p(a) in Z_2 be such signs and suppose their product is conserved by every nonzero literal move, with g,t_d,lambda nonzero. Vacuum moves require

    p(R_a)=p(L_a).

The switch requires p(D_s)=p(X_s), and terminal moves require

    p(M_s)=p(D_s)+p(E)+p(L_0).

The port E R_0 -> M_A D_B then gives

    p(R_0)=p(D_A)+p(D_B)+p(L_0),

so p(D_A)=p(D_B). The conversion D_A R_1 -> D_B R_0 now gives

    p(R_1)=p(R_0).

Therefore

    p(L_0)=p(L_1)=p(R_0)=p(R_1).                          (21)

Such a grading is blind to every stored payload bit. The two allowed vacancy colors in a hop across R_0 also require p(V_0)=p(V_1). The result is independent of the magnitude of the nonzero couplings and of their complex phases.

On the frozen vacuum code the marker and head counts are fixed, and (21) makes the product sign independent of all thirteen payload occupations. It cannot equal (-1)^(N_payload). The exact physical parity exp(i pi Q) is one such bit-blind symmetry. The quantity exp(i pi Q/2) gives the desired sign on the code but squares to -1 on odd-Q carriers.

[Numerical, exact finite enumeration in integer arithmetic] checks.py examines all 8192 sign assignments against the 37 oriented rules, omitting their nonzero amplitudes. Sixteen assignments survive; all satisfy (21). This is a redundant finite control of the algebraic proof.

[Proved, limitation of the obstruction] This excludes diagonal product gradings at the microscopic letter level. It does not exclude a non-diagonal block grading, a dressed representation, a graded auxiliary enlargement with a separately justified equivalence, or elimination of odd-Q carriers before fermionizing the logical code. A claim that no local representation of any kind exists would go beyond (21). C48 requires one of these interfaces to be supplied, not just the code's charge selection rule.

## 6. Finite K and the order of choices

[Proved, imported finite-order identity] C88 constructs a local integer-penalty conjugation at kappa=JK,

    U_K H U_K^* = kappa N_bad+D_2+R_K,
    [D_2,N_bad]=0,
    ||R_K||_a <= C_3(a) nu^3/(JK)^2.                     (22)

Its legal second-order correction is

    D_color+ D_vacuum+ D_port
       =-[t^2 C+J^2 A_counter+2J^2 lambda^2 Z]/(JK).       (23)

The counter C counts R_0 or E immediately followed by a physical vacancy. The mixed tJ correction vanishes. C_3(a) is not numerically specified by C88, and (22) is a finite-order identity retaining a remainder, not an exact invariant-sector projection. Its physical support norm must also be transported to any chosen rooted coefficient representation.

[Proved, explicit nonuniform field in the decoupled reference] In the product of the frozen vacuum clock and the homogeneous free gap reference, a stored bit r is in the R tail when the clock index l<=r. Its contribution to C is then (1-n_r) times the nonempty-gap indicator. Define

    c_r=sum_(l=0)^r (2/16)sin^2(pi(l+1)/16),
    w_r=(1-alpha_c)c_r,                  r=0,...,12.

Projecting the displayed order-1/K term onto that reference gives, apart from its scalar,

    D_color,ref=+(t^2/(JK)) sum_(cells c,r) w_r n_(c,r).   (24)

Every increment w_(r+1)-w_r is strictly positive. This is not a common chemical potential or a scalar and cannot be removed by changing the occupied-number source. Nor can one common payload field cancel its thirteen unequal coefficients. The centered gap and clock factors remain before taking this reference projection. Formula (24) is a coefficient of the stated reference calculation, not an all-order finite-J self-energy theorem.

[Numerical] At b=15 the smallest and largest w_r are 0.00334284159219 and 0.68643707665319, with spread 0.68309423506100. These numbers check the strict nonuniformity. No finite-volume eigenvalue inference is made from them.

[Definition] Use W_ref=J alpha_c lambda^2 E_unit as the reference payload energy unit, where E_unit is C88's alpha_end^4 beta^2 g^2. This is a normalization unit, not the smallest velocity. The latter can be much smaller. Let epsilon_band be a positive tolerance in this energy unit.

[Proved, constraint on this coefficient] Requiring the largest coefficient in (24) to be at most epsilon_band W_ref is equivalent to

    lambda >= lambda_min,
    lambda_min^2 = w_max/(epsilon_band alpha_c E_unit K s^2).
                                                                  (25)

For fixed finite J,K this is a positive lower bound. Thus taking lambda arbitrarily small last does not keep the same certified reference perturbatively close. The finite-g certificate does not suppress this term. At K=infinity (25) disappears, which is a useful consistency check. At t=0 the color term also disappears, in agreement with the frozen material result.

This is a requirement for treating the specific leading coefficient as a small perturbation of the stated quadratic reference. It is not a theorem that every finite-J chain violating (25) is gapped or has fewer modes. Large corrections could instead define a different reference needing a new certificate. No cancellation of (24) against uncomputed higher terms is proved here.

[Proved, the explicit ultraviolet and color budgets can coexist] The objections above are not an algebraic incompatibility of every finite parameter choice. For example, set t=1 and consider

    J=2^-5000,  K=2^14000,  lambda=2^-608,  delta=2^-250.

This is only an illustration for the legal ultraviolet remainder and the leading color coefficient. Since C_R<357, q<2^-4224. Also C_UV<2^80 and 2^15 M0 C_R^2<2^45. Equation (13) gives

    ||R_UV||_(1,0)/J
       <2^-2420+2^1295 exp(-2^498)<2^-2419.

The recorded lower bound alpha_end^4>=2^-18, together with d_15<1004^15<2^150, gives E_unit>2^-862. Since alpha_c>1/4,

    W_ref/J>2^-2080.

Thus the displayed ultraviolet first-moment remainder is below 2^-35 W_ref. The largest color coefficient divided by W_ref is less than

    2^(10000-14000+2080)=2^-1920.

This comparison leaves out the finite-K remainder, the strong coefficient conversion, the low-frequency filtered operator, sources and the dressed gap. In particular this polynomial cutoff choice does not meet (20) at a fixed positive spatial rate. It is not a certified mobile member. It shows why failure of the existing proof cannot be restated as a contradiction among all finite scale choices.

[Numerical, scale illustration] For C88's g=2^-272 and lambda=2^-608,

    log_2 E_unit = -858.57741008...,
    log_2(W_ref/J) = -2076.32713907....

At the allowed floor K=2^2560 and s=2^-250, the largest reference color field divided by W_ref is about 2^15.7843. To put this coefficient alone below 2^-35 W_ref at that K requires s>2^-224.608 approximately. These are arithmetic comparisons of known coefficients. They do not certify admissibility or inadmissibility of the full interacting point. Strong-norm conversion and the smallest-velocity tolerance require a stricter budget.

[Proved, how a lambda-last construction could be arranged conditionally] Suppose a new theorem provided a positive upper bound Lambda(J) for port matching, uniform in all K above a threshold. Fix J, and choose

    K > 4 w_max/[epsilon_band alpha_c E_unit s^2 Lambda(J)^2].

Then lambda_min<Lambda(J)/2 and one may choose lambda=Lambda(J)/2 last. This is a genuine interval construction, not a limit lambda->0 at fixed K. It also needs the other grammar and matching budgets. Uniformity of Lambda in K is an assumption here; if the theorem's constants depend adversely on K, joint feasibility must be proved instead.

[Proved, further budget restriction] A crude ultraviolet bound of order J^2/t without a port factor, divided by W_ref, has order s/(alpha_c lambda^2 E_unit). The power-ten bound (13) has the still larger factor s delta^-10/(alpha_c lambda^2 E_unit). It cannot be treated as a lambda^3 remainder. Some of its zero-port contributions may be bit-independent and absorbable into a dressed clock and charge reference, but that cancellation requires a construction. A total operator norm does not prove it. The same issue applies to the finite-K remainder in (22).

## 7. The matched logical action that is justified, and the one requested

[Proved, imported first-order reference matching] The free charge sea gives the mean contact alpha_c=1/2-2/pi^2. Its first-order material Hamiltonian is J alpha_c H_mat. Within a fixed cyclic sector the coefficient follows by projecting the total number of adjacent occupied pairs and dividing by N. This is the finite-ring projection identity in the collaborator's doping_check/README.md, with the same thermodynamic limit for the unique legal coloring. No finite-ring charge gap is used as a uniform denominator.

[Proved, imported frozen logical coefficient] Let H_quad^(13) be the dimensionless quadratic reference specified by C88's two-carrier resolvent, including its normalized field 2^-120. Let V_g be the static nonquadratic part of that frozen second-order port coefficient. The formally averaged retained Hamiltonian has the reference terms

    H_ref,ret=tH0 + E_scalar
                +J alpha_c lambda^2 E_unit
                         (H_quad^(13)+V_g)
                +(t^2/(JK)) sum_(c,r) w_r n_(c,r)
                +other retained grammar terms.                    (26)

The material quadratic reference has 26 simple Fermi points and the imported margins delta_2>0.0156742281 and delta_4>3.2513114e-6. Equation (26) records the coefficients that are known in the combined reference expansions. It is not asserted to be an exact finite-J Hamiltonian, or a matched action with a controlled remainder. In particular averaging the material term first and then eliminating carriers is justified as this reference calculation, not as commuting exact matching operations.

[Heuristic, expected low-energy form] In cell coordinate y=j/15, a possible initialized charge action would have reference part

    S_charge=(1/(2 pi K_gap)) integral dy d tau
                 [v_c^-1 (partial_tau phi)^2
                      +v_c (partial_y phi)^2],
    K_gap=4,              v_c=t/15.                       (27)

A mixed logical action would then take the form

    S_match = beta M e_match + S_charge
              +sum_(a=1)^13 sum_(eta=R,L)
                 integral dy d tau
                 psi_bar_(a,eta)(Z_(a,eta) partial_tau
                         -i v_signed,(a,eta) partial_y)psi_(a,eta)
              +S_forward+S_nonforward+S_retarded.          (28)

The forward part includes charge-charge, payload-payload and charge-payload currents. The color field and all self-energies shift the quadratic reference and must be retained in the matching. All scalar terms and physical sources must also be matched. There are thirteen proposed payload channels and one charge channel; the vacancy colors add none. The symbols in (28) do not stand for coefficients calculated in this report.

[Assumed, quantitative contract needed to upgrade (28)] One needs a valid local graded or replacement representation, an explicit dressed threshold Delta_*>0 and heat bound, an exact selected-trace matching with its source representatives, and coefficient estimates for the remainder in a fixed strong norm. For example, with spatial rate alpha>0 and temporal rate kappa_t>0, a sufficient all-degree bound is

    sum_p H^p sup_(one external anchor)
        sum_(species,relative sites) integral_(relative times)
        exp(alpha diam_y+kappa_t diam_tau)
        |K_p| <= epsilon W_min,                           (29)

together with scalar-density control, the required degree-order hierarchy or an alternative proved summability theorem, and bounds after two time moments. W_min must account for the slowest branch, not just W_ref. Fixed source derivatives must satisfy the same estimates on one common source interval. A fixed-volume cutoff cannot enter the constants.

If those data satisfy (18), C48 supplies its conditional leakage-linear endpoint bound. In its notation

    ||delta S||_LE <= A_dec r,
    A_dec=2^22 exp(L_w)[2+16(16Theta/L_w)^2],
    ||delta S(omega)-delta S(0)||_LE
        <=A_dec r |omega|/Delta_*.

Its massive thermal winding term must be retained. This implication is imported from C48. It is not verification of any of the missing input data for the mobile chain.

[Proved, conditional harmonic dictionary only] If (27)-(28) and the physical number-field identification are established, a cell vertex with momentum q and payload dimension Delta_pay has branches

    Q_h=q+2pi l,
    k_phys=Q_h/30 mod 2pi,
    Delta_ref=Delta_pay+[Q_h/(30pi)]^2.

Uniform branches have Q_h=60pi h and charge dimension 4h^2. This is C88/C90's conditional resonance rule. It neither constructs (28) nor bounds (29). A central-charge statement would additionally require the infrared and thermal/completeness inputs distinguished by C90. None is concluded here.

## 8. Verification and final boundary

[Numerical] checks.py uses only local source, NumPy, SciPy and the standard library. Its largest dynamical matrix is the 35-dimensional four-gap, four-vacancy phase Hamiltonian. It checks the Gaussian filter's sign and first-order identity, anti-Hermiticity, and convergence of the second-order Taylor coefficient for a single exponential. It checks the actual length-fifteen vacuum clock leakage, the color-field weights and parameter arithmetic, the exact sign-assignment enumeration, and the geometric-sum constants. It does not diagonalize the d=13 many-body chain, determine its dressed gap, or measure a central charge.

[Proved, outcome] Equations (1)-(14) provide exact coordinates for the legal compression and an explicit volume-uniform ultraviolet unitary and remainder. Equations (16)-(17) show why these do not isolate a massive material code. Equations (18)-(20) separate the available spatial estimates from C48's strong coefficient input. Equation (21) excludes a simple microscopic grading repair. Equations (24)-(25) quantify the finite-K obstruction to an unrestricted lambda-last limit.

[Assumed, still open] The task's central requested object, a finite-J mobile material code with a verified dressed gap and source-uniform all-degree bounds at one jointly admissible finite schedule, has not been supplied. The precise next work is to construct that code at zero port while retaining the charge bath, then prove the coefficient and grading interface before invoking C48. A successful proof would turn (28) into a matched action. The present results do not do so.
