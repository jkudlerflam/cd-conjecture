# C85. Three colored vacancies for the eleven-letter grammar

[Scope] The finite result is a d=14 lift with n=13 active flavors and physical material-cell length 15. Its occupied quadratic reference has a regenerated integer/rational certificate. F3 holds for full frozen marker sectors. F4 holds as a finite-order local penalty normal form. The mobile central-charge conclusion remains conditional on the infrared inputs stated below.

[Tags] Proved denotes an algebraic argument or an integer/rational acceptance test. Numerical denotes finite floating-point calculations, with dimensions stated. Assumed denotes an imported infrared input. Heuristic denotes the resulting physical interpretation. Inherited proved results are attributed. None of the ED calculations measures a central charge.

## 1. The grammar has three states

[Definition] The material alphabet is C84's

    L_0,L_1,R_0,R_1,D_A,D_B,X_A,X_B,E,M_A,M_B.

Use p for data sites, n=p-1 for active bits, and b=p+1=n+2 for cell length. The extra data endpoint carries the fixed zero bit and does not convert. It cannot be removed from the carrier propagator.

[Definition] The outgoing state o and incoming set I are

| Letter | o | I |
| :-- | :-- | :-- |
| L_0 | 0 | {0,2} |
| L_1 | 2 | {0,2} |
| R_0 | 1 | {0,1,2} |
| R_1 | 1 | {1} |
| D_A,D_B,X_A,X_B | 1 | {0,2} |
| E | 0 | {0,1} |
| M_A,M_B | 0 | {0} |

A pair ab is allowed if o(a) is in I(b). State 2 records the L_1 boundary pin. In particular L_1 E and L_1 M_s are forbidden, while L_0 E and L_0 M_s are allowed. Merging states 0 and 2 would erase this pin. This table reproduces all 121 ordered pairs of C84.

[Definition] Add V_0,V_1,V_2 with o(V_q)=q and I(V_q)={q}. Let p_i be the mismatch projector on physical bond i, N_bad=sum_i p_i, and P=1_(N_bad=0). The dimension is

    d=11+3=14.                                                (1)

[Proved] Suppose consecutive material letters a,b have a nonempty vacancy string between them. The first physical bond fixes its first color to o(a). Every vacancy-vacancy bond copies that color. The last bond is legal if and only if o(a) belongs to I(b). Therefore every legal material word, at every fixed nonempty occupancy pattern, has exactly one legal coloring. The all-vacancy sector has three constant legal words and is outside this statement.

Each illegal squeezed material pair forces at least one physical mismatch in its intervening interval. These intervals are disjoint on the ring, including the seam. Hence

    N_material_bad <= N_bad                                  (2)

as diagonal operators in every nonempty occupancy sector. The argument includes the one-material-letter ring.

[Proved, restricted minimum] Three colors are necessary for any nonnegative diagonal nearest-neighbor repair that admits all legal words at arbitrary vacancy placements and rejects every illegal squeezed word. Consider the three legal diagonal pairs

    (L_1,L_1), (R_1,R_1), (M_A,M_A).

For a single intervening vacancy color q, allowed endpoints form a rectangle A_q times B_q. Two of these three diagonal pairs cannot share a color, since at least one of their two crossed pairs is illegal. L_1 R_1 and R_1 L_1 are both illegal; L_1 M_A is illegal; R_1 M_A and M_A R_1 are both illegal. Each diagonal pair occurs in a complete legal material ring, so requiring all legal embeddings really requires its rectangle coverage. Thus three distinct colors are needed even if one permits nondeterministic color assignments. The table attains this bound with unique coloring for arbitrary gap lengths. This is a restricted minimum for this unchanged grammar and pair-penalty mechanism, not a no-go for redesigned alphabets.

[Numerical, exhaustive] `check_lift.py` compares all material permissions with the local C84 copy and checks (2) on every nonempty four-site word, 14^4-3^4=38335 words. It also checks uniqueness among all legal colored four-site words.

## 2. The full bond

[Definition] Let z_A(a)=(-1)^a and z_B(a)=1. The following amplitude multiplies output ket times input bra. Add the Hermitian conjugate of every row, sum a=0,1 and s=A,B, and sum over all three markers m where shown.

| Input | Output | Amplitude |
| :-- | :-- | :-- |
| R_0 R_a | L_a R_0 | -1 |
| R_0 E | L_0 E | -1 |
| D_s R_a | L_a D_s | -z_s(a) |
| X_s R_a | L_a X_s | -exp(i phi_s) z_s(a) |
| X_s E | L_0 M_s | -exp(i phi_s) |
| D_s E | L_0 M_s | -t_d exp(i psi_s) |
| m D_s | m X_s | -1 |
| D_A R_1 | D_B R_0 | g |
| E R_0 | M_A D_B | lambda exp(-i pi/4) |
| E R_0 | M_B D_A | lambda exp(+i pi/4) |

There is no conversion at the last data endpoint. The one-site diagonal o_site assigns Delta to D_s,X_s,M_s, an additional D to X_s,M_s, mu_sel to every marker, and f a to L_a,R_a. Let w denote the transition table, and set

    b_11=w+(o_site tensor I+I tensor o_site)/2,
    h_11=K p_material+b_11.                                  (3)

The complete colored bond is

    h_14 = JK p_ext + J embed(b_11)
       - t sum_(a material) sum_(q in I(a))
           (|a,V_o(a)><V_q,a|+h.c.)
       - (mu/2)(n_occ tensor I+I tensor n_occ).               (4)

Here embed extends by zero outside the material-material block. In particular the half-site energies stay inside occupied contacts. Moving them onto physical sites would change this Hamiltonian. Its occupied block is J h_11. Equation (4) is a single Hermitian translation-invariant nearest-neighbor bond, not an occupancy-dependent choice of bonds.

[Proved] The exact charges are

    q(L_a)=q(R_a)=2a, q(E)=q(V_q)=0,
    q(D_A)=q(X_A)=q(M_A)=-1,
    q(D_B)=q(X_B)=q(M_B)=+1.

All rows preserve total q. Conversion preserves -1+2=1+0; each port output has total zero. Vacancy hopping transports the same material letter. Consequently Q=sum q, occupied number N, and marker number M commute with the full chain, including every illegal state. Vacancy color is not conserved. Material order changes only by the displayed internal transitions. Marker order and the number of material sites between consecutive markers are conserved up to cyclic relabeling.

[Numerical] `bond14.py` constructs (3)-(4) from the printed table. Its material block agrees entry for entry with C84's `bond11.py` at two generic parameter points, including a nonzero field and generic phase. Hermiticity and the N,M,Q commutators vanish entry by entry. The local bond dimension is 196 by 196.

## 3. Projection and illegal-weight suppression

[Proved] In a legal word a reverse hop a V_o(a) -> V_q a has one legal choice of q, fixed by the preceding site's outgoing state. A forward hop V_q a -> a V_o(a) is legal whenever its input is legal. Thus every projected swap has amplitude -t, without a multiplicity factor. All material moves preserve the outgoing state of their rightmost site. Any change of incoming permission at their leftmost site is tested by the preceding physical bond. The coloring map therefore intertwines their legal projections.

Write U_col for the coloring isometry and g_j for the vacancy number after material letter j. Including the physical root and cyclic word rotation on a ring,

    U_col^dagger P H P U_col
      = H_hop + J sum_j 1_(g_j=0) P_mat b_11,j P_mat - mu N. (5)

The coloring is rotation-equivariant and adds no legal internal degree of freedom. Three colors do not give three charge bosons.

[Proved] An illegal squeezed junction cannot be concealed by a vacancy string. For example R_0 V_1 R_0 R_a can undergo a vacuum move to R_0 V_1 L_a R_0. The new V_1 L_a bond costs JK. The new boundary-pin obstruction is L_1 V_2 E, which costs JK at V_2 E. Recoloring that vacancy to 0 moves the mismatch to L_1 V_0. Neither obstruction can be removed by lengthening the vacancy interval.

[Proved, fixed finite volume] Set kappa=JK and H=kappa N_bad+B. Let P be the legal projector, Q_bad=1-P, and assume the invariant subspace has a legal trial state. A ground energy E obeys E<=||B||. For kappa>2||B||,

    Q_bad(H-E)Q_bad >= (kappa-2||B||)Q_bad,
    Q_bad psi = -[Q_bad(H-E)Q_bad]^(-1) Q_bad B P psi,
    ||Q_bad psi|| <= ||B||/(kappa-2||B||).                    (6)

The illegal-material projector is bounded by Q_bad using (2), so its weight is bounded by the square of (6). This is K^-2 suppression at fixed positive J and fixed volume. It is not a bound uniform in volume on total leakage. Section 7 supplies the local statement.

## 4. Every leading grammar energy, including colored hopping

[Definition] Write B=t T+J W-mu N, where T has unit hopping amplitude and W is the unpenalized material contact Hamiltonian with its port amplitude lambda included. On P, every illegal intermediate reached by T or W has exactly one mismatch. The second-order penalty correction is

    D_2|P = -P B Q_bad N_bad^(-1) Q_bad B P/(JK).            (7)

The following diagonal counters give its entire compression. All predicates refer to physical sites, with cyclic indices.

    C = sum_i 1_(a_i material, a_(i+1) vacancy) c(a_i),
    c(a)=|I(a)|-1,
    A = sum_i 1_(o(a_(i-1))=1, a_i=R_0,
                  a_(i+1) in {R_0,R_1,E}),
    Z = sum_i 1_(o(a_(i-1))=1, a_i=E, a_(i+1)=R_0).

[Proved] Their coefficients are

    D_2|P = -P [ t^2 C + J^2 A + 2J^2 lambda^2 Z ] P/(JK). (8)

This equality includes nonzero bits and all legal auxiliary states. There are no other terms at order 1/K, and in particular the mixed tJ/(JK) term vanishes.

Here is the channel-by-channel proof. Starting from a legal a V_o(a) pair, the reverse hop has |I(a)|-1 wrong colors. Each creates one mismatch immediately before the new vacancy. Its only one-step return to a legal word is its inverse. Thus T^dagger Q_bad N_bad^(-1)Q_bad T is the diagonal counter C. The full table is

| Letter | Wrong reverse channels c(a) |
| :-- | :-- |
| L_0,L_1 | 1 |
| R_0 | 2 |
| R_1 | 0 |
| D_A,D_B,X_A,X_B | 1 |
| E | 1 |
| M_A,M_B | 0 |

A material move can leave P only when it narrows the incoming set of its first letter. There are just two such moves. The vacuum move R_0 R_a -> L_a R_0, or R_0 E -> L_0 E, rejects predecessor state 1 and contributes unit squared amplitude, counted by A. The port E R_0 -> M_s D_other also rejects predecessor state 1. It has two orthogonal outputs of squared amplitude lambda^2, counted by 2 lambda^2 Z. Every such illegal state has only the inverse move leading back to P. Rail moves, switches, conversion, and their inverses either keep the incoming set or enlarge it. Diagonal terms have no P-to-Q_bad matrix elements.

For a hopping-generated illegal state the second site of its sole mismatch is a vacancy. For a material-generated illegal state that second site is material. These sets of intermediate basis vectors are disjoint. This proves P T Q_bad N_bad^(-1) Q_bad W P=0, not just cancellation of its expectation. It also proves that there are no g^2/K, t_d^2/K, or mixed vacuum-port terms at this order. The port phases do not enter the squared amplitudes.

[Proved, physical names and payload content] The three terms in (8) are respectively the color-ambiguity contact self-energy, the vacuum-tail self-energy, and the wrong-end port self-energy. Their scales are t^2/(JK), J/K, and 2J lambda^2/K. The first two have no lambda^2 factor.

The color self-energy is bit-dependent because c(R_a)=2(1-a). At a fixed legal zero-excitation clock configuration and fixed gaps,

    C = sum_(L sites) 1_(vacancy after site)
        +2 1_(vacancy after vacuum head, if present)
        +2 sum_(tail bits r) (1-a_r)1_(vacancy after R_(a_r))
        +1_(vacancy after E).                              (9)

In the terminal all-L configuration the head is absent. Equation (9) names the vacuum-head contact energy, the L-prefix contact energy, the empty-marker contact energy, and the tail-bit contact field. In auxiliary configurations the D/X contact energy has coefficient one, while the M-marker contact coefficient is zero. L_0 and L_1 have the same coefficient despite their different outgoing colors.

After projecting a bit-independent clock/charge reference, the tail part of -t^2 C/(JK) becomes a scalar plus

    +(2t^2/(JK)) sum_r w_r a_r,
    w_r=<1_(bit r in R tail and vacancy immediately after it)>_reference.

Before averaging, the difference from this mean is a payload-density coupling to the vacancy and clock variables at the same scale. This is the t^2/(JK) effect that must not be omitted from the mobile budget. It is not an extra flavor or a harmless common scalar.

Similarly, at full occupancy the vacuum-tail term projected onto C84's vacuum clock is

    -(J/K) sum_(r=0)^(n-1) (1-a_r) sum_(j=0)^r |v_j|^2,
    v_j=sqrt(2/(p+2)) sin(pi(j+1)/(p+2)).                   (10)

Bit r sits at data position r+1 when j<=r, and its R_0 hop then has a state-1 predecessor. This proves (10). With vacancies the same term carries the physical contact indicator in A and becomes another density-clock-charge coupling. The counter Z is independent of the bit values on the zero-excitation legal clock basis, but depends on clock endpoints and contacts. Its projections through auxiliary dressing need not be scalar, so its full operator is retained in the error budget.

[Numerical, exhaustive] `check_virtual.py` builds P-to-illegal matrices directly from the literal bond on every legal nonempty ring of lengths 3,4,5, with respective legal dimensions 313,2127,14497 and 144,1256,10690 distinct illegal intermediates. It verifies all three counters, diagonality, unit denominators, and the vanishing mixed operator. The maximum discrepancy in (8) is below 4.5e-16. This includes auxiliary letters, nonzero g,t_d,lambda, and a nonzero field.

[Numerical, energy test] `check_energy.py` uses a p=3 cell plus one vacancy, all four active-bit strings, J=t=1, lambda=0. The constrained energies agree at -3.29813082968. The color and vacuum coefficients of the negative 1/K energy shift are

| Bits | Color coefficient | Vacuum coefficient |
| :-- | --: | --: |
| 00 | -1.3308967188 | -0.4480522957 |
| 01 | -1.0220888489 | -0.1024562307 |
| 10 | -1.2719499702 | -0.3455960650 |
| 11 | -0.9631421003 | 0 |

The invariant dimensions are 440,310,310,190, with legal dimension 80 each. At K=10000, K(E_K-E_infinity) differs from the sum of these coefficients by less than 3.1e-4. The reference clock amplitudes are independent of bits after applying the coloring isometry. The color coefficients consequently have zero two-bit discrete second difference to numerical precision. Finite K is not asserted to preserve this one-body form.

[Proved, higher-order bookkeeping] Equation (8) is the complete leading grammar correction, not an assertion that all higher terms vanish. At the next order the possible homogeneous strengths are

    t^3/(JK)^2,
    t^2 J B_material/(JK)^2,
    t J^2 B_material^2/(JK)^2,
    J^3 B_material^3/(JK)^2.                               (11)

They represent corrected vacancy hopping, material-energy and clock-dependent color denominators, mixed contact-dependent hopping, and higher material grammar corrections. They can generate bit-dependent hopping, bit fields, and multi-bit interactions. Insertions of the diagonal Delta,D,mu_sel,f terms, rail amplitudes 1,t_d, conversion g, or port lambda are included in B_material, not dropped. The chemical potential -mu N is constant in each conserved-N sector and cancels from internal virtual denominators. In particular there are t^2 J Delta/(JK)^2 and t^2 J D/(JK)^2 auxiliary-energy corrections, and field-dependent corrections with f in place of Delta. Terms with one port may be off-diagonal in auxiliary excitation number and can generate further payload terms when that sector is eliminated.

A complete finite-volume definition of every higher correction is the Feshbach operator

    H_eff(E)=PHP-PBQ_bad
        [JK N_bad|Q_bad+Q_bad B Q_bad-E]^(-1)Q_bad BP.       (12)

Expanding its inverse generates (7) and all insertions of Q_bad B Q_bad-E, with their signs and integer penalty denominators. The local unitary recursion in Section 7 gives the energy-independent fixed-order version. All terms beyond (8) are kept in a named local remainder. Neither (11) nor (12) is used to claim exact Gaussianity, and no unnamed bit-dependent scalar is silently discarded.

## 5. The regenerated n=13 finite certificate

[Inherited proved input, C84 Sections 2-3] A legal p-data cell has (5p+3)2^(p-1) states, with p+1 vacuum-clock positions and 2^(p-1) vacuum code states. The clock energy and endpoint amplitude are

    E_v=-2 cos(pi/(p+2)), alpha=sqrt(2/(p+2)) sin(pi/(p+2)).

The two-species CAR reduction remains valid after removing the last conversion. Its scalar carrier has ordered vertices

    M,X_(p-1),...,X_0,D_0,...,D_(p-1),

diagonal 1004 on M,X and 4 on D, path entries -1 and closing entry H_(M,D_(p-1))=-i t_d. For p=14 it has 29 vertices. Only D_0,...,D_12 convert. The inactive D_13 and its X partner remain in the full 841-dimensional product resolvent.

[Definition] Let R=(H_c tensor I+I tensor H_c)^(-1), x=R|M,D_0>. U_i retains row D_i of x and V_i column D_i, viewed as product-graph vectors. With i,j=0,...,12,

    C_ij=U_i^dagger R V_j,
    A_ij=U_i^dagger R U_j+V_i^dagger R V_j,
    T_quad=-2 alpha^4 g^2 Re C,
    mu_quad=-2i alpha^4 g^2 Im A.                           (13)

The g^2 port coefficient is scalar plus this number-conserving quadratic operator by C84's CAR calculation. Higher even g orders can contain interactions.

[Definition] Fix a single microscopic material member by

    n=13, p=14, b=15, D=1000, Delta=4-2 cos(pi/16),
    phi_A=pi/28, phi_B=0, psi_A=pi, psi_B=pi/2,
    d_0=1, d_1=1004, d_(j+1)=1004 d_j-d_(j-1), beta=1/d_15,
    t_d=beta 2^-66, g=2^-272, lambda=2^-608, K>=2^2560,
    E_unit=alpha^4 beta^2 g^2, f=lambda^2 E_unit 2^-120.      (14)

The selector below fixes mu_sel. These couplings stay fixed in cells of every other length. `MEMBER.json` records their exact definitions. For the mobile problem K may need to be larger as specified in Section 9.

[Proved, computer-assisted] `run_certificate.py` reruns C84's scalar pipeline at n=13 rather than reusing C83's n=13 outputs. Tree and derivative solves use 256-bit fixed point; the finite-loop solve uses 512 bits. Here the express block has 15 vertices and the data block 14. The Jacobi contraction is at most 1/2+t_d/4. The inherited ten-unit real/imaginary component radii cover truncation and iteration error with bit_length(d_15)+512+20 steps in the loop solve. Product radii are 10(|x_r|+|x_i|+|y_r|+|y_i|)+200 integer units. Integer LDL identities and rational congruences transport these errors to the dyadic metric W.

The Hermitian hopping margin is greater than 0.5560750. The normalized matrix error is less than 6.339e-92. Numerical generalized eigenvalues propose brackets only. For rational t=tan(k/2), use cos k=(1-t^2)/(1+t^2) and sin k=2t/(1+t^2). A rational approximation to a proposed eigenbasis is proved invertible by its Gram matrix. Strict rational diagonal dominance after congruence, including all entry errors, proves the endpoint signs and inertia.

There are 26 disjoint intervals. Inertia increases from 0 to 13 and decreases back to 0. Each interval contains a sign-changing determinant zero. Multiplication of det h(k) by exp(13ik) gives a polynomial of degree at most 26, so these roots exhaust the degree and all are simple for this nearest-neighbor quadratic pencil.

Exact rational unit phases are formed from bracket midpoints. A bracket's chord error is bounded by its t-width. Product phases for unordered pairs have error at most the sum of the two chord errors. Subtracting these radii from exact adjacent phase distances proves

    delta_2 > 0.0156742281,
    delta_4 > 3.2513114e-6,                                (15)

including all 26*27/2=351 unordered pairs with repeats. This also separates all sea centers modulo pi. The permitted additional Hermitian Bloch norm in the dyadic metric is 2^-35. This is low-order nonresonance, not rational independence against every integer momentum relation. The exact brackets are in `root_certificate_exact_n13.json`.

[Proved, finite material error budget] The source norm is at most sqrt(2) alpha^2, pair conversion norm at most 2g, and source gap at least three. Odd conversion orders do not return the source species. The inherited even Neumann tail gives

    ||H_2(g)-H_2(0)-g^2 H_22|| <= alpha^4 g^4/4.            (16)

The selected cell has fifteen sites, onsite norm below 1007, and transition norm sum below 64 per bond, so B_0=2^16 remains a valid cell norm bound. Its vacuum gap exceeds 1/16. The new boundary pin is part of N_bad, not B_0. The high-block inverse at |z|<=2B_0 is at most 2/K for K>=8B_0, giving one-cell Schur correction 2B_0^2/K. The two-cell spectral projection and polar-overlap estimates inherited from C84 give energy shift at most 2^36/K and code rotation at most 2^44/K.

The second-order coefficient uses the quasidegenerate formula

    (H_2)_ab=-1/2 sum_q V_aq V_qb
        [(E_q-E_a)^(-1)+(E_q-E_b)^(-1)].                   (17)

A common denominator would erase the finite-K bit splitting. The source-resolvent and rotation estimates bound each column error by 2^48/K. There are 4^13 columns, so the column-to-operator factor is at most 2^13. Allowing symmetrization and source contractions gives the enlarged bound 2^68/K. Including the physical zero-port energies, field compression, and code rotation is covered by 2^84/K, without a lambda^2 prefactor. A contour of radius 1/64 and analytic port disk of radius 2^-12, with the two-cell excitation parity making the effective operator even in lambda, give the conservative isolated-pair remainder 2^52 lambda^4. These are finite-cell bounds, not a many-cell infrared theorem.

Write tau_W=tr(W W^dagger). A normalized quadratic coefficient is bounded by four times operator norm; summing onsite and hopping terms is bounded by the factor 16n. Since alpha^4>=64/16^6=2^-18, exact rational arithmetic in `finite_budget.py` proves

    err_g = 4n beta^-2 g^2 tau_W < 1.727e-37,
    err_K = 16n tau_W (2^84/K)/(lambda^2 E_unit) < 3.024e-83,
    err_lambda,local = 16n tau_W 2^52 lambda^2/E_unit < 2.396e-54.

Their sum is below (2^-35)/8. Thus the finite-g, finite-K quadratic projection through second order in lambda retains the certificate. The root routine's scalar threshold g<=2^-228 is not the final two-species choice 2^-272. The latter is used in every final budget. No finite-K degeneracy, exact finite-g Gaussianity, or full interacting root-degree bound is assumed.

## 6. F3 for full frozen marker sectors

[Inherited proved bound, C84 Section 6] Omitting E_1 is compression of the C80 material space and the new K penalty on L_1 followed by a marker is positive. Therefore for every cell length l, on the complete cell Hilbert space including illegal states and arbitrary auxiliary excitations,

    E_l^infinity-(l-2)_+/(K-4) <= E_l(K) <= E_l^infinity,
    E_l^infinity=mu_sel-2 cos(pi/(l+1)).                    (18)

The upper trial state fixes the last bit to zero. The hypotheses Delta>2+g and |t_d|<=1 hold at (14). A nonnegative field preserves the lower bound and vanishes on the all-zero trial state. The marker-free ring density is at least -1/(K-4), by the same compression argument.

[Definition] Set a_l=2[1-cos(pi/(l+1))], delta_l=(l+1)a_l-l a_(l+1), and

    mu_sel=2-(delta_15+delta_14)/2.                         (19)

[Proved] For x>=1,

    a''(x)=4pi sin(pi/(x+1))/(x+1)^3
          +2pi^2 cos(pi/(x+1))/(x+1)^4 >0.

Hence delta_l decreases strictly, while e_l=E_l^infinity/l obeys

    e_(l+1)-e_l=(2-mu_sel-delta_l)/(l(l+1)).

The midpoint (19) uniquely selects b=15; the closest competing density is at 14 or 16. Rational alternating-series enclosures for pi and cosine in `selector_certificate.py` prove

    mu_sel=1.8891880765550189...,
    delta=min_(l!=15)(e_l-e_15)>2.8116274e-5>2^-16,
    e_15<-0.0048254989<-2^-8,
    vacuum gap>0.11381149>1/16.                            (20)

[Proved, cell bookkeeping] All three markers have outgoing state 0, so a marker-to-first-letter mismatch belongs to the right cell. The new last-L_1 penalty depends only on the last data letter and the fact that the next site is a marker, not its type, so it belongs to the left cell. The switch is the identity on the preceding marker and also belongs to the right cell. Thus at zero port the full cells are tensor factors even off grammar. The sole coupling between cells is the port, of norm B_port=sqrt(2)|lambda|, the singular value of one input with two orthogonal outputs.

Let A_m be the exact ground energy of m selected complete cells in an open chain. Cutting one port gives

    |A_(m+r)-A_m-A_r|<=B_port.

Almost additivity gives u=lim A_m/m. Applying the bound to q copies of the m-cell chain and letting q grow gives |A_m-m*u|<=B_port. Thus

    A_m >= m*u-B_port,
    S_m <= m*u+2B_port,
    u <= E_15(K) <= E_15^infinity.                         (21)

The last variational inequality uses a product of zero-excitation cell ground states, for which every port expectation is zero. Put gamma=delta-1/(K-4). For a wrong cell,

    E_l(K)>=l e_15+l gamma>=l*u/15+l gamma.

On L=15m with k wrong cells, cut at most 2k ports touching them. The remaining selected cells form at most k open runs. Applying (21) to these runs and the preceding inequality to the wrong cells gives

    E_wrong-S_m >= gamma sum_wrong l-3k B_port-2B_port
                 >= k(gamma-5B_port).                     (22)

This also covers all-wrong rings. Marker-free rings are excluded by (18) and e_15<0. At (14), gamma>2^-17>5sqrt(2)lambda. Thus F3 transfers to C84 at nonzero port, with all illegal components retained in every full marker sector.

[Numerical] `check_frozen.py` compares full 576-dimensional sectors with lengths (2,2) and (1,3), without a charge cut. It uses the b=2 midpoint selector, K=100, D=4, g=.03, t_d=.001. For lambda=0,.001,.003 the energy differences are .5833105881,.5833106366,.5833110243, exceeding respective bounds .1744288125,.1602866769,.1320024057. Maximum eigenpair residual is below 2.0e-11. This tests the bookkeeping at resolved parameters, not the tiny member (14).

[Scope] Equation (22) is frozen selection. Mobile contact factors share a gapless charge fluid, so mobile selection still needs a uniform common-bulk-subtracted dressed-defect comparison. F3 is not promoted to that stronger theorem.

## 7. F4 and a local grammar budget

[Proved, applying C82/C83] Write H=kappa N_bad+B, with kappa=K frozen and kappa=JK mobile. The bond mismatch projectors commute and have integer eigenvalues. For local A supported on X, resolve

    A_m=(1/(2pi)) integral_0^(2pi)
        exp(i theta N_bad) A exp(-i theta N_bad) exp(-im theta) dtheta.

Only penalty projectors meeting X enter. Thus A_m is supported on X enlarged by one site at each end, ||A_m||<=||A||, and [N_bad,A_m]=m A_m with |m|<=|X|+1. For a two-site term there are at most three participating projectors. The extra color does not change this support bound.

Let B_0 denote the m=0 component in this section, B_off=B-B_0, and nu=sup_i sum_(X contains i)||B_X||. Set

    S_1=sum_(m!=0) B_m/(kappa m),
    F_2=[S_1,B_0]+[S_1,B_off]/2,
    S_2=sum_(m!=0) (F_2)_m/(kappa m).

S_1,S_2 are anti-Hermitian and [S_1,kappa N_bad]=-B_off. Expanding exp(ad(S_1+S_2)) through second order gives

    U H U^dagger=kappa N_bad+B_0+(F_2)_0+R_2,
    [N_bad,B_0+(F_2)_0]=0,
    ||R_2||_local,exponential <= C_3 nu^3/kappa^2.          (23)

On P, B_m P=0 for m<0 and P[S_1,B_0]P=0, so

    P(B_0+(F_2)_0)P
        =PBP-sum_(m>0) PB_(-m)B_mP/(kappa m),              (24)

which is (7). The explicit compression (8) is stronger than a generic norm bound at this order.

For completeness, nonzero nested commutators involve connected overlapping intervals. At each next step, the sum of generator norms meeting an interval I is bounded by |I| times its interaction norm, and the commutator contributes a factor two. The support grows linearly with the order, giving a bound k! C^k. BCH's 1/k! cancels the factorial. For sufficiently small nu/kappa, the exact finite-generator conjugation has a geometrically bounded tail in a weaker exponential interaction norm. This proves (23) with constants independent of volume. Dividing nonzero integer Fourier components recursively proves for each fixed r

    U_r H U_r^dagger=kappa N_bad+D_r+R_r,
    [D_r,N_bad]=0,
    ||R_r||_local,exponential<=C_r nu(nu/kappa)^r.           (25)

No convergence as r tends to infinity and no global many-body band gap are asserted. For a legal vector phi, P[S,p_i]P=0, so its dressed local mismatch expectation starts at second order and is O((nu/kappa)^2). Its summed mismatch can still be extensive. Identifying actual ground states with this dressed space needs stability under R_r.

[Proved] In the ordinary local interaction norm the explicit diagonal representative (8) has strength at most

    (4t^2+3J^2+6J^2 lambda^2)/(JK).                        (26)

The first counter has range two and coefficient at most two; A and Z have range three and coefficient at most one. Counting terms meeting a site gives the displayed factors. This controls the full counters, including their scalar and bit-dependent pieces, without pretending the first term carries a port factor.

[Numerical] `check_sw.py` constructs the literal 78-dimensional invariant closure of L_0 E V_0, including nonzero port. For K=40,80,160,320, the values of K^2||R_2|| are 8.36368,8.33570,8.32186,8.31497. The second-order compression agrees with (24) below 6.2e-21; the unitary residual is below 6.7e-16. These are sign and scaling checks of the finite-order construction.

## 8. Small exact diagonalizations through vacancies

[Numerical] `check_lift.py` exhausts each invariant closure of the literal local rules, without discarding states by mismatch count. All main tests use J=t=1, D=4, Delta=4-2cos(pi/4), g=.07, t_d=.003, zero selector/field/charge chemical potential, and phi_A=pi/4. These are diagnostic parameters, not the member (14). The first case has two p=3 cells, each with active bits (1,0), and one vacancy, with lambda=0. The last two have a p=2 cell, nonzero bit and total charge Q=2, lambda=.1, and respectively one or two vacancies.

| Test | Dimension | K=100 illegal weight | K=1000 | K=10000 |
| :-- | --: | --: | --: | --: |
| Two cells, uncolored | 1296 | .0239352 | .0230743 | .0229959 |
| Two cells, three colors | 3348 | 8.30676e-5 | 8.39127e-7 | 8.39976e-9 |
| Nonzero port, one vacancy | 368 | 2.30364e-7 | 2.67569e-9 | 2.71784e-11 |
| Nonzero port, two vacancies | 1945 | 1.31616e-7 | 1.52114e-9 | 1.54426e-11 |

The uncolored nonzero-port counterparts approach 2.51578e-5 and 2.64428e-5 at K=10000. In all three cases the legal colored matrix agrees entry for entry with the uncolored matrix restricted to squeezed-legal words. Maximum residual over these runs is below 2.7e-11.

[Numerical] A separate one-vacancy Q=2 test uses J=.01,t=1 and K=10^4,10^5,10^6, so kappa=100,1000,10000. The illegal weights are 2.70365e-11,3.26038e-13,3.31954e-15. The physical mismatch expectation is 1.14906e-4,1.16418e-6,1.16557e-8. This distinguishes illegal material weight from wrong-color weight and uses the actual penalty JK. These ground states belong to the stated invariant closures, not an exhaustive claim about arbitrary rings or charge sectors.

## 9. What is certified, and the conditional infrared count

[Proved, sufficient local hierarchy] Fix t>0, a positive J/t, and a material-energy tolerance tau>0. Let A2=4t^2+3J^2+6J^2 lambda^2. The unpenalized local strength admits, for example, the bound nu<=2^13 J+4t+2|mu| at (14). Choose K finite so that

    K>=2^2560,
    nu/(JK) is within the fixed-order conjugation domain,
    K>2 A2/(J^2 tau),
    K>sqrt(2 C_3 nu^3/(J^3 tau)).                          (27)

Then (26) and (23) are each below J tau/2. Increasing K improves the frozen certificate. The material quadratic tolerance can be chosen as

    tau_band=lambda^2 E_unit 2^-35/(128 n tau_W).            (28)

For a half-occupancy reference, include the positive mean-contact factor alpha_contact=1/2-2/pi^2 by taking tau smaller than alpha_contact tau_band. The selector, massive-sector and infrared stability tolerances impose additional upper limits. In particular K=2^2560 alone is not a certificate of the mobile model for arbitrarily small J. Nor is JK much greater than t enough. The color field requires control of t^2/(JK) relative to J times the very small logical scale.

[Assumed, inherited reference dictionary] C83's charge reference and C81's harmonic dictionary apply to (5) because the coloring is a rotation-equivariant bijection. The reference charge fluid has one compact boson and mean contact alpha_contact at half occupancy. Positive compressibility permits tuning the density at finite J within the matching assumption. The physical root and cyclic zero-mode gluing are retained.

[Proved conditional on that dictionary] With b=15, a payload vertex of cell momentum q and reference dimension Delta_pay has harmonic branches

    Q_h=q+2pi ell,
    k_phys=Q_h/30 modulo 2pi,
    Delta_ref=Delta_pay+[Q_h/(30pi)]^2.

A uniform term requires Q_h=60pi h. For h!=0 its charge dimension is at least four. For h=0, (15) excludes nonforward neutral quadratic and quartic resonances. Higher neutral fermion vertices have reference dimension at least three. Exact Q and N forbid anomalous payload pairing and charge-number-changing terms. For vertices below dimension three the inherited mismatch estimate gives

    delta_phys >= min(delta_internal/30, pi(2-sqrt(3)))
               >1.0837704e-7,
    delta_internal=min(delta_2,delta_4).

The first marker harmonic has wavevector pi/15 and dimension 1/225. It is an oscillating observable, not a uniform relevant Hamiltonian term. All integer branches are retained. The new density-contact field in (9) must be included in the same matching and stability problem.

[Assumed] A microscopic c=14 conclusion still requires massive-sector matching in local coefficient norms, low-frequency matching of the gapless vacancy fluid at fixed positive J, a uniform mobile dressed-defect comparison, multichannel stability of the interacting retarded theory, and enough current/stress/thermal control to identify its conformal modes. No positive numerical window for all these inputs is computed here. exp(i pi Q/2) acts as payload parity on the code but can be +/-i on odd-Q auxiliary states, so a microscopic Z2 grading cannot be imported unchanged.

[Heuristic, conditional on those inputs] The finite member has thirteen payload modes and one vacancy mode,

    n+1=14 > d-1=13, d=14.                                 (29)

Thus the eleven-letter compression does not reduce the colored lift below C83's d=14. Saving E_1 requires an additional vacancy state and an additional material site per fixed flavor count. It does give a separately certified length-15 member. There is no measured central charge or unconditional counterexample here.

## 10. Reproduction and provenance

[Definition] Run from work_C85 with Python, NumPy and SciPy only.

    python3 -B run_certificate.py
    python3 -B selector_certificate.py
    python3 -B check_lift.py
    python3 -B check_virtual.py
    python3 -B check_energy.py
    python3 -B check_frozen.py
    python3 -B check_sw.py

The certificate sources are local copies of C84's C80/C75/C33 pipeline, with its active/inactive endpoint convention unchanged. The driver changes n to 13 and loop exponent to 66. The final finite budget enlarges the column and port constants and uses (14). The independent literal lift starts from C83's table implementation, removes E_1 and endpoint conversion, and introduces the three-state table. `bond11.py`, `parent_shift.py`, and `bond18.py` are local C84 sources used for the material comparison. No parent output was relabeled as a new certificate.

[Scope] All files and computations stay in this folder. No network commands, cluster jobs, or subagents were used. No file outside this folder was modified. Local `STATE.md` records the handoff because the project-level memory lies outside the authorized write scope.
