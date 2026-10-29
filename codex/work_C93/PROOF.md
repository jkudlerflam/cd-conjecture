# C93. A block grading and fixed-cutoff bounds for the d=13 chain

## 1. Outcome and conventions

[Proved] There is an integer block charge implementing payload parity for the full finite-K Hamiltonian lifted to the marker-root cover in sectors with an even number of markers. It also works on an open material chain. This is stronger than a grading of the legally compressed operator. Section 2 gives the construction and its boundary qualifications.

[Proved] A Gaussian filter at a fixed cutoff has a remainder controlled at any prescribed spatial interaction exponent. Section 4 supplies finite constants and an explicit positive J/t schedule. It corrects the shrinking-exponent problem without replacing the filtered contact by a scalar.

[Proved] The leading color field can be absorbed into an onsite logical chemical-potential matrix. Absorption is a change of reference, not an extra microscopic tuning parameter. Section 5 bounds the leading coefficient and identifies the terms that absorption leaves behind.

[Assumed, unresolved] These constructions do not establish the dressed mobile material threshold or the massive influence action requested in the question. Controlled free-sea statements about contact fluctuations are proved below, but their transfer to the material-label contact coupled to the clock is missing. The action in section 7 is consequently a conditional reference action, not a certified matched action for the literal finite-J chain. No value of the dressed gap or of its matching error is asserted.

[Definition] A tag applies until the next tag or section. Proved includes explicitly attributed algebraic and finite-certificate inputs with their original scope. Numerical means a finite floating-point check. Heuristic identifies a continuum interpretation. Assumed identifies an unproved input. Logarithms are natural. The material cell has thirteen active bits, fourteen data letters and its terminal marker, hence b=15. The symbol K denotes the grammar penalty, and K_gap denotes the charge stiffness.

[Proved, provenance] The inputs are C91, C48, C88 and SEPARATE_JW.md. C91 was read in full, including its check source and small output. C48 and C88's notes and source files were read. Their large certificates were not rerun. One exact scalar, the metric trace in C88's finite_budget_exact.json, was extracted for the bound in section 5. All new execution and writes stay in work_C93. No network commands or delegated agents were used.

## 2. An integer charge on rooted blocks

[Definition] Label the markers in their conserved order. Cell j contains the data letters after marker j-1 and the terminal marker j. Vacancies carry none of the charges below. Their assignment to gap factors has no effect on the definitions. Put

    eta_j=(-1)^j,
    E_j=sum_(letters a in cell j) e(a),
    Q_j=sum_(letters a in cell j) q(a),
    S=sum_j eta_j E_j,
    Nhat_j=(Q_j-eta_j E_j)/2,       Nhat=sum_j Nhat_j.       (1)

Here e(a)=1 for D_A,D_B,X_A,X_B,M_A,M_B and zero otherwise. This is an actual count, not a projector onto the single-carrier legal subspace. Multiple carriers in an illegal word are included. The charge q is C88's microscopic charge.

[Proved, integrality] For every microscopic letter, q(a)-e(a) is even. Thus Nhat_j is an integer on every cell word, legal or illegal. More explicitly, its single-letter contributions are

| Letter | eta=+1 | eta=-1 |
| :-- | --: | --: |
| L_0,R_0,E,V_0,V_1 | 0 | 0 |
| L_1,R_1 | 1 | 1 |
| D_A,X_A,M_A | -1 | 0 |
| D_B,X_B,M_B | 0 | 1 |

Negative integer auxiliary charges are allowed. Neither a square root of microscopic parity nor a half-integer auxiliary number is used.

[Proved, conservation on the full rooted sector] Every term preserves Q. To verify S, classify the literal off-diagonal rules.

1. Vacuum head motion changes no e count. Rail motion, switches and A/B conversion preserve E_j in the cell on which they act. The conversion D_A R_1 to D_B R_0 also preserves Q_j because -1+2=1+0.
2. Terminal moves D_s E or X_s E to L_0 M_s move a charged letter to the terminal marker of the same cell. They preserve E_j and Q_j, even with additional illegal charged letters elsewhere in that cell.
3. A port E R_0 to M_A D_B, or to M_B D_A, increases E_j and E_(j+1) by one each. Its change in S is eta_j+eta_(j+1)=0. The inverse has the opposite changes and also preserves S. The left and right microscopic charges change by opposite units.
4. A physical vacancy hop preserves the squeezed sequence of material letters. Each letter retains its cell membership, including a marker moving across vacancies. Therefore every E_j and Q_j is unchanged. The color chosen by the hop is immaterial to (1).
5. All diagonal terms, including every grammar penalty, commute with (1).

It follows term by term that

    [H,Q]=[H,S]=[H,Nhat]=0.                               (2)

The argument uses no legal projection, no perturbation expansion, and no restriction to one charged letter per cell. It holds at finite K, nonzero g,t_d,lambda and arbitrary vacancy placements. It also applies to the local penalty Fourier components used by C88. Those components, and all their nested commutators, conserve Nhat.

[Proved, payload parity] On the frozen vacuum code E_j=0 and Q_j=2 N_payload,j. Hence

    Nhat_j P_code=N_payload,j P_code,
    P_j=(-1)^Nhat_j,
    product_j P_j|code=(-1)^N_payload.                    (3)

This does not contradict C91. The table depends on the parity of the cell label, which cannot be determined from a microscopic letter alone. Grouping two neighboring cells makes the convention periodic on the rooted lattice. It is a block representation, not a finite-depth transformation local in physical distance.

[Proved, boundary and quotient qualifications] On an even-cell ring the alternating signs are consistent across the seam. On an odd-cell ring the seam joins equal signs. A port there changes Nhat by one, so (2) fails for that one seam term. Cutting its port repairs the grading at an operator cost at most sqrt(2) J |lambda|. This is a bounded boundary modification, not equality of the original finite-ring action. The marker-free sector is outside the construction.

Changing the root by one cell sends S to -S and exchanges Nhat=(Q-S)/2 with Nhat'=(Q+S)/2. Both charges are conserved in the even rooted sector. Their parities agree on the code, where S=0, but need not agree on arbitrary auxiliary states. The physical cyclic quotient must retain this frame change. Neither a periodic fermion trace with an arbitrary twist nor a projection onto one frame is automatically the physical trace. At the level of pressure, removing a bounded seam changes log Z by at most beta times the seam norm. This does not prove equality for charged sources or finite-size levels.

[Proved, code-connected sector] Every vacuum-code state has S=0, and (2) implies that every virtual trajectory starting from it remains in S=0. On this entire invariant sector, Nhat=Nhat'=Q/2 is integer and the total parity is independent of the root frame. Thus the total payload parity does descend to the physical quotient for these code-connected amplitudes. The assertion that the full thermal or ground-state problem is exhausted by that sector requires an energetic comparison. It does not follow from conservation. Local Jordan-Wigner representatives still require their boundary twists.

## 3. The graded embedding and the exact zero-port subspace

[Proved] Choose the cell code as the frozen clock vectors labelled by the thirteen bits. This subspace reduces Nhat_j. Its orthogonal complement therefore carries integer Nhat_j eigenvalues as well. At fixed material length 15 the unrestricted material cell dimension is

    D_mat=3*8^14=13194139533312,                          (4)

because each of the fourteen nonmarker sites has eight choices and the terminal marker has three. Formula (4) counts illegal material configurations too. The gaps, and their color words away from legality, are separate infinite factors and are not counted in D_mat.

Define the algebraic isometry on the material factor by

    J_j |u,code> = |u>_log tensor Omega,
    J_j |s,complement> = |0>_log tensor |s>_aux.

Assign charge zero to Omega and charge Nhat_j(s) to |s>_aux. Then

    Nhat_j = N_log,j + N_aux,j

on the image. There is no spectral assertion in this isometry. In particular it does not make the frozen code invariant under the mobile Hamiltonian.

[Proved, ordinary CAR] With p_j=(-1)^N_log,j, pi_j=(-1)^N_aux,j and P_j=p_j pi_j, put

    S_j=product_(k<j) P_k,
    psi_(j,a)=S_j c_(j,a),
    Xi_j(O_r)=S_j^r p_j^r O_r,
    pi_j O_r pi_j=(-1)^r O_r.                            (5)

This is the single-family specialization of SEPARATE_JW.md. The psi obey ordinary CAR. Xi_j is an isometric star homomorphism into a graded auxiliary algebra. Odd Xi coefficients anticommute with every psi, and odd coefficients on distinct auxiliary sites anticommute. To check isometry, decompose into eigenspaces of the commuting involution S_j p_j in the commutant of the original auxiliary algebra. Xi is either the identity on that algebra or conjugation by pi_j in each block. Each block has the original operator norm.

For a term neutral under Nhat, substituting (5) leaves only parity strings between its endpoints. Every string strictly to the left occurs to even power and cancels. Interval support and operator norm are preserved. The conserved gap count supplies a second integer charge if a separate charge-field representation is constructed; (1) supplies the missing payload charge in that two-family construction. It does not supply a length qubit for the unbounded gap space.

[Proved, reference extension] A logical hopping must be extended as psi_j^dagger T psi_(j+1). In the original tensor representation it is

    c_j^dagger p_j pi_j T c_(j+1).                       (6)

The auxiliary parity pi_j cannot be omitted. For any neutral code operator O_X, the extension defined by (5) is unitarily gauge-equivalent to O_X in every fixed auxiliary parity sector, with gauge sign product_(a<=k<j) pi_k at site j. It is therefore isometric. For a locally code-reducing B_X,

    K_X=J_X B_X J_X^dagger-E_JW(P_X B_X P_X),
    ||K_X||<=2||B_X||,
    K_X Omega_X=Omega_X^dagger K_X=0.                    (7)

Equation (7) is available if local block reduction has been proved. Applying it to the undressed mobile Hamiltonian would assume the desired reduction.

[Proved, coefficient bound with a retained bath] On r cells, logical matrix units and their normal-order expansions give

    sum_I 4^degree(I) ||O_(X,I)|| <=128^(13r)||O_X||.     (8)

The coefficients in (8) may act on the auxiliary material space and the retained gaps. Their norm has no dimension factor from either. This is a coefficient bound on operators at one time. It does not turn noncommuting retained-gap operators into scalar Grassmann kernels or prove a connected expansion with their dynamics.

[Proved, an exact zero-port factorization with a larger retained bath] In the legal pin-free length-15 sector at K=infinity, lambda=f=0, the subspace with no charged letters is invariant. Decode each material word into its thirteen bits and one of fifteen vacuum-clock positions. Every vacuum head move has amplitude -1 independent of the decoded bits. Thus, retaining all clock positions and all gaps,

    H_no-carrier = H_gap-clock tensor I_payload.          (9)

The conserved bits are not entangled with this bath by the Hamiltonian in (9). However H_gap-clock is an interacting gap and clock Hamiltonian, not a product of on-cell massive propagators. Retaining it all removes frozen-code leakage only by retaining the clocks themselves. It does not integrate them out. Diagonalizing this entire bath in a finite volume gives no uniform local dressed representation.

[Proved, why the bare C48 substitution fails] C91's single-contact vacuum-clock matrix element is

    ||(1-P_v) k_0 v||
       =sqrt(v_0^2+v_1^2-4v_0^2 v_1^2)
       =0.1507149408065...,
    v_l=sqrt(1/8) sin(pi(l+1)/16).                       (10)

A retained gap configuration realizes that contact. The leakage is order J, independently of lambda and J/t. Meanwhile the averaged vacuum-clock splitting is

    Delta_clock,ref=2 J alpha_c [cos(pi/16)-cos(2pi/16)],
    alpha_c=1/2-2/pi^2.                                  (11)

The coefficient ratio in (10) divided by Delta_clock,ref/J is order one, not a small parameter. A stronger ultraviolet remainder bound leaves this leading leakage in the filtered operator. There is no positive operator lower bound for the bare contact-weighted excitation energy on the entire gap space, as shown by C91's all-positive-gap vector. A threshold relative to the interacting ground state could still exist.

## 4. A fixed cutoff with a prescribed strong spatial norm

[Proved, scope] This section concerns C91's bounded legal compression

    H=t H0+J V,
    H0=-sum_j(E_j^+ E_(j+1)^-+h.c.),
    V=sum_j v_j,
    ||v_j||<=M0=8192,       supp(v_j) subset [j-2,j+2].   (12)

It does not eliminate the finite-K grammar remainder. All distances in this section are material-site distances. Set s=J/t. The interaction norm is C91's ||Phi||_(1,a), with weight (1+diam X) exp(a diam X).

[Definition] Fix 0<eta<1 and a>=0. The following constants are independent of volume, s and gap occupation. Fix the energy-transfer cutoff delta=1/2. Define

    mu=4a+4,
    v_mu=64 exp(mu),
    B=4v_mu/mu+17,
    T=ceil(max(1,
        16(a B+1)/delta^2,
        sqrt[(8/delta^2) log(98 M0 B^2 exp(2a)/eta)],
        log[25088 M0 B^2/(delta eta)]/(v_mu-2))),
    R=ceil(4v_mu T/mu+16).                               (13)

The letter B in (13) is a numerical constant, not C48's interaction. These constants can be extremely large. Their size is not suppressed in the conclusions.

[Proved, weighted locality] The overlapping-bond path estimate in C91, now weighted by exp(mu distance), gives

    ||[tau_u^I(v_j),H_boundary]||
       <=128 M0 exp[-mu(R-4)+v_mu |u|].                 (14)

Indeed the path sum is bounded by sum_(n>=d)(12|u|)^n/n! <= exp[-mu d+12 exp(mu)|u|]. The two boundaries and the initial five-site support fit in the displayed prefactor and the larger velocity 64 exp(mu). All local hopping operators are bounded. Thus the proof is independent of the infinite gap-factor dimension.

Use C91's Gaussian g_delta and anti-Hermitian primitive F_delta, but truncate at the T in (13):

    g_delta(u)=delta exp(-delta^2 u^2/4)/(2sqrt(pi)),
    F_delta(u)=-(i/2)sign(u) erfc(delta|u|/2),
    S_j=integral_(-T)^T F_delta(u) tau_u^I(v_j) du,
    w_j=integral_(-T)^T g_delta(u) tau_u^I(v_j) du.

Then ||S_j||<=2M0/delta, ||w_j||<=M0, and integration by parts gives

    v_j+[S_j,H0]=w_j+e_j,
    ||e_j||<=M0 exp(-delta^2 T^2/4)
       +(256M0/delta)exp[-mu(R-4)+v_mu T].              (15)

The support of e_j has diameter at most 2R+2. Summing incidence and first moment yields

    ||sum_j e_j||_(1,a)
      <=49 M0 R^2 exp[a(2R+2)]
           {exp(-delta^2 T^2/4)
             +(256/delta) exp[-mu(R-4)+v_mu T]}
      <=eta.                                           (16)

For completeness, R<=BT. The endpoint part is bounded by

    49 M0 B^2 exp(2a)
       exp[2T+2aBT-delta^2 T^2/4]
      <=49 M0 B^2 exp(2a) exp(-delta^2 T^2/8)<=eta/2.

Here T^2<=exp(2T), and the second entry in (13) bounds the positive linear exponent. For the boundary part,

    -(mu-2a)R+v_mu T+4mu+2a <=-v_mu T,

using mu-2a>=mu/2 and R>=4v_mu T/mu+16. It is at most

    (12544 M0 B^2/delta) exp[-(v_mu-2)T]<=eta/2.

This establishes (16) with all errors retained. Increasing R alone at C91's original unweighted velocity would not establish (16).

[Proved, nonlinear tail] Define

    D_a=(20 M0 R/delta) exp(6aR),
    C_N=49 R^2(44 D_a^2+12 M0 D_a).

The same connected commutator count as C91 gives, if 0<s<=1/(2D_a),

    U=exp(s sum_j S_j),
    U H U^dagger=t H0+J W_delta,T+R_UV,
    ||R_UV||_(1,a)/J <=eta+C_N s.                      (17)

To see the power of s, substitute q_a=D_a s in C91's bound
49R^2[44q_a^2+12sM0q_a] for the nonlinear remainder in units of t, then divide by J/t=s. No extensive operator norm or inverse charge gap enters. Taking

    0<s<=min(1/2,1/(2D_a),eta/C_N)                     (18)

gives ||R_UV||_(1,a)<=2eta J. The cutoff, time range and spatial range are fixed before choosing s. Unlike delta=s^z, this schedule controls a fixed strong spatial exponent.

[Proved, explicit coefficient conversion] For output spatial exponent alpha_out=1, put

    L=8(1+log32+10)=115.7258872224...,
    a_13=L+13log128+1=179.8022806534...,
    C_blk=15 exp(a_13),      a=256.                     (19)

A material interval of diameter d meets at most d+1 cells. For the per-cell incidence sum, sum over the fifteen material sites in the cell. Equation (8) then implies

    v_L,operator-coefficients(R_UV)
       <=C_blk ||R_UV||_(1,256).                        (20)

The extra +1 in a_13 supplies a support reserve. The same choice a=256 also exceeds the hypothetical fourteen-mode exponent L+14log128+1=184.6543109173..., if a valid finite logical charge-mode representation were supplied.

Let W_ref=J alpha_c lambda^2 E_unit, with E_unit from section 5. For any 0<epsilon<1 choose

    eta=epsilon alpha_c lambda^2 E_unit/(4 C_blk),       (21)

then (13),(18) give v_L(R_UV)<=epsilon W_ref/2. This is a finite, volume-independent bound in the requested strength of spatial coefficient norm for this ultraviolet remainder, with retained gap operators still in the coefficients. It is not a bound on the retained W_delta,T, on a temporal influence kernel, or on a full matched action. No numerical representation of exp(6aR) is needed to specify the finite schedule.

## 5. Absorbing the finite-K color fields

[Proved, imported coefficient and notation] Use C88's constants

    alpha_end=sqrt(1/8)sin(pi/16),
    d_0=1, d_1=1004, d_(j+1)=1004d_j-d_(j-1),
    beta_c=1/d_15,
    g=2^-272, lambda=2^-608, t_d=beta_c 2^-66,
    E_unit=alpha_end^4 beta_c^2 g^2,
    f=lambda^2 E_unit 2^-120.                            (22)

The subscript on beta_c distinguishes the carrier constant from inverse temperature. The leading projected color field in the frozen-clock and free-sea reference is

    D_col=diag(d_r),
    d_r=t^2 w_r/(JK),
    c_r=sum_(l=0)^r (1/8)sin^2(pi(l+1)/16),
    w_r=(1-alpha_c)c_r,             0<=r<=12.            (23)

This is a reference projection of C88's order-1/K operator, not the full mobile self-energy. Its sign is positive. Its extreme weights are 0.00334284159219 and 0.68643707665319.

[Proved, absorption and its meaning] In the convention H_onsite=-psi^dagger mu_log psi, define

    mu_log,new=mu_log,old-D_col.                         (24)

Then the displayed color field is part of the quadratic Hamiltonian, with no error incurred by this reassignment. The matrix need not commute with the hopping matrix, and its diagonal orbital entries are not thirteen separately conserved microscopic charges. Equation (24) does not authorize adding thirteen counterterms to the original bond. In a band basis it is generally not diagonal at all momenta. Its effect on roots and velocities must be included in that basis change.

[Proved, a finite bound preserving the reference coefficient budget] Let tau_W be the dyadic metric trace in C88's exact certificate. Direct integer comparison of that scalar gives tau_W<2^118. The inherited conservative conversion factor for a physical onsite perturbation is at most 16n tau_W, with n=13. Thus

    16n tau_W ||D_col||/W_ref
       <2^126/[K s^2 alpha_c lambda^2 E_unit].           (25)

C91's elementary bound alpha_c lambda^2 E_unit>2^-2080 follows from alpha_c>1/4, alpha_end^4>=2^-18 and d_15<2^150. Therefore the explicit choice

    K>=max(2^2560,2^2243 s^-2)                           (26)

makes (25) less than 2^-37. This is one quarter of C88's 2^-35 additional coefficient allowance. It permits the leading color field to enter the reference without exhausting that allowance. Other mobile corrections, longer-range terms and the dressed velocities are not certified by (26). If the field is absorbed without a smallness bound, the old root certificate must be replaced by a new one. Merely renaming it a chemical potential does not preserve thirteen Fermi channels.

[Proved, what remains after absorption] Before taking the charge-sea average, the bare vacuum-clock projection of the bit-dependent part is

    (t^2/(JK)) sum_(c,r) c_r (1-A_(c,r+1)) n_(c,r),
    A_j=|g_j=0><g_j=0|.                                 (27)

The position r+1 is in zero-based data coordinates: when the clock lies at l<=r, bit r occupies that physical material position in the R tail. Subtracting (23) leaves

    R_col,bit=(t^2/(JK)) sum_(c,r)
                   c_r(alpha_c-A_(c,r+1)) n_(c,r).       (28)

The per-cell operator norm of (28) is at most 13t^2/(JK). With field weight 4 its per-cell logical coefficient sum is at most 208t^2/(JK). The centered scalar contact terms have a conservative bound 15t^2/(JK). Hence the whole centered bare-clock color projection has v_L at most

    223 exp(L) t^2/(JK).                                (29)

For example, adding K>=892 exp(L)/(epsilon alpha_c lambda^2 E_unit s^2) to (26) puts (29) below epsilon W_ref/4. This bound keeps the gap operators external; no irrelevance assertion is needed for this order-1/K contribution. Clock off-diagonal terms, the other grammar counters and the higher grammar remainder remain. C88's unspecified strong-norm C_3 and its coordinate transfer must be supplied to bound those at the same schedule.

## 6. What can be proved about replacing contact by its mean

[Proved, free-sea mean and derivative] For physical free fermions at density rho, the equal-time two-point function is G_r=sin(pi rho r)/(pi r), with G_0=rho. Wick's rule gives

    b(rho)=<n_i n_(i+1)>=rho^2-sin^2(pi rho)/pi^2,
    alpha(rho)=b(rho)/rho
              =rho-sin^2(pi rho)/(pi^2 rho).            (30)

The second formula is the contact count per occupied particle. At half filling it agrees with alpha_c. For mean vacancies per material letter u=rho^-1-1,

    alpha(1)=alpha_c,
    d alpha/d rho|_(rho=1/2)=1+4/pi^2,
    d alpha/d u|_(u=1)=-1/4-1/pi^2.                     (31)

Here alpha(1) means alpha evaluated at u=1. Equation (31) is an exact derivative of a reference expectation; it is not an operator Taylor expansion of A_j in the interacting chain.

[Proved, a controlled block quantity in the sea] Let B_i=n_i n_(i+1), and define the physical block observable

    C_ell=(1/(rho ell)) sum_(i=1)^ell B_i.

Its mean is alpha(rho). At rho=1/2,

    <(C_ell-alpha_c)^2> <=44/ell.                        (32)

To prove it, write the four-site expectation at separation r>=2 as the determinant of the four-by-four matrix G on {0,1,r,r+1}. Subtracting <B_0><B_r> cancels the permutations that stay within each two-site set. Every remaining permutation has at least two cross-set entries. Each cross entry is bounded by 1/[pi(r-1)], and each other entry is bounded by one. There are at most 24 permutations, so

    |Cov(B_0,B_r)|<=24/[pi^2(r-1)^2],       r>=2.

For r=0,1 use |Cov|<=1. Expanding the variance of the sum and extending its tail gives

    Var(ell^-1 sum_i B_i)
       <=ell^-1[1+2+48 pi^-2 sum_(s>=1)s^-2]
       =11/ell.

Division by rho^2=1/4 proves (32). This is a state-dependent mean-square bound. There is no full retained-space operator-norm estimate tending to zero with ell.

[Proved, low-patch compression] A different controlled replacement is available in the free charge Hilbert space. On a finite physical ring choose a set of low momentum modes around the two Fermi points, and freeze every other mode in its sea occupation. Let P_L denote this subspace, nu_L its number of low modes divided by the physical ring length, and b_L=<B_i> its finite-ring reference expectation. Then

    ||P_L(B_i-b_L)P_L||<=14 nu_L,                       (33)

provided the reference sea is one of the low-mode configurations with the fixed high occupations. The same bound divided by rho holds for C_ell-b_L/rho, uniformly in ell. In the thermodynamic limit b_L approaches b(rho) in (30).

For a proof decompose c_i=d_i+h_i into low and high momenta. Then ||d_i||=sqrt(nu_L), ||h_i||<=1. Since B_i=c_i^dagger c_(i+1)^dagger c_(i+1)c_i, its expansion has sixteen terms. Terms with an odd number of high operators vanish on the fixed high Slater state. The zero-low term is a scalar. There are six terms with two low operators, each bounded by nu_L, and one four-low term bounded by nu_L^2. The nonscalar part D has norm at most 6nu_L+nu_L^2<=7nu_L. Subtracting its expectation in the low sea changes its norm by at most another ||D||, proving (33). For two patches of half-width Lambda, nu_L tends to 2Lambda/pi in the thermodynamic limit.

[Proved, why (32) and (33) do not finish matching] The exact global identity is that the total number of material contacts equals sum_i B_i. A particular A_j is the contact after the j-th material particle. Its physical position is an operator depending on the cumulative particle number. Thus (32) and (33) control an unweighted physical contact count, not A_j h_j with a dynamical clock and material-label-dependent h_j.

Moreover P_L is not invariant under the interaction. Projecting onto it is not an exact change of representation of the original chain. A sharp momentum cutoff also has algebraic spatial tails, so (33) alone cannot establish C48's exponential spatial coefficient norm. These facts prevent promotion of the patch bound to a dressed mobile theorem. The finite-cutoff unitary of section 4 and the projection P_L are different operations.

[Heuristic, the non-irrelevant term] If the mobile charge field admits the local hydrodynamic expansion delta u=-partial_x phi/pi in material coordinate x, (31) predicts

    A_smooth=alpha_c+[(1/4+1/pi^2)/pi] partial_x phi+... . (34)

A smooth payload density has scaling dimension one, and partial_x phi also has dimension one. Their product is dimension two. Thus a contact-weighted payload chemical potential can generate a marginal charge-payload forward term. A scalar times partial_x phi is a density source; a normal-ordered payload kinetic energy density times partial_x phi has dimension three. The dimensions differ, so the whole contact deviation cannot be declared irrelevant. In cell coordinate y=x/15, the corresponding linear expansion coefficient is

    g_A=(1/4+1/pi^2)/(15pi).                             (35)

Equations (34)-(35) assume the field/operator dictionary and a controlled local derivative expansion. They are not bounds on its remainder. Oscillatory harmonics and the root quotient also require the momentum dictionary, not just this smooth term.

## 7. Conditional action and the missing quantitative input

[Proved, explicit reference matrices imported from C88] The carrier matrix H_c has ordered vertices M,X_13,...,X_0,D_0,...,D_13. Its diagonal is 1004 on M,X and 4 on D, its path entries are -1 and H_(M,D_13)=-i t_d, with the conjugate reverse entry. Set

    R_c=(H_c tensor I+I tensor H_c)^-1,
    x=R_c|M,D_0>,
    U_r=projection_(first carrier at D_r) x,
    V_r=projection_(second carrier at D_r) x,
    C_rs=U_r^dagger R_c V_s,
    A_rs=U_r^dagger R_c U_s+V_r^dagger R_c V_s,
    T_f=-2 alpha_end^4 g^2 Re C,
    M_f=-2i alpha_end^4 g^2 Im A.                        (36)

U_r,V_r are vectors on the full product-carrier space, with the unwanted entries set to zero. The reference physical Bloch matrix, including (23), is

    h_ref(k)=J alpha_c lambda^2
           [T_f exp(ik)+T_f^dagger exp(-ik)
                        +M_f+E_unit 2^-120 I_13]
              +D_col.                                  (37)

Equations (22),(23),(36) specify every entry of this reference. The unshifted C88 coefficient has 26 simple roots and the recorded separation margins. Under the leading color allowance in (25)-(26), that perturbation fits inside the imported coefficient budget. No assertion is made that (37) is the full finite-J quadratic self-energy.

[Heuristic, reference action with the required terms visible] A possible local low-energy form in cell coordinates is

    S_ref=beta M e_ref
       +(1/(8pi)) integral dy d tau
          [(15/t)(partial_tau phi)^2+(t/15)(partial_y phi)^2]
       +sum_(a=1)^13 sum_(sigma=R,L) integral dy d tau
          psi_bar_(a,sigma)
            [Z_(a,sigma) partial_tau
                    -i v_signed,(a,sigma) partial_y]psi_(a,sigma)
       +S_forward+S_oscillatory+S_retarded.               (38)

The charge coefficients in (38) are the imported free reference values K_gap=4 and v_c=t/15. At quadratic reference level Z=1 and v_signed is the derivative of the eigenvalue of (37) at the corresponding root. The chemical matrix has already been shifted by (24). S_forward must allow charge-payload density products as well as intra-payload forward terms. Equations (34)-(35) describe one possible contribution before further matching. At finite J the stiffness, velocities, residue, scalar and forward matrix can all be renormalized. Their exact values are not calculated here.

[Assumed, the theorem needed to turn (38) into a matched action] One needs a volume-uniform representation of the zero-port mobile gap-clock bath with a massive material complement, while keeping the gapless charge degree of freedom external. It must provide a threshold Delta_*>0, a product auxiliary heat bound, separately neutral local vertices, and a justified scalar/external-field or mixed-bath connected expansion. It must also provide the actual dressed observables, scalar term, and bounds for the retained contact vertices on a common density-source interval. The unbounded gap source and the physical root quotient cannot be omitted. Neither the exact symmetry (2) nor the static factorization (9) supplies these dynamical data.

[Proved, conditional C48 constants rather than an application] Suppose such a representation has a finite material auxiliary dimension at most D_mat in (4), with heat bound Tr_Q exp(-u A)<=D_mat exp(-u Delta_*) for u>=0, and suppose its retained-field interface really satisfies C48. Then at alpha_out=1 one may set

    C0=D_mat,
    Theta=max(2L,log C0+2)=2L=231.4517744448...,
    tau=Theta/Delta_*,
    kappa=Delta_*/32,
    beta>=4L/Delta_*,
    A_dec=2^22 exp(L)(16386),
    log_10 A_dec=61.0962470892... .                       (39)

These follow by substituting Theta=2L in C48, not by assigning the unknown mobile Delta_* the averaged value (11). If the dressing needs a larger or infinite auxiliary space, its actual heat constants must replace C0; the finite count (4) cannot then be imported unchanged.

Writing k and r for C48's auxiliary coefficient norms, its sufficient hypothesis and conclusion are

    k+r<=Delta_*/(64Theta),
    ||delta S||_LE<=A_dec r,
    ||delta S(omega)-delta S(0)||_LE
                      <=A_dec r |omega|/Delta_* .        (40)

Before the time-moment conversion the thermal baseline is bounded by

    B_* exp[-beta Delta_*/(4Theta)],
    B_*=2^16 exp(L) Delta_*/Theta.                       (41)

Its two-moment bound uses the same factor 16386 as in (39). Scalar terms must be kept. These constants are explicit, but Delta_*, k and r for the actual mobile model are not known. In particular r cannot be replaced by the ultraviolet remainder alone, since that would discard the contact-clock leakage in (10).

## 8. Checks and final status

[Proved, exact finite checks] checks.py exhausts two full 576-state material sectors, with cell lengths (2,2) and (1,3), without legality restrictions. All 2070 directed literal material moves preserve Q,S,Nhat. It checks all 8192 code bit strings and verifies that the odd-cell seam is exceptional. These finite controls supplement the termwise proof, which applies to any fixed material cell lengths. Vacancy hopping conservation follows separately because it leaves the squeezed material word unchanged.

[Numerical] The two-cell, sixteen-dimensional payload/auxiliary Jordan-Wigner check has zero recorded residual for CAR, graded commutation and the parity-twisted hopping extension. Determinantal free-sea variances at block lengths 8,32,128 obey (32). A finite-difference check of (31) has error below 1.7e-11. The logical color weights and C48 constant arithmetic agree with the displayed values. These are small checks, not many-body spectroscopy or a central-charge measurement.

[Proved, limits of the result] C91's parity obstruction is resolved in the rooted bulk by (1). Its shrinking spatial-rate obstruction is resolved for the ultraviolet Hamiltonian remainder by (13)-(21). The leading finite-K color term can be kept in the logical chemical matrix, with the explicit preservation budget (26). C91's contact and dressed-threshold obstruction remains. Equations (32)-(33) give controlled replacement observables in restricted free-sea senses, but neither justifies the needed weighted material-contact replacement. Consequently there is no proved matched action for the finite-J d=13 chain in this report.

[Assumed, next step] The specific missing estimate is for the zero-port mobile gap-clock problem, including its low-frequency material-label contact insertions. A construction proving Delta_*, a compatible auxiliary heat trace and the norms in (40) would make the block grading and fixed-cutoff bounds usable in C48. Declaring all contact fluctuations irrelevant, or declaring all color shifts removable by physical chemical potentials, would not supply that estimate.
