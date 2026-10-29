# C88. A thirteen-letter mobile candidate

[Scope] The local dimension is d=13. The eleven material letters are C84's, with its boundary pin removed. Two vacancy colors implement the pinless grammar. C85's n=13 material certificate applies without changing a carrier vertex, matrix, coupling, or error budget. The conditional count is c=14>12. A ground-state statement for the full mobile ring still requires the inherited mobile comparison and infrared inputs, with pin defects explicitly included.

[Tags] Proved denotes an algebraic argument, an exact finite enumeration, or an integer/rational certificate. Numerical denotes floating-point calculations at the stated finite sizes. Assumed denotes an input not established here. Heuristic denotes the resulting infrared interpretation. Imported results are identified as inherited. Numerical tests below do not measure central charge.

## 1. Literal alphabet and bond

[Definition] The material letters, in matrix order, are

    L_0,L_1,R_0,R_1,D_A,D_B,X_A,X_B,E,M_A,M_B.

Append V_0,V_1. The markers are E,M_A,M_B. A cell has p data sites, n=p-1 active bits and b=p+1=n+2 material sites. The certified member has n=13, p=14, b=15. The code's last data bit is zero. There is no E_1 letter or conversion at the last data endpoint.

[Definition] The outgoing state o and incoming set I are

| Letter | o | I |
| :-- | :-- | :-- |
| L_0,L_1 | 0 | {0} |
| R_0 | 1 | {0,1} |
| R_1 | 1 | {1} |
| D_A,D_B,X_A,X_B | 1 | {0} |
| E | 0 | {0,1} |
| M_A,M_B | 0 | {0} |
| V_q | q | {q} |

A pair ab is legal when o(a) belongs to I(b). Let p_ext be the diagonal projector on the other pairs, and p_mat its material restriction. Unlike C84, L_1 followed by any marker is allowed. This is C83's automaton with E_1 deleted.

[Definition] Set z_A(a)=(-1)^a and z_B(a)=1. Each row denotes output ket times input bra with the displayed coefficient. Include the Hermitian conjugate. Sum a=0,1, s=A,B, and m over markers where indicated.

| Input | Output | Coefficient |
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

Write w for this matrix. The diagonal u assigns Delta to D_s,X_s,M_s, an additional D to X_s,M_s, mu_sel to every marker, and f a to L_a,R_a. Define

    b_11=w+(u tensor I+I tensor u)/2,
    h_11^0=K p_mat+b_11,
    h_13=JK p_ext+J embed(b_11)
         -t sum_(a material) sum_(q in I(a))
           (|a,V_o(a)><V_q,a|+h.c.)
         -(mu/2)(n_occ tensor I+I tensor n_occ).             (1)

The embedding is zero outside the material-material block. In particular the material half-site energies are contact terms, not independent physical onsite terms. The 169 by 169 matrix is implemented in bond13.py. Its parameter named n counts data sites p, a convention retained from the parent diagnostic scripts; the certificate drivers use n for active bits. Thus the certified literal bond uses params(n=14), not params(n=13).

[Proved] Equation (1) is a single Hermitian nearest-neighbor bond on C^13, repeated at every physical site. The charges are

    q(L_a)=q(R_a)=2a, q(E)=q(V_q)=0,
    q(D_A)=q(X_A)=q(M_A)=-1,
    q(D_B)=q(X_B)=q(M_B)=+1.

Conversion preserves -1+2=1+0, both port outputs have total charge zero, and all other moves transport the same charges. Hence Q, occupied number N, and marker number M commute with H on the full physical Hilbert space. Marker order and the material lengths between markers are conserved, up to cyclic rotation. Neither vacancy color number is conserved.

## 2. What selecting the pin sector means

[Definition] For a physical basis word w, remove vacancies and call the cyclic material word s(w). For N>0 put

    F(w)=sum_j 1_(s_j=L_1, s_(j+1) in {E,M_A,M_B}).        (2)

Set F=0 when N=0. This is a diagonal, translation-invariant operator on the original physical ring. It is generally nonlocal in physical coordinates, because the adjacent material letters can have an arbitrarily long vacancy interval between them. Let P_pin=1_(F=0). Restricting to it is selecting an invariant sector, not adding a local pin penalty or changing the onsite dimension.

[Proved] [H,F]=0 at every K,J,t,g,t_d,lambda, including on illegal words. Hopping preserves the entire squeezed word. A material rule acts on two consecutive material letters. An L_1 output occurs only on the left of a nonmarker partner, or is unchanged; every rule with L_1 input has a nonmarker right partner. A marker changes type only while remaining a marker. A rule whose right letter changes between marker types has left letter X_s,D_s or L_0, never L_1. No rule changes a right letter into L_1. These observations check the internal junction and both exterior junctions of a move. Thus no squeezed L_1-marker junction is created or removed. In the frozen chain each boundary pin projector is separately conserved.

[Proved] On full occupancy,

    H_C84=H_pinless+K F,
    H_C84|_(F=0)=H_pinless|_(F=0).                         (3)

This equality holds at finite K and nonzero port, not only after legal projection. All virtual paths starting in F=0 remain there. The physical mobile bond also preserves F=0 at finite K. However C85's three-color mobile bond is not identical to (1) at finite K, even on F=0. Its virtual wrong-color channels differ. Only the material certificate and the legally projected squeezed Hamiltonian are unchanged.

[Proved] At zero mismatch, starting after a marker gives an L prefix, followed either by one head T in {R_0,D_A,D_B,X_A,X_B} and an R tail ending in E, or by an all-L word ending in any marker. In F=0, the all-L word ends in L_0. These are C84's (5p+3)2^(p-1) cell states. In F=1 for a given cell, its legal states are

    L_(a_0)...L_(a_(p-2)) L_1 m, m in {E,M_A,M_B}.         (4)

There are 3*2^(p-1) of them for p>=1. Their zero-port intracell kinetic operator vanishes. For m=E their minimum cell energy is mu_sel+f, and the charged markers cost an additional Delta+D. Thus at f=0 a pinned cell loses the vacuum clock energy -2cos(pi/(b+1)). It is a physical defect, not an extra payload bit. At finite K the pinned sector can include illegal fluctuating words; the statement that every such word is frozen would be too strong.

[Proved, frozen ground-state selection] Section 6 proves that on compatible full-occupancy rings the ground energy is attained in F=0 and all cells have b=15, even at nonzero port. This uses an energy bound in addition to (2).

[Assumed, mobile ground-state selection, explicit hypothesis S_pin] For the full mobile chain, a claimed c=14 ground state must lie in the sector with F=0 and all material cells of length 15. We require the inherited common-bulk-subtracted mobile defect comparison to include pinned cells of length 15 as well as wrong-length cells and marker-free competitors. Concretely, at fixed compatible N=15M and positive vacancy density, let E_good(L,N) minimize energy in F=0 with all lengths 15, without imposing an artificial Q restriction. For every other sector sigma containing k defects, including pinned length-15 cells, require

    E_sigma(L,N)-E_good(L,N) >= J m_def k-C J,
    m_def>0, C independent of L,N,k,                      (5)

together with positivity of the finite-number defect energies that the boundary constant in (5) does not settle. A stronger bound without C on compatible rings is sufficient. Marker-free and all-vacancy competitors and the grand-canonical choice of N must also lie above the selected half-density ground state after tuning mu. All these comparisons are part of ground-state selection, not consequences of sector conservation. No theorem about the full ring is obtained by discarding other sectors by definition.

[Proved reference cost; assumed transfer] In the half-density projected charge reference, alpha_contact=1/2-2/pi^2>0 multiplies the material Hamiltonian. At lambda=f=0, replacing one b=15 clock cell by a pinned one costs

    J alpha_contact * 2cos(pi/16)>0.                       (6)

This is an order-J defect cost with no small port factor. Small frozen port and grammar corrections fit well below it. The transfer of this positive reference cost to a mobile common-bulk-subtracted comparison is precisely the pin part of S_pin. A gapless vacancy fluid precludes inferring (5) just from a finite-volume perturbation expansion. The candidate is conditional on S_pin, along with the other infrared inputs. We do not claim an unconditional mobile sector-selection proof.

## 3. Unique coloring, projection, and suppression

[Proved] Between consecutive material letters a,b, a nonempty vacancy interval must start with color o(a). Every following vacancy has that same color. The final bond is legal if and only if o(a) belongs to I(b). Thus at every fixed nonempty occupancy pattern, every legal squeezed word has exactly one legal coloring. The result includes pin defects and one-material-letter rings. The all-vacancy sector has two constant legal colorings and is excluded from this uniqueness statement.

Each illegal squeezed junction forces a physical mismatch in its own intervening interval. These intervals partition the physical bonds, so on N>0

    N_matbad <= N_bad=sum_i p_ext,i.                       (7)

Pin count F is not counted by N_matbad. A pin is allowed in the pinless grammar, and its exclusion is the separate energy issue in Section 2.

[Proved, restricted minimum] Two colors are needed even within the pin-free selected-cell sector, if every legal word and vacancy placement must be retained using a nonnegative diagonal pair penalty. That sector contains L_0 L_0 and R_0 R_0 junctions. A single vacancy would have to accept L_0 V, V L_0, R_0 V, and V R_0. Then it admits R_0 V L_0 although the squeezed junction is illegal. Two colors attain the bound. This is not a bound for arbitrary redesigned mechanisms.

[Proved] A projected forward swap V_q a -> a V_o(a) has amplitude -t. In the reverse direction the predecessor fixes exactly one q. All material moves preserve the outgoing state of their rightmost letter. The left exterior bond tests any narrowed incoming set of their first letter. Therefore the unique-color isometry U satisfies

    U^dagger P H P U
      =H_hop+J sum_j 1_(g_j=0) P_mat b_11,j P_mat-mu N,    (8)

where P=1_(N_bad=0), g_j is the vacancy count after material letter j, and the physical root and cyclic rotation are retained. The isometry commutes with rotation and with F. Colors add no independent legal degree of freedom or extra boson.

[Proved, fixed finite volume] In an invariant subspace containing a legal trial state, write H=kappa N_bad+B, kappa=JK, and Q_bad=1-P. A ground energy obeys E<=||B||. For kappa>2||B||,

    Q_bad(H-E)Q_bad >= (kappa-2||B||)Q_bad,
    Q_bad psi=-[Q_bad(H-E)Q_bad]^-1 Q_bad BP psi,
    ||Q_bad psi||^2 <= ||B||^2/(kappa-2||B||)^2.           (9)

By (7), illegal material weight is no larger than this. At fixed positive J this is K^-2 suppression. It applies to nonzero bits and port, with any finite number of vacancies. It is not volume-uniform total leakage. Section 7 gives the dressed local statement.

## 4. The complete leading grammar correction

[Definition] Write B=t T+J W-mu N, with lambda included in W. On P define the physical cyclic counters

    C=sum_i 1_(a_i in {R_0,E}, a_(i+1) vacancy),
    A=sum_i 1_(o(a_(i-1))=1, a_i=R_0,
                a_(i+1) in {R_0,R_1,E}),
    Z=sum_i 1_(o(a_(i-1))=1, a_i=E, a_(i+1)=R_0).

[Proved] Every one-step illegal intermediate has one mismatch. The full order-1/K correction is

    D_2|P=-P B Q_bad N_bad^-1 Q_bad BP/(JK)
          =-P(t^2 C+J^2 A+2J^2 lambda^2 Z)P/(JK).          (10)

For a reverse hop a V_o(a) -> V_q a, exactly |I(a)|-1 colors are wrong. This equals one for R_0,E and zero otherwise. The sole mismatch is immediately before the new vacancy. Only the inverse transition returns to P in one step. This gives the color-ambiguity contact self-energy -t^2 C/(JK).

The vacuum moves R_0 R_a -> L_a R_0 and R_0 E -> L_0 E narrow the incoming set from {0,1} to {0}. They fail only when the predecessor has state 1. Each has unit squared amplitude, giving the vacuum-tail self-energy -J A/K. The two port outputs similarly narrow E to M_s and each has squared amplitude lambda^2, giving the wrong-end port self-energy -2J lambda^2 Z/K. Their illegal outputs are distinct and have only their inverse one-step legal return.

Rail moves, switches, conversion, and their inverses either keep or enlarge incoming permissions and preserve the right outgoing state. They create no other channel. A hopping-created mismatch ends on a vacancy; a material-created mismatch ends on a material letter. These intermediate sets are disjoint, so the mixed tJ operator vanishes. There are no g^2/K, t_d^2/K, or mixed vacuum-port terms in this compression. Diagonal terms have no P-to-Q_bad element.

[Proved, bit dependence] In a zero-excitation clock word,

    C=1_(vacancy after head R_0, if present)
      +sum_r (1-a_r)1_(bit r in R tail, vacancy after it)
      +1_(vacancy after E).                               (11)

There is no L-prefix contact term and no D/X contact term. The terminal word has no head. Averaging a bit-independent clock/charge reference gives a scalar plus the positive field

    +(t^2/(JK)) sum_r w_r a_r,
    w_r=<1_(bit r in R tail, vacancy after it)>_ref.         (12)

Before averaging, the centered contact indicator gives a payload-density coupling to clock and vacancy fluctuations at the same scale. This term has no lambda^2 factor and is not a common scalar. Relative to C85 its tail-bit coefficient is halved and its L/D/X coefficients disappear.

[Proved, inherited vacuum projection] At full occupancy the A term reduces to

    -(J/K) sum_(r=0)^(n-1) (1-a_r) sum_(j=0)^r |v_j|^2,
    v_j=sqrt(2/(p+2))sin(pi(j+1)/(p+2)).                  (13)

The tail position r+1 is occupied by R_(a_r) whenever j<=r. Only a_r=0 allows the offending vacuum hop, yielding (13). With vacancies A retains its physical contact factors. Z retains its full clock/contact dependence. Auxiliary dressing of these operators need not be scalar.

[Proved, higher-order definition] Every further finite-volume correction is retained by

    H_eff(E)=PHP-PBQ_bad
      [JK N_bad|Q_bad+Q_bad B Q_bad-E]^-1 Q_bad BP.         (14)

The next-order strengths include t^3/(JK)^2, t^2 J B_mat/(JK)^2, t J^2 B_mat^2/(JK)^2, and J^3 B_mat^3/(JK)^2. These include corrected hopping, fields, contact-dependent denominators, and multibit interactions. B_mat contains Delta,D,mu_sel,f, all rail terms, g and lambda. The chemical potential cancels within a conserved-N denominator. No Gaussianity or scalar form is asserted for this remainder. Section 7 provides an energy-independent local version.

## 5. The material certificate is unchanged

[Proved by (3); inherited C84/C85 reduction] On F=0 the finite-K material operator, its code, and every virtual port path equal C85's material operator. Consequently its n=13 certificate applies without modification. The folded scalar carrier is

    M,X_13,...,X_0,D_0,...,D_13.

Its diagonal is 1004 on M,X and 4 on D, path entries are -1, and the closing entry H_(M,D_13)=-i t_d. Only D_0,...,D_12 convert. There are 29 vertices and an 841-dimensional product resolvent. D_13 and its X partner cannot be dropped merely because their conversion is absent.

Let R=(H_c tensor I+I tensor H_c)^-1, x=R|M,D_0>, and U_i,V_i retain respectively row D_i and column D_i of x. For i,j=0,...,12 set

    C_ij=U_i^dagger R V_j,
    A_ij=U_i^dagger R U_j+V_i^dagger R V_j,
    T_quad=-2 alpha^4 g^2 Re C,
    mu_quad=-2i alpha^4 g^2 Im A.                          (15)

C84's two-species CAR calculation makes the g^2 port coefficient scalar plus this number-conserving quadratic operator. Higher even orders need not be quadratic.

[Definition] The same material schedule as C85 is

    n=13, p=14, b=15, D=1000, Delta=4-2cos(pi/16),
    phi_A=pi/28, phi_B=0, psi_A=pi, psi_B=pi/2,
    d_0=1, d_1=1004, d_(j+1)=1004d_j-d_(j-1), beta=1/d_15,
    alpha=sqrt(2/16)sin(pi/16), E_unit=alpha^4 beta^2 g^2,
    t_d=beta 2^-66, g=2^-272, lambda=2^-608, K>=2^2560,
    f=lambda^2 E_unit 2^-120.                             (16)

The selector is defined in Section 6. These are fixed parameters for every cell length in the same Hamiltonian. Extra mobile conditions on K are in Section 7. No finite numerical infrared window is claimed.

[Proved, regenerated integer/rational certificate] The local unchanged C85 pipeline uses 256-bit tree and derivative solves and a 512-bit finite-loop solve. Integer component radii, dyadic congruence, rational endpoint inertia tests and rational phase-distance bounds give

    26 simple Fermi points,
    351 distinct unordered pair sums, including repeats,
    delta_2>0.0156742281,
    delta_4>3.2513114e-6,
    additional normalized Bloch allowance=2^-35.           (17)

Endpoint inertia increases 0 to 13 and decreases back to 0. Multiplying the determinant of the nearest-neighbor n=13 pencil by exp(13ik) gives a polynomial of degree at most 26. The 26 disjoint sign-changing brackets exhaust that degree, proving simplicity and completeness for this quadratic reference. Pair sums use exact unit phases at rational tan(k/2) bracket midpoints, with chord errors subtracted. The statement is low-order nonresonance, not all-order rational independence. The certified Hermitian hopping margin exceeds 0.5560750 and the normalized entry error is below 6.339e-92.

[Inherited proved finite material bounds, independently regenerated arithmetic] C85 bounds the finite-g coefficient tail by alpha^4 g^4/4, the full physical finite-K correction including zero-port splitting by 2^84/K, and the isolated two-cell port remainder by 2^52 lambda^4. Equality (3) transfers all three. With tau_W=tr(W W^dagger), the normalized bounds are

    err_g=4n beta^-2 g^2 tau_W <1.727e-37,
    err_K=16n tau_W (2^84/K)/(lambda^2 E_unit) <3.024e-83,
    err_lambda=16n tau_W 2^52 lambda^2/E_unit <2.396e-54.   (18)

The exact Fraction test gives sum<(2^-35)/8, using alpha^4>=64/16^6=2^-18. The denominator lambda^2 E_unit in err_K is necessary because the zero-port splitting has no lambda^2 factor. The inherited root routine's threshold g<=2^-228 is an intermediate scalar threshold; the final schedule and (18) use g=2^-272. These are finite material bounds. They do not certify the full interacting mobile chain or turn its longer-range effective operator into a degree-26 pencil.

## 6. F3, including frozen pin defects

[Inherited proved all-string bound; proved application] The pinless material bond is the compression of C80 obtained by deleting E_1. The positive pin is unnecessary for this compression bound. For every complete cell of material length l, including illegal words and all auxiliary excitations,

    E_l^infinity-(l-2)_+/(K-4) <= E_l(K) <= E_l^infinity,
    E_l^infinity=mu_sel-2cos(pi/(l+1)).                    (19)

The lower bound follows by compressing C80's all-string form inequality. The upper trial state is its vacuum with last bit zero. Delta>2+g, |t_d|<=1 and f>=0 hold at (16). The field is nonnegative and vanishes on the empty-bit trial state. The marker-free density is at least -1/(K-4). This bound is uniform in l, not an extension of a small-cell matrix norm estimate.

[Definition] Put a_l=2[1-cos(pi/(l+1))], delta_l=(l+1)a_l-l a_(l+1), and

    mu_sel=2-(delta_15+delta_14)/2,
    e_l=E_l^infinity/l,
    gamma=min_(l!=15)(e_l-e_15)-1/(K-4).                  (20)

[Proved] The continuous extension has

    a''(x)=4pi sin(pi/(x+1))/(x+1)^3
           +2pi^2 cos(pi/(x+1))/(x+1)^4>0, x>=1.

Hence delta_l strictly decreases and

    e_(l+1)-e_l=[2-mu_sel-delta_l]/[l(l+1)].

The midpoint selects l=15; the nearest competing densities are at 14 or 16. Rational alternating-series enclosures regenerated by selector_certificate.py give

    mu_sel=1.8891880765550189...,
    min_(l!=15)(e_l-e_15)>2.8116274e-5,
    e_15<-0.0048254989,
    vacuum gap>0.11381149,
    gamma>2^-17.                                         (21)

[Proved, pinned length-15 cells at finite K] In the F=1 cell block write H_cell=K N_bad+B_cell. The zero-mismatch subspace has B_cell>=mu_sel by (4). The conservative bound ||B_cell||<=B_*=2^16 from C85 applies unchanged. For any normalized vector with legal and illegal norms x,y,

    <H_cell> >= mu_sel x^2-2 B_* xy+(K-B_*)y^2.

Subtract mu_sel, use 0<mu_sel<B_* and x<=1, and complete the square in y. For K>2 B_* this gives

    E_pin,15 >= mu_sel-epsilon_pin,
    epsilon_pin=B_*^2/(K-2 B_*).                          (22)

The pin-free cell has trial energy E_15^infinity=mu_sel-2cos(pi/16). Thus the pin defect excess over that trial is at least

    r_pin=2cos(pi/16)-epsilon_pin>1                       (23)

at (16). The field only improves (22). This avoids assuming that the entire finite-K pin block is kinetically inert.

[Proved, nonzero-port F3] All markers have outgoing state 0. Assign the left-marker-to-first-data penalty to the cell on its right; the switch is the identity on that left marker. All zero-port terms then belong to individual complete cells. The port alone couples cells, with norm B_port=sqrt(2)|lambda|.

Define A_m as the exact ground energy of m open good cells, each length 15 with pin count zero, and u=lim_m A_m/m. Cutting one port gives |A_(m+r)-A_m-A_r|<=B_port. Applying this to q copies and dividing by q proves |A_m-m*u|<=B_port. Consequently,

    A_m>=m*u-B_port,
    S_m<=m*u+2B_port,
    u<=E_15^infinity.                                    (24)

S_m is the periodic good-sector energy. For the last inequality use a product of empty-bit zero-excitation clock states, whose port expectation vanishes. No degeneracy of the finite-K code is needed.

For each wrong-length cell, (19) gives E_l>=l*u/15+l*gamma. For each pinned cell of length 15, (22) gives E_pin,15>=u+r_pin. On a compatible frozen ring L=15m with k defective cells, cut every port touching a defect, at most 2k. There are at most k good open runs. Applying (24) and the cell bounds yields

    E_defect-S_m >= gamma sum_(wrong lengths) l
                   +r_pin k_(pinned length 15)
                   -3k B_port-2B_port
                >= k[min(gamma,r_pin)-5B_port]>0.         (25)

The all-defect case obeys the same weaker bound without good runs. Marker-free rings are excluded by (19) and e_15<0. Thus on compatible full-occupancy rings F3 selects length 15 and F=0. This explicitly proves frozen ground-state pin selection, not just wrong-length exclusion. Mobile contact factors invalidate the cell-factorization step, so (25) is not a proof of (5).

## 7. F4 and the mobile hierarchy

[Proved, inherited local construction applied to the new bond] Write H=kappa N_bad+B with kappa=JK and local strength nu=sup_i sum_(X contains i)||B_X||. The p_ext,i commute and have integer eigenvalues. For a local operator O_X define

    (O_X)_m=(1/(2pi)) integral_0^(2pi)
       exp(i theta N_bad) O_X exp(-i theta N_bad)exp(-im theta)dtheta.

Only penalties meeting X enter. Support grows by one site at each end, the norm does not increase, [N_bad,(O_X)_m]=m(O_X)_m, and |m|<=|X|+1. Set

    S_1=sum_(m!=0) B_m/(kappa m),
    F_2=[S_1,B_0]+[S_1,B_off]/2,
    S_2=sum_(m!=0)(F_2)_m/(kappa m).

The generators are anti-Hermitian and [S_1,kappa N_bad]=-B_off. With U=exp(S_1+S_2), the BCH expansion gives

    U H U^dagger=kappa N_bad+B_0+(F_2)_0+R_2,
    [N_bad,B_0+(F_2)_0]=0,
    ||R_2||_local,exponential<=C_3 nu^3/kappa^2.           (26)

Since B_mP=0 for m<0,

    P(B_0+(F_2)_0)P=PBP-sum_(m>0)PB_(-m)B_mP/(kappa m), (27)

which gives the negative sign and denominator JK in (10). Local Fourier components preserve F because F commutes with N_bad and with each literal term. Thus this layer does not mix pin sectors.

[Proved, volume-independent fixed-order scope] Nonzero nested commutators have connected overlapping interval supports. Summing generator norms meeting an interval I costs at most |I| times the local strength. Each commutator contributes a factor two and support grows linearly with the order. The resulting k! C^k bound is canceled by BCH's 1/k!, leaving a convergent geometric tail for the fixed finite generator at sufficiently small nu/kappa in a weaker exponential norm. Recursively dividing nonzero integer transfers by kappa m gives, for each fixed r,

    U_r H U_r^dagger=kappa N_bad+D_r+R_r,
    [D_r,N_bad]=0,
    ||R_r||_local,exponential<=C_r nu(nu/kappa)^r.          (28)

This is not convergence as r tends to infinity or a global many-body band gap. For legal phi, P[S,p_i]P=0, so a dressed local mismatch is O((nu/kappa)^2). Identification with actual low-energy states needs stability under the remainder. Summed outside-P weight can be extensive.

[Proved sufficient grammar budget] The local interaction norm of the explicit diagonal representative (10) is at most

    A2/(JK), A2=2t^2+3J^2+6J^2 lambda^2.                 (29)

C has range two and coefficient at most one; A,Z have range three. Fix t>0, J>0 and a material tolerance tau>0. A safe unpenalized bound is nu<=2^13 J+4t+2|mu| at (16). Choose finite K with

    K>=2^2560,
    nu/(JK) in the fixed-order conjugation domain,
    K>2 A2/(J^2 tau),
    K>sqrt(2 C_3 nu^3/(J^3 tau)).                         (30)

Then both (29) and the remainder (26) are below J tau/2. For the material bands alone use

    tau<alpha_contact lambda^2 E_unit 2^-35/(128 n tau_W).

Selector, pin, auxiliary, and infrared tolerances can require smaller tau. Increasing K preserves (18). The existence of a sufficiently large finite K for any fixed positive tolerance follows from (30), but C_3 and the full infrared tolerance have not been numerically specified. Thus the mobile member is an implicitly bounded family, not a numerically certified finite-J infrared point. JK much greater than t alone is insufficient because of the t^2/(JK) contact field.

## 8. Finite checks and their limits

[Numerical] check_lift.py compares the pinless material matrix with the C84 matrix minus its pin at two generic parameter points, with error below 1e-13. Hermiticity and N,M,Q commutators vanish entrywise. Exhaustion of all 13^4-2^4=28545 nonempty four-site words verifies (7), and every legal uncolored pattern has one coloring. Legal colored and squeezed uncolored matrices agree entrywise in all four invariant-closure tests.

[Proved, exact enumeration] check_sector.py checks all 28561 four-site physical words and 50024 directed bond transitions, including illegal colors, for conservation of the integer F. Its legal two-cell graph used for letter deletion has 3528 states and 19376 directed transitions. It supplies explicit input/output witnesses for all thirteen letters, as explained in Section 9.

[Numerical] check_virtual.py builds the entire legal-to-illegal map on every nonempty legal ring of lengths 3,4,5. The legal dimensions are 376,2685,19267 and the intermediate counts 96,840,7370. Every intermediate has one mismatch, the mixed operator vanishes entrywise, and the largest error in (10) is below 4.5e-16. These calculations include nonzero g,t_d,lambda, field and selector, and include legal pinned words.

[Numerical] The following ED uses full invariant closures, without a mismatch-number cutoff. J=t=1, D=4, g=.07, t_d=.003, Delta=4-2cos(pi/4), phi_A=pi/4, mu=mu_sel=f=0. The first seed has two p=3 cells with bits (1,0), lambda=0, and one vacancy. The other seeds have one p=2 cell, Q=2, lambda=.1 and one or two vacancies. Parameters are diagnostic, not (16).

| Closure | Dimension | Illegal weight K=100 | K=1000 | K=10000 |
| :-- | --: | --: | --: | --: |
| Two cells, two colors | 2484 | 8.31442e-5 | 8.39207e-7 | 8.39984e-9 |
| Nonzero port, one vacancy | 268 | 2.31273e-7 | 2.67675e-9 | 2.71795e-11 |
| Nonzero port, two vacancies | 1045 | 1.33164e-7 | 1.52294e-9 | 1.54444e-11 |

The uncolored controls have dimensions 1296,180,450 and approach respectively .0229959,2.51578e-5,2.64428e-5 at K=10000. The largest residual is below 1.5e-11. These are ground states of their specified invariant closures. The one-cell port tests retain auxiliary excursions through illegal states; they do not constitute a multi-cell interacting infrared test. The legal counter test separately checks allowed intercell port channels.

[Numerical] A fourth check has two material cells of lengths 3 and 1, one active bit equal to one, two vacancies, and lambda=.1. The length-1 cell has no data site, but its marker supports a legal intercell port. The full colored closure has dimension 3078, legal dimension 1080, and illegal weights 5.66133e-14,6.28110e-16,6.34945e-18 at K=100,1000,10000. Its physical mismatch expectations are 1.27153e-4,1.31873e-6,1.32358e-8. The uncolored closure has dimension 1200 and illegal weight approaching 6.24315e-13. These small material weights are an additional check, not the main resolved leakage evidence in the table. The largest residual across all four closures is below 5.1e-11. Neither the unequal diagnostic cell lengths nor this test replaces selection at b=15.

[Numerical] check_energy.py resolves all four p=3 two-bit strings, lambda=0, one vacancy. Their constrained ground energy is -3.29813082968. The coefficients of the negative 1/K shift are

| Bits | Color coefficient | Vacuum coefficient | Full dimension |
| :-- | --: | --: | --: |
| 00 | -.6691032812 | -.4480522957 | 320 |
| 01 | -.5146993463 | -.1024562307 | 230 |
| 10 | -.6396299070 | -.3455960650 | 230 |
| 11 | -.4852259720 | 0 | 145 |

All legal dimensions are 80. At K=10000 the rescaled shift agrees with their sum within .0002. This checks the non-scalar color field in (12). At J=.01 and kappa=100,1000,10000, a Q=2 nonzero-port closure of dimension 268 gives illegal weights 3.99344e-11,3.40562e-13,3.33422e-15. Its mismatch expectations are 5.90198e-5,5.85049e-7,5.83795e-9. Color mismatch and illegal squeezed material are distinct observables.

[Numerical] Frozen F3 uses complete 576-dimensional marker sectors with lengths (2,2) and (1,3), at the b=2 midpoint selector and K=100. At lambda=0,.001,.003 the energy differences are .5833105881,.5833106366,.5833110243, above the respective bounds .1744288125,.1602866769,.1320024057. Separate pin blocks in the (2,2) sector have dimensions 441,126,9 for F=0,1,2. At mu_sel=1.1 and lambda=.003 their energies are .1999991604,1.1999991768,2.2. Pin-free energy is lowest. Residuals in these frozen checks are below 2.4e-13.

[Numerical] F4 uses the full 60-dimensional closure of L_0 E V_0 at nonzero port. At K=40,80,160,320, K^2||R_2|| is 7.95144,7.93017,7.91948,7.91413. The second-order compression agrees with (27) below 6.2e-21 and unitarity error is below 4.5e-16.

## 9. Can another letter be removed the same way?

[Proved, restricted answer] No remaining single letter can be deleted by discarding only components disjoint from the retained clock/port sector. check_sector.py starts at the legal pin-free word

    L_1 R_0 E L_0 R_0 E V_0,

with two p=2 cells and Q=2. Its legal connected component contains every material letter and both vacancy colors. More strongly, for each of the thirteen letters a it exhibits an allowed transition from a word with no a to a word containing a. Thus the subspace avoiding a is not invariant even within the legal pin-free selected-length component. Exact witness words and nonzero amplitudes are in sector_checks.json. These paths use only the literal local table, so their existence persists at any nonzero g,t_d,lambda in (16).

For the material letters the mechanism is visible directly. Vacuum motion exchanges L_a storage with R_a storage; R_0 also supplies the moving head. The port creates M_A,D_B and M_B,D_A from E,R_0. Switches exchange D_s and X_s. The terminal moves join X_s,E and L_0,M_s. All these states occur within the retained dynamics, unlike the isolated legal terminal-L_1 configurations. E_1 was different because its deletion also changed the conversion and endpoint rules; it was not removed as an invariant unused letter of the original parent.

Vacancy colors can change when moving across R_0 or E, and the two-color obstruction in Section 3 survives in pin-free cells. Thus neither V_0 nor V_1 is a superselection label. This excludes a further saving by this specific deletion argument. It is not a universal lower bound d>=13, nor does it exclude a redesigned ten-letter material mechanism or a different vacancy encoding.

## 10. Conditional infrared conclusion

[Assumed, inherited C81/C83/C85 dictionary] The projected squeezed reference has one compact charge boson at half occupancy and contact alpha_contact. Positive compressibility permits the density tuning at finite J. Physical cyclic gluing and the root seam are retained. The rotation-equivariant coloring introduces no new branch or twist.

[Proved conditional on that dictionary] For b=15, a payload vertex of cell momentum q and reference dimension Delta_pay has branches

    Q_h=q+2pi ell,
    k_phys=Q_h/30 mod 2pi,
    Delta_ref=Delta_pay+[Q_h/(30pi)]^2.

A uniform term needs Q_h=60pi h. For h!=0 the charge contribution is at least four. For h=0 the certified root and pair-sum margins exclude nonforward neutral quadratic and quartic resonances. Higher neutral fermion vertices have reference dimension at least three. Exact Q and N exclude anomalous pairing and number-changing charge terms. Below dimension three the inherited mismatch bound is

    delta_phys>=min(delta_internal/30,pi(2-sqrt(3)))
              >1.0837704e-7.

The first marker harmonic has wavevector pi/15 and dimension 1/225. It is an oscillating observable, not a uniform relevant perturbation. The bit/contact operator (12) remains part of the matching and stability problem.

[Assumed] The microscopic ground-state inference requires massive-sector matching in local coefficient norms, low-frequency vacancy matching at fixed positive J, S_pin and the full mobile defect comparison, multichannel stability of the interacting retarded action, and current/stress/thermal control sufficient to identify conformal modes. exp(i pi Q/2) has auxiliary eigenvalues +/-i, so an undressed microscopic Z2 fermion grading cannot be imported unchanged.

[Heuristic, conditional on those inputs] Thirteen payload modes and one vacancy mode give

    d=13, c_candidate=13+1=14>12=d-1.                     (31)

This improves the conditional dimension from C85's 14 to 13. It is not an established counterexample or a measured central charge. Ground-state pin selection is proved for the frozen ring and explicitly assumed as an energetic comparison for the full mobile ring.

## 11. Reproduction and provenance

[Definition] Run python3 -B run_all.py from work_C88 with NumPy and SciPy. It rebuilds every certificate and diagnostic locally. Source copies of C85's scalar certificate pipeline retain their C84/C80/C75/C33 attribution. They are unchanged. bond13.py and the diagnostic scripts are local adaptations of C85; check_sector.py is new. No parent output is relabeled as a new certificate. The unchanged numerical material margins are regenerated.

[Proved, exact comparison] certificate_transfer.json records byte equality of the ten certificate source files with C85 and structural equality of the regenerated exact root, finite-budget and selector JSON certificates with C85. Timing fields in diagnostic summaries are not part of these comparisons.

[Scope] All Python and writes stay in work_C88. No network commands, subagents or cluster jobs are used. Project-level STATE.md is outside the allowed write scope; local STATE.md supplies the handoff. The unresolved items are those explicitly tagged assumed above.
