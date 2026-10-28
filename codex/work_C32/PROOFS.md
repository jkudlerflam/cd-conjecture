# C32 finite constructions and bounds

[Proved, scope] The construction below gives a mobile positive folded-bus
core, an arbitrarily accurate second-order equalization theorem, and a
finite-parameter all-n crossing theorem. It does not prove an all-n
distinct-center or nonresonance theorem. The latter limitation is separate
from all the bounds proved here.

## 1. A positive initialized core

[Definition] Write A for the new left marker and B for the old right marker.
Use the following 28 letters.

| Letters | Number | Onsite charge |
| :-- | --: | :-- |
| L_a,R_a, a=0,1 | 4 | a |
| d_(s,a), s=v,0,+,- | 8 | a+q_s |
| x_(s,a), s=0,+,- | 6 | a+q_s |
| A_empty,B_empty | 2 | 0 |
| A_s,B_s, s=v,0,+,- | 8 | q_s |

[Definition] Here q_v=q_0=0, q_+=1, q_-=-1. The initialized occupied word
has k=n+2 sites, one A, n data sites and one B. Its legal configurations
contain one token, on A, B, a data rail d, or an express rail x. A vacuum
token v cannot occupy x. Before the token the data letters are L; after
it they are R. The exact legal-sector dimension is (7n+8)2^n.

[Definition] All rules are two-site rules. Data-rail motion, express
motion, the express-to-B link, the weak data-to-B link, and conversion
are the corresponding rules of routeA_spatial, with its M replaced by B.
The junction switch is enabled by A immediately to the left of d_0.
There is an additional left link

    A_s R_a <-> A_empty d_(s,a),       s=v,0,+,-.

[Definition] Its amplitude is -1. The vacuum right link is
d_(v,a) B_empty <-> L_a B_v, also amplitude -1. All rail amplitudes are
one. A rightward neutral rail hop carries (-1)^a on the departure data.
Each rightward neutral express hop also carries exp(i phi). There are n
such phase-carrying hops between d_0 and B. Charged rails are real. The
switch has amplitude -1. The optional loop has amplitude
-t_d exp(i psi_s)(-1)^(a [s=0]), with psi_0=0.9, psi_+=psi_-=0 in the checks.

[Definition] Replace every kinetic matrix element c |v><u|+h.c. by

    E(c;u,v)=|c|(|u><u|+|v><v|)+c|v><u|+conj(c)|u><v|.       (1)

[Proved] This matrix has eigenvalues 0 and 2|c|. Thus every kinetic term
is positive. A vacuum edge has c=-1 and annihilates the equal-amplitude
sum of its endpoints. This replacement is a genuine change of the
23-letter Hamiltonian, not a claim that its original sine vacuum was
frustration free.

[Definition] On each letter add m P_exc, with an additional D P_exc on
A, B and x. On a data token add

    g (|d_(+,0)><d_(0,1)| + |d_(-,1)><d_(0,0)| + h.c.).   (2)

[Proved] The onsite blocks are positive for m>|g|. Each nontrivial
conversion block has eigenvalues m+g,m-g. All other excited data entries
are m, and high-rail or marker entries are m+D. Vacuum entries vanish.
Split onsite terms equally between adjacent bonds and put the missing
endpoint halves in the positive boundary operators a_pay,b_pay.

[Definition] Positive diagonal penalties impose the legal arrows. Bulk
L,R,d/x have arrows L->L, R->R, L->R. A_empty has outgoing L, A_s
outgoing R; B_empty has incoming R, B_s incoming L. No A may be the
right member of an internal bond and no B its left member. Boundary
penalties require A first and B last. The remaining pairs must match
arrows. All illegal-pair and boundary penalties have coefficient K>0.

[Proved] Every kinetic move preserves the outer arrows and marker
positions. Every violation projector therefore commutes with the
Hamiltonian. This includes the extra A link: its internal matching
changes on both sites together, and its outer arrows do not change.
Conversion preserves each arrow. Illegal sectors have energy at least K.

[Proved] In a legal zero vector, (2) and m>|g| first exclude every excited
token. The remaining graph, for each fixed n-bit data string, is the path
A,d_0,...,d_(n-1),B. Equation (1) forces a constant vector on this connected
path. Thus the complete core kernel is

    G|a_0...a_(n-1)> = (n+2)^(-1/2)
                      sum_(j=A,d_0,...,d_(n-1),B) |v at j; a>.       (3)

[Proved] It has dimension 2^n. The legal vacuum excitation gap is
2[1-cos(pi/(n+2))], the legal excited sector costs at least m-|g|,
and illegal sectors cost at least K. Consequently

    gamma_core >= min(K,m-|g|,2[1-cos(pi/(n+2))]) > 0.               (4)

[Proved] The source operators F_0,F_+,F_- on A or B are the maps in
mobile.FMAP. They only change v to an excited species or the reverse,
so G^dag F G=0. The port bond is

    V=F_0^B F_0^A+i F_+^B F_-^A-i F_-^B F_+^A.                    (5)

[Proved] Each term preserves total Q=sum_site(a+q_s). This is an onsite
integer U(1), with parity exp(i pi Q). On (3) it is the logical number
and parity. All claims of exact degeneracy refer to the unperturbed
core; a later small uniform field splits that code as intended.

## 2. The fixed-payload mobile chain

[Assumed input, explicit attribution] Use the proved local wrapper and
selector of C8's length_filter.md, arrow_wrapper.md, selector_wrapper.md
and helix_equalizer.md with the positive initialized core above. Their
global exclusion arguments are inherited rather than independently
enumerated here. The local positivity and initialization hypotheses
needed by those arguments were established in Section 1.

[Definition] The complete onsite space is

    C^4_selector tensor
      [ (C^3_tag tensor C^2_filter tensor C^28_core) direct-sum C^3_hole ]
      tensor C^2_helixA tensor C^2_helixB.

[Proved count] This gives

    d0=4(3*2*28+3)*4=2736.                                       (6)

[Definition] The occupied core length is k=n+2, and the two physical
cell lengths are b=k+1=n+3 and b+1=n+4. Both contain exactly n logical
fermions. All added tracks are neutral under Q. Vacancy swaps move the
entire charged core label, so they conserve Q exactly. The helix,
occupation filter, selector and marker terms also conserve Q.

[Proved using the inherited wrapper] At p=1,2 holes the exact code is

    W_p=chi_p tensor G tensor g_filter tensor G_helix,k+p.          (7)

[Definition] In ordered occupied coordinates, let o_1<...<o_k be the
particle positions, D_j=1_(o_(j+1)=o_j+1), and O_L,O_R the physical
endpoint occupation projectors. The charge Hamiltonian is

    H_ch=U*(number of occupied-hole interfaces)-w*(adjacent swaps).

[Proved] Its graph is connected when w>0. Shifting a multiple of the
identity and applying the elementary Perron eigenvector theorem gives
a unique strictly positive ground vector chi_p. Reflection gives equal
endpoint weights r_p=<O_L>=<O_R>>0. Subtract its ground energy. The full
charge-core operator in this sector is

    H_p=(H_ch-e_p) tensor I
        +s[sum_j D_j tensor B_j+I tensor(a_pay+b_pay)].             (8)

[Proved] Here the B_j and a_pay,b_pay are positive. Equation (7) is
annihilated term by term after the charge shift. Conversely a zero
vector first lies in chi_p tensor core; strict positivity of every
<D_j> then forces every core term to annihilate it. This proves the
complete code in the ordered sector. C8 supplies the occupation,
arrow, marker-free and wrong-length exclusions on the unrestricted ring.

[Proved] The endpoint helix compensator with expectations m_p=1/r_p
equalizes the source projected onto the fast wrapper ground space.
It does not equalize the complete resolvent. In particular, after that
projection the core Hamiltonian is

    B_p=sum_j d_(p,j) B_j+a_pay+b_pay,
    d_(p,j)=<chi_p|D_j|chi_p>.                                  (9)

[Proved] These adjacency factors generally depend on p. This is the
same obstruction identified by C25. Source equalization alone cannot
remove it, even though every first-order port compression is zero.

[Proved] The local converted marker move transfers a hole and leaves
the ordered occupied core states unchanged. Its code matrix element is

    -u beta_n I_payload,
    beta_n=A_b^2 beta_v cos(phi_h/2)^(2b)>0,
    beta_v=[sum_(j=0)^k chi_1({j})chi_2({j,k+1})]^2.               (10)

[Definition] A_b is the sine-head overlap in C8; take phi_h=1/b.
The required hole conversion is H_0->H_2 and the transferred-site helix
rotation is the one specified in C8. There is no transport of an
occupied data site through the cut.

[Proved] At fixed cell number M, the leading length operator is the
nearest-neighbor hard-core hopping chain on the cell index,

    -J_l sum_j (s_j^+ s_(j+1)^-+h.c.),  J_l=u beta_n.               (11)

[Proved] Jordan-Wigner turns (11) into free fermions. At density nu in
(0,1) its cell-coordinate Fermi velocity is 2J_l sin(pi nu). On a
physical ring the cyclic cell-label constraint still glues the length
and payload momenta, as in C8 and C13. Equation (11) is a statement
about the leading local operator, not an unrestricted tensor product
of finite-ring spectra.

## 3. A noncommuting second-order equalization bound

[Definition] In a two-cell port source sector write the complete shifted
operator as A+sW. A is the fast charge and helix Hamiltonian, P its
ground projector, Q=I-P. Let

    QAQ>=delta Q,  B=PWP>=a P,  ||W||<=L,
    C=QWP,  D=Q(A+sW)Q,  t=s/delta.

[Proved] One may take L=2||H_core|| and
a=2 d_min(m-|g|), where d_min is the least adjacency weight over the
two cells. W is positive even though its different core terms need
not commute. The charge-core excitations, including the neutral ones,
are not projected out of B.

[Definition] Let Z be the entire port-source map from the two codes,
including all three terms in (5) and both compensators. Write X=PZ,
Y=QZ and x=||X||, y=||Y||. A conservative bound is
||Z||<=3||M_L||||M_R||. Equalization makes X the bare source map.

[Proved] Block inversion gives the operator identity

    s Z^dag(A+sW)^(-1)Z
      =(X-s C^dag D^(-1)Y)^dag
        (B-s C^dag D^(-1)C)^(-1)
        (X-s C^dag D^(-1)Y)+s Y^dag D^(-1)Y.                    (12)

[Proved] Since D>=delta and ||C||<=L, the condition t L^2<=a/2 gives
||(B-s C^dag D^(-1)C)^(-1)||<=2/a. The resolvent identity then gives

    ||s Z^dag(A+sW)^(-1)Z-X^dag B^(-1)X||
       <= 2t L^2 x^2/a^2
          +(2/a)(2tLxy+t^2 L^2 y^2)+t y^2.                    (13)

[Proved] To derive (13), the inverse difference contributes the first
term; replacing X by X-s C^dag D^(-1)Y contributes the two terms in
parentheses; the last term of (12) contributes t y^2. This argument
is in operator norm, not just for a diagonal logical matrix element.
Thus it applies to the noncommuting hopping core.

[Proved] At w=0 the charge minima are the two end-hole blocks, of
energy U. Other configurations cost at least 2U. Let P_end project
onto those two minima. The eigenvector equation on its complement gives,
for N=k+p and w(N-1)<U,

    ||(I-P_end)chi_p||<=w(N-1)/(U-w(N-1)).                     (14)

[Proved] Every occupied bond is present in both end configurations.
Therefore

    0<=1-d_(p,j)<=zeta^2,
    zeta=max_(p=1,2) w(k+p-1)/(U-w(k+p-1)).                   (15)

[Proved] Put L_bond=2 sum_j ||B_j|| for a pair of cores. Equations
(9),(15) give ||B-B_bare||<=zeta^2 L_bond. Both inverses have a positive
lower bound a, decreased if necessary. A second resolvent identity gives

    ||X^dag B^(-1)X-X^dag B_bare^(-1)X||
        <= x^2 zeta^2 L_bond/a^2.                            (16)

[Proved] Combining (13),(16), for all four ordered pairs p,q,

    s K2_(p,q)=K2_bare+E_(p,q),
    ||E_(p,q)||<=E_fast(t)+x^2 zeta^2 L_bond/a^2.             (17)

[Proved] The bound includes the scalar. It proves arbitrarily accurate
equalization by first choosing w/U small and nonzero, then s/delta small
and nonzero. It does not prove equality at finite values. A scalar
coefficient c_(p,q) contains a length-density interaction
c_22-c_21-c_12+c_11, which cannot be removed by a chemical-potential shift.

[Proved] The shortest charge path between the two end minima uses kp
swaps. The off-diagonal Schur complement first appears at that order.
All shortest path contributions have the same sign and nonzero positive
denominators, while reflection equates the diagonal entries. Thus

    delta_ch,p = C_(k,p) w^(kp)/U^(kp-1)+O(w^(kp+1)/U^(kp)),
    C_(k,p)>0.                                                (18)

[Proved] In particular the slow p=2 charge splitting is of order
U(w/U)^(2n+4). Treating U as the wrapper gap would give a false hierarchy.

[Proved explicit coefficient] If a shortest path has intermediate
occupation configurations S_1,...,S_(kp-1), let I(S) be their interface
counts. Then

    C_(k,p)=2 sum_(shortest paths) product_(a=1)^(kp-1) 1/[I(S_a)-1].

[Proved] This is a finite positive rational function of the integer
parameters k,p, defined by the ordinary adjacent-swap graph. All intermediate
interface counts exceed one. It supplies an explicit coefficient in (18)
without assuming an n-independent prefactor.

## 4. Selection rules for higher terms

[Proved] The unperturbed core preserves vacuum versus excited species.
Vacancy transport and the marker move also preserve this distinction
for each ordered core. A port insertion changes it in both neighboring
cells. Hence PVP=0. Two different boundary insertions leave an excited
outer cell, so no three-cell term occurs at second order in lambda.

[Proved] A word with just one V and any number of marker insertions
cannot return all tokens to vacuum. In particular the mixed u lambda
coefficient vanishes. A pure marker excursion stays in the vacuum-token
sector, where every core matrix is the identity on its data string.
Its effective operator is consequently scalar on the payload, although
it can give length interactions at order u^2.

[Proved] On an open cell chain, multiplication by (-1) for an excited
token on alternate cells sends V to -V while fixing the unperturbed
Hamiltonian and marker motion. Odd powers of lambda vanish. This
argument is not asserted on an odd periodic cell ring.

## 5. The leading Gram matrix for arbitrary n

[Definition] Set t_d=0, phi=pi/(2n), m=1 for the explicit schedule below.
Temporarily expand the second-order bond in g. Its first nonconstant
quadratic hopping is g^2 T2(D). Write z=1/D and p_n=2n+4. Eliminating
the high A, x and B vertices leaves the real stub matrix

    A_n=m I+diag(3,2,...,2,1)-adjacency_path, n>=2,
    A_1=[m+2].                                               (19)

[Proved] A_n>=m I. The diagonal 3 at the first stub site counts its
data, A and switch edges. The diagonal 1 at the last site counts its
single data edge. For n=1 there are just the A and switch edges.

[Definition] Let K(t)=exp(-t A_n), and number data sites j=0,...,n-1.
Define

    f_j(t)=integral_0^t K_(0,j)(u)K_(j,0)(t-u) du,
    Gamma_ij=integral_0^infinity f_i(t)f_j(t) dt.             (20)

[Proved] The norm bound ||K(t)||<=exp(-mt) proves convergence. The
unique shortest path from 0 to j gives

    K_(0,j)(t)=t^j/j!+O(t^(j+1)),
    f_j(t)=t^(2j+1)/(2j+1)!+O(t^(2j+2)).                    (21)

[Proved] In a linear combination of the f_j, the least index with a
nonzero coefficient has a nonzero lowest Taylor coefficient. Thus the
f_j are linearly independent. The integral of the square of any nonzero
real combination is positive, and Gamma is positive definite.

[Proved, source calculation] A first-order conversion in the left
source heat kernel has, after eliminating A, the factor
-g z^2 f_j(t)/(n+2). In the right source it has
-g z^(2n+2) f_j(t)/(n+2), multiplied by the neutral rail phase and the
right-port parity string. The two arms both reach the stub at site 0.
For charged rails the transfer is real; the neutral right transfer has
phase exp(+-i n phi). The opposite charged conversions have opposite
parity signs. In (5) the two allowed interference terms therefore sum
to -2 sin(n phi), rather than canceling as they do at zero express phase.
After the ordinary Jordan-Wigner identification of logical fermions,

    T2(D)=-c_n z^p_n Gamma+O(z^(p_n+1)),
    c_n=2 sin(n phi)/(n+2)^2=2/(n+2)^2>0.                  (22)

[Proved, fermion signs] The neutral-token gauge at data position j
multiplies its data vector by Z_0...Z_(j-1), on either rail. At A this
product is empty; at B it is the full data parity. Each departure-Z
hop becomes a scalar under this gauge. Conversion at j becomes a
Jordan-Wigner annihilator or creator, since its local bit flip is
preceded by precisely that string. Express phases can then be removed
on the tree and retained at the B source as exp(+-i n phi). Thus the
two-cell term with one conversion in each cell is a fermion bilinear,
with the right-cell grading supplied by the B parity. This establishes
the Gaussian nature of the g^2 hopping and its relative particle/hole
signs, rather than assuming it from a one-particle fit.

[Proved, justification of the expansion] The coefficient of g^2 in the
pair resolvent is R0 C R0 C R0, where C is the sum of the two unit
conversion matrices. For an intercell one-particle matrix element one
conversion must occur in each cell. Each source must reach a data vertex
and return. The A arm requires two high-vertex factors; the B arm requires
2n+2. Any excursion that never reaches the low-low pair sector has extra
high resolvents and starts at a higher order. Schur elimination of the
high vertices thus gives no power below z^(2n+4). Its leading factors
are the real stub propagators in (20). This also derives (22) directly
from the finite pair resolvent, without interchanging an uncontrolled
large-D limit with an infinite-time integral.

[Proved] At this order and t_d=0 the one-cell quadratic matrix vanishes
at half filling. In a one-cell term the other cell has no conversion;
the two conversions in the active cell have opposite particle and hole
contributions. The charged propagators coincide, and a phase on a tree
cancels in a return amplitude. Their difference is zero. Terms with
four or more conversions are included in the norm remainder below.

[Numerical, independent extraction] check_asymptotic.py computes
-Z^dag R0 C R0 C R0 Z by three sparse solves in dimensions 576 and 5184
at n=1,2. The normalized matrix -D^p_n T2/c_n approaches Gamma.
This calculation takes the coefficient directly and does not divide a
small difference of two finite-g Hamiltonians by g^2.

## 6. An explicit finite express-energy schedule

[Proved] Gamma can be defined by finite rational linear algebra. Use
column vectorization, E=|0><0|, and

    Bcal = [[I tensor A_n, -I tensor E],
            [0,            I tensor A_n]],
    v=(0,vec I),
    Bcal P+P Bcal^T=v v^T.                                  (23)

[Proved] Bcal has positive eigenvalues, so this Lyapunov equation has
a unique solution. Its first block evolves as
Y(t)=integral_0^t exp(-(t-u)A_n) E exp(-u A_n)du, whose diagonal is f_j.
Taking the corresponding n diagonal-coordinate rows and columns of P
therefore gives Gamma. For rational m, Gamma is a positive rational
matrix. Define the explicit positive number

    rho_n=det(Gamma)/(tr Gamma)^(n-1) <= lambda_min(Gamma).    (24)

[Proved] The inequality follows by bounding the other n-1 eigenvalues
by tr Gamma in their product. For n=1 the empty power is one. Equation
(24), not an empirical exponential fit, will control the narrowest band.

[Definition] A deliberately conservative analytic bound is obtained with

    L=2m+12,   r=min(1,m)/(16 L^2),
    B0=288/[(n+2)^2 m^3].                                   (25)

[Proved] At g=0 the pair source Hamiltonian has the form z^(-1)N+H_f,
where N counts tokens on high vertices, N>=1 on the high subspace,
||H_f||<=L, and the low block is at least 2m. For complex |z|<=r, the
high inverse has norm at most 2r. Its Schur correction is at most
2r L^2<=m/8. Completing the block inverse gives ||R0(z)||<=2/m.
The source norm is at most 3/(n+2), and ||C||<=2. Hence each matrix
element of T2(z) is analytic in this disk and bounded by B0.

[Proved] Cauchy's coefficient bound and the vanishing coefficients
below p_n from (22) give, for 0<z<=r/2,

    ||z^(-p_n) T2(z)+c_n Gamma||
        <=2n B0 z/r^(p_n+1).                               (26)

[Definition] It is sufficient to choose the fixed finite number

    D_n >= max(2/r, 8n B0/[c_n rho_n r^(p_n+1)]).            (27)

[Proved] Then the Hermitian part of -D_n^p_n T2 is at least
(3/4)c_n rho_n I. Thus a finite energy, rather than D=infinity, suffices.
The bound is intentionally large; it is not a practical simulation point.

[Definition] Put a_n=c_n rho_n D_n^(-p_n). Choose

    0<g<=m/4,
    g^2 <= a_n m^5/(10^5 n),
    eta_n=g^2 a_n,
    epsilon=eta_n/4,
    0<t_d<=eta_n m^2/(10^5 n).                              (28)

[Proved] At fixed D the convergent Neumann expansion of the pair
resolvent has ||R0||<=1/(2m) and ||C||<=2. Charge neutrality permits
only even conversion counts in a code-to-code process. The norm of
the part beginning at g^4 is bounded, for g<=m/4, by

    9(2g)^4/[(2m)^5(1-g/m)] <= 6g^4/m^5.                  (29)

[Proved] Each loop edge is positive with norm 2t_d; the two-cell
change is at most 4t_d. The resolvent identity bounds the change of
the full second-order bond by 9t_d/m^2. Extracting its quadratic
Bloch matrix increases this bound by at most 16n. One can obtain that
constant by using the orthogonal traceless fermion-bilinear basis:
each coefficient is bounded by four times the input operator norm,
and the onsite plus two hopping matrices have at most 4n such row sums.

[Proved] Equations (27)-(29) consequently keep the total quadratic
error well below eta_n/8. At the two momenta where the phase-rotated
hopping enters as plus or minus its Hermitian part, every eigenvalue
has opposite signs, with a strict margin exceeding eta_n. Ordered
Hermitian eigenvalues are continuous functions of momentum. Each of
the n ordered bands therefore crosses the common field on both arcs.
All parameters in (27),(28) are finite and nonzero at each fixed n.

[Proved limitation] This sign argument does not guarantee distinct roots,
simple velocities, distinct sea centers, or a quartic nonresonance margin.
These conclusions cannot be appended to (26) on the basis of positivity.

## 7. Why the all-n Sidon step remains open

[Proved counterexample to an inference] Take T=-diag(t_1,...,t_n), with
distinct t_a>0 and a common field 0<epsilon<2 min t_a. All bands cross,
their velocities are nonzero, and T is definite. Nevertheless each sea
is centered at zero. Every pair of opposite Fermi points sums to zero.
For n>=2 there is a number-conserving interband Cooper-transfer resonance.
Exact total U(1) does not forbid it.

[Proved] Distinct centers alone also do not imply that every sum of two
Fermi points is distinct modulo 2pi. A finite list of further relations
must be excluded. C20's volume argument varies the Fermi points as
independent coordinates. A few-parameter spatial core traces a much
smaller subset of that space, so that argument cannot be pulled back
without proving that no required resonance is an identity on this family.

[Assumed, missing lemma] A completion needs, for every n, a point or a
finite parameter search with a proved termination bound at which the
2n Fermi zeros are simple and every nontrivial charge-conserving pair-sum
mismatch is nonzero. Analyticity alone gives such a search only after
each possible identity has been ruled out. The loop coupling breaks
the obvious tree symmetry, but that is not a proof of this missing
nonidentity statement for arbitrary n.

[Numerical] The n=1,2 witness in band_checks.json passes those tests.
No witness or nonidentity proof for every n is supplied here.

## 8. Scale closure, and its exact conditional scope

[Proved] All microscopic occupation, filter, selector and helix gaps
are positive at fixed n and fixed w>0. The selector has b=n+3 and
margin of order t_sel b^(-4). The helix gap is at least
J_h/[(b+1)(2b+1)]. The charge gap is (18). The core gap is (4).
The full wrapper gap is bounded below by C8's expression

    delta_ch d_min gamma /
       [(sqrt(delta_ch)+sqrt(M_core))^2+d_min gamma],           (30)

with the filter included in gamma and M_core. Additional occupation and
selector exclusion bounds are minimized with (30). Denote the resulting
positive unperturbed many-cell exclusion gap by Delta_n.

[Definition] Suppose a particular finite core has a certified minimum
quadratic velocity v_*>0 and a finite payload mismatch delta_*>0. These
are dimensionless coefficients in K2_bare before lambda^2/s. Let

    tau_n = min(eta_n/100, v_* delta_*/100, v_*/100).             (31)

[Assumed] The v_*,delta_* premise is verified numerically for the small
examples, but it is not established for all n. The unconditional band
crossing hierarchy uses tau_n=eta_n/100 and makes no stable-phase claim.

[Proved] Choose 0<w/U<1/[2(k+1)] small enough that the right side of
(16) is below tau_n/4. This is an explicit scalar inequality in w/U;
all its other constants are finite. It requires roughly
w/U << sqrt(tau_n)/poly(n), not merely a small fixed vacancy hopping.
Next choose 0<s/delta small enough that (13) is below tau_n/4.
The compensator norms and the tiny tunneling gap are included in that
inequality. This order of choices avoids an illicit uniform limit.

[Proved] After this step every scaled second-order bond differs from
K2_bare by less than tau_n/2, including its scalar part. The physical
logical scale and narrowest crossing scale are

    J_payload=(lambda^2/s)||T_bare||,
    J_narrow >= (lambda^2/s) eta_n/2.                         (32)

[Definition] Set u=(lambda^2/s) j_l/beta_n, with a fixed j_l>0 chosen
in the desired length-fluid regime. Thus J_l=(lambda^2/s)j_l. The large
length-independent scalar in K2 changes the cell chemical potential
and is retuned there; the residual scalar length interaction is covered
by (17). The uniform field is implemented by
(lambda^2/s) epsilon Q_total. It acts equally on both lengths.

[Proved, finite-volume perturbation estimate] Let v=||V_bond|| and
t_move=||T_move,bond||, including helix factors. On any fixed cluster
of r cells put q_r=(r-1)(|lambda|v+|u|t_move). For q_r<Delta_n/4,
the contour-resolvent expansion for the isolated unperturbed code band
converges. The discarded terms beyond second order are bounded by
C_r q_r^3/Delta_n^2, with a finite C_r obtained from the contour of
radius Delta_n/2. This follows by expanding the resolvent in a norm
geometric series on that contour and integrating term by term.

[Proved] Since u is proportional to lambda^2, the ratio of that bound
to lambda^2/s tends to zero as lambda tends to zero. More refined
selection rules remove odd lambda and u lambda and leave terms of
orders u^2/Delta_n, lambda^2 u/Delta_n^2 and lambda^4/Delta_n^3.
Either bound suffices: choose a fixed nonzero lambda after all previous
parameters so that these errors are below (lambda^2/s)tau_n/4.
There is no requirement to send any coupling to zero with ring length.

[Assumed thermodynamic input] A local massive-sector elimination and
the weak-coupling bosonization stability input of C13/C20 are required
to turn these local scales into an infrared phase. A global finite-ring
code-band perturbation radius does not stay uniform with the volume,
so the preceding contour estimate is not a thermodynamic proof.

[Proved, length-density choice conditional on payload points] For a
given nonresonant payload witness, only finitely many values of nu are
forbidden by mixed length backscattering and payload bilinears. Remove
their points from [1/4,3/4], take the midpoint of a largest component,
and then a sufficiently small irrational displacement. The resulting
mixed mismatch is strictly positive. C13 gives the explicit bound
pi/[4(K+1)] with K=4 binomial(2n,2), which can be used conservatively
also in the U(1)-preserving subfamily. No new all-n payload witness is
created by this choice of nu.

[Heuristic, conditional conclusion] Given the missing payload witness
and the thermodynamic inputs, this closed hierarchy would yield the
length mode plus n payload modes. Without them, the proved conclusion
is the mobile finite construction and band crossings, not c=n+1.

[Proved scale interpretation] The observed loss of a factor 50 to 100
per flavor is not an asymptotic bound. Equations (23)-(28) give an actual
n-dependent schedule, much smaller than that empirical extrapolation.
No lower bound uniform in n or experimentally accessible energy window
is claimed.
