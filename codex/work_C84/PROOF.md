# C84. Eleven letters with one unused boundary bit

[Proved finite ingredients; assumed infrared transfer] This construction
has alphabet a=11 and a certified member with n=11 payload flavors.
It retains the folded carrier. The saving comes from eliminating the
marker's data bit and losing only one payload flavor at each cell length.
The finite result is not a microscopic central-charge theorem.

## 1. Alphabet, bond, and conserved quantities

[Definition] Write p for the number of nonmarker sites in a cell,
b=p+1 for its physical length, and n=p-1 for its payload flavors when
p>=1. This distinction matters. The certified member has p=12, b=13,
n=11. The alphabet is

    L_0, L_1, R_0, R_1, D_A, D_B, X_A, X_B, E, M_A, M_B.

The marker set is M={E,M_A,M_B}. The charges are

    q(L_a)=q(R_a)=2a, q(E)=0,
    q(D_A)=q(X_A)=q(M_A)=-1,
    q(D_B)=q(X_B)=q(M_B)=+1.                            (1)

Here E is C80's E_0. E_1 is absent.

[Definition] The complete allowed ordered-pair table is

| Left letter | Allowed right letters |
| :-- | :-- |
| L_0 or any marker | L_0,L_1,R_0,D_A,D_B,X_A,X_B,E,M_A,M_B |
| L_1 | L_0,L_1,R_0,D_A,D_B,X_A,X_B |
| R_0,R_1,D_A,D_B,X_A,X_B | R_0,R_1,E |

Put K on every other ordered pair. Relative to C80 restricted to the
remaining letters, the only new forbidden pairs are L_1 followed by a
marker. They fix the unused boundary bit without adding a letter.

[Definition] Set z_A(a)=(-1)^a and z_B(a)=1. The table lists the output
ket times input bra coefficient. Include its Hermitian conjugate and
sum over the displayed indices a=0,1 and s=A,B.

| Input | Output | Coefficient |
| :-- | :-- | :-- |
| R_0 R_a | L_a R_0 | -1 |
| R_0 E | L_0 E | -1 |
| D_s R_a | L_a D_s | -z_s(a) |
| X_s R_a | L_a X_s | -exp(i phi_s) z_s(a) |
| X_s E | L_0 M_s | -exp(i phi_s) |
| D_s E | L_0 M_s | -t_d exp(i psi_s) |
| m D_s, m in M | m X_s | -1 |
| D_A R_1 | D_B R_0 | g |
| E R_0 | M_A D_B | lambda exp(-i pi/4) |
| E R_0 | M_B D_A | lambda exp(+i pi/4) |

There is no conversion on a token followed by E. The last two rows are
lambda V and are the whole port. The one-site diagonal o assigns Delta
to every D_s,X_s,M_s, an additional D to X_s,M_s, mu_sel to every marker,
and lambda^2 epsilon times a to L_a,R_a. Thus

    h = h_table + K Pi_bad + (o tensor I+I tensor o)/2.   (2)

[Proved] Equation (2) is a one-site translation-invariant Hermitian
nearest-neighbor bond on C^11. Every rule preserves the two-site charge
(1). For the conversion both sides have charge one; both port outputs
have charge zero. Every rule also preserves whether each of its sites
is a marker. Hence each on-site marker projector commutes with H, and
the cells are frozen. H_0 preserves excited-token number separately in
each cell. A port changes that number by one in each incident cell.
These statements hold on illegal words as well as legal words.

[Numerical] `bond11.py` implements the complete 121 by 121 bond. At a
generic diagnostic point it has 166 nonzero entries. Hermiticity, the
charge commutator, and both marker commutators vanish in the literal
matrix check. `alphabet_audit.py` also checks the charge equations using
exact rational row reduction, independently of the phase coefficients.

## 2. Legal words and the saved letter

[Definition] For bits a_0,...,a_(n-1), p=n+1, the words with a head at
j=0,...,p-1 are

    L_(a_0)...L_(a_(j-1)) T
    R_(a_j)...R_(a_(n-1)) E,                          (3)

where T is R_0,D_A,D_B,X_A,X_B. At j=n the R tail is empty.
The three terminal words are

    L_(a_0)...L_(a_(n-1)) L_0 m,  m=E,M_A,M_B.        (4)

The first token choice and m=E describe the vacuum walker. The other
choices describe the two excited species and their two rails.

[Proved] The pair table characterizes (3)-(4). Starting after a marker,
one sees an L prefix. The first non-L letter is R_0 or an excited token.
After it every nonmarker letter is R. The right marker must then be E.
If all data letters are L, the final one must be L_0, and the marker
can be any of the three. The first non-L letter locates the head even
when it is R_0. No bit or clock position is identified. For p>=1 the
legal cell dimension is

    (5p+3) 2^(p-1).                                    (5)

For p=0 the three marker states remain.

[Proved] Take C80's p-bit legal words, fix their last logical bit to
zero, and remove its conversion. Equations (3)-(4) are their images.
Every surviving move agrees with the parent, including its phases and
Jordan-Wigner sign. The last token hop crosses a fixed zero. Conversion
at that position is absent. This is a unitary equivalence of the legal
finite-cell graphs for every p>=1. It is not restriction to an invariant
subspace of the original parent with its last conversion still present.
The literal bond has been changed by deleting E_1 and its transitions.

The clock has p+1 positions, although only p-1 payload bits remain.
Therefore its constrained vacuum is

    E_v(p)=-2 cos(pi/(p+2)),
    v_j=sqrt(2/(p+2)) sin(pi(j+1)/(p+2)), 0<=j<=p,
    G|a>=sum_j v_j |j;a>,  alpha=v_0=v_p.               (6)

On this code Q=2 sum_j a_j. The vacuum gap is
2[cos(pi/(p+2))-cos(2pi/(p+2))]. Each bit string gives one ground state.
At finite K this band can split. It is not claimed to remain degenerate.

[Proved] At fixed n the construction encodes n qubits in n+2 physical
sites. Losing a constant number of payload bits has no asymptotic cost
for the fixed-alphabet mechanism. This is why marker bit storage is
not logically necessary for unbounded flavor count.

[Numerical and exact finite enumeration] Legal words, completeness, and
injectivity are checked through p=7, with 2432 legal states at that size.
At n=1,2 the independent C75 bit-carrying parent, with its last bit fixed
and its last conversion compressed away, gives identical H_0 and V
matrices. The pair dimensions are 676 and 5184. These comparisons do
not assume the shifted-bit rules when assembling the reference.

## 3. Scalar reduction with an inactive endpoint

[Definition] At p=12 choose

    D=1000, Delta=4+E_v(12),
    phi_A=pi/24, phi_B=0, psi_A=pi, psi_B=pi/2.          (7)

For diagnostic members, replace 12 by p and phi_A by pi/(2p).
Within the final microscopic Hamiltonian all couplings stay fixed at (7),
including on cells of other lengths. The length selector does not retune
the bond cell by cell.
Gauge the A carrier by the parity of the crossed bits and remove its
express phases. Both species then have the same scalar graph

    M, X_(p-1),...,X_0, D_0,...,D_(p-1).

In energies above E_v, its diagonal is 4+D on M,X and 4 on D. Path
edges are -1. The closing entry H_(M,D_(p-1)) is -i t_d. The conjugate
entry closes the reverse edge. There are 2p+1 vertices. Only D_0 through
D_(p-2) convert, so there are n=p-1 active sites.

[Proved] Let

    R=(H tensor I+I tensor H)^(-1), x=R|M,D_0>.

Define U_i by retaining row D_i of x, and V_i by retaining column D_i,
viewing each as a vector on the product graph. For i,j=0,...,n-1 set

    C_ij=U_i^dagger R V_j,
    A_ij=U_i^dagger R U_j+V_i^dagger R V_j.              (8)

The gauged code source, up to a common phase and alpha^2, is
|B,A>-P_left|A,B>. One conversion gives the two orthogonal sectors

    Y_AA=sum_i U_i f_left,i^dagger-sum_i V_i f_right,i^dagger,
    Y_BB=P_left(sum_i U_i f_left,i+sum_i V_i f_right,i).

The g^2 coefficient of the port effective bond is
-alpha^4 g^2 (Y_AA^dagger R Y_AA+Y_BB^dagger R Y_BB).
In the BB term the two adjacent P_left factors cancel. In the AA term
use f_i f_j^dagger=delta_ij-f_j^dagger f_i. The result is scalar plus
number-conserving quadratic, with no quartic term. Its hopping and
combined onsite matrices are

    T=-2 alpha^4 g^2 Re C,
    mu=-2i alpha^4 g^2 Im A.                            (9)

These are the active principal blocks of the p-conversion parent
coefficient, calculated with the full 2p+1-vertex propagator. One must
not shorten the scalar carrier to 2n+1 vertices. The unused D endpoint
and its X partner still propagate and still close the flux loop.

[Numerical] Full pair Fock-space calculations at n=1,2 give a relative
non-Gaussian residual at most 1.83e-16 for this coefficient. An independent
scalar resolvent calculation reproduces T to relative error below
1.2e-15. The combined onsite cancellation has relative error below
2.9e-11 at the diagnostic couplings, where it is a subtraction of larger
numbers. This numerical cancellation is not used as the proof of (9).

[Proved] In the selected cell the source norm is at most sqrt(2) alpha^2,
the pair conversion norm is at most 2g, and the free excited pair gap is
at least three. Odd conversion orders do not return the source species.
Consequently the even Neumann tail satisfies, for g<=1/2,

    ||H_2(g)-H_2(0)-g^2 H_22||
      <= (2 alpha^4/3)(2g/3)^4/[1-(2g/3)^2]
      <= alpha^4 g^4/4.                                (10)

Removing one conversion cannot increase the conversion norm. The exact
finite-g operator generally has interactions.

## 4. Exact arithmetic certificate and finite schedule

[Definition] Define

    d_0=1, d_1=1004, d_(j+1)=1004 d_j-d_(j-1),
    beta=1/d_13, alpha=sqrt(2/14) sin(pi/14),
    t_d=beta 2^-55, g=2^-216, K=2^2048, lambda=2^-480,
    E_unit=alpha^4 beta^2 g^2, epsilon=E_unit 2^-101.     (11)

The physical data field is lambda^2 epsilon. These exact definitions,
not a floating-point microscopic array, specify the small couplings.
The marker chemical potential is specified in Section 6.

[Proved, computer-assisted] `run_certificate.py` regenerates the complete
certificate. The C80 arithmetic is copied locally, with the following
structural change. Its argument n now counts active flavors. The express
block has n+2 vertices, the data block has n+1 vertices, and only its
first n rows and columns receive conversion sources. Thus n=11 uses
the 25^2=625-dimensional product graph. No parent output file was used
as the new certificate.

The tree solve uses 256-bit fixed-point integers and the finite loop
uses 512 bits. The Jacobi contraction is at most 1/2+t_d/4. With
bit_length(d_13)+512+20 iterations, both source and resolved columns
are enclosed by ten units per real or imaginary component. The product
radius is 10(|x_r|+|x_i|+|y_r|+|y_i|)+200. Integer LDL congruences are
checked as identities. Dyadic rescaling gives a metric with diagonal
entries between one and four. Its propagated row error is below
1.144e-100. This bounds a Hermitian matrix norm.

Rational diagonal dominance proves a Hermitian hopping margin greater
than 0.56526 in that metric. Numerical eigenvalues only propose root
brackets. At every endpoint t=tan(k/2) is rational and

    cos k=(1-t^2)/(1+t^2), sin k=2t/(1+t^2).

An integer approximation to a numerical eigenbasis is proved invertible.
The rational transformed matrix is strictly diagonally dominant with
certified diagonal signs after subtracting all entry errors. Its inertia
is therefore certified without trusting the numerical eigenvalue signs.

[Proved, computer-assisted] There are 22 disjoint brackets. Their inertia
changes are 0 to 1 through 10 to 11, then 11 to 10 through 1 to 0.
Each bracket contains a sign-changing determinant root. Multiplying
det h(k) by exp(11ik) gives a polynomial of degree at most 22. These
22 distinct roots exhaust its degree, so all roots are simple and no
additional roots occur for this nearest-neighbor quadratic pencil.

Rational phase midpoints are sorted by exact orientation. A root phase's
chord error is bounded by its rational t-bracket width. A product of two
unit phases has error at most the sum of the individual errors. Subtract
these errors from adjacent max-component chord distances. The resulting
lower bounds are

    delta_2 > 0.0164202883,
    delta_4 > 9.11916e-7.                               (12)

All 22*23/2=253 unordered pairs with replacement are included. In
particular, the eleven sea centers are distinct modulo pi. The endpoint
inertia certificates allow an additional Hermitian Bloch error of 2^-35
in the dyadic metric. `root_certificate_exact_n11.json` retains all
rational brackets and margins. This is finite low-order nonresonance,
not independence of arbitrary integer momentum relations.

[Numerical rejected parameter] The first run reused C80's loop exponent
66. Its smallest numerical pair margin was about 4.46e-10; brackets of
half-width 1e-9 did not certify it. The rational comparison rejected two
overlapping pair enclosures, with lower bound -6.15e-9. Exponent 55
passes. `rejected_parameters.json` preserves the rejected parameter and
failure. Nothing in the final result uses that failed acceptance test.

## 5. Finite K and the physical error budget

[Proved] For one selected cell write H_K=K N_bad+B with ||B||<=B_0=2^16.
There are 13 sites, each onsite norm is below 1007, and the move norms
per bond sum to less than 64. The field is below one. The new boundary
pin belongs to N_bad, not B. With P=ker N_bad and Q=1-P, for |z|<=2B_0
and K>=8B_0,

    QH_KQ-z >= K/2,
    ||PBQ(QH_KQ-z)^(-1)QBP|| <= 2B_0^2/K.              (13)

The legal vacuum code has dimension 2^11 and gap greater than 1/16.
Use the polar overlap isometry, separately in charge sectors, to compare
its finite-K spectral subspace with the constrained code. For two cells,
replace B_0 by 2B_0 in (13), include the energy shift, and divide by the
vacuum gap for the spectral projection. Conservative bounds are

    max_a |E_a-2(E_v+mu_sel)| <= 2^36/K,
    ||G_K-G_infinity|| <= 2^44/K.                       (14)

These bounds control a band, not one exactly degenerate eigenvalue.

[Proved] PVP=0 because the code has zero excited tokens. Two different
adjacent ports cannot restore the outer cells' excitation numbers, so
second order acts on at most two adjacent cells. On two cells the parity
of the left excited-token number negates V and preserves H_0. The
canonical effective operator is even in lambda.

For finite-K vacuum eigenvectors with eigenvalues E_a, the Hermitian
second-order coefficient is

    (H_2)_ab=-1/2 sum_q V_aq V_qb
               [(E_q-E_a)^(-1)+(E_q-E_b)^(-1)].         (15)

The sum includes the full (1,1) source sector. A single common denominator
would omit the finite-K splitting.

[Proved] The source-sector Schur correction and the shifts (14) are
less than 2^39/K. The resolvent identity, including the high-high block
bounded by 2/K and both mixed blocks, gives inverse error at most 2^42/K.
The ideal inverse is extended by zero outside the legal space. Thus the
full source, including every unintended port process, is retained.
Combining the code rotation with ||V||=sqrt(2), each column in (15)
changes by less than 2^48/K. There are 4^11 columns, so passage from
column norms to operator norm costs at most 2^11. Symmetrizing and
contracting with the source is covered by

    ||H_2(K,g)-H_2(infinity,g)|| <= 2^64/K.              (16)

The zero-port effective onsite term, its two copies, the field correction,
and lambda^2 times (16) are together less than 2^80/K in physical units
for lambda<=1. The bound is not multiplied by lambda^2, since it includes
the zero-port splitting.

[Proved] For completeness, an isolated two-cell effective operator has
a remainder at most 2^50 lambda^4. A contour of radius 1/64 around the
vacuum band and an analytic port disk of radius 2^-12 give a convergent
resolvent series. The spectral projector and its intertwiner have norm
below two there. The effective operator relative to its central energy
has norm below one. Cauchy's coefficient estimate, evenness, and a
geometric tail bound give 2^49 lambda^4 for lambda<=2^-13; the factor
two left in the displayed bound covers finite K. This estimate concerns
two cells. It is not a uniform many-cell remainder theorem.

[Proved] Let tau=tr(W W^dagger) for the certificate's dyadic congruence.
A quadratic coefficient in normalized Hilbert-Schmidt projection is
bounded by four times the operator norm. Summing the onsite and hopping
blocks gives a conservative Bloch bound 16n times that norm. Equations
(10), (14)-(16) therefore give

    error_g = 4n beta^-2 g^2 tau,
    error_K = 16n tau (2^80/K)/(lambda^2 E_unit),
    error_lambda,local = 16n tau 2^50 lambda^2/E_unit.    (17)

The denominator lambda^2 E_unit in error_K is essential. Dividing only
by E_unit would fail to budget the zero-port code splitting.

[Proved, computer-assisted] The inequality sin(pi/14)>=2/14 gives
alpha^4>=64/14^6>2^-17. `finite_budget.py` uses this rational bound, the
integer d_13, and exact Fractions for every acceptance. It obtains

| Contribution | Normalized upper bound |
| :-- | :-- |
| Finite conversion strength | 3.380e-21 |
| Finite K, including zero-port splitting | 7.946e-59 |
| Local two-cell port remainder | 2.518e-29 |

Their sum is less than one eighth of 2^-35. The exact fractions are
retained in `finite_budget_exact.json`. The root routine's
finite_g_power=-199 is an inherited scalar threshold, not the final
choice g=2^-216; the budget here uses the two-species prefactor in (17).
In particular, all root brackets
and both margins survive for the quadratic projection of the finite-K,
finite-g effective Hamiltonian through second order in lambda. The local
port remainder is an additional finite check, not a claim about the
range or root degree of the full many-cell effective Hamiltonian.

## 6. Uniform zero-port selection

[Definition] For integer b>=1 let

    a_b=2[1-cos(pi/(b+1))],
    delta_b=(b+1)a_b-b a_(b+1),
    mu_sel=2-(delta_13+delta_12)/2.                      (18)

[Proved application of C80] For any fixed marker configuration let J
embed the eleven-letter Hilbert space into C80 by omitting E_1. At zero
port, including all illegal strings and all excited-token numbers,

    H_11(K)=J^dagger H_12(K) J+K sum_i Pi_(L_1 M),i.    (19)

The added term is positive. C80's uniform all-string bound for k data
sites is

    H_12(K) >= mu_sel+E_v(k)-(k-1)_+/(K-4).

Compress this form inequality by J and use (19). For the upper bound
choose the C80 vacuum with its last logical bit zero, which lies in the
new legal space and has expectation mu_sel+E_v(k). Hence, also for
k=0 by direct inspection,

    mu_sel+E_v(k)-(k-1)_+/(K-4)
       <= E_11,K(k) <= mu_sel+E_v(k).                  (20)

This is a bound uniform in k, not a large-cell use of the finite-cell
norm B_0. C80's underlying argument splits the offending vacuum moves
into partial isometries that create one predecessor bad pair. Cutting
m bad pairs leaves at most m+1 fragments of kinetic row sum two.
Young's inequality with eta=K-4 absorbs these moves into K N_bad and
leaves (k-1)_+/(K-4). Equation (19) transfers that already established
bound without a new fragment assumption.

[Proved] At K=infinity the density is e_b=(mu_sel-2+a_b)/b, and

    e_(b+1)-e_b=[2-mu_sel-delta_b]/[b(b+1)].             (21)

The function a(x)=2[1-cos(pi/(x+1))] has

    a''(x)=4pi sin(pi/(x+1))/(x+1)^3
           +2pi^2 cos(pi/(x+1))/(x+1)^4 >0, x>=1.

Thus delta_b-delta_(b+1)=(b+1)(a_b-2a_(b+1)+a_(b+2))>0.
The midpoint (18) makes the density decrease to b=13 and increase
thereafter. The nearest competing density is one of its two neighbors.

[Proved, computer-assisted] Rational alternating-series bounds for pi
and cos establish

    min(e_12-e_13,e_14-e_13)>5.4743300e-5,
    e_13<-0.0072057833,
    vacuum gap>0.14791808,
    mu_sel=1.8561806413628503... .

Equation (20) changes any density by at most 1/(K-4), below the selector
gap at (11). The same chemical potential therefore selects length 13
uniquely at finite K. On marker-free rings the bond is unchanged from
C80, and its density is at least -1/(K-4), above the selected-cell trial
density. Compatible rings have the selected frozen cells and their
translation-related sectors at zero port.

[Numerical] The selector form bound is checked in full two-cell invariant
spaces at p=1,2,3, dimensions 94, 2221, 42424. The lowest form eigenvalues
are no smaller than -8.0e-15. Actual energies exceed (20). Exact local
context enumeration also confirms that only a vacuum move changes a
bad-pair projector, and only at its predecessor.

## 7. Diagnostic finite-size checks

[Numerical] `checks.py` closes all legal two-cell words under every
literal rule. It includes all reachable illegal words, with no cutoff
on their number of violations. The invariant dimensions at n=1,2 are
2221 and 42424. These are small diagnostic calculations, not ED of
the extreme schedule (11).

For D=4, Delta=0.5, g=0.04, t_d=1e-4 and K=100, all four or sixteen
code levels agree with (15) with the following energy error/lambda^4.

| n | lambda=0.2 | lambda=0.1 | lambda=0.05 |
| :-- | :-- | :-- | :-- |
| 1 | 0.0017438122 | 0.0017393831 | 0.0017382762 |
| 2 | 0.0005318533 | 0.0005304506 | 0.0005301001 |

The largest eigenpair residual is below 1.6e-13. At n=2 the zero-port
pair-band width is 0.012741366 at K=100 and 0.0012761712 at K=1000.
The corresponding H_2 errors, after polar overlap alignment, are
0.0024119506 and 0.0002498303. Both decrease at the expected 1/K scale.
The code gap remains positive. Finite-g residuals at g=0.08,0.04,0.02
also have the expected fourth-order scaling. `small_checks.json`
retains the complete outputs and sizes.

## 8. Lower bounds on the alphabet

[Proved, restricted] Eleven is optimal among faithful on-site encodings
of this shifted-bit folded graph with a frozen on-site marker projector,
two locally distinct rail states per charged species, and both left and
right backgrounds independently storing a binary payload. The necessary
disjoint groups have sizes

    backgrounds: 2 orientations times 2 bit values = 4,
    carriers: 2 species times 2 rail states = 4,
    markers: vacuum and the two source species = 3.       (22)

This is a theorem for a specified encoding class, not for every possible
route-A mechanism. Here is why the distinctions in (22) are required.

First, changing a payload bit at a fixed head position gives orthogonal
code states, so its two symbols must differ. The orientation labels
locate the vacuum head and enforce a single sweep with a two-site
grammar. Merging L_0 and R_0 identifies different legal clock states.
Merging L_1 and R_1 admits extra zero-penalty paths and loses the scalar
vacuum clock. These are not errors removable by taking K larger.

Second, at a fixed head site and fixed payload the D_s and X_s states
are different vertices with different diagonal masses and hopping
phases. All their other sites agree in the faithful product encoding.
They must therefore have distinct on-site letters. The two species
also differ in additive charge, as shown below.

Third, an adjacent intercell port can change only its marker and the
first data site. Its neutral input and its two oppositely charged
outputs must be distinct states on that marker. Moving the distinction
to a neighboring data site would change the port mechanism. Markers
cannot be identified with nonmarkers while retaining the frozen on-site
marker projector. No bit-storage marker is needed, which is the saving
from twelve to eleven.

[Proved] The general additive charge of the literal nonzero rules has
three parameters. Write q(L_0)=q(R_0)=u, q(L_1)=q(R_1)=u+delta, and
q(E)=v. The switch makes q(D_s)=q(X_s). Conversion and the two port
rows imply

    q(D_B)-q(D_A)=delta,
    q(D_A)+q(D_B)=2u.

The terminal hops then give

    q(D_A)=q(X_A)=u-delta/2,
    q(D_B)=q(X_B)=u+delta/2,
    q(M_A)=v-delta/2, q(M_B)=v+delta/2.                 (23)

For a nontrivial logical charge delta is nonzero. Thus the two charged
carriers cannot collapse, and the three markers cannot collapse. The
four distinct nonmarker charge values u,u+delta,u-delta/2,u+delta/2
also prevent identifying any carrier with a payload background. Exact
rational row reduction finds rank eight on the eleven charge variables,
so the three-dimensional solution space is spanned by identity, marker
number, and (1). Equation (23) is not just an arbitrary chosen charge
assignment that a literal quotient could evade.

[Proved, restricted] If one retains the on-site charge roles, terminal hops, conversion
and neutral two-channel port equations (23), but does not require two rails or
two payload orientations, the charge and marker argument alone gives
only a>=7. It counts two payload charge values, two carrier charge
values, and three distinct marker states. Retaining both independent
payload orientations strengthens this to a>=9. The faithful two-rail
encoding strengthens it to a>=11. A replacement fold must be assessed
under the weaker bound, not assigned two rail letters by definition.

[Proved, more general but weak] For an unbounded dense-payload family
with one frozen on-site marker per cell and b=n+O(1), the code dimension
alone implies at least two nonmarker states. If d_data=1 and the fixed
marker space has dimension d_M, the whole fixed-cell Hilbert space has
dimension d_M, which cannot contain 2^n states for unbounded n. Since
a marker is disjoint from the data set, a>=3. Without the conversion,
port, and charge assumptions, this argument does not give seven or
more. It is not asserted for arbitrary encodings with a different
relationship between cell length and flavor count.

[Numerical, with exact word counts] All ten charge-preserving single
identifications for (1) were checked at p=3, n=2. D_s=X_s collapses
72 legal words to 60. L_0=R_0 collapses them to 64; its lowest four
vacuum energies at K=100 are -3,-2.61803,-2,-1.80194, instead of a
fourfold -1.61803. L_1=R_1 is injective on those words but has lowest
energies -2.020008,-1.624405,-1.618034,-1.549903. The other six merges
identify a marker with data and remove the required marker projector.
`alphabet_checks.json` records the list, spectra and exact charge rank.

[Proved, computer-assisted] Both background quotients have a zero-penalty
integer Rayleigh witness strictly below -2. For L_0=R_0 its quotient is
-300000000/100000000. For L_1=R_1 it is -201539920/100004126, which is
less than -2 because -201539920 < -200008252. The target vacuum energy
is -2cos(pi/5)>-2. Thus each quotient is undercut for every K, without
trusting a numerical spectrum. The vectors and zero-penalty basis order
are retained in the audit and its source. In each of the six marker/data
quotients, the graph connecting on-site letters changed by local rules
is connected. A diagonal on-site projector conserved by all these rules
must be constant on this graph, hence trivial. This excludes any
replacement diagonal frozen-marker label in those literal quotients.

[Heuristic] A ten-letter construction appears to require a different
spatial encoding, selector, port, or charge representation. This is a
useful design constraint, not a proof of impossibility. The alphabet
stores functions, not one compulsory letter for every descriptive role.
This construction itself combines the vacuum token with R_0 and the
vacuum marker with E. A lower bound that treats every function as a
disjoint on-site state would already overcount these examples.

## 9. Scope of the infrared claim

[Assumed] The conditional frozen-chain prediction is c=11>10. It uses
volume-uniform elimination of massive sectors and stability of the
weakly interacting multichannel Luttinger liquid. The exact U(1), the
finite low-order momentum margins and the controlled finite-cell errors
are ingredients for that inference. They are not the inference itself.

[Open] A common-bulk-subtracted dressed-defect comparison at nonzero
lambda and the many-cell infrared transfer remain necessary, as in
C75 and C80. Zero-port selection and a two-cell analytic remainder do
not prove them. No central charge was measured and no completed
microscopic counterexample is claimed. No independent referee was used.
The smallest alphabet certified here is eleven. Global minimality is
open.
