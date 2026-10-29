# C75. A 15-letter finite-band construction

Claims in Sections 1 through 7 concern finite cells and the exact coefficient
of the port coupling squared. They do not prove the thermodynamic infrared
limit. The finite-cell arguments have not had an independent referee.

## 1. The explicit bond

[Definition] The alphabet is

    L_a, R_a, D_(A,a), D_(B,a), X_(A,a), X_(B,a),   a=0,1,
    M_0, M_A, M_B.

There are 12 data letters and three marker letters. The integers carried by
the on-site charge are 2a on L,R, 2a-1 on D_A,X_A, 2a+1 on D_B,X_B,
and 0,-1,+1 on M_0,M_A,M_B. Every term below conserves their sum Q.
On the code Q=2 N_data, so exp(i pi Q/2) is logical fermion parity.
The normalization of Q is immaterial for U(1).

[Definition] The allowed ordered pairs are given by this complete table.

| Left letter | Allowed right letters |
| :-- | :-- |
| Any marker | Every letter |
| L_a | Every letter |
| R_a,D_(s,a),X_(s,a) | R_0,R_1,M_0 |

Put K on each disallowed pair. The grammar penalty is diagonal. It need
not commute with the remaining bond at finite K.

[Definition] Set z_A(a)=(-1)^a, z_B(a)=1. Include these transitions and
their Hermitian conjugates, summing all displayed bit and species indices.
The amplitude is the coefficient of the output ket times the input bra.

| Input | Output | Amplitude |
| :-- | :-- | :-- |
| R_a R_b | L_a R_b | -1 |
| R_a M_0 | L_a M_0 | -1 |
| D_(s,a) R_b | L_a D_(s,b) | -z_s(a) |
| X_(s,a) R_b | L_a X_(s,b) | -exp(i phi_s) z_s(a) |
| X_(s,a) M_0 | L_a M_s | -exp(i phi_s) z_s(a) |
| D_(s,a) M_0 | L_a M_s | -t_d exp(i psi_s) z_s(a) |
| m D_(s,a) | m X_(s,a) | -1, for each of the three markers m |
| M_0 R_a | M_A D_(B,a) | lambda exp(-i pi/4) |
| M_0 R_a | M_B D_(A,a) | lambda exp(+i pi/4) |

The last two rows define lambda V. They are the entire port. In particular
there is no additional product of one-site raising and lowering operators
that would add unspecified transitions.

[Definition] The on-site operator o gives energy Delta to every D,X or
excited marker letter, an additional D to X and excited markers, and
mu_sel to every marker. Add

    g ( |D_(B,0)><D_(A,1)| + h.c. ).

The logical field is lambda^2 epsilon times the data bit a on every data
letter. The repeated Hermitian bond is

    h = h_moves + K Pi_bad + lambda V + (o tensor I + I tensor o)/2.   (1)

Here h_moves excludes the last two rows, and o includes the conversion
and field. Formula (1) has one-site translation invariance and range one.
`bond15.py` constructs its 225 by 225 matrix. `bond18.py` supplies the
unmerged parent used to implement the rule identification. Duplicate
rules after identification are included once, not added twice.

## 2. Exact constrained code and the alphabet reduction

[Proved] In a cell of n data sites ending at a marker, the complete legal
space consists of

    L^j R^(n-j) M_0,                         0 <= j <= n,
    L^j D_s R^(n-j-1) M_0,
    L^j X_s R^(n-j-1) M_0,                  0 <= j < n, s=A,B,
    L^n M_s,                                s=A,B.

All n data bits are independent. To verify completeness, once R or an
excited data letter occurs, every later data letter must be R. An excited
marker requires an all-L prefix. Thus the dimension is (5n+3)2^n.

[Proved] The first line is a path of n+1 vacuum positions, even though
there is no vacuum-token letter. Its head is the L/R boundary. Projecting
the first two moves of the table onto legal words advances or retracts
that boundary by one site. Consequently, at K=infinity,

    E_v = -2 cos(pi/(n+2)),
    G|a> = sum_(j=0)^n v_j |L^j R^(n-j) M_0; a>,
    v_j = sqrt(2/(n+2)) sin(pi(j+1)/(n+2)),
    alpha = v_0 = v_n.

The marker energy adds mu_sel. The vacuum gap is
2[cos(pi/(n+2))-cos(2pi/(n+2))].

[Proved] This legal Hamiltonian is unitarily equivalent to an 18-letter
parent. Its extra letters are a vacuum token D_(v,a) and a distinct
vacuum marker M_v. Identify

    D_(v,a) -> R_a,       M_v -> M_0.

On legal words this is injective, including at j=n. At j=n-1 the image
ends in R M_0, whereas the vacuum-marker image is all L. The excited
words are unchanged. The projected move table, conversion, energies,
and port agree term by term with the parent. Before identifying D_v,
identifying only M_v gives the intermediate 17-letter design.

[Proved limitation] The identification is not a unitary on the complete
on-site Hilbert space. For example, an RR move in the middle of an R tail
creates an illegal RL pair. Likewise the new port can excite a marker
while the left cell's domain wall is elsewhere. Discarding those moves
at finite K would be incorrect. Sections 6 and 7 retain and bound them.

## 3. Two excited species suffice

[Definition] Choose

    phi_A = pi/(2n), phi_B=0,
    psi_A = pi, psi_B=pi/2,
    Delta = 4 + E_v.

Gauging the A rail by its prefix of Z operators turns the conversion into
a Jordan-Wigner annihilation operator A -> B, and a creation operator
B -> A. The right A source acquires -i P_A, where P_A is the parity of
the left cell. Gauging the express phase away makes both species have
the same scalar graph H. It is the folded path

    M, X_(n-1), ..., X_0, D_0, ..., D_(n-1),

with diagonal 4+D on M,X, diagonal 4 on D, nearest-neighbor entries -1,
and a loop entry -i t_d from D_(n-1) to M. The reverse entry is its
complex conjugate. All energies in H are measured above E_v.

[Proved] The port source, up to a common phase and alpha^2, is

    S = |B,A> - P_A |A,B>.                                    (2)

The two terms have net excited charge zero. Relative charges -1/2,+1/2
in units where a data particle has charge one are sufficient. The old
neutral,+,- triplet is not necessary.

Let R=(H tensor I+I tensor H)^(-1), x=R|M,D_0>, and let U_i be x
restricted to row D_i, V_j be x restricted to column D_j. Define

    C_ij = U_i^dagger R V_j,
    A_ij = U_i^dagger R U_j + V_i^dagger R V_j.

[Proved] To second order in the conversion, the pair source after one
conversion has two orthogonal species sectors. In global fermions f,
their operator parts are

    Y_AA = sum_i U_i f_Ai^dagger - sum_j V_j f_Bj^dagger,
    Y_BB = P_A [sum_i U_i f_Ai + sum_j V_j f_Bj].

The second-order port coefficient at order g^2 is
-alpha^4 g^2 (Y_AA^dagger R Y_AA+Y_BB^dagger R Y_BB).
In the BB term the two adjacent P_A factors cancel. Using the CAR in
the AA term gives a scalar, quadratic hopping, and onsite terms only.
The results are

    T/g^2  = -2 alpha^4 Re C,
    mu/g^2 = -2i alpha^4 Im A.                                (3)

The onsite matrix in (3) includes the contributions of both neighboring
bonds. There is no quartic term at order g^2. These matrices are one
half of C33's common-loop matrices, with the same sign of both T and mu.
This is an algebraic identity, not a numerical extrapolation to n=15.

[Numerical] The independent literal-bond calculation and the species-source
enumeration agree on the full 4^n-dimensional pair code at n=1,2 with
relative Frobenius error below 1.1e-15. Reconstructing the full operator
from its vacuum and one-particle matrix elements leaves a relative
non-Gaussian residual below 3.2e-16.

## 4. Parameters and finite root certificate

[Definition] For the 15-letter member take n=15, cell length b=16, D=1000,
and define the positive integers

    d_0=1, d_1=1004, d_(j+1)=1004 d_j-d_(j-1), beta=1/d_16.

Choose

    t_d = beta 2^-84,       g=2^-256,       K=2^1152,
    E_unit = alpha^4 beta^2 g^2,
    epsilon = E_unit 2^-140,                lambda=2^-640.     (4)

Use the phases and Delta of Section 3. The physical common field is
lambda^2 epsilon. At zero port and infinite K, the inherited marker
selector is

    a_b = 2[1-cos(pi/(b+1))],
    delta_b = (b+1)a_b-b a_(b+1),
    mu_sel = 2-(delta_16+delta_15)/2.                          (5)

The finite-K all-length selector is proved in Section 8.
For computations inside fixed cells mu_sel is a common scalar and can
be omitted.

[Proved, computer-assisted] `certify_tree.py`, `rational_bands.py`,
`full_loop_integer.py`, and `certify_roots.py` reuse C33's scalar routines.
They were copied into this folder, with attribution in their first lines.
The change of energy unit from 2 alpha^4 beta^2 g^2 to E_unit absorbs the
factor 1/2 in (3). No C33 output at n=23 was substituted for n=15.

The scalar product graph has dimension 31^2=961. The tree solves use
256-bit fixed-point integers, and the full finite-loop solves use 512.
The integer Jacobi recurrence encloses the source and resolved columns
by 10 last-place units. The loop denominator is d_16 2^84. The infinity
norm contraction is at most 1/2+t_d/4, so the recurrence length
bit_length(d_16)+512+20 controls its tail. Propagating the entry enclosures
and performing fraction-free congruences gives a normalized row error
less than 7.873e-86.

For clarity, a minor conservative correction was made to the inherited
product-radius formula. If both real and imaginary component errors are
at most 10 units, the quadratic product-error term is 200, not 100.
The full-loop implementation here uses 200. This change is included in
the final certificates.

[Proved, computer-assisted] A rational diagonal-dominance test gives
Herm(-T/E_unit)>0 with margin greater than 0.5510 after the stated
congruence. Numerical generalized eigenvalues propose root brackets,
but all acceptances use rational endpoint inertias. At each endpoint
t=tan(k/2) is rational and

    cos k=(1-t^2)/(1+t^2), sin k=2t/(1+t^2).

The endpoint inertia changes by one in each of 30 disjoint brackets.
The counts increase from 0 to 15 and then decrease to 0. Multiplying
det h(k) by exp(15 i k) produces a polynomial of degree at most 30.
The 30 sign-changing roots therefore exhaust the degree and are simple.

[Proved, computer-assisted] All 465 unordered Fermi-point pairs with
replacement are included in the pair-sum test. Their rational unit-circle
phase products are sorted by exact orientation tests. Subtracting root
bracket errors from neighboring chord distances gives

    delta_2 > 0.0145540,
    delta_4 > 1.1444e-7.                                      (6)

The second test includes pairs from the same band and repeated points.
It consequently excludes coincident sea centers modulo pi. The exact
rational brackets and lower bounds are in root_certificate_exact_n15.json.
The numerical root centers are in full_loop_bands_n15.json.

[Scope] Equation (6) is finite low-order nonresonance. It does not prove
rational independence against every integer linear combination of the
momenta. Higher-degree operators are handled by the separate RG argument,
not by interpreting a finite decimal expansion as an irrationality proof.

## 5. Finite conversion strength

[Proved] The free pair inverse obeys A_0>=3 I for (4). The conversion
operator on one legal cell has norm g, so the pair conversion has norm
at most 2g. The port source norm is at most sqrt(2) alpha^2. Odd conversion
orders cannot return a source to its species sector. For g<=1/2 the
Neumann tail starting at g^4 therefore obeys

    ||H_2(g)-H_2(0)-g^2 H_22||
      <= (2 alpha^4/3) (2g/3)^4/[1-(2g/3)^2]
      <= alpha^4 g^4/4.                                    (7)

This bounds the complete operator remainder, not just its quadratic
projection. A matrix element of the Hilbert-Schmidt quadratic projection
is bounded by 4 times the operator norm. Summing two onsite blocks and
two hopping blocks gives a Bloch norm bound of 16n times that norm.
Thus (7), in E_unit units and after the certificate congruence W, gives

    error_g <= 4n g^2 beta^-2 tr(W W^dagger).                 (8)

The prefactor is twice that in the inherited three-species energy unit.
The final g=2^-256 is smaller than the inherited threshold 2^-247 and
accounts for this change. The exact rational value of (8) is less than
1.137e-15. The root certificate allows a total additional norm 2^-32.

## 6. A finite-K comparison for a fixed cell

[Proved] Write a single cell, with frozen marker positions, as

    H_K = K N_bad + B,       ||B|| <= B_0=2^16.               (9)

This deliberately loose bound holds at n=15, D=1000 and (4), including
the marker energy. One can obtain it by summing the norms of the listed
local moves and onsite terms over 16 sites. There are at most 12 unit
move blocks per bond and onsite norm below 1005. The bound also holds
on sectors with more than one excited token. On P=ker N_bad, PBP is
the 18-letter parent. On its complement N_bad>=1.

For |z|<=2B_0 and K>=8B_0, the Schur complement is

    PBP-z - PBQ [QH_KQ-z]^-1 QBP,

and the norm of the last term is at most 2B_0^2/K. This follows from
QH_KQ-z>=K-3B_0>=K/2. It is an ordinary finite-matrix estimate.

[Proved] In the zero-excited-token sector the data bits are spectators
at every K. Its scalar Hamiltonian acts on L/R orientations only.
The ideal ground is simple and has gap greater than delta_0=1/16.
The Schur estimate and the spectral gap imply

    |E_K-E_v-mu_sel| <= 2B_0^2/K,
    ||G_K-G_infinity|| <= 16 B_0^2/(K delta_0).               (10)

Here the isometries use the same data basis, and G_infinity is embedded
in the full orientation space. To see the second estimate, the high
component is at most 2B_0/K. The component in the orthogonal legal
eigenspaces is at most (4B_0^2/K)/delta_0. Normalization and the choice
of the positive scalar ground-vector phase fit inside the factor 16.
The code is exactly 2^15 dimensional, with no data-dependent splitting.

[Proved] The source changes each cell's excited-token number from zero
to one. H_0 preserves that number separately in each cell, even outside
the grammar. Hence the complete second-order resolvent is the (1,1)
excitation block, with no truncation in the number of illegal pairs.
Its ideal legal part has gap at least 3 above the pair vacuum energy.
Repeating the Schur estimate with pair norm 2B_0, and using (10) to
shift the energy, bounds its inverse difference from the embedded ideal
inverse by 2^40/K. The high-high block contributes O(1/K) and the
off-diagonal blocks O(B_0/K), both covered by this bound.

More explicitly, the pair Schur correction is at most 8B_0^2/K, and
the vacuum-energy correction is at most 4B_0^2/K. Their sum is below
2^36/K. The resolvent identity between legal inverses with gaps at least
2 and 3 costs at most 1/6 of this. The high block inverse is at most
2/K; each mixed block is bounded by 2B_0/K times a legal inverse.
These contributions are smaller than 2^40/K.

The two-cell code-isometry difference is at most 2^41/K by (10).
Since ||V||=sqrt(2), replacing both code vectors and the resolvent in
-G_K^dagger V(H_0-2E_K)^-1 V G_K gives, conservatively,

    ||H_2(K,g)-H_2(infinity,g)|| <= 2^64/K.                   (11)

In applying the resolvent identity, keep the full source V G_infinity,
including its illegal component, and extend the ideal inverse by zero
outside P. This is why (11) includes the unwanted marker-port process.
It is not assumed away.

[Proved, computer-assisted] A lower bound sufficient to convert (11)
to certificate units is

    alpha^4 >=64/17^6 >2^-19,
    E_unit >= 2^-19 d_16^-2 2^-512.

The first inequality follows from sin(pi/17)>=2/17. Therefore

    error_K <=16n tr(W W^dagger) (2^64/K)/E_unit.              (12)

`finite_budget.py` compares (8)+(12) with the allowed budget using
Fractions, not decimal eigenvalues. The log10 upper bound from (12)
is below -27.8885. The sum uses less than one eighth of 2^-32. Every root
bracket and both pair margins in (6) therefore hold for the exact
finite-K, finite-g quadratic projection of H_2.

[Scope] Bounds (9)-(12) are for n=15 cells. They have not been promoted
to a uniform estimate over all cell lengths or to a many-cell infrared
matching theorem.

## 7. Exact second order and finite-size diagonalization

[Proved] Every port changes the excited-token number of both incident
cells by one. Its code compression is zero. Two different adjacent
ports cannot restore the excitation number in their outer cells, so
there are no three-cell terms at order lambda^2. On two cells,
(-1)^(N_exc,left) commutes with H_0, negates V, and is identity on the
code. The analytic effective Hamiltonian is even in lambda.

[Numerical] `checks.py` constructs the full invariant closure of the
literal move table, including illegal states. At n=1,2 these spaces have
dimensions 436 and 15404. Their legal subspaces have dimensions 256
and 2704. The exact H_2 uses sparse linear solves in the entire source
sector. Its energies are compared with shift-invert diagonalization of
H_0+lambda V. All four or sixteen code levels are retained.

At D=4, Delta=0.5, g=0.04, t_d=1e-4, K=100, the energy error divided
by lambda^4 approaches 0.00581 at n=1 and 0.00173 at n=2 as lambda
decreases through 0.2,0.1,0.05. The raw values, residuals and dimensions
are in small_checks.json. Separate K=100,1000 comparisons show the
expected O(1/K) correction to H_2. These are checks at resolved diagnostic
parameters, not direct ED at the exponentially small parameters (4).

[Rejected numerical check] Unshifted Lanczos with a small Krylov space
missed multiplicities in the almost degenerate n=2 code band. Its
purported low 16 levels included an excited level about 0.776 above
the code. That failed comparison is preserved in
small_checks_unshifted_rejected.json. Shift-invert with 80 Krylov
vectors resolves all 16 and passes the fourth-order scaling test.

## 8. The zero-port selector at finite K, uniformly in cell length

[Proved] There is a stronger estimate for selection than the finite-cell
Schur bound. Write k for an arbitrary number of data sites. Every excited
token hop, switch and conversion preserves each individual bad-pair
projector. A vacuum move changes a bad-pair projector only when its
predecessor is R,D or X. Split off these vacuum moves as -(a_j+a_j^dagger),
with a_j creating a bad pair immediately to the left of site j. Then

    a_j a_j^dagger <= Pi_bad,(j-1),
    sum_j a_j^dagger a_j <= max(k-1,0) I.

There is no such term at the first data site because its predecessor is
a marker. The last-site rule uses M_0 in place of the right R projector.
The rest, H_good, commutes with every bad-pair projector.

[Proved] In a block with m bad pairs, cutting at these pairs gives at
most m+1 legal fragments. Each fragment has at most one active head,
either a vacuum domain wall or an excited token. Each head has total
kinetic row sum at most two. This includes both express and data paths,
the switch and the weak link, whose weight is below one. The onsite
conversion costs at most g per excited token and is absorbed by Delta.
Thus in blocks with m>=1,

    H_good >= mu_sel-2(m+1) >= mu_sel+E_v(k)-4m.

For m=0 the vacuum block has lower energy mu_sel+E_v(k), and the excited
block has lower energy at least mu_sel+Delta-2-g, which is larger because
Delta>2+g at (4). The k=0 cell is immediate. Therefore, on the entire
fixed-cell Hilbert space,

    H_good >= mu_sel+E_v(k)-4 N_bad.

For eta>0, completing the square gives

    -(a_j+a_j^dagger) >= -eta a_j a_j^dagger-eta^-1 a_j^dagger a_j.

Choose eta=K-4 and restore K N_bad. A constrained vacuum trial state
gives the upper bound, so the complete cell ground energy obeys

    mu_sel+E_v(k)-max(k-1,0)/(K-4)
       <= E_K(k) <= mu_sel+E_v(k).                           (13)

This estimate includes illegal strings and arbitrary excited-token number.
It is uniform in k. It avoids applying the k-dependent Schur norm bound
outside its range of use.

[Proved] Let e_b(infinity)=(mu_sel-2 cos(pi/(b+1)))/b. The inherited
C22 selector lemma says this sequence decreases to b=16 and increases
thereafter for (5). Thus the smallest density difference from another
integer length is the smaller of the two neighboring differences.
`selector_checks.py` uses rational alternating-series enclosures of pi
and cos to prove

    min(e_15-e_16,e_17-e_16) > 2^-16,
    e_16 < -2^-8,
    omega_16 > 1/16.

Equation (13), divided by b=k+1, changes any cell energy density by at
most 1/(K-4). Since 1/(2^1152-4)<2^-17, the SAME marker chemical potential
still selects b=16 uniquely at finite K. All excited and illegal sectors
are covered by the lower bound; the finite-cell gap in Section 6 resolves
the selected cell itself.

[Proved] A marker-free ring with m bad pairs has at most m fragments.
The same argument gives energy density at least -1/(K-4). Its zero-bad
sectors are the all-L and all-R strings. This lower density exceeds the
selected-cell trial density. Marker positions remain exactly conserved,
so the zero-port ground on compatible rings consists of length-16 cells
and their translations. The code is exactly 2^15 per cell.

[Numerical] Full invariant two-cell spaces at n=1,2 check both the
H_good form inequality and (13). Their dimensions are 436 and 15404.
The exact selector inequalities and these checks are saved in
selector_certificate.json.

## 9. What the finite calculation does not establish

[Assumed for the central-charge interpretation] The n=15 quadratic
reference has 15 massless Dirac pairs. Exact U(1), (6), and sufficiently
small residual interactions are the usual inputs to multichannel
Luttinger-liquid stability. The conversion remainder (7) and the
constraint remainder (11) are far below the smallest reference band
scale at (4). A volume-uniform integration of the massive sectors and
the infrared RG are separate inputs, as in the parent candidate.

[Open] Section 8 concerns lambda=0. At nonzero lambda, comparing a defect
configuration with the interacting selected bulk requires subtraction
of the common bulk energy. A bound on the extensive total port norm is
not sufficient. The uniform dressed-defect comparison and the final
infrared transfer remain separate from the finite band certificate.

[Conditional] If those selection and infrared steps hold, the frozen
15-letter chain has c=15>14. Adding one ordinary vacancy produces a
16-letter candidate with 15 payload modes and one charge mode, c=16>15,
subject also to the previously identified vacancy matching and dressed
defect comparison. Neither c was measured. This report establishes the
finite band ingredients, not a completed counterexample.
