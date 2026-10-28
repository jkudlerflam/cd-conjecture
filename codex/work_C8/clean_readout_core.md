# An autonomous initialized high-channel core

Status: proved finite construction, with a literal collision table and
finite checks. This closes the clean-readout program left open in run 7.
It does not by itself prove a microscopic CFT at nonzero coupling, and it
does not realize C19's second-order spectral filters.

## 1. The theorem

Let q be divisible by four, s=q/2, m=q^2/8, and n=m+6. There is an
explicit positive, uniform nearest-neighbor Hamiltonian on an open
length-n segment, with one-site positive boundary terms, whose site
dimension is

    d_core = 446400(q+1).

Its complete zero space has dimension 2^m. Its product grading acts on
this space as logical fermion parity. There are 2m specified odd
single-site operators at each endpoint and m specified even two-site
operators at each end with exact compressions

    V† A_(2f-1) V = alpha_Lf Gamma_(2f-1),
    V† A_(2f)   V = alpha_Lf Gamma_(2f),
    V† B_(2f-1) V = alpha_Rf i Gamma_(2f-1) P,
    V† B_(2f)   V = alpha_Rf i Gamma_(2f) P,
    V† D_uf V = beta_uf Z_f.

Every alpha and beta is a positive, explicitly computable rational
number. The gate table, initialization projectors and readout matrices
are supplied below and in clean_readout_machine.py. They do not require
a preloaded program, a universal-machine theorem, position-dependent
instructions or a control state for each flavor.

For positive propagation scale J, initialization strength eta and clock
and role penalties kappa_clock,kappa_role, a valid gap lower bound is

    min{kappa_clock,kappa_role,
        J[1-cos(pi/C_n)], eta/(2C_n), J/(4C_n^2)},
    C_n=(2n+4)2^n+n-4.

The bound is for this finite core. It is not a gap above an interacting
many-cell low-energy band.

## 2. Registers, grammar and initialized input

Each ordinary work site contains a static role, one cursor bit c, one
logical-data bit d, one READY bit, and one record v in

    {e_0,...,e_(s-1),o_0,...,o_(s-1),*}.

The ordinary factor is 8*2*2*2*(q+1)=64(q+1). A site can instead carry a
head overlay consisting of one of 58 direction-fixed controls, a fresh
or updated phase, a neutral parity bit p, and a graded carrier bit x.
The overlay holds the same site record, not an additional q-state record.
Thus

    w = 64[1+2*58*2*2](q+1) = 29760(q+1).

The data bit d, the head's x bit and each odd record o_j are odd. All
other labels are even. The site's grading is the product of these signs.
For an ordinary site omit the head's x factor.

Allowed neighboring roles are

    L A, A B, B S, S D, S E, D D, D E, E F, F R.

The left endpoint is pinned to L and the right endpoint to R. For each
n>=7 these constraints have the unique solution

    L A B S D^(n-7) E F R.

Indeed L forces A and then B and S. A string of D's can only terminate
through E,F,R, and cannot enter a second S. Every other role sector has
a positive conserved diagonal penalty. At n=m+6 the S,D sites number m.
They are the logical input bits. Roles A,B hold the left row and column
counters; F,E hold the right row and column counters.

During the clock's I initialization sweep impose the following one-site
pins, conditioned on the local immutable role.

| Role | Cursor | Data | READY | Record | Head |
| :-- | :-- | :-- | :-- | :-- | :-- |
| L | 0 | 0 | 0 | * | fresh L0, p=x=0 |
| A,B,E,F | 0 | 0 | 0 | e_0 | absent |
| S | 1 | arbitrary | 0 | * | absent |
| D | 0 | arbitrary | 0 | * | absent |
| R | 0 | 0 | 0 | * | absent |

For each role the pin is I minus the projector onto the displayed local
subspace. Every entry outside it is penalized, including additional heads,
unused record values, unexpected carrier bits, cursor bits, READY bits,
and surplus data bits on the six endpoint/counter sites. No supposedly
blank untrusted input is tested using uninitialized scratch. Only the m
S,D data bits remain free. The initialization sweep acts as identity on
work, including its activation at the far end; the initial head is already
pinned at L. This is allowed by the inherited clock theorem.

This regular layout avoids the program-generation problem. It differs
from a compiler that accepts an arbitrary instruction tape and later has
to certify its contents.

## 3. Local quantum permutations

The complete record alphabet is represented in code by integers
0,...,s-1,s,...,2s-1,2s, with *=2s. All modular arithmetic below is modulo s.
The finite menu is defined on the entire local space, not only on valid
inputs.

The onsite operations are as follows.

* PARITY sends p to p xor d.
* TAKE exchanges d and x.
* LOAD exchanges (x,v)=(0,*) with (0,e_0), and exchanges (1,*) with
  (0,o_0), fixing all other pairs.
* ODD_P sends o_j to o_(j+p), fixing even records and *.
* EVEN_P sends e_j to e_(j+p), fixing odd records and *.
* FLAG toggles READY. Identity is the seventh operation.

Their inverses are literal inverses; LOAD, PARITY, TAKE and FLAG are
involutions. In particular LOAD removes x from the head. Copying x into
a port while retaining x elsewhere would not work.

The two-site record gates use departure record a and arrival record b.
They leave all static roles, cursor and data bits unchanged.

* SWAP exchanges the records.
* EVEN_ADD sends (e_i,e_j) to (e_(i+j),e_j).
* ODD_ADD sends (o_i,e_j) to (o_(i+j),e_j).
* SUM sends (e_i,e_j) to (e_i,e_(j+i)).
* INC sends (e_a,e_b) to (e_(a+carry),e_((b+2) mod s)), where carry=1
  if b+2>=s and zero otherwise.

Each gate is identity on every type pair not displayed. Every displayed
type sector is invariant and the map on it is bijective, so these are
whole-space permutations. INC has inverse b_old=(b-2) mod s,
a_old=a-1_(b_old+2>=s). Starting at (0,0) with even s, its orbit has
length s^2/2 and lists

    (a,b)=(0,0),(0,2),...,(0,s-2),(1,0),...,(s-1,s-2).

The loop counters are therefore on two sites of the existing record
track. The controller does not store a or b. Arbitrary q-dependent
matrix entries in a fixed gate family are allowed by the microscopic
problem; their count does not introduce a q-dependent control alphabet.
All listed onsite and bond gates preserve total product grading.

There is a stronger exact symmetry. Define the integer onsite charge
as data d plus the indicator of an odd record, plus x if a head is present.
Every listed operation conserves the sum of this charge, not merely its
parity. TAKE and LOAD transfer a unit between registers; the other gates
preserve each charged occupation. Neutral-block completion uses identity
quantum gates and therefore preserves the charge as well. Initialization
is diagonal in it. On the core code the charge is sum_f x_f. Thus the core
has an exact on-site U(1) whose restriction is total logical fermion number.
It does not give U(1)^m flavor symmetry. This stronger conservation was
checked on all primitive columns at q=4,8,12 and on 12000 full-space G,F
samples, in clean_readout_number_checks.json.

## 4. One readout and its inverse

Use zero-based flavor f. The cursor is on physical site j=f+3. During
the left readout the counters are a=f div(s/2), b=2(f mod(s/2)). A head
starting at L with p=x=0 follows these eleven preparatory segments.
An endpoint operation occurs after the move, on arrival. A scan operation
occurs at each intermediate S,D site, but not at the endpoint cursor.

| Control | Segment | Record gate on each move | Scan | Endpoint operation |
| :-- | :-- | :-- | :-- | :-- |
| L0 | L to cursor | I | PARITY | TAKE |
| L1 | cursor to L | I | I | LOAD |
| L2 | L to A | EVEN_ADD | I | I |
| L3 | A to B | SWAP | I | I |
| L4 | B to L | I | I | I |
| L5 | L to A | ODD_ADD | I | I |
| L6 | A to B | SWAP | I | I |
| L7 | B to L | I | I | ODD_P |
| L8 | L to B | I | I | I |
| L9 | B to A | SUM | I | EVEN_P |
| L10 | A to L | I | I | Enable READY |

The only nontrivial moving gates in this table have adjacent start and
end sites; the longer segments have identity moving gates. On L0 the
head accumulates the prefix parity, then swaps x_f into x, leaving the
cursor's data bit zero. L1 moves it to the physical endpoint and LOAD
sets x=0. The endpoint record is now e_0 or o_0.

L2 adds the row index to an even endpoint. L3 temporarily swaps the
counter records so L5 can add the column index to an odd endpoint.
L6 restores their positions. L7 adds p to the odd endpoint. L9 changes
A from e_a to e_(a+b+p). B remains e_b. Thus at activation,

    p = sum_(i<f) x_i mod 2,
    v_L = e_a if x_f=0, o_(b+p) if x_f=1,
    v_A = e_(a+b+p),    v_B=e_b,    head x=0,
    data_f=0,          data_i=x_i for i!=f.

Every record outside v_L is independent of x_f. The stored p is permitted
because it depends only on the other bits.

With READY enabled, the head walks L to R and back using identity gates
on all data and record registers. On its return it clears READY and
executes the literal inverse of all eleven preparatory segments. It
restores the logical input, both counters, all port and color records,
and p=x=0. No active readout occurs during uncomputation.

The right readout is the reflection of this construction. Replace
(L,A,B,R) by (R,F,E,L), reverse every direction, and compute the suffix
parity instead. Its active endpoint is e_a or o_(b+p) with
p=sum_(i>f)x_i mod 2. Its adjacent color is again e_(a+b+p).

For clarity about the inverse implementation, each reversed control
moves opposite its corresponding forward control and applies the adjoint
moving gate in the reversed departure/arrival convention. A reversed
collision is the inverse of the forward collision at the arrived-at site.
The endpoint READY clearing switches to reverse L10 or R10 before any
record changes. This is the move-then-collide normal form, not an
informal instruction to reverse a computation while keeping its gates.

## 5. The complete autonomous cycle and collision completion

After the left readout and its inverse, the head travels neutrally to R
and performs the right readout and its inverse. It then executes INC on
F,E, travels to L and executes INC on A,B. It finds the cursor, clears
its mark and moves one site right. If that site is D, it places the mark
there and returns to L. If it is E, it returns to the first data role S,
marks S, and returns to L. The directions of these scans depend only on
roles and cursor, never on quantum input bits or counter values.

There are 44 preparatory/reverse controls, four readout-walk controls,
and ten travel/counter/cursor controls, giving Q=58. The local neutral
collision labels are (control,role,cursor), 58*8*2=928 possibilities.
clean_readout_machine.py supplies 799 prescribed transitions, each with
one of the seven onsite permutations above. All other entries are
completed explicitly within each immutable-role sector, first fixing
common unused domain/range entries and then pairing the remaining lists
in increasing order. clean_readout_table.json contains all 928 blocks.

The prescribed neutral map is an injection. A phase's scan loop excludes
its starting mark and its terminal mark. Entry from the preceding phase
has the starting mark, and exit to the next has the terminal mark. A
cursor-clearing exit has c=0 in its new NEXT phase. Its subsequent marked
entry to RETURN has c=1, excluded from RETURN's scan loop. The last-bit
branch enters RESET at E, excluded from RESET's loop; RESET reaches S,
places c=1, and enters RETURN on S rather than D. These ranges do not
collide. Forward and reversed controls have disjoint labels. The generator
checks every domain and range, including malformed neutral labels.

Since every decoration is itself a whole-space graded permutation, the
completed block map C is a whole-space graded permutation as well. This
supplies a finite collision table with no undefined quantum operation.

The work head makes one move and one collision per full clock sweep.
Use the fresh/updated G,F construction proved for endpoint_shuttle.md.
On a bond with one fresh head in its outgoing orientation, G applies the
specified moving gate and transfers the full head, including p,x, to the
neighbor in updated phase. On the reverse subspace use the adjoint. All
other two-site sectors, including two-head sectors, are identity. F maps
updated through C to fresh, and fresh through C^{-1} to updated. These
are Hermitian unitaries on their full spaces. The fixed direction and
phase prevent a head from moving twice in an ordered outward sweep.

Each flavor iteration takes 10n+46+2f moves, including both clean
readouts, both counter updates and cursor advance. To verify this,
the left preparation has length 2j+12, its inverse the same, and its
readout walk length 2(n-1), totaling 4j+2n+22. The reflected right total
is 6n-4j+18. Their sum is 8n+40. Travel to the right costs n-1; the
remaining update/advance/return route costs n+2j+1. Adding gives the
stated expression. On the last flavor the return via S has the same
length as the ordinary return route.

Consequently the complete cycle has

    T_m = sum_(f=0)^(m-1) (10n+46+2f) = 11m^2+105m

moves. After m=s^2/2 iterations, INC has returned both counter pairs to
zero and the cursor returns to S. The whole quantum and neutral state
returns to its initialized value. At m=2, T_m=254<256=2^(m+6).
For m>=2, T_(m+1)<2T_m, since 2T_m-T_(m+1)=11m^2+83m-116>0.
Induction gives T_m<2^(m+6). Every history therefore contains a full
cycle, including both readouts of every flavor. Truncating a later
cycle cannot introduce an unclean active time.

## 6. The microscopic positive core

Tensor this work site with the inherited 15-letter sweeping clock.
For a clock edge u to v decorated with work unitary U use the positive
term

    (J/2)[ |u><u| tensor I + |v><v| tensor I
           - |v><u| tensor U - |u><v| tensor U† ].

Decorate outward bond rules with G, departing-site return rules with F,
and the final left return with F. All other clock rules carry identity,
including initialization, activation, binary carry and right bounce.
clean_core_bond.py evaluates every column of the resulting complete
bulk and boundary matrices without restricting the work sector.
The complete neutral clock rules and legality projectors are given in
sweeping_clock.md and implemented in sweeping_clock_checks.py. No clock
rules are added or removed.

Add the forbidden-role-pair projectors, the role endpoint pins, the
clock legality projectors and the I-sweep initialization pins of Section
2. All terms are positive. Bulk one-site terms can be divided equally
between neighboring bonds, with the missing half restored in the one-site
boundary term. Thus this is one repeated Hermitian NN bulk matrix and
one fixed one-site matrix at either end, in the sense required by the
existing cell wrapper.

There is no dimension-dependent norm hidden in the controller. The 13
bond rewrite terms each have norm J. The shared onsite turn and input
pin contribute at most J+eta to one bulk bond. Hence

    ||h_pay|| <= 14J+eta+kappa_clock+kappa_role.

The left and right boundary norms are bounded respectively by
(9J+eta)/2+kappa_clock+kappa_role and
(7J+eta)/2+kappa_clock+kappa_role. At unit scales these bounds are
17,7,6. They are uniform in q and n, although the finite-core gap and
normalized readout strengths are not.

Role defects are conserved since neither G nor F changes a role.
Clock defects are conserved by the inherited clock rules. A zero vector
therefore has legal roles and lies on the unique legal clock path.
On that path gauge away the work permutations. Propagation becomes the
scalar path Laplacian tensored with the entire work space. Its zero
vectors are constant in clock time. The I-sweep gates are identity, so
its initialization pins restrict their common input to exactly the
product subspace in Section 2, of dimension 2^m. This proves completeness
of the ground code, including all malformed head-number and work sectors.
No work-machine-sector assumption remains.

Every rejected computational-basis input violates at least one initial
projector of weight eta. After gauging, each such input is a path with a
positive pin. The run-7 path-pin bound gives the stated gap. Positivity
and the conserved defects add the minima kappa_clock and kappa_role.
On the accepted inputs the product grading is (-1)^(sum x_f), and all
gates preserve it. Hence the full history isometry intertwines grading
with logical parity.

## 7. Physical readouts and all-time correctness

For each f define the two signed edges by

    (a_(f,p),b_(f,p)) = (floor(2f/s), (2f mod s)+p), p=0,1,

joining e_a to o_b. Define record matrices

    X_f=sum_(p=0,1) (-1)^p (|e_a><o_b|+|o_b><e_a|),
    Y_f=sum_(p=0,1) (-1)^p i(|o_b><e_a|-|e_a><o_b|).

They annihilate *. At the left endpoint use X_f,Y_f; at the right use
Y_f,-X_f. Multiply each by the local READY projector and the no-head
projector, and, if desired, by its endpoint-role projector. They are
one-site Hermitian odd operators, identity on the clock factor.

During every active clock vertex, the head carrier x is zero, the
selected input site is zero and all other registers depend only on
z_f, the bit string with x_f removed. The head's p is prefix or suffix
parity of z_f. This statement follows from Section 4 during the readout
walk. At partial sweeps, work changes only at the unique G head move
and the unique F collision. All these moves during READY use identity
record gates. At activation or clearing the head is at the endpoint,
where the readout projector is zero. Thus no partial update is measured.

The rest configurations for x_f=0 and x_f=1 are equal at each active
time. Conversely the data and neutral registers distinguish all z_f's.
The edge matrix elements are therefore the Jordan-Wigner prefix/suffix
signs. Distinct edge pairs give zero for a different flavor operator.
Averaging over the orthogonal clock states proves the four port identities
in Section 1 with alpha_uf=N_port,uf/C_n.

For the density readout use the proper edge color c=(a+b) mod s. On the
endpoint and its neighbor define D_f to have value +1 on (e_a,e_c)
and -1 on (o_b,e_c) for the two edges of f, and zero elsewhere. A fixed
vertex and color determine a unique edge, so these assignments never
conflict and different f have disjoint supports. Multiply by READY and
no-head projectors on both sites. Its compression is beta_uf Z_f,
where beta_uf=N_density,uf/C_n. These even NN fields permit independently
chosen chemical potentials. In general beta is not alpha.

The schedule is independent of all logical values. All visibility counts
are therefore input independent. Each readout lasts two full traversals
of a segment with n>=8, so every port and density count is positive.
Both left and right readouts are proved, including the right parity
dressing, and both are implemented by the same repeated bond.

## 8. Counts without exponentially long simulation

For a full work sweep k, let a_k count its active vertices including the
begin, all n-1 outward moves, the bounce, and all n return/finish moves.
There are 2n+1 such vertices. Let b_k be the active indicator after the
completed sweep. Both are periodic in k with period T_m.

There are 2n inactive initialization/return vertices before the first
begin. Between sweep k and k+1, 1<=k<2^n, the binary clock adds
2 nu_2(k)+1 vertices, all with work configuration k. The final overflow
adds n-1 vertices. Therefore, componentwise for every port and density,

    N = sum_(k=1)^(2^n) a_k
        +sum_(k=1)^(2^n-1) [2 nu_2(k)+1] b_k
        +(n-1)b_(2^n).

This reproduces C_n when all indicators are one and the initialization
is included. The only exponentially large integers have O(n) bits.
Use nu_2(k)=sum_(t=1)^(n-1) 1_(2^t divides k). At fixed t, b_(2^t j)
has period T_m/gcd(T_m,2^t), so complete periods and one residual prefix
suffice. clean_readout_counts.py implements this finite formula with
integer arithmetic and simulates only one polynomial-length machine
cycle. It is an explicit definition of every rational normalization,
without a circuit list of length 2^n.

## 9. Alphabet accounting and scope

The direct-sum realization gives the advertised K_direct=446400 in
q_core=K_direct(q+1). If the literal tensor shape of the old contract is
desired, take s0=32 for role,cursor,data, the separate READY bit, and two
lanes each of dimension r=465=1+8*58. Embed the head overlay in the first
lane with the other empty. Extend all operations by identity outside this
local subspace. Initialization excludes those additional states. Then

    K_padded=30 s0 r^2=207576000.

This is a concrete but unnecessary padding. Its outward gate is the
ordered G above, not the simple lane SWAP of the old packet toy. The
clock works for an ordered product of arbitrary G's; commutation of
neighboring work gates is not required. The smaller direct-sum count
is the implemented one and should be used for the microscopic dimension.

This construction supplies first-order Clifford channels. It does not
supply C19's port-accessible excited energies or spectral weights. Its
port-touching gates also do not satisfy the hypotheses of the C18/C19
static-Bell, port-idle obstruction, so that obstruction is respected.
The stronger fixed-d unbounded-c program remains separate.

## 10. Reproduction and actual checks

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 clean_readout_machine.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 clean_readout_checks.py --clock-q8
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 clean_readout_physical_checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 clean_readout_counts.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 clean_core_bond.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 clean_core_kernel_checks.py

The tests include all local permutation/inverse/grade columns at q=4,8,12;
all 4 and 256 initialized inputs for complete cycles at q=4 and q=8;
24 q=12 inputs; 51200 paired clean environments at q=8; all 5124 clock
vertices for all q=4 inputs; all 524298 vertices for one q=8 input;
and comparison with the machine after every full sweep. Separate tests
apply the actual physical endpoint matrices to all four q=4 histories,
recovering the four Clifford operators and both sets of densities with
zero integer-matrix discrepancy. They do not insert the target matrix
as the measured operator.

The full immutable-role grammar is counted through length 30. The entire
onsite I-pin alphabet at q=4 rejects 148790 columns and retains exactly
one or two columns according to the role, including every head phase
and carrier state. Another 12000 random full-work-space G,F involution
and grading tests include inputs outside the promised machine sector.

At q=4, the port counts are (254,254) on the left and (260,260) on the
right; the density counts are (218,218) and (220,220). At q=8 the counts
are flavor dependent; their full lists are in clean_readout_checks.json.
The polynomial counting formula agrees with both complete-clock checks
and also evaluates the q=12 history, whose clock has over 800 million
vertices, without traversing that history.

The sparse core-bond oracle passed Hermiticity and exact number-conservation
checks on 1600 arbitrary bulk columns at each q=4,8, together with both
boundary matrices. Applying that actual full-space Hamiltonian to every
one of the four q=4 history columns (5124 configurations per column)
gives zero residual, using 21002 local matrix entries per column. This
checks the core Hamiltonian independently of reading a program trajectory.
It is a zero-state test, not a numerical full-spectrum completeness test;
completeness and the gap follow from the positive initialization proof.

These are finite algebraic and controller tests. No microscopic central
charge was measured and no thermodynamic extrapolation was performed.
