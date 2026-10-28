# An initialized binary clock with a constant alphabet

Status: proved finite-chain construction and analytic spectrum, with
small-system numerical checks. It supplies a clock and workspace
initialization, not the high-channel program or a c-d counterexample.
The clock has no independent tape-length or program sector. Every
penalty-free configuration belongs to one explicitly specified path.

## 1. Alphabet and complete configuration graph

Use nine clock letters

    L, 0, 1, C0, C1, R0, R1, I, Z.

The head letters are C0,C1,R0,R1,I. The symbols 0 and 1 to the right of
a counter head are ordinary binary digits, with the least significant
digit at the left. L means a zero digit to the left of the head. I is
the initialization head, and Z is an unvisited initialization site.
These are orthogonal physical states, not external time labels.

The allowed internal ordered pairs are

    LL, Lh (all five heads h),
    hb (h=C0,C1,R0,R1 and b=0,1),
    ab (a,b=0,1), IZ, ZZ.                                (1)

At the left boundary allow L or a head. At the right allow 0,1,Z or
a head. All other pairs and boundary letters have penalty K>0.
Let Q_j,Q_left,Q_right denote the corresponding diagonal projectors.

Every zero-penalty word of length n>=1 has one of the forms

    L^(i-1) h s,  h=C0,C1,R0,R1,  s in {0,1}^{n-i},
    L^(i-1) I Z^(n-i),                                  (2)

with 1<=i<=n. To prove this, the left endpoint permits only L or a
head. L can be followed only by L or a head. A string of only L's
fails the right condition. After a head only its specified suffix
alphabet is allowed, and neither suffix alphabet permits another
head. This proves both necessity and sufficiency, including n=1.

Use the following directed transitions and include their adjoints
in the Hamiltonian:

    C1 b -> L Cb                  b=0,1,
    C0   -> R1,
    L Rb -> R0 b                  b=0,1,
    left-boundary Rb -> Cb        b=0,1,
    I Z  -> L I,
    right-boundary I -> R0.                              (3)

The second rule is one-site and is included once per site. The last
two boundary rules are boundary operators, or two-site rules involving
a marker when the clock is used on a ring.

Theorem 1. The complete zero-penalty graph is one directed path with

    C_n = n + 4(2^n-1)                                  (4)

vertices. Its first word is I Z^(n-1), its last is L^(n-1) C1.
Every Q_j, including the two boundary projectors, commutes with every
propagation operator. This statement holds on the full clock space.

Proof. The initialization head moves from left to right, making n
initialization configurations. Activation at the right gives
L^(n-1) R0. The return head moves left through the zero digits and
restarts as C0 at site 1. Thus the counter starts at binary zero.

From a C head at site 1 with number x, carry flips each trailing one
to L and moves right. At the first zero, C0 -> R1 changes that zero
to one. The return head moves left through the zeros just created,
then restarts at site 1. The number is now x+1. The all-ones number
carries to C1 at site n and has no outgoing transition.

Every C configuration in (2) occurs once: its i-1 erased digits were
ones immediately before that carry, and its unchanged suffix fixes
the number being incremented. Every R1 configuration is the unique
turnaround at its site. For an R0 configuration, the first one to its
right fixes the turnaround from which it is returning; if there is
no one to its right, it belongs to the initial all-zero return sweep.
Hence every R configuration occurs once as well. This also proves
there are no extra components, cycles, or missing configurations.
Counting (2) gives n+4 sum_i 2^(n-i), which is (4).

For defect preservation, the two sides of every bulk rule have the
same allowed predecessor set and the same allowed successor set,
and the internal pair is allowed on both sides. C0 and R1 have
identical allowed neighbors. The boundary substitutions preserve
the adjacent internal constraint and the boundary constraint. This
checks each Q separately, even if surrounding letters are illegal.
QED.

This is a concrete version of the "longest penalty-free path"
property for a clock. There is only one such maximal path. A quantum
program must still be attached without altering this property.

## 2. A positive initialized core

Attach a work space W=C^w to every clock letter. For each directed
rule a->b choose a unitary U on the work sites touched by that rule.
The same rule has the same U at every position. One-site rules have
one-site U's. The initialization move IZ->LI is assigned identity. The activation
I->R0 can carry a fixed seed unitary after the last initialization
check. It is identity in binary_clock_checks.py. No gates run while
the head is I.

For a rule with support on one or two sites put

    P(a,b;U) = (J/2)[ |a><a| tensor 1 + |b><b| tensor 1
                   - |b><a| tensor U - |a><b| tensor U^dag ], (5)

where J>0. This operator is positive because its quadratic form is
(J/2)||v_b-U v_a||^2. Add K times each legality projector.

Let p be a projector on W and add

    eta sum_j |I><I|_j tensor p_j,       eta>0.           (6)

The accepted input space is A=(ker p)^{tensor n}. Endpoint-specific
positive pins may be added during the endpoint I configurations.
For the statements below, take those pins to be commuting projectors
on the corresponding work factors, so they select a tensor-product
input space. The numerical example uses a logical qubit and an
ancilla at each site, p=1_logical tensor |1><1|_ancilla.

This has the uniform-core form required by C8's wrapper,

    H_n = sum_(j=1)^(n-1) h0_(j,j+1) + a_1 + b_n,       (7)

with h0,a,b positive and independent of n. Explicitly, h0 contains
K Q, all two-site terms (5), and half of each one-site turn term and
pin in each slot. The left boundary a contains K Q_left, the restart
terms, and the missing half of the one-site terms. The right boundary
b contains K Q_right, the activation term, and the missing halves.
This formula also works at n=1, where a+b supplies both halves.
The core local dimension is 9w. It is not the final ring dimension.

Let |t>, t=0,...,C_n-1, be the path (3). Let W_t be its cumulative
work circuit, with W_0=1. The complete ground isometry is

    V |psi> = C_n^(-1/2) sum_t |t> W_t |psi>,
    |psi> in A.                                        (8)

Theorem 2. Equation (8) is the complete zero space of (7). There
are no additional ground states in malformed clock sectors. A lower
bound on the spectral gap is

    min{ K,
         J[1-cos(pi/C_n)],
         eta/(2 C_n), J/(4 C_n^2) }.                    (9)

Proof. In the legal sector, conjugate by sum_t |t><t| tensor W_t.
Propagation becomes the path Laplacian

    L_C=(J/2) sum_(t=0)^(C-2) (|t>-|t+1>)(<t|-<t+1|).

Its spectrum is J[1-cos(pi k/C)], k=0,...,C-1. During all n I
configurations W_t=1. The pins therefore become |j-1><j-1| tensor
p_j and commute on the input space. Decompose that space into their
joint eigenspaces. In the accepted eigenspace all pins vanish and
the unique clock zero mode is constant. In an invalid eigenspace
at least one time s has potential eta |s><s|.

For any scalar path vector z and any vertex s, the telescoping
identity z_t=z_s+sum of increments gives

    sum_t |z_t|^2
        <= 2 C |z_s|^2 + 2 C^2 sum_edges |z_(t+1)-z_t|^2.

Comparing with eta |z_s|^2+(J/2) sum_edges |difference|^2 gives the
last two entries of (9). The same inequality holds for vector-valued
z. In an illegal sector at least one conserved Q equals one. The
propagation and pins remain positive, so energy is at least K.
QED.

This proof allows arbitrary unitary decorations after initialization,
including entangling gates. Their circuit may create nontrivial codes,
but it is the circuit generated by the fixed rules, not an arbitrary
externally scheduled circuit. There is no program register selected
by assumption.

For a physical left or right endpoint operator F, (8) gives the exact
compression

    V^dag F V = (1/C_n) sum_t P_A W_t^dag F_t W_t P_A,   (10)

where F_t is its diagonal clock block at time t. Off-diagonal clock
blocks require their corresponding two-time matrix elements instead.
No useful Clifford identity follows from (10) without a suitable
choice of rule gates and initialized work space.

## 3. Hopping clock and cell-length selection

Replace (5) by its off-diagonal part,

    -J[ |b><a| tensor U + |a><b| tensor U^dag ],          (11)

while retaining legality and initialization penalties. In the accepted
input sector its spectrum is

    E_k(n)=-2J cos(pi k/(C_n+1)),  1<=k<=C_n.            (12)

Its ground history has amplitudes

    phi_t=sqrt(2/(C_n+1)) sin(pi(t+1)/(C_n+1)).           (13)

The invalid input sectors have the same adjacency matrix plus
nonnegative potentials. Their lowest energy is strictly greater
than E_1: equality in the variational inequality would require the
strictly positive path ground vector to vanish at a penalized time.
For a quantitative, conservative lower bound on that increase, put

    delta=2J[cos(pi/(C+1))-cos(2pi/(C+1))],
    r_min=2 sin^2(pi/(C+1))/(C+1).

The increase is at least delta eta r_min/(delta+eta). Indeed,
A=H_path-E_1>=delta(1-|phi><phi|). Adding eta |s><s| and diagonalizing
the two-dimensional span of phi and |s> gives a smaller eigenvalue

    [delta+eta-sqrt((delta+eta)^2-4 delta eta |phi_s|^2)]/2
        >= delta eta |phi_s|^2/(delta+eta).

Additional pins only raise the energy. This bound need not be sharp.

There is also a uniform illegal-clock estimate. A word with v defects
has at most v+1 heads, since cutting at each forbidden internal pair
leaves words with at most one head. Each head admits at most two
incident moves. The weighted version of 2 Re<u,Uv> <= ||u||^2+||v||^2
therefore gives H_hop>=-2J N_head. Since v is conserved,

    H_hop+K v >= -2J+(K-2J)v.                            (14)

For K>=4J every illegal sector has nonnegative energy, while the
accepted history has negative energy. Thus it wins at every n.
Initialization pins are nonnegative and do not weaken (14).

The same clock can select cell lengths. With a scalar marker # and
nine work-carrying letters, d=9w+1. Replace the boundary conditions
by the allowed pairs

    #L, #h, 0#,1#,Z#,h#,

along with (1). The pair ## is forbidden. The restart and activation
moves now use # as their second support site. The marker is static.
Include onsite turn terms once and a chemical potential mu per marker.
This is one repeated Hermitian nearest-neighbor bond.

A legal cell with n interior sites occupies b=n+1 ring sites when
its left marker is counted. Its accepted energy is mu+E_1(n).
Define

    a_n=2J[1-cos(pi/(C_n+1))], delta_mu=2J-mu,
    e_n=(a_n-delta_mu)/(n+1),
    d_n=(n+2)a_n-(n+1)a_(n+1).                          (15)

Theorem 3. For any n>=2, choosing d_n<delta_mu<d_(n-1) gives n as
the unique minimum-density legal cell length. The choice has
0<delta_mu<2J. With K>4J, on rings divisible by n+1, illegal sectors
and marker-free sectors are above the tiling by these selected cells.
For w=1 the ground multiplicity is n+1, the translates of the tiling.
For nontrivial accepted work inputs, each cell also carries its code.
This is a frozen selector, not a conformal ring model.

Proof of convexity. Extend n to x>=1 and put

    u=4*2^x, y=u+x-3, a(x)=2J[1-cos(pi/y)].

Here y>=6, y'>0, y''=(log 2)^2 u. Direct differentiation gives

    a''(x)=2J*pi/y^4 *
       [(2y(y')^2-y^2 y'') sin(pi/y)
          + pi(y')^2 cos(pi/y)].                        (16)

The two trigonometric factors are positive. Also

    2(y')^2-y y''
      =(log 2)^2 u(u-x+3)+4(log 2)u+2 > 0.

The last inequality follows from u-x+3>0 for x>=1. Hence a is
strictly convex and decreasing. Therefore

    d_n-d_(n+1)=(n+2)(a_n-2a_(n+1)+a_(n+2))>0,
    e_n-e_(n+1)=(d_n-delta_mu)/[(n+1)(n+2)].             (17)

These identities prove the global minimum and e_n<0. For n>=2,
delta_mu<d_1<3a_1=3J(2-sqrt(3))<2J, so mu>0.

For the full ring, defect positions and markers are conserved. Cut
at v forbidden bonds. Every resulting legal string decomposes into
complete marker-bounded cells and at most two marker-free end
fragments. A fragment has at most one head and hopping energy at
least -2J. Any marker chemical potentials left in those fragments
are nonnegative. Complete cells have energy at least e_n times their
length, also for invalid initialized inputs by the preceding argument.
Since e_n<0, replacing the total complete-cell length by the full
ring length only lowers that bound. Thus the illegal sector obeys

    H_ring >= e_n L + (K-4J)v.                           (18)

A legal marker-free ring contains no head. Its possible words use
only L, only Z, or only binary suffix letters, so its energy is zero.
The selected tiling attains e_n L. This proves the claim. QED.

There is no claim here of the target d0 near 45 or 130. The compiler
and payload alphabet must still be counted. Nor does replacing the
C7 walker by this clock automatically give C8's mobile-cell interface.
The positive core in Section 2 can, however, be used as an input to
that interface once its boundary code is useful.

## 4. Verification and the remaining issue

Run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 binary_clock_checks.py

The script uses the local rule list, not the analytic path enumeration,
to enumerate the complete clock space at n=1,...,5. It checks every
edge against every defect projector, counts all connected components,
and checks the head and degree bounds. At n=5 there are 59049 words,
14509 components, and exactly one legal component of 129 vertices.
Full illegal-component spectra were computed through n=4, including
multi-head components and cycles; their FF energies are >=8 for K=8,
and their hopping energies obey (14).

Both complete legal spectra were compared with (12) and the Laplacian
formula through n=8, C_n=1028, with maximum error below 4e-15.
Cell selection was checked for selected n=2,...,39 against all lengths
1,...,99 using a stable sine formula for a_n.

A noncommuting work-circuit example uses n=2, w=4, one logical qubit
and one initialized ancilla per site. Every counter edge is dressed
with the fixed rule unitary, with identity during initialization.
The legal physical matrix has dimension 224. Its complete ground
space is four-dimensional, with history residual below 1.7e-15.
The gap is 0.0062877901, above the bound 1/784=0.0012755102.
These are initialization and clock checks, not CFT estimators or a
high-channel port test. Output is binary_clock_checks.json.

The next problem is now more specific. Find rule unitaries and a
constant-overhead work alphabet whose generated circuit implements
the clean boundary readouts of RESULT Section 25, with all workspace
and any program selected by the I-sweep pins. Alternatively, compile
a reversible interpreter into this deterministic clock without adding
unchecked input sectors. Fixed arbitrary gates on the rules give only
one recursively structured circuit family. They are not a proof of
universality or of the required high-channel compression.


A separate full periodic-bond check is

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 clock_ring_checks.py

It enumerates every word and every component for the ten-letter ring
clock (nine letters plus #) at L=3,4,5, of dimensions 1000,10000,100000.
Choosing the selected interior length n=L-1 gives ground multiplicities
3,4,5 and energies -0.3874779687, -0.0827385044, -0.0226164370, agreeing
with (12) and the marker potential within 3e-16. Every illegal component
obeys (18). Their minimum energies exceed 6.08 for K=8. Some illegal
components are larger than the legal one; the conserved penalties, not
a bound on their unpenalized path length, exclude them. These are complete
small-ring selector spectra, not conformal data. Output is
clock_ring_checks.json.
