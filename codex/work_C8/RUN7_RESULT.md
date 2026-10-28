# C8 seventh-run result snapshot

## Verdict

Seventh-run update. The verdict remains nonconclusive under the task's
criteria. No verified microscopic counterexample, full-class bound, or
fixed-d construction with unbounded c has been obtained. Sections 36-39
contain the new work and supersede earlier next-step rankings.

The clock part of the initialized-core problem now has explicit local
rules. A nine-letter positive NN clock has exactly one legal history at
every length, including an initialization sweep. A fifteen-letter extension
makes 2^n full work-tape sweeps at length n, still with one legal component.
All malformed clock sectors are excluded by conserved local penalties.
Both constructions have explicit finite-length gaps. Their hopping forms
can select any cell length by a chemical potential, but remain frozen
selectors when used directly on a ring.

Each full sweep executes a fixed local work update. An injective graded
collision table has an explicit completion to a local permutation, and
input pins select its complete quantum history code. A small packet machine
checks this interface. It gives a nonlocal diagonal parity readout, while
its off-diagonal readout vanishes because the bit remains recorded on the
tape. The high-channel program must clear every such record before its
port becomes active. That program and its alphabet count are still missing.

Four new scripts passed, including exhaustive open-clock spaces through
59049 words, complete periodic-ring spaces through dimension 100000,
quantum initialization, and the full local-update matrix. None is a CFT
measurement. All new proofs, limitations, reproduction commands, and the
remaining finite-program contract appear below.

The earlier two-length wrapper, scalar marker transport, exact leading
port equalization at d=96q+48 for a q-state core, and combined projected
length/payload spectrum remain available. They still need an actual
high-channel core and finite-coupling infrared control. C13 handles the
latter question separately; its final progress was read, not rederived.

No external files were changed, and no network, cluster job, subagent,
independent referee, commit, or push was used. Earlier restricted results
are not promoted to statements about the full class.


## 36. An initialized constant-alphabet binary clock

The numbered subsections and equation references in this section are local.

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


## 37. Full sweeps and a local reversible work interface

The numbered subsections and equation references in this section are local.

Status: proved clock and gate interface, with explicit local rules and
checks. The final clean-readout machine is not supplied. This note does
not assert a universal compiler or a high-channel code.

## 1. Fifteen letters and one complete path

Start with initialized_binary_clock.md. Add six letters

    A0,A1,D0,D1,E0,E1.

The old counter and initialization pair rules remain, except that L
may precede only C0,C1,R0,R1,I. Add exactly the pairs

    Aa Ab, Aa Db, Aa Eb, Db c, Eb c,    a,b,c in {0,1}.  (1)

At the left boundary also allow A0,A1 and all new heads; at the right
allow all new heads. A new sweep configuration is

    A_(x1) ... A_(x_(i-1)) H_(xi) x_(i+1) ... x_n,
    H=D or E.                                           (2)

Replace the left-boundary restart Rb->Cb by

    left Rb -> Db,
    Da b -> Aa Db,
    right Db -> Eb,
    Aa Eb -> Ea b,
    left Eb -> Cb.                                      (3)

The first step enters an outward sweep. The right boundary reverses
its direction. The return restores the binary tape and exits at the
left. All other counter rules are unchanged.

Theorem 1. At length n>=1 the complete legal clock graph is one path
with

    C_n=(2n+4)2^n+n-4                                   (4)

vertices. It makes exactly 2^n full outward/return sweeps. Every
local defect projector, including boundary defects, is conserved.
All positive-core and illegal-hopping-sector estimates of the
nine-letter clock apply with this C_n.

Proof. Legal configurations are exactly the old counter and I words,
together with (2). The boundary and pair constraints force one head
and the appropriate prefix/suffix alphabets. For each of the 2^n
binary values, the old restart edge is replaced by one detour through
2n sweep configurations. Each word (2) occurs on the unique detour
with its tape value. Thus the single old path acquires 2n2^n vertices
and no additional component. This proves (4) and the sweep count.
The predecessor and successor sets of each rule match, just as in
the nine-letter proof. Every defect is therefore conserved. A legal
fragment still has at most one head, and every head has at most two
incident transitions. The positive and hopping comparisons follow.
QED.

The shifted hopping energy is again decreasing and strictly convex
in n. To check this without relying on numerics, use the earlier
formula for a'' with

    u=(2x+4)2^x, y=u+x-3, k=log 2+1/(x+2),
    beta=(log 2)^2+2 log 2/(x+2).

Then y'=ku+1, y''=beta u, and

    2(y')^2-y y''
      = beta u(u-x+3)+2u^2/(x+2)^2+4ku+2 > 0.           (5)

Here u-x+3>0 for x>=1. Both sine and cosine in the a'' formula are
positive, since y>=10. Thus the same chemical-potential construction
selects every interior length n>=2. It is still a frozen selector.

## 2. One sweep executes a partitioned local update

Let the work space be

    W = S tensor B_left tensor B_right,

with dim B_left=dim B_right=r. S is a static tape register. Assign
the same two-site unitary G, acting only on B_right at the left site
and B_left at the right site, to every outward rule in (3). Distinct
G bonds have disjoint subregister supports, so they commute even
though the physical bonds overlap.

Assign a one-site work unitary F to the departing right site's work
space in every return rule Aa Eb->Ea b. Assign F at site 1 to the
left exit Eb->Cb. Assign identities to the entry, bounce, old counter
rules, and initialization moves. Then each full sweep executes

    U_n = (tensor_(j=1)^n F_j) product_(j=1)^(n-1) G_(j,j+1),
    F_j=F.                                              (6)

This is an equality of matrices on the full work space. The outward
pass first applies every G once. The return then applies F_n through
F_2, and the exit applies F_1. The F's act on separate sites.

Boundary collisions can differ. To implement F_left at site 1, use
that unitary on the left exit. To implement F_right at site n, use
F^{-1} F_right on the right bounce. The ordinary return subsequently
applies F, so the net right collision is F_right. For n>=2, (6) then
holds with these two boundary replacements. The n=1 special case
can be defined separately and is not needed by the wrapper.

A seed unitary on the activation I->R0 can create the initial work
head after the last initialization check. The checks are still
conjugated by identity at every I time. The assumption of identity
activation in the first note is therefore unnecessary for the
initialization theorem; identity on the IZ sweep is sufficient.

This is a fixed alphabet, repeated-bond implementation of 2^n
iterations of one spatially uniform local update. It does not hard
code 2^n time labels in the local dimension. Its complete legal
Hamiltonian is gauge equivalent to a path tensored with the full
work space, regardless of whether U_n itself has cycles or extra
invariant sectors. Input pins, rather than a claim about U_n's
entire dynamics, select the work input.

## 3. Mechanical interface for a reversible local machine

Here is a precise finite-state input that this clock can execute.
The static tape alphabet is S. Each transport lane is either empty
or holds a packet state from a finite set Q. A local collision table
specifies an injective partial map on configurations consisting of
one packet and one tape symbol. Its output also has one packet,
possibly in the other lane, and a new tape symbol. Boundary collision
tables may differ. Require the specified map to preserve the chosen
on-site grading, if any.

An injective partial map between finite basis sets always extends
to a permutation. For each grading sector separately, list the unused
domain labels and unused range labels in lexicographic order and
pair them. Injectivity and grade preservation imply that those two
lists have the same size. This gives a definite matrix, not a
spectral existence assumption. The routine complete_permutation
in sweep_machine_checks.py implements this rule.

Choose G to swap B_right(j) and B_left(j+1). Choose F to swap the
two on-site lanes and then apply the completed collision permutation.
Away from the ends, a right-moving packet travels one site to the
right and a left-moving packet one site to the left before the new
collision. On the ends the unused outer lane reflects under the
on-site swap; a boundary collision may use a separate packet label
to signal this event. Equivalently use the separate F_left,F_right
of Section 2. Thus a specified reversible packet-and-tape machine
has a literal nearest-neighbor unitary implementation.

If each lane has r=|Q|+1 states and dim S=s, the work dimension is

    w=s r^2,              q_core=15 s r^2.              (7)

A port register of dimension q and a constant number of control bits
multiply s by O(q), not by the number of clock steps. If the collision
table and ordinary tape alphabet can be kept fixed as q changes,
this would give q_core=K q+O(1) for a fixed K. The conditional "if"
is the remaining program issue, not a property established by the
clock alone.

The complete initialized-space theorem is useful here. Pin the input
tape and empty packet lanes during the I sweep, retaining the desired
logical input factors. Seed one packet by the activation gate. No
malformed packet arrangement, extra head, or arbitrary initial program
can then be a zero mode. The circuit is a bijection on the entire
work space, and the path is connected. Any zero-energy history is
therefore determined by one input vector, which all initialization
checks constrain. A work-machine cycle introduces no additional
clock component.

These statements prove execution of a given finite collision table.
They do not prove that a particular high-channel readout algorithm
has already been converted into such a table. Appealing to a generic
universal machine without accounting for tape, parity, boundary
operations, and initialization would leave that step unproved.

## 4. Additional verification checks during the computation

Initialization need not select the final logical input subspace by
itself. Positive one-site checks can be attached to a work control
state encountered later. Let V_t>=0 be the resulting work operator
at clock slice t after the clock has been restricted to its unique
path. Assume 0<=V_t<=M_t and define the input-space operator

    S_check = sum_t W_t^dag V_t W_t.                    (8)

Its complete kernel A is precisely the set of allowed inputs. A
zero mode of propagation is a constant path after gauging; positivity
then says every check annihilates it if and only if its input lies
in A. Consequently the complete history zero space is (8) of the
previous note with this A, even when the checks do not commute.

There is a quantitative gap statement. Suppose S_check>=gamma(1-P_A)
with gamma>0. Let M>=max_t ||V_t||, and put

    Delta=J[1-cos(pi/C)], lambda=gamma/C.

The legal-sector gap is at least

    Delta lambda / [ (sqrt(Delta)+sqrt(M))^2+lambda ].   (9)

If M=0 and A is the full input space, the gap is just Delta; the
formula is intended for nontrivial checks. Include K as a minimum
to cover illegal clock sectors.

Proof. Gauge the circuit away. Let P project onto the constant
clock vector and write psi=x+eta, x=P psi. For psi orthogonal to
the full zero code, its work vector in x is orthogonal to A. If
E=<psi|H|psi>, then ||eta||^2<=E/Delta and

    <x|V|x> >= lambda ||x||^2,
    sqrt(<x|V|x>) <= sqrt(E)+sqrt(M)||eta||.

Eliminate eta and use ||x||^2>=1-E/Delta. The resulting inequality
is (9). QED.

This allows, for example, an oblivious program to test that all
input bits outside a computed logical region are zero. It is not
necessary to build a separate physical letter for every input
position. But the tests must actually be supplied and (8) analyzed;
merely declaring a chosen program sector does not meet this condition.

For the intended classical readout program there is a sharper simplification.
If every W_t is a basis permutation and every V_t is diagonal, all
W_t^dag V_t W_t commute. If each nonzero check has weight at least eta,
then each rejected basis input encounters a potential at least eta at
some time. The path-pin inequality from the nine-letter construction
therefore gives directly

    gap >= min{K, J[1-cos(pi/C)], eta/(2C), J/(4C^2)}.   (10)

This does not require the check times to be the same for rejected
inputs. Each diagonal input sector is a scalar path with nonnegative
potentials, at least one of which is >=eta. Accepted inputs have no
potential. Thus arbitrary superpositions of accepted classical
inputs form the exact quantum code, with the stated gap. No simulation
error or limit in the clock coupling is involved in this assertion.

## 5. Explicit packet machine and its limitation

sweep_machine_checks.py supplies a literal small machine. Its tape
symbol is a logical bit x, and each lane is empty or carries an even
affine accumulator bit p. The collision toggles p by x when exactly
one packet is present. The specified map is completed to a permutation
on the remaining basis states by the prescription above. There are
18 local work states. Activation seeds a right-moving packet with
p=0 at the rightmost site; all other lanes were pinned empty.

The packet bounces between boundaries. At its first left-moving
arrival at site 1, its accumulator is the parity of all input bits.
On the return trip the same updates uncompute it. This repeats. The
physical site-1 observable

    O_left = |L,p=0;R=empty><L,p=0;R=empty|
           - |L,p=1;R=empty><L,p=1;R=empty|              (11)

therefore compresses, in the uniform FF history, to

    alpha_n Z_1 Z_2 ... Z_n,   alpha_n>0.              (12)

It is zero when the packet is absent or right-moving. During a
partial outward sweep an arriving leftward packet occupies the right
lane until the local collision, so (11) also avoids reading an
unfinished parity. The binary carry intervals only idle a completed
configuration. This proves (12) at every clock slice, not just at
completed whole-sweep times.

For n=2,3,4,5, the active clock counts are respectively 8,8,28,36 out
of 30,79,192,449. Full initialized basis histories verify (12) with
zero numerical residual. The code dimensions are 4,8,16,32.

The accumulator X observable has zero compression. Every original
logical bit remains in its static physical tape, so changing only
the accumulator never connects two code basis inputs. This is an
explicit failure of the naive "compute the parity and read the head"
route to anticommuting boundary maps. Clean readout must transfer a
logical bit into the port and remove every other record of that bit,
as required in RESULT Section 25. This toy is a test of the machinery,
not a new claim about boundary complexity or high central charge.

## 6. Reproduction

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 sweeping_clock_checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 sweep_machine_checks.py

The first script exhausts the full clock alphabet at n=1,...,4, up
to 50625 words and 16477 connected components. Only one component is
legal. At n=4, 3855 illegal components contain cycles; each remains
penalized because every defect projector is conserved. Complete legal
hopping spectra through C_6=1026 agree below 5.3e-15. Length selection
is checked for n=2,...,34 against lengths through 79.

A separate quantum update test uses a four-dimensional work site,
a noncommuting onsite F, and transport swaps on disjoint subregisters.
At n=3 the full work dimension is 64. Every completed sweep agrees
with the corresponding power of (6), through all eight sweeps, with
error below 8e-15. This tests the update on the full work space, not
only on one initialized state.

The second script verifies the literal collision permutation, its
grade preservation, all initialized packet trajectories through n=5,
the nonlocal diagonal readout, and the vanishing off-diagonal readout.
Outputs are sweeping_clock_checks.json and sweep_machine_checks.json.


## 38. The remaining clean-readout machine contract

The numbered subsections and equation references in this section are local.

Status: a proved conditional reduction and a specification for the
unfinished readout machine. No collision table meeting the complete
specification is supplied. This is progress, not a conclusive answer.

## 1. The remaining finite program

The new sweeping clock removes the need for a separate program-dependent
clock alphabet. A work machine is now one fixed finite collision table,
executed 2^n times on a length-n tape, with all input work space pinned
at the beginning. Its counter is external to the work tape and has only
one legal history component. Malformed work configurations cannot be
additional zero modes because they are rejected at initialization.

The missing program must do more than run a conventional reversible
computation. At every active boundary readout it must remove all other
records of the logical bit being exposed. A history tape that retains
that bit defeats the required off-diagonal matrix elements.

Here is a precise sufficient interface for the first-order route. Take
q divisible by four, split its letters into q/2 even and q/2 odd letters,
and set m=q^2/8. Allocate two distinct even-odd edges to each flavor,
using a fixed lexicographic enumeration of all q^2/4 bipartite edges.
Use the readout operators and left/right sign conventions of RESULT
Section 25. Their physical definition is unchanged by this note.

The input work tape has m logical bits x=(x_1,...,x_m) and otherwise
initialized work space. The program can first verify that any surplus
input bit slots are zero. It must execute the following finite loop,
with all timing independent of the logical input:

    for boundary u in {left,right}:
        for f=1,...,m:
            prepare the clean (u,f) readout with port flag off
            turn the port flag on
            keep the readout clean for at least one clock vertex
            turn the port flag off
            undo the readout preparation
    return to the same baseline and repeat

Only the selected port is active. For left readout define
p=parity(x_1,...,x_(f-1)); for right define
p=parity(x_(f+1),...,x_m). The clean active configurations have the form

    |ell_f(x)>_left tensor |g_(left,f,t)(z_f(x))>,
    |r_f(x)>_right tensor |g_(right,f,t)(z_f(x))>.        (1)

Here ell_f,r_f are exactly the graded edge readouts of Section 25,
and z_f is the input string with bit f deleted. The rest state g is
normalized, records z_f injectively, and has no dependence on x_f.
It includes every other physical register, head state and scratch bit.
At a given active time the boundary and flavor are independent of x.
The record of z_f may depend on time. These conditions must hold at
partial sweep configurations as well as at completed machine updates.

To retain the independently tunable flavor densities needed by the
incommensurate target, also write the proper edge color
(a+b) mod(q/2) in the adjacent site's port register during each readout,
as in RESULT Section 25. This record depends on the chosen edge and
z_f, not on x_f. The associated diagonal two-site readout then compresses
to the flavor Z operator with the same active-time counting. It uses the
(q+1)-state register already allocated at that neighboring site. The
final machine must implement and clear this color record as well; the
conditional Clifford identities alone do not establish the detuned IR
target.

All other times have both flags off. The physical endpoint operators
are zero in the inactive flag sector. Preparation may then use arbitrary
scratch records, including records of x_f, provided it clears them
before activating the flag. Reversal must turn the flag off first.

The finite primitive operations already available in the earlier
scheduled circuit suffice at the algorithmic level: adjacent data
swaps, reversible even-control arithmetic, parity accumulation, transfer
of a graded logical bit to an even/odd port reference letter, cycles
within each grade, and an active-flag toggle. A transfer acts by the
transposition (carrier=1,port=e0) <-> (carrier=0,port=o0), fixing
(carrier=0,port=e0). All these preserve total grading.

There are two points that the final collision table must address.
First, its controller and ordinary tape alphabet must remain bounded
as q changes. Flavor indices and loop counters must live in tape
registers, not in one distinct head state per flavor. A possible source
of q is the order of a within-grade cyclic port permutation, which the
controller can measure by a reversible loop while the flag is off.
Second, controlled operations depending on p must take the same number
of steps for both values of p. Otherwise binary-clock dwell times can
produce state-dependent readout weights, even if each branch looks
correct at its endpoint. The ordinary history clock is uniform in
clock vertices, not in completed work-machine updates.

These are proposed finite-program implementation choices, not a supplied
finite transition table or an invoked universality theorem.

## 2. Conditional high-channel core theorem

Assume such a collision table, boundary tables, input tests, and a finite
length n are given, and that the complete first readout loop finishes
before the 2^n available sweeps. Assume its number r of lane states and
ordinary tape dimension s0 are independent of q. Give each site one
(q+1)-state port register and an active flag, even at sites where that
register stays inactive. The work dimension is

    w=2 s0 r^2 (q+1),
    q_core=15 w=K(q+1),   K=30 s0 r^2.                  (2)

This count deliberately uses a tensor product rather than claiming a
smaller unverified direct-sum count. It is linear in q if the stated
machine condition holds.

Then the sweeping-clock Hamiltonian is a uniform positive NN core with
one-site positive boundary terms and a complete 2^m-dimensional zero
code. Its two endpoint compressions are

    V^dag A_(2f-1,2f) V = alpha_(left,f) Gamma_(2f-1,2f),
    V^dag B_(2f-1,2f) V = alpha_(right,f) R_(2f-1,2f),
    alpha_(u,f)=N_(u,f)/C_n > 0,                        (3)

where N_(u,f) counts active clock vertices for that readout. Different
flavors and boundaries may have different counts. Rescaling their
physical operators separately compensates these positive scalars.
The norm cost can be large but is finite for fixed q,n.

Proof. The full clock and work-space kernel is determined by one input
vector. Initial pins reject arbitrary workspace, head, and program data.
The supplied later tests leave exactly the m logical bits unconstrained.
The initialization and check-kernel theorems in sweeping_clock.md imply
that there are exactly 2^m zero states, with a positive finite-length gap.
For a classical permutation program
with diagonal checks of weight eta, the explicit lower bound is

    min{K_penalty, J[1-cos(pi/C_n)], eta/(2C_n), J/(4C_n^2)}.

Here K_penalty is the clock-legality penalty, distinct from the
alphabet constant K in (2). The bound follows by decomposing into
input basis strings after gauging away the permutations. Every rejected
input meets at least one positive path pin; every accepted input meets
none. This does not require a globally separated many-cell low band.

A diagonal-in-clock physical endpoint operator compresses to the sum of
its same-time matrix elements divided by C_n. It vanishes at every
inactive time. At a clean active time (1), only input strings differing
in bit f can be connected by the corresponding edge operator, because
the rest records z_f. Their rest states agree. The two edge labels and
the prefix or suffix parity give the Gamma or R matrix element, as in
Section 25. An operator belonging to another flavor has zero matrix
element on that edge. Counting the active times proves (3).

There is no extra work-machine-sector assumption in this argument.
The collision gates are genuine permutations on the full local space,
and all propagation components over the unique legal clock are related
by the gauge transformation. Work-machine reversibility is used to
supply local unitaries, not to assume that every machine orbit is valid.
QED, conditional on the explicitly stated machine and tests.

For comparison with the existing wrapper, its physical ring dimension
would be

    d=96 K(q+1)+48,        c_target=m=q^2/8.             (4)

For K>=1, choosing q=1536 K, which is divisible by four, makes
c_target>d-1. This is conditional resource accounting. K is not supplied
by this run because the collision table is not supplied. The optional
length mode would add one to the target central charge but is unnecessary
for this inequality. Equations (2)-(4) do not prove a microscopic CFT.
The finite-coupling analysis remains separate, as in C13.

## 3. Why a generic history-preserving simulator is insufficient

Suppose an active basis history has

    |x> -> |t> |ell_f(x)> |g_t(x)>.

The boundary matrix element between x and x^f is proportional to

    <g_t(x^f)|g_t(x)> .                                 (5)

If a reversible simulator has retained an orthogonal copy of x_f in
its work or history tape, (5) is zero. Averaging over time does not
restore it because the clock states are orthogonal. Idling the copied
readout also leaves it zero. The packet-parity example in
sweep_machine_checks.py verifies this failure directly: its diagonal
readout is nonzero, while its accumulator X compression is zero.

A generic reversible computation theorem allows workspace that is
clean only at the end of the entire computation. That is insufficient
here unless the active flag is postponed until a clean checkpoint.
Conversely a genuinely clean, parity-preserving implementation of the
scheduled readout permutation would meet this part of the interface.
It is not necessary to reject all malformed intermediate work strings
locally, since the positive history construction already rejects their
initial inputs. This is the useful simplification supplied by the clock.

The next run should write and test that clean reversible finite program,
including surplus-input tests and its local collision permutation. It
should not repeat clock construction, cell selection, or the toy parity
readout, and should not report (4) as an explicit d0 until K is counted.


## 39. Seventh-run audit and ranked continuation

The verdict remains nonconclusive. No high-channel collision table, explicit
microscopic counterexample, full-class upper bound, or fixed-d unbounded-c
construction was completed. The new contribution is the explicit clock
and initialization mechanism, together with a precise remaining program
interface. The numerical work concerns these finite constructions only.

Four new scripts passed with single-thread BLAS and compiled:

- binary_clock_checks.py exhausts open clock graphs through 59049 words,
  checks every conserved defect, diagonalizes complete small illegal
  components, checks both legal spectra through 1028 vertices, and verifies
  a 224-dimensional entangling initialized core.
- sweeping_clock_checks.py exhausts the fifteen-letter graph through
  50625 words and checks a full 64-dimensional work update through eight
  sweeps. The largest gate-composition error is below 8e-15.
- sweep_machine_checks.py specifies an actual eighteen-state collision
  permutation and checks all initialized inputs through five sites. Its
  parity readout succeeds; its off-diagonal accumulator readout vanishes.
  That failure is recorded as a limitation, not omitted from the results.
- clock_ring_checks.py diagonalizes every component of full periodic clock
  spaces through dimension 100000. Ground energies and multiplicities
  agree with the selector theorem. These frozen rings are not CFTs.

The source audit checked rule orientation, the initialization seed, all
boundary constraints, individual defect conservation, positive boundary
terms, the two convexity calculations, head-degree and ring-fragment
bounds, local update order, and cleanup at active readout times. No
independent referee was used. Previous numerical suites were not rerun.
There was no network access, cluster job, subagent, commit, or push by
this run. Only this folder was modified.

C13 finished while this run was in progress. Its final PROGRESS reports
conditional Gaussian-phase control and calibrated projected-model checks,
but leaves the high-channel microscopic core and generic microscopic
central-charge estimates open. It also distinguishes unique finite-ring
ground states from convergence to one fixed vacuum sector at every
integer length. Those are source-reported results, not new proofs here.

Remaining work, ranked:

1. Supply the actual clean-readout collision table for Section 38. Use
   the fifteen-letter clock and its input pins. Store counters on the tape,
   keep the control alphabet independent of q, preserve grading, test
   surplus inputs, include colored density readouts, and verify both ports at every active partial-sweep
   time. Count s0 and r to make K=30 s0 r^2 concrete. A generic history
   simulation that retains copies of the exposed bit does not suffice.
2. Attach that proved core to the existing two-length wrapper and equalizer,
   and calculate the actual local finite-coupling corrections. The clock
   gap is exponentially small in cell length but positive for each fixed
   construction. Couplings must be chosen below the actual cell scales;
   no size-uniform low-band isolation follows just from a cell gap.
3. Use the C13 resonance and sector analysis to assess the resulting
   moving-cell Hamiltonian. Obtain calibrated microscopic estimators or
   a complete reduction to stated infrared inputs. This run did not
   duplicate that separate project.
4. Only after the first-order program is explicit, consider a universal
   quantum work machine for route A. The full-sweep clock supplies an
   initialized executor, but it does not by itself preserve prescribed
   excited-state spectral weights or supply a fixed-d unbounded-c family.

Do not redo the nine-letter counter, fifteen-letter full-sweep clock,
conserved-defect proof, initialization-gap proof, or bare-clock selector.
Those parts now have explicit local rules and all-sector checks. The
clean-readout work machine is the next construction task.
