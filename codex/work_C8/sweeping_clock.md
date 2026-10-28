# A full-sweep clock and a local work-machine interface

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
