# What remains after the clock construction

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
