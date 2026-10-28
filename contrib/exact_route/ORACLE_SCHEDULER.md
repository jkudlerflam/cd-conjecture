# A reversible neutral scheduler for the charge conserving circuit

## Precise claim

[proved, using the classical reversible-computation theorem cited below]
Fix a classical algorithm which, from a finite binary parameter string
\(x\), produces a nearest-neighbour circuit on a tape alphabet
\(\Sigma_q\) of dimension \(D(q)\). Suppose:

1. The circuit uses a fixed number of gate families \(U_k(q)\).
   Their matrices may depend on \(q\), but the number of families does not.
2. Every \(U_k(q)\) preserves each total species charge.
3. The circuit description, width and generation time are bounded by a
   polynomial \(P(r)\), where \(r\) is the numerical size parameter
   encoded by \(x\). This is a polynomial in \(r\), not a polynomial in
   the available physical cell length.
4. A valid initialized tape provides the classical parameter and blank
   neutral work tracks, the specified quantum input, and exactly one
   neutral interpreter head.

There is a single-head reversible interpreter for this family, with a
constant neutral tape alphabet and a constant number of control states,
which implements the circuit on every quantum input. The classical
schedule never branches on a quantum tape symbol. If its neutral tape
alphabet has size \(B\) and its finite control has size \(K\), its
ordinary tape dimension is \(B D(q)\), and the sweep construction uses

\[
D_{\mathrm{interpreter}}=(1+4K)B D(q).
\tag{1}
\]

Here \(B,K\) depend on the fixed generating/interpreting algorithm and the
fixed number of opcode families, not on \(q,r\), the circuit length or
the quantum input. They are not estimated to be small.

This is an existential constructive compiler lemma. It certifies the
constant-factor controller overhead; it does not certify a particular
minimal transition table, an initialized history Hamiltonian or its
illegal sectors. Hypothesis 4 is an input to this lemma, not its
conclusion.

An [explicit dispatcher table](../controller_instantiation/controller/README.md)
now implements a preloaded instance of this interface. Its exact count is
88 neutral tape symbols and 150 direction-fixed controls for the 17-family
preparation menu. It does not generate the program from blank: its program
symbols are conserved. The length-to-program computation and energetic
verification remain separate parts of the full compiler claim.

The [autonomous continuation](../autonomous_controller/README.md) now has
one finite table composing a generated binary-digit demonstration program
with quantum execution, endpoint readout and reset. It also freezes the
full sparse generator as 3,115 word instructions. The latter's arithmetic,
addressing and journal have not yet been expanded into this local interface.
An explicit [parameter query](../autonomous_controller/audit/PARAMETER_ORACLE.md)
can supply p from the already p-dependent local matrices, eliminating the
need to infer p from physical length; capacity and full-input checks remain.

## 1. Keep classical control separate from quantum scratch

An ordinary site is

\[
|c\rangle_{\mathrm{classical}}\otimes
|a\rangle_{\Sigma_q},\qquad c\in\Gamma,\quad |\Gamma|=B.
\tag{2}
\]

The classical factor has a fixed number of tracks, each with a fixed
finite alphabet. It stores program bits, integer addresses, classical
work, delimiters and a fixed number of head/position marks. All its
states have zero species charge.

The quantum factor includes the charged data, neutral *quantum* work
bits and loaded symbols \(T_{j,c,s}\) of CIRCUIT. These neutral quantum
bits must not be confused with the classical metadata. A branch label,
a queried species bit, or a routing control can be in superposition;
the interpreter never reads such a value into its classical control.
Instead it executes the corresponding controlled quantum gate.

Integer parameters \(p,q,N,j\), circuit addresses and counters are
binary strings on the classical track. They are not control-state
labels. In particular storing \(q\) does not introduce \(q\) classical
states. The allowed data alphabet itself has indexed \(H,T\) symbols,
but those remain in the single \(\Sigma_q\) factor in (2).

## 2. The classical compilation is a legitimate theorem input

[read] Bennett's 1973 construction converts any fixed ordinary machine
into a fixed reversible machine which preserves its input, writes its
output and clears its history. It separates read/write and movement
steps and arranges disjoint transition domains and ranges. See the
theorem and proof on pp. 527–529 of
[Logical Reversibility of Computation, original paper](https://worrydream.com/refs/Bennett_1973_-_Logical_Reversibility_of_Computation.pdf);
the [Princeton transcription](https://www.cs.princeton.edu/courses/archive/fall06/cos576/papers/bennett73.html)
is another source for the same theorem.

Apply this result to the fixed algorithm which generates the circuit
description, and to the fixed classical routines for addresses, loop
indices and navigation. This produces fixed finite control and a fixed
number of constant-alphabet neutral tracks. For the polynomial-size
circuit at issue, even retaining a polynomial amount of classical
history fits a sufficiently large polynomial cell.

The fact that Bennett states his main construction for several tapes
does not create several quantum payload lanes. Combine the *classical*
tapes into a fixed number of tracks in \(\Gamma\). Each simulated tape
head is a mark on its own classical track. A single physical head scans
between those marks, collecting the finitely many scanned classical
symbols, performing the simulated reversible transition, and updating
the marks. There are a fixed number of simulated tapes and a fixed
alphabet on each; the tuple temporarily held in finite control therefore
has constant size. No species label is in that tuple.

One can retain the old tuple through this operation and erase it using
the inverse simulated transition. Travel and return phases and marked
positions make the scans reversible. This changes time by a polynomial
factor but uses a constant number of tracks and no additional charged
record. The one-tape reversible reduction can also be used directly:
[Morita, Shirasaki and Gono, original publication](https://globals.ieice.org/en_transactions/transactions/10.1587/e72-e_3_223/_p)
explicitly converts an ordinary machine to a one-tape reversible machine
with a fixed alphabet.

Stationary read/write transitions are finite partial bijections on
neutral symbol and control pairs. Add finitely many phase tags so that
their domains and ranges are distinguished, and extend these partial
bijections to permutations on the unused pairs. Shift steps are
reversible head moves. Equivalently, use the read/write/shift normal
form of Bennett's proof. A required stationary step can be surrounded
by an identity move and its return, so the move-then-collision interface
of SWEEP_INTERFACE imposes no restriction. All these changes multiply
\(B,K\) by constants only.

## 3. Insert quantum oracle calls explicitly

The classical generator first creates a list of instructions
\((k,i)\): apply opcode \(k\) to physical sites \(i,i+1\).
One-site gates are treated as \(U\otimes I\). The list is classical and
independent of the quantum input.

For each instruction:

1. Reversibly compute its opcode and operand address, retaining the
   instruction counter and program.
2. Reversibly move the neutral interpreter head to the left operand.
   Only classical metadata controls this path.
3. Enter a distinct finite control state \(\mathrm{EXEC}_k\).
   Move right once and apply \(U_k(q)\) on the two quantum records.
4. Move left once with the identity quantum opcode. The head is back
   at its pre-call site, in a distinct post-call classical phase.
5. Reverse the classical navigation and instruction computation, using
   identity quantum moves, then reversibly increment the instruction
   counter.

Step 5 does not apply the inverse quantum gate. It is the inverse of
the *classical preparation/navigation*, in a separately tagged control
phase. Quantum execution occurs only at the single designated step 3.
This is the usual compute–use–uncompute structure with a quantum oracle
between the two classical computations.

In the sweep interface, the execution block is explicitly

\[
|H^0_{R,\mathrm{EXEC}_k,(c,a)},(d,b)\rangle
\longmapsto
\sum_{a',b'}(U_k)_{a'b',ab}
|(c,a'),H^1_{R,\mathrm{EXEC}_k,(d,b')}\rangle.
\tag{3}
\]

Use its adjoint on the reverse block, and the identity elsewhere.
The return call uses the same formula with \(U=I\) and a distinct
opcode. Different instruction families, direction and phase tags
give orthogonal blocks. The local collision changes the classical
phase from execute to return and subsequently to the reversing
classical subroutine.

Classical computation blocks use \(U=I\) on \(\Sigma_q^{\otimes2}\).
Their neutral collisions act on \(\Gamma\) and finite control only.
Thus every constituent gate is defined on the entire local Hilbert
space and conserves charge there, including outside initialized inputs.
There is no appeal to an arbitrary isometry synthesis.

If \(J_c\) denotes placement of one head and classical configuration
\(c\) over the quantum tape, induction over machine steps gives

\[
\mathcal U^t J_{c_0}|\psi\rangle
=J_{c_t}V_t|\psi\rangle.
\tag{4}
\]

The configuration \(c_t\) is independent of \(|\psi\rangle\), and
\(V_t\) is exactly the prefix of executed circuit gates. Equation (4)
holds for all superpositions, including superposed branch labels and
species-query results. In particular all logical inputs have identical
head trajectories, running times and selected readout times.

## 4. The gate menu for this application really is finite

One possible menu consists of identity, whole-record SWAP, QUERY, LOAD,
EXCHANGE, SHIFT, neutral type conversion, neutral bit gates \(X,H,T\),
neutral CNOT and a head-bit/ordinary-bit CNOT. Include adjoints as needed.
Neutral Toffoli is implemented by its fixed Clifford+T decomposition;
it is not a new growing family of opcodes.

For the input verifier also include a constant number of TYPE-query
families. The record types are ordinary charged data, blank, ordinary
neutral bit, unloaded \(H\), and loaded \(T\). For each bit of this
five-valued type code, a diagonal control XORs that bit into an adjacent
neutral query bit. The number of type bits is constant, independently
of \(q\). An invalid type can therefore be rejected before interpreting
species-query output. This operation never moves a species label into
classical finite control.

There is one initialization detail if all ordinary \(\Sigma_q\) records
are initially untrusted. Trusted binary scratch alone does not supply
a trusted indexed \(H_{j,c}\) query register. A minimal interface is to
pin **one boundary** \(\Sigma_q\) record to \(H_{0,0}\) during the clock's
initialization sweep. This is a fixed neutral reference, not an extra
logical degree of freedom or an extra tensor factor in the alphabet.
Such a boundary pin is allowed in the marker-bounded core.

For a helper at the right boundary, move it left by a compact SWAP walk
until it is immediately to the right of the record being checked. Do
not cross that target record. Shift its index to \(j\), apply QUERY or
TYPE, copy its result bit into a separately initialized trusted neutral
bit track, apply the same query again, and shift its index back to zero.
Then reverse the SWAP walk. The target and every displaced
\(\Sigma_q\) record return to their original positions, the helper
returns to \(H_{0,0}\), and the trusted bit holds the queried function.
The complete map on a target superposition is

\[
\sum_a\psi_a|a\rangle\,|H_{0,0}\rangle\,|0\rangle
\longmapsto
\sum_a\psi_a|a\rangle\,|H_{0,0}\rangle\,|f_j(a)\rangle.
\tag{5}
\]

Malformed records may themselves have \(H\) or \(T\) type. This does
not confuse the scheduler: it follows classical position marks, and
the known SWAP sequence identifies the transported tensor factor.
An untrusted head-like quantum symbol is never treated as the neutral
interpreter head. TYPE can act on every first-slot record type, while
species QUERY may act trivially outside the ordinary charged-data type;
the acceptance predicate rejects those types separately.

The fixed boundary helper is omitted from the untrusted-input count.
All other accepted untrusted neutral records must be the unique blank:
TYPE plus species queries do not distinguish ordinary bit 0 from bit 1.
Convert designated accepted blanks to zero bits with an explicit
blank-to-bit0 involution after verification; include its inverse in
cleanup. Trusted verifier Q scratch is separately pinned. The
[format audit](../controller_instantiation/audit/README.md) proves and
checks this exact input contract, including malformed H/T records.

After the format check, the other quantum work records are known to
have their specified blank/reference states on the accepted subspace.
The helper may then be moved into the auxiliary workspace and a charged
output port routed to its former boundary position. It contributes
zero to every charge offset. This avoids an extra \(O(\log q)\)
trusted-index tensor factor and preserves (1).

For the exact variant add the single indexed controlled-rotation
family in [EXACT_AMPLITUDES](EXACT_AMPLITUDES.md). Its \(O(\log q)\)
index values live in \(\Sigma_q\); SHIFT selects an index by a fixed
sequence. Its block matrices, and the species-query matrix, depend
on \(q\). Neither dependence increases the number of EXEC labels.
All gate schedules, including selection of an angle index, are
computed from the classical description.

To see explicitly why the loaded-head operation does not cost \(q^2\),
consider its nontrivial action

\[
(T_{j,1,s},D_t)\longmapsto(T_{j,1,t},D_s).
\tag{6}
\]

When (3) executes it, the outer head initially overlays the one tape
symbol \(T_{j,1,s}\). Afterwards the left site contains
\(T_{j,1,t}\), and the outer head overlays the one right-site symbol
\(D_s\). Neither site ever contains independent \(s,t\) registers.
The pair Hilbert space necessarily has a \(q^2\) sector, as any two
ordinary \(q\)-species sites do, but the **single-site** alphabet is
only a constant times \(D(q)\).

## 5. Space, physical length and boundary counting

There are two distinct size parameters. The sparse circuit width,
gate count, graph arithmetic and instruction list are polynomial in
\(p\) (and in a precision parameter if the approximate finite numerical
gate set is chosen). They need not scale with the unused physical
padding. Keep this computation in a work region of polynomial length
in \(p\); its single-head simulation, classical histories and cleanup
then also use polynomial length in \(p\).

Choose the physical cell length \(n\ge p^K\), for a fixed sufficiently
large exponent \(K\). The work region fits, with additional room for
the output transport and readout schedule. The local alphabet is
unchanged when \(K\) is increased. The exact amplitude variant has
no precision parameter.

Long transport to the physical endpoint must also use a compact loop,
not a list of \(O(n)\) addressed gates. A right-moving SWAP opcode gives

\[
|H^0_{R,\mathrm{SWAP},a},b\rangle
\longmapsto |b,H^1_{R,\mathrm{SWAP},a}\rangle.
\tag{7}
\]

Thus the outer head carries the same single quantum record \(a\)
rightwards, exchanging it with each encountered record until the
immutable RIGHT mark is reached. Repeating the corresponding left
SWAP walk reverses the permutation. Remember the starting position
with a neutral mark or a binary address; no step history is needed.
The control decisions depend only on those marks and the immutable
endpoints. The quantum data do not determine the stopping point.

For example, arrange the prepared left port at the first site and the
right port at a known interior site, then transport the latter to the
last site without visiting the first. Every intervening record is
shifted by the same predetermined permutation. This changes neither
endpoint compression: it only permutes the other tensor factors.
The inverse preparation later reverses the same walk. These SWAPs
exchange the \(\Sigma_q\) records, not static endpoint marks or classical
program tracks. Dedicated per-site READY flags can be set only after
transport is finished.

Consequently the gate-list part remains polynomial in \(p\), whereas
transport and idle walks may take polynomial time in \(n\) with only
constant mark tracks and \(O(\log n)\) counters. There is no circular
requirement to fit an \(O(n\log n)\) instruction list into \(n\) sites.

It would be incorrect to apply a history-storing simulation to an
arbitrary \(\operatorname{poly}(n)\)-time length-counting routine and
then assume its history fits in \(n\) sites. Here is a separate
reversible length count which avoids that mistake.

Give the classical tape static endpoint marks \(L I^{n-2}R\), and one
movable frontier mark on another neutral track. Initially the frontier
is at \(L\), and a binary counter near \(L\) is 1. Repeat:

* Shift the frontier one site to the right, with a tagged reversible
  move, unless it is already at \(R\).
* Walk to \(L\), increment the binary counter reversibly, and walk
  back to the frontier.

Each outward and return walk has its own direction/phase tag and stops
at a uniquely marked site; it is undone by the corresponding reverse
walk. The counter increment is a bijection on canonical nonnegative
binary integers. It can be implemented by computing \(x+1\), using the
inverse \(y-1\) to clear the old copy, and uncomputing the two short
Boolean histories. This uses \(O(\log n)\) counter/work symbols on a
fixed number of tracks. Arithmetic subroutines work only near \(L\);
long travel steps accumulate no history. The frontier occupies one
mark on a track already provided at every site, rather than a list of
the visited positions.

At the right endpoint the counter equals \(n\). The position and
counter together determine the predecessor at every macrostep.
Finite start/end phases distinguish entry, completion and reverse
traversal. A counter delimiter may grow as the integer grows, so
knowing \(\log n\) in advance is unnecessary. For all sufficiently
large \(n\), the \(O(\log n)\) arithmetic region fits; any finite set of
smaller cases may be handled by separate fixed control states.
This construction takes polynomial time (the repeated travel costs
\(O(n^2)\)) and \(O(\log n)\) nontrivial counter/work cells, plus the
constant neutral mark tracks.

[read] The general reversible time/space trade-off of
[Bennett 1989, author institution](https://research.ibm.com/publications/timespace-trade-offs-for-reversible-computation)
also replaces a time-length history by \(O(S\log T)\) work for
ordinary time \(T\), space \(S\), at polynomial time overhead.
Our explicit marked walks make it unnecessary to store or checkpoint
a length-\(n\) traversal history.

The static endpoint-mark pattern itself is compatible with positive
nearest-neighbour format checks: permit \(LI,II,IR,LR\) and pin the
two endpoints to \(L,R\). The scheduler preserves these marks. This
observation only addresses that one classical track; it does not
prove all input-format or history checks.

## 6. What the lemma closes and what remains

Combining this lemma, EXACT_AMPLITUDES and SWEEP_INTERFACE proves the
existence of fixed on-site and nearest-neighbour gates which run the
exact sparse encoding, from a specified initialized input, with local
dimension \(Cq\log q\) for a constant \(C\). All gates commute exactly
with every species charge. There is no approximation, no exponential
matrix synthesis and no extra \(q\)-valued control register.

The interpreter may end a finite computation prefix in a distinguished
classical state. An absorbing halt cannot be part of a reversible
unitary evolution. For a longer prefix append a reversible idle counter
or a separately tagged walk, of prescribed duration; its duration is
independent of the quantum input.

The Hamiltonian requires an additional initialized-core contract:
exactly the accepted input sectors, enough legal sweeps, penalties for
every unwanted component, and the desired endpoint readout.
The companion [EXISTENTIAL_CORE](../torus_bottleneck/EXISTENTIAL_CORE.md)
supplies that composition using initialized trusted scratch, an all-time
FAIL penalty and a READY readout protocol. The existence of the
interpreter alone does not establish those spectral or sector claims.

This proof was developed with GPT-6-family assistance in the current
research session. The algebraic opcode blocks and short ordered sweeps
are exhaustively checked by sweep_checks.py; there is no executable
full universal interpreter transition table in this folder.
The independent torus/core worker reviewed the classical/quantum
separation, fixed opcode count, reverse navigation and resource argument.
The trusted-helper interface was added during adversarial cross-review;
helper_checks.py exhausts 29,359 small-alphabet query round trips,
including malformed records of the same types as the helper.
