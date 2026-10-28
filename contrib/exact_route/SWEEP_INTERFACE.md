# A single moving head without quadratic species overhead

## Statement and scope

[proved] An ordered left-to-right sweep of one fixed two-site unitary
can execute one move and one selected two-site gate of a single
interpreter head. A second fixed on-site unitary supplies a local
collision or control update. If the underlying tape alphabet has
dimension \(D\) and there are \(K\) gate/control labels, one explicit
alphabet has dimension

\[
D_{\rm head}=(1+4K)D.
\tag{1}
\]

The head holds only its current tape symbol, not two independent
charged records. Consequently this construction preserves the
\(O(q\log q)\) alphabet bound of [CIRCUIT](CIRCUIT.md) whenever \(K\)
is independent of \(q\). It requires no pair of \(q\)-dimensional
transport lanes.

This note proves the local execution interface. The companion
[ORACLE_SCHEDULER](ORACLE_SCHEDULER.md) supplies an existential
constant-control interpreter construction for the specified circuit
family. A valid initialized input format and its endpoint-readout
schedule remain additional Hamiltonian obligations. No executable full
universal-machine transition table is claimed here.

## 1. Direct sum head states

Let \(\Sigma=\{|a\rangle\}_{a=1}^D\) be the underlying tape alphabet.
Each symbol has species charge \(\chi(a)\). It may include charged
data, neutral program symbols, and the temporary loaded-head symbols
already used by CIRCUIT.

For every \(a\), include the ordinary symbol \(|a\rangle\) and the
head symbols

\[
|H^{0}_{d,k,a}\rangle,\quad |H^{1}_{d,k,a}\rangle,
\qquad d\in\{R,L\},\quad 1\le k\le K.
\tag{2}
\]

The superscript is a fresh/updated phase. Every head symbol in (2)
has the same charge \(\chi(a)\) as its underlying tape symbol.
The direction, phase and instruction are neutral. Equation (1)
counts exactly these states.

Let \(U_{d,k}\) be any specified unitary on \(\Sigma\otimes\Sigma\)
that conserves the total species charge. The finite list may contain
identity moves, neutral gates, species queries, and the controlled
load/exchange operations from CIRCUIT. The matrix of one such gate may
depend on \(q\); the number of opcode families need not.

## 2. One fixed two-site unitary

On the following mutually orthogonal pairs of subspaces define \(G\):

\[
\begin{aligned}
|H^0_{R,k,a},b\rangle
&\longmapsto
\sum_{a',b'}(U_{R,k})_{a'b',ab}
|a',H^1_{R,k,b'}\rangle,\\
|a,H^0_{L,k,b}\rangle
&\longmapsto
\sum_{a',b'}(U_{L,k})_{a'b',ab}
|H^1_{L,k,a'},b'\rangle .
\end{aligned}
\tag{3}
\]

On the reverse subspaces use the adjoint maps; everywhere else use
the identity. Each block has the form

\[
\begin{pmatrix}0&U^\dagger\\U&0\end{pmatrix}.
\tag{4}
\]

Thus \(G^\dagger=G\), \(G^2=I\), and \(G\) is a well-defined unitary
on the *whole* two-site Hilbert space, including configurations with
zero or two interpreter heads. All its matrix elements preserve
the sum of species charges.

The block domains and ranges are disjoint: a fresh right head occurs
only on the left of its forward block, an updated right head only on
the right of its reverse block, and the left-head blocks have the
opposite orientation. The fresh/updated and direction tags prevent
any overlap between these four cases.

## 3. Exactly one move in an ordered sweep

Suppose a chain contains exactly one fresh interpreter head, and its
specified move stays inside the chain. Apply

\[
G_{n-1,n}\cdots G_{2,3}G_{1,2},
\tag{5}
\]

with the rightmost gate acting first.

For a right-moving head at site \(j\), the bond \((j-1,j)\) does
nothing: it has a fresh right head in the wrong slot. The bond
\((j,j+1)\) applies (3), producing an updated right head at \(j+1\).
The subsequent bond \((j+1,j+2)\) also does nothing: the inverse
right rule would require that updated head in the right slot,
while the inverse left rule requires an updated *left* head.
All later bonds miss the head.

For a left-moving head, the move occurs on \((j-1,j)\), leaving
an updated left head at \(j-1\), strictly behind every remaining
bond in the sweep. Therefore there is exactly one move in either
case, with exactly the specified \(U_{d,k}\) on the underlying
ordered pair.

Without the direction tag, an updated right-moving head could be
mistaken for the inverse output of a left move at the next bond.
Without the phase tag, it could keep moving through the entire sweep.
Both tags are needed for this implementation.

## 4. Fixed local collision and interpreter requirements

At minimum, take the on-site gate \(F\) to exchange fresh and updated
phase for every head, fixing ordinary tape symbols. Then (5) followed
by \(F\) on every site executes the same instruction repeatedly, one
move per sweep.

More generally choose a charge-preserving unitary collision \(C\)
from the updated-head subspace to the fresh-head subspace. It may
change instruction, direction, neutral program symbol and neutral
head memory. Extend it to an on-site unitary
\(F=\left(\begin{smallmatrix}0&C\\C^\dagger&0\end{smallmatrix}\right)\)
in fresh/updated ordering. After a valid sweep the head is updated,
so this supplies the desired next control update. Endpoint-specific
collisions may reverse a head before its next move leaves the chain.

A local collision cannot change a charged species into a different
species on its own. Such manipulations must use the two-site
charge-preserving \(U_{d,k}\), as in CIRCUIT's routing construction.
This is a real restriction on a purported controller, not a choice
of notation.

To use the interface as a compiler, the finite controller must:

1. Runs a reversible neutral interpreter for the circuit-generating
   algorithm and its scratch cleanup.
2. Moves to each operand and emits the appropriate opcode without
   merging distinct computation histories.
3. Handles every operand gate using the finite charge-preserving
   \(U_{d,k}\) list, with at most one charged symbol under the head.
4. Finishes within the available number of clock sweeps and gives
   both endpoint records a specified nonzero readout window.

The direction/phase construction proves that doing so will not itself
require a \(q^2\)-dimensional transport site. The existential construction
in ORACLE_SCHEDULER establishes items 1–3 with constant control overhead;
it assumes a valid initialized input and leaves the available time and
readout contract in item 4 to the Hamiltonian construction. In particular,
an absorbing halt state cannot simply be appended to a reversible
unitary machine; idle time must be realized by a reversible clock
or an explicitly accounted cycle.

To preserve the intended port visibility, the execution schedule and
readout weights must also be independent of logical occupation and
the coherent branch register. An oblivious scan can make conditional
gate actions depend on those values while keeping the schedule fixed.
An input-dependent halting time cannot silently be treated as the
same history-state encoding.

## 5. Connection to the existing clock

The 15-letter C8 clock visits every bond once on its outward sweep,
and every site once on return. Decorate those moves with \(G\) and
\(F\), respectively. No disjoint-lane assumption is required:
the executed gate is the *ordered* product (5), followed by the
on-site collisions. The original neutral clock path and legality
projectors are unchanged.

This introduces a factor 15 on top of (1), still \(O(q\log q)\) for
fixed \(K\). Selecting a unique interpreter head and the right program
requires the initialized-input or terminal-verifier construction;
it is not implied by the outer clock's unique head. These are two
different heads on two different registers.

The later [initialized-core composition](../torus_bottleneck/EXISTENTIAL_CORE.md)
uses an initial sweep with no interpreter heads and a single boundary
activation seed, together with initialized trusted scratch and an
all-time acceptance penalty.

The tests in sweep_checks.py construct full two-site matrices for
two quantum opcode blocks, verify unitarity and exact charge
conservation on all two-site states, and test the complete ordered
sweep on every eligible one-head basis input for chains of length
two through four. The opcode list includes a nontrivial complex
particle-number-preserving gate; the tests do not merely follow
classical head positions.
