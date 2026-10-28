# Exact full-sweep interface for a prospective interpreter

## Claim and provenance

The existing 15-letter extension in
`codex/work_C8/sweeping_clock_checks.py` remedies the sparse schedule
of the nine-letter clock. This note states and proves its general
interface, including an extension to computational input checks.
It is not a completed interpreter or a new discovery of that clock.

[proved] At length n the complete legal clock graph is one path of
C_n=(2n+4)2^n+n-4 vertices. It executes exactly 2^n complete outward
and return sweeps. Arbitrary fixed two-site work unitaries can decorate
the outward sweep; disjoint transport lanes are unnecessary for this
claim. A repeated sweep therefore realizes an ordered local circuit.
Commuting pulled-back initialization/acceptance checks give an exact
complete code and a quantitative gap.

Author: OpenAI Codex, GPT-6 family, 2026-10-28. The C8 code and PROGRESS
were read. The graph and general circuit/acceptance statements below
were independently checked. No human review is claimed.

## 1. Clock alphabet, legality and transition rules

Use the fifteen letters

    L, 0, 1, C0, C1, R0, R1, I, Z, A0, A1, D0, D1, E0, E1.

The counter heads are C0,C1,R0,R1 and sweep heads D0,D1,E0,E1.
Allowed internal pairs are LL, IZ, ZZ, Lh for counter heads or I,
hb for any counter or sweep head and b in {0,1}, arbitrary ordinary
bit pairs, arbitrary A-bit pairs, and (Aa,h) for sweep heads h.
At the left boundary allow L,A0,A1 or any head, including I.
At the right allow 0,1,Z or any head. Every forbidden pattern has a
positive diagonal penalty K.

The zero-penalty words are exactly

    L^(i-1) h bit^(n-i), h in {C0,C1,R0,R1};
    L^(i-1) I Z^(n-i);
    A-bit^(i-1) h bit^(n-i), h in {D0,D1,E0,E1}.         (1)

There is one head: the prefixes cannot terminate at the right boundary,
and after a head only its appropriate suffix is allowed.

Retain the binary carry/return rules

    C1 b -> L Cb, C0 -> R1, L Rb -> R0 b,
    IZ -> LI, right-boundary I -> R0,

and replace the left restart by a full sweep:

    left-boundary Rb -> Db;
    Da b -> Aa Db;
    right-boundary Db -> Eb;
    Aa Eb -> Ea b;
    left-boundary Eb -> Cb.                             (2)

Here a,b are ordinary bits. The sweep preserves the counter value:
D moves right, changing the already visited digits to A-digits; E
moves left, restoring ordinary digits. Every bit string has exactly
one such sweep interlude before its increment.

Deleting these interludes recovers the proved nine-letter counter path.
Conversely every sweep word in (1) specifies one counter value and
one unique sweep location/direction, so lies on its interlude. There
are no other legal components. Counting (1) gives

    C_n = 4(2^n-1)+n+2n*2^n = (2n+4)2^n+n-4.            (3)

Each rule preserves every legality projector, even with illegal
neighbours: its two sides have identical allowed predecessor and
successor sets, and their internal pairs are both allowed. The boundary
rules preserve the corresponding endpoint conditions. The positive
clock Hamiltonian thus separates illegal sectors by at least K.

## 2. An arbitrary ordered sweep circuit

Attach work dimension w to each clock letter. Decorate outward rule
Da b -> Aa Db by a fixed arbitrary unitary G on its two work sites.
Decorate each return rule Aa Eb -> Ea b by I on the left work site
and a fixed F on the departed right site. Decorate the left finish
Eb -> Cb by F on site 1. All other counter transitions may be identity.
The activation I->R0 can prepare a fixed right-end seed by a unitary B.

Every complete sweep then implements exactly

    U_n = (F_1 F_2 ... F_n) G_(n-1,n) ... G_(1,2).       (4)

The two-site gates are in the displayed order and need not commute.
Special endpoint F operators may be inserted using the existing
boundary rules. After r sweep finishes, the work state is U_n^r B psi.
At terminal overflow it is U_n^(2^n) B psi.

This is stronger than the restricted partitioned-CA example in the C8
test, but follows from the same literal transition sequence. Requiring
G to swap disjoint lanes is only one way to interpret (4) as a standard
partitioned cellular automaton. It is not a condition for history-state
legality, positivity, completeness or the spectral gap.

This distinction matters for alphabet accounting. Two independent
q-dimensional charged lanes cost q^2 states per site. The clock itself
only multiplies the chosen work dimension by fifteen. A work alphabet
built as a direct sum of neutral, unloaded and one-species-loaded head
states could cost O(q log q), if its actual reversible rules fit (4).
That efficient interpreter is a separate obligation.

## 3. Positive propagation and computational acceptance checks

For each decorated transition a->b with unitary U use

    (J/2)[|a><a|+|b><b|-|b><a| tensor U-|a><b| tensor U^dagger].

Let W_t be the cumulative circuit on the unique path. Add local positive
pins eta_i |clock-pattern_i><clock-pattern_i| tensor Q_i. A local
clock pattern may occur at multiple times; treat each occurrence
separately. Suppose its pulled-back projectors

    P_i = W_(t_i)^dagger Q_i W_(t_i)                     (5)

commute, and let A be their common kernel. For eta_i>=eta_min>0 the
complete ground code is

    V psi = C_n^(-1/2) sum_t |t> W_t psi, psi in A,       (6)

and the gap is at least

    min{K,J[1-cos(pi/C_n)],eta_min/(2C_n),J/(4C_n^2)}.    (7)

The proof is the same controlled-circuit conjugation and scalar
path-potential estimate as for initialization. In a common eigenspace
outside A at least one vertex has a pin of strength eta_min.

This includes diagonal local format checks during the initial I sweep.
It also includes a classical reversible format verifier whose terminal
FAIL projector pulls back to a diagonal predicate of the input. If the
initial checks are diagonal in that basis, the pulled-back checks
commute. Thus the accepted input need not be recognizable by static
nearest-neighbour tags alone. It can be selected by a verified classical
computation, without changing the history-space proof.

This does not supply the verifier, its runtime or its local gate table.
It states precisely how one can use it without assuming a preferred
program or padding sector. If quantum branch preparation follows the
check, a saved classical acceptance flag can be kept untouched; one
must show that its final projector still pulls back to that predicate.

## 4. Readiness at both physical endpoints

The code compression of an endpoint operator is the time average

    V^dagger F V = (1/C_n) sum_t P_A W_t^dagger F_t W_t P_A. (8)

A sufficient clean-readout contract is: a local endpoint flag gives an
operator F that vanishes at every inactive time, and at every active
time its compressed action is the same desired logical A, independent
of the input. Then (8) is |S| A/C_n, with S the active times. Left and
right sets S may differ; their nonzero known weights need not be equal.

Terminal clock position alone is insufficient: only the right endpoint
uniquely sees the final C1. The left endpoint must acquire its readiness
from the work computation, or use a separately proved readout protocol.
Any claimed interpreter must give both active-time sets and show that
no intermediate computation contaminates the port compression.

Reversible idling cannot use an absorbing READY state of one repeated
unitary. Such a finite invariant subspace has no incoming transition
from its complement. Idling can instead move a timer or other garbage
while preserving the output, with a proved no-wraparound interval.
The existing exponentially long clock allows such time, but does not
construct the timer.

## 5. Error bounds and what remains unproved

If the actual work unitaries differ from intended ones by operator
norm errors delta_s, telescoping gives

    ||W_t-W_t^ideal|| <= sum_(s<=t) delta_s.

For a norm-bounded readout, its compressed error is consequently at most
2||F|| max_t sum_(s<=t)delta_s, before any separate input-subspace error.
Repeated use of a fixed imperfect rule must be counted repeatedly.
For any exact unitary decorations the legal propagation spectrum is
the exact path-Laplacian spectrum. The gap bound (7) holds whenever
the stated commuting-pin hypotheses hold. An approximate algorithm
must preserve those acceptance hypotheses or receive a separate
analysis of the resulting code splitting and subspace error.

A universal reversible machine is not automatically a solution: its
embedding into the concrete sweep, initialization from blank and length,
charged alphabet count, input acceptance, clean two-endpoint readouts,
and bounded reversible idling still need an explicit proof. Neither
(4) nor a generic universality assertion supplies those components.

## Reproducing

`sweeping_interface_checks.py` independently enumerates the local rules
on complete clock spaces through n=4, verifies legality-defect
preservation, and tests (4) for arbitrary noncommuting two-site gates
on a three-qubit work chain. NumPy is its only external dependency.
