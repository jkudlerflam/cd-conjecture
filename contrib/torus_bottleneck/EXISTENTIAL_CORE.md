# An initialized species-preserving core with a superlinear logical payload

## Statement

For every sufficiently large prime `p`, put

\[
q=2(p^2+p+1),\qquad N=(p+1)(p^2+p+1).
\]

There is a computably specified family of finite, marker-bounded,
one-dimensional history Hamiltonians, with the same Hermitian
nearest-neighbour interaction repeated throughout each segment and
fixed local boundary instructions, with the following properties.

1. The full local dimension is at most `C q(1+log₂q)` for a constant
   `C` independent of `p` and segment length.
2. A selected segment length `n(p)=C₀p^K`, with fixed sufficiently large
   integer constants `C₀,K`, has an exactly `2^N`-dimensional ground
   space, ground energy zero, and a strictly positive finite-segment gap.
3. Every local term preserves every species count exactly. On the code,
   the species torus acts with the incidence-graph logical charges
   `w_(αβ)=e_(1α)-e_(0β)`. Physical upper-species parity restricts to
   logical fermion parity.
4. Local, flagged endpoint operators compress exactly to positive
   computable time weights times the port matrices of the sparse
   two-port isometry `V_p`. In particular, the first-order intercell
   hopping has exactly `N` independently weighted channels.

The constant `C` is not optimized or numerically evaluated. The
construction specifies a finite algorithm for its local matrices;
an expanded list of their entries is not needed for this existential
dimension bound.

**Status: [proved, by composition of the interfaces cited below].** This
is a finite-core theorem. The mobile wrapper, local finite-coupling
reduction and infrared identification lie outside its statement and
are treated in the companion route audit. The original conjecture is
about the full periodic chain, so this theorem alone is not a
microscopic counterexample to c–d.

Author: Codex, GPT-6 family, 2026-10-27, with the parallel project agents.
The clock, circuit, interpreter and acceptance arguments were exchanged
and checked between those agents. No independent human review is
claimed. The construction is analytic; no numerical matrix diagonalization
at these very large selected lengths is claimed.

## 1. Ingredients and their precise use

The proof uses the following results, with distinct responsibilities.

- [The sparse two-port isometry](../exact_route/README.md) specifies the
  finite code and its exact charged port matrices.
- [The nearest-neighbour circuit](../exact_route/CIRCUIT.md) prepares
  that isometry by coherent branch preparation, a reversible basis-word
  map, and species-preserving routing. Its work alphabet has dimension
  `D(q)=q+3+2(q+1)ceil(log₂q)` and its time and space bounds are polynomial
  in `N+q`.
- [Exact branch amplitudes](../exact_route/EXACT_AMPLITUDES.md) enlarge
  the index range to `m=O(log q)` and replace approximate synthesis by
  a fixed menu of algebraic rotation matrices. The resulting exact
  alphabet dimension is `D_exact=q+3+2(q+1)m=O(q log q)`.
- [The oracle interpreter](../exact_route/ORACLE_SCHEDULER.md) uses a
  constant number of neutral tracks and a fixed finite-state controller
  to schedule a constant-cardinality menu of such local gates. It
  does not store a second charged species in its control state.
- [The moving-head sweep interface](../exact_route/SWEEP_INTERFACE.md)
  turns each interpreter move into one fixed ordered nearest-neighbour
  sweep followed by one fixed on-site collision. Its overhead is
  `(1+4K_h)D` for a fixed number `K_h` of control/opcode states.
- [The fifteen-letter clock](../graded_route_a/sweeping_clock_interface.md)
  realizes exactly `2^n` such sweeps on a segment of length `n`, has a
  complete known legal path, and energetically excludes its other clock
  sectors.
- [The acceptance bound](COMPILER_INTERFACE.md#2-acceptance-need-only-be-controlled-inside-the-initialized-space)
  separates initialized but incorrectly formatted quantum data without
  assuming that acceptance projectors commute on arbitrary malformed
  controller inputs.

The oracle-controller result is used as a constructive finite-machine
compilation lemma, not as an assumption that a universal machine
somehow supplies initialized inputs or valid readout. Those two steps
are implemented below.

## 2. Fixed geometry and initialized workspace

Before the outer clock and interpreter-head overlay, a site has three
distinct work factors: constant-dimensional neutral classical metadata
`Γ`, a fixed number of neutral quantum flag/scratch qubits `Q`, and the
single charged/neutral circuit record `Σ_q`. FAIL, READY and trusted
quantum scratch belong to `Q`; static marks belong to `Γ`. The data
SWAP and endpoint transport operations act on `Σ_q` only, leaving
physical flags and static marks at their original sites. The finite
opcode menu includes gates on `Q` and controlled queries between
`Σ_q` and `Q`. These constant neutral factors multiply the alphabet
by a constant, and are distinct from the classical scheduler state.
Only `Σ_q` contains a species-valued register.

Attach a neutral static mark with values LEFT, BULK, RIGHT. Boundary
pins require LEFT at the first site and RIGHT at the last. Allowed
interior pairs are

\[
({\rm LEFT},{\rm BULK}),\quad({\rm BULK},{\rm BULK}),\quad
({\rm BULK},{\rm RIGHT}),\quad({\rm LEFT},{\rm RIGHT}).
\]

All other pairs have a positive diagonal penalty. For `n≥2`, the
unique zero-penalty mark string is `LEFT BULK^(n-2) RIGHT`. These marks
are never moved by a payload gate. The interpreter can read them to
recognize the endpoints; the factor three is included in the alphabet
constant.

During the clock's initial I sweep, local pins initialize the neutral
classical controller and its program tracks, dedicated quantum scratch
tracks, READY bits and FAIL bit. They require the absence of every
interpreter-head overlay at every site. At the unique initial
right-boundary I-clock vertex, one additional positive boundary pin
fixes its underlying `Σ_q` record to the known neutral query helper
`H_(0,0)`. This is an initialization pin, not a permanent restriction
on the right output port. All I-sweep data gates are identity, so it
commutes with the other pins in the initial clock gauge.
At the end of this sweep,
the clock's right-boundary activation applies a seed unitary `B`:
for every underlying record `σ`, exchange ordinary `σ` with the
start-head overlay of `σ`, and fix the remaining orthogonal states.
This creates exactly one interpreter head without a position-dependent
bulk pin. It preserves charge because the overlay carries the charge
of `σ`.

Every other underlying quantum-data symbol can initially range over
the entire work alphabet `Σ`, including its unloaded and loaded
record types. Their position-dependent format will be verified
computationally. The right-end helper is an additional fixed reference
input, so it contributes no free logical degree of freedom. Classical scheduling tracks are distinct from coherent
quantum scratch, FAIL and READY: scheduler movement never branches on
any of those quantum values.

This distinction avoids a circular initialization claim. The format
verifier has trusted, zero-initialized neutral scratch from local pins.
It is not asked to certify an untrusted work register using that same
register. The initial I sweep applies no quantum-data gate, so these
pins are ordinary commuting initialization constraints in the clock
gauge. Mark and clock legality penalties are also positive.

Let `H_init` denote this positive initialized history Hamiltonian, `P`
its exact kernel projector, and `γ_n>0` the gap supplied by the clock
interface. At this stage `P` includes all quantum-data inputs compatible
with the deliberately permissive initial pins. It does not yet have
the desired code dimension.

## 3. The uniform program and its resources

The controller first counts the distance between the immutable endpoint
marks and decodes `p` from the designated length `n=C₀p^K`. Integer-root
decoding and the prime-field incidence test are fixed algorithms; their
descriptions do not grow with `p`. The prime need not be stored in an
onsite control-state label.

The reversible length count does not require storing a polynomial-time
history in only `n` sites. Ordinary endpoint counting takes polynomial
time and logarithmic working space. A reversible implementation uses
`O(log²n)` space by the standard time/space trade-off, or can use a
visited-marker track and an explicit reversible counter. The relevant
general result preserves the original input while simulating time `T`,
space `S` in reversible space `O(S log T)` and polynomial time.
[Bennett, 1989, primary-source abstract](https://research.ibm.com/publications/timespace-trade-offs-for-reversible-computation).

The program next generates the incidence graph, the reference word,
the circuit instructions and the classical format checks. The sparse
circuit uses polynomial space and time in `p`; its exact rotation menu
has only `O(log q)` indexed blocks. Reversible computation of its classical schedule,
including the needed input/output addresses, consequently takes at
most `p^A` times a fixed power of `log p` in workspace for some fixed
`A`. Choose `K>A+1` and enlarge fixed `C₀` if needed. The workspace then
fits in the segment with a constant number of neutral tracks per site.
Physical head travel adds polynomial time and counter bits; it does
not require a record of every travel step.

Operations which traverse the full physical length use compact
reversible loops, not a stored `n`-instruction or `n log n`-bit gate
list. In particular, an endpoint SWAP walk carries one underlying
`Σ` record under the interpreter head, applies the same adjacent SWAP
opcode at every step, and stops at the immutable endpoint mark.
Other records undergo a fixed permutation; the inverse walk restores
them after readout. No second charged species register is introduced.
The two output walks can be chosen successively without crossing the
already placed opposite endpoint. They also move the initial right-end
query helper into the auxiliary interior. Their permutations merely
relocate unobserved auxiliary records, so they do not alter either
port matrix.
The complete-input verifier below likewise uses compact scans and a
short counter. Thus increasing `K` does not create a new workspace
requirement proportional to a higher power of `n`.

The resulting computation cycle has time polynomial in `n` and `p`.
For all sufficiently large `p`, at least one entire cycle fits before
the outer clock completes its `2^n` sweeps. This is the only lower
bound on the number of completed cycles needed here.

The q-dependent species queries, indexed exact rotations, load
operations and adjacent exchanges form a fixed-cardinality gate menu: the species and bit-index
dependence is in the matrices and in the single work record, not in
a second q-valued controller. A constant number of neutral control
tracks tensored with this record therefore costs `C₁D(q)`. The head
and outer clock multiply this by constants. All markers, program,
scratch and readout flags are included, giving

\[
d_{\rm core}\le Cq(1+\log_2q).
\tag{1}
\]

No bound on time, cell length or coefficient magnitude is asserted to
be uniform in `p`. Only the alphabet constant is uniform.

## 4. Exactly the desired accepted input space

The classical format predicate requires:

1. One two-valued occupation site for every incidence edge, in exactly
   its two specified species states.
2. The prescribed charged reference records needed by the two-stage
   sparse isometry, in their prescribed places.
3. Every other quantum-data/work slot in its prescribed neutral or
   charged reference state.

A TYPE query checks the complete underlying alphabet before species
queries are interpreted: an unexpected loaded, unloaded or other
record type must fail, even if some of its species bits agree with a
valid record. TYPE and label queries have constant-cardinality opcode
families, with their finite matrices defined on every alphabet state.
The indexed species query has a trusted target: transport the pinned
`H_(0,0)` helper by `Σ`-only SWAPs until it lies immediately to the
right of the record being tested. SHIFT selects a known index `j`;
QUERY writes the corresponding species bit into the helper's neutral
bit, and an exact copy gate transfers it to trusted `Q` scratch.
Undo QUERY and SHIFT, then reverse the SWAP walk. All input records
and the helper return to their original positions. TYPE handles
malformed records, on which the species query is defined by its
unitary extension. The helper's trajectory is a known tensor-factor
permutation even when an untrusted record happens to have H or T type.
Thus no untrusted work slot is borrowed as a fresh query register,
and no additional `O(log q)` tensor factor is needed.
The predicate is computed coherently by these diagonal queries and
reversible neutral Boolean gates. Its trajectory is independent of
the values it reads: it executes a fixed classical schedule, rather
than branching its head motion on logical occupations.

The scan checks every untrusted `Σ` record, including the unused
padding; the trusted boundary helper is already fixed by its initial
pin. This does not require a Boolean history of length `n log n`.
For each position `i`, compute its violation predicate into trusted
quantum scratch, conditionally increment a trusted quantum counter
with `ceil(log₂(n+1))` bits, and uncompute that local predicate.
The reference/type predicate is a function of classical `p,i` and the
queried quantum record. After the full scan, compute whether the
counter is nonzero, copy that one bit to FAIL, undo the comparison,
and reverse the counter scan. This clears all counting scratch.
The counter cannot overflow because there are at most `n` violations.
Compact reversible position loops need only `O(log n)` classical bits;
quantum queries and counter arithmetic use `O(poly(p)+log n)` work,
with a constant number of trusted neutral tracks per physical site.
The polynomial preparation program is kept only for the `poly(p)`
work region. Full-length scan and travel loops have no stored
step-by-step gate list.

Use a single exact CNOT opcode to copy the finished violation bit to
FAIL. FAIL is not an intermediate target in a synthesized Boolean
subroutine: it is touched only by that copy and its later inverse.
Thus valid inputs never have even a transient FAIL component at an
elementary history vertex. Keep that bit untouched by the subsequent
encoding circuit. The scheduler executes the same
complete preparation word on every input, including inputs with FAIL
set. There is no failure-dependent halt or branch. Each primitive has
a charge-preserving unitary extension on every malformed record type;
that extension cannot remove the already recorded violation.

The accepted input subspace is exactly the tensor product of the `N`
two-valued occupation spaces with all reference states. Its dimension
is therefore exactly `2^N`. In particular, unused neutral tape cells
do not remain as free low-energy witnesses.

Add an onsite penalty `η` for every occurrence of the dedicated FAIL
symbol. This is the same local term everywhere; its total norm is at
most `B_n=ηn`. Every incorrect initialized basis word produces FAIL at
least once in the first completed verification. On the normalized
history space,

\[
PAP\ge a_n(P-P_{\rm good}),\qquad
a_n\ge\frac\eta{C_n},\qquad
C_n=(2n+4)2^n+n-4.
\tag{2}
\]

The checker is a reversible classical circuit before preparation, and
FAIL is preserved during preparation, so its restricted history
expectation is diagonal in the initialized input basis. This proves
(2) also for arbitrary coherent superpositions of accepted and rejected
inputs. The invalid control/scratch sectors excluded by `H_init` need
no further commutativity claim.

The combined core has exact zero kernel `P_good` and gap

\[
\Delta_n\ge
\frac{a_n\gamma_n}{a_n+B_n+\gamma_n}>0.
\tag{3}
\]

Additional positive mark and clock penalties retain their stated
lower bounds. Equation (3) is intentionally conservative. It suffices
that the gap is positive for each fixed selected segment.

## 5. Repeated preparation, clean readout, and exact reset

After checking the format, execute the exact charge-preserving sparse
encoding, using the indexed algebraic rotations. Route its left and
right output records to the corresponding physical endpoints while
READY is zero. The whole preparation is a unitary word `U_p` on the
larger work space. Restricted to valid inputs, it defines exactly the
sparse isometry `V_p`, including its clean auxiliary reference records.

Use the following data-independent schedule.

1. After every operation on the quantum encoding has finished, move
   the interpreter head neutrally to the left endpoint and set its
   READY bit.
2. Travel with identity quantum-data gates to the right endpoint and
   set its READY bit.
3. Travel through the interior and back, or execute another explicitly
   reversible neutral control loop. This supplies a nonzero interval
   of clean endpoint readout.
4. Clear the right READY bit, travel neutrally to the left, and clear
   the left READY bit.
5. Execute the literal inverse word `U_p†`, then undo the format
   computation and controller bookkeeping, returning to the initial
   logical input and the initial neutral work configuration.

The interpreter phase distinguishes preparation, neutral readout,
and inverse preparation. Its finite reversible cycle has no absorbing
halt state. The controller then repeats it. FAIL is penalized whenever
it appears, so resetting it during uncomputation cannot hide an invalid
input from (2).

Define each local readout operator to include both its READY projector
and the projector onto the absence of the interpreter-head overlay.
When the head is at an endpoint, the readout vanishes. At any head-absent
active time, the same prepared endpoint record is read out. The
intervening operations act only on neutral controller tracks and have
identity action on the quantum data. Thus their conjugation cannot
change the desired endpoint matrix element.

This argument applies at every elementary clock vertex, including
within a sweep. Setting and clearing READY while the head is present
cannot contribute an intermediate unwanted readout. The two endpoint
flags need not become active simultaneously; after preparation,
neutral motion leaves each endpoint observable unchanged.

Every complete cycle restores the logical input exactly because the
actual gate word is followed by its actual inverse. Truncation of the
final cycle by the finite outer clock introduces no contamination:
each included time is either inactive or has the same clean readout.
There is at least one active time at each end because a full cycle fits.

If `S_L,S_R` are the resulting active time sets, their sizes and the
history path length are computable from the classical controller
schedule. Write `p_L=|S_L|/C_n>0`, `p_R=|S_R|/C_n>0`. For every graph
edge `a`, the full history isometry `G_p` obeys

\[
\begin{aligned}
G_p^\dagger\widehat F_{L,a}G_p
 &=p_LV_p^\dagger F_{L,a}V_p,\\
G_p^\dagger\widehat F_{R,a}G_p
 &=p_RV_p^\dagger F_{R,a}V_p.
\end{aligned}
\tag{4}
\]

The weights need not agree. Their smallness affects the effective
energy scale, not the number of logical channels.

## 6. Exact symmetry, port amplitudes, and the order of limits

All species queries are diagonal on charged labels. Load/exchange
gates conserve the combined charges of their two records. Neutral
gates, control collisions, static marks, initialization and FAIL
penalties carry no species charge. The moving-head overlay carries
exactly the charge of its underlying tape record. Thus every term of
the assembled Hamiltonian commutes with all total species counts.

The charge action on the initialized occupation code is exact. The
fixed reference word contributes a constant, and its upper-species
count is even in the two-stage sparse construction. Logical parity is
therefore exactly physical upper-species parity throughout the history.

Let `u=p²+p+1`, `b²=1/(2N)` and `ε=1/(2sqrt(uN))`. The ideal sparse
isometry has the port matrices

\[
V_p^\dagger F_{L,a}V_p=\epsilon f_a^\dagger,\qquad
V_p^\dagger F_{R,a}V_p=
\epsilon(1-b^2)f_a^\dagger P_{\rm log}.
\tag{5}
\]

Together with (4), these identities give the exact first-order
hopping scale

\[
\lambda\,p_Lp_R\,
\frac{1-1/(2N)}{4uN}.
\tag{6}
\]

There is no gate-approximation error in this statement. The coefficient
matrices depend on `p`, as the Hamiltonian family is allowed to do;
they are computably specified algebraic matrices. No matrix depends
on the total ring length.

An optional finite-gate-set variant may replace `V_p` by `V'_p` with
`||V'_p-V_p||≤δ_p`. Each unit-norm endpoint matrix then changes by
at most `2p_end δ_p`. Time weights multiply both the signal and this
error, so long histories do not require exponentially accurate gate
synthesis. All target graph and port parameters are polynomial in
`p`; inverse-polynomial accuracy suffices for any specified
inverse-polynomial first-order tolerance. The exact construction
avoids this additional choice entirely.

For each fixed `p`, first fix the exact gate word, the finite
controller, the segment length and all core couplings. Then choose a
strictly positive intercell coupling small compared with the finite
core gap and any selector gaps. Only after these choices is the
thermodynamic limit contemplated. No parameter is adjusted with the
total ring length.

## 7. Consequence and remaining scope

The finite initialized high-channel core can be constructed with

\[
\frac{N}{d_{\rm core}}
\ge\frac{p+1}{2C[1+\log_2(2(p^2+p+1))]}
\longrightarrow\infty.
\]

This removes the existential alphabet/controller obstruction for the
species-preserving first-order route. It does not establish unbounded
central charge at a fixed local dimension; `d_core` grows with `p`.

Applying this core to c–d additionally uses spontaneous-cell placement
with all physical states counted, the C8 mobile transport wrapper,
a local finite-coupling reduction, and identification of the periodic
chain with its stable Gaussian theory. Those are separate composition
steps; no assertion that they fail follows from the scope of this
finite-core theorem. The companion route audit treats them explicitly.
The sparse graph supplies the protected infrared target, including
the mobile-mode corollary in the companion README.
