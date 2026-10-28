# Charge-preserving clock decoration, acceptance, and repeated readout

## Scope and status

This note isolates the precise interfaces between the new
[sparse circuit](../exact_route/CIRCUIT.md), the
[full-sweep clock](../graded_route_a/sweeping_clock_interface.md), and
a prospective uniform cell. The mathematical lemmas below are [proved].
The final compilation implication is conditional on the stated local
moving-head interpreter interface; this note does not substitute a
generic universality assertion for that interface.

Author: Codex, GPT-6 family, 2026-10-27. The proofs were derived directly
and communicated to the coordinating agent. No independent human
review is claimed. No new thermodynamic or microscopic simulation is
claimed here.

## 1. Unitary decoration preserves scalar exclusion bounds

Let `X` be a finite graph of classical control configurations. Each
vertex has a copy of the same payload Hilbert space `D`, and every
oriented edge `x→y` carries a unitary `U_yx` on `D`, with the opposite
orientation carrying its adjoint. For nonnegative edge weights `J_e`
and scalar penalties `v_x≥0`, define the quadratic form

\[
\langle\psi,H\psi\rangle=
\sum_{e=(x,y)}J_e\|\psi_y-U_{yx}\psi_x\|^2+
\sum_x v_x\|\psi_x\|^2.
\tag{1}
\]

Put `r_x=||ψ_x||`. The reverse triangle inequality gives

\[
\langle\psi,H\psi\rangle\ge
\sum_{e=(x,y)}J_e(r_y-r_x)^2+\sum_xv_xr_x^2.
\tag{2}
\]

Since `Σ r_x²=||ψ||²`, the lowest eigenvalue is no smaller than that of
the scalar weighted graph Laplacian plus its scalar potential. This
holds on graphs with branching and cycles, not just legal histories.
Thus a scalar lower bound excluding an unwanted control component
survives arbitrary unitary gates on its edges. It cannot repair a
component that already has an unpenalized scalar zero mode.

For a legal path, write `W_t=U_t⋯U_1`. The block diagonal unitary
`Σ_t |t><t|⊗W_t` conjugates propagation to the scalar path Laplacian
tensored with the payload identity. Its complete spectrum and its
history kernel follow exactly, with no approximation in the gates.

If a local control rewrite acts on a pair of sites and `U_e` acts on
the same pair's work registers, its decorated propagation term is still
two-site. Repeating the same rewrite everywhere preserves one-site
translation invariance. If every `U_e` commutes with every two-site
species count, all propagation terms preserve the species torus
exactly. Diagonal initialization and legality penalties also do so.

This statement avoids a misleading alphabet estimate: a two-site
species SWAP is directly an admissible edge unitary. It need not first
be encoded in a Turing-machine head that simultaneously stores another
arbitrary species label.

## 2. Acceptance need only be controlled inside the initialized space

Let an initialized clock Hamiltonian `H_init≥0` have kernel projection
`P` and gap `γ>0`. Add a positive operator `A`, for example the sum of
penalties for a dedicated FAIL symbol. Suppose a subspace projection
`P_good≤P` obeys

\[
A P_{\rm good}=0,\qquad
PAP\ge a(P-P_{\rm good}),\qquad \|A\|\le B,
\quad a>0.
\tag{3}
\]

Then the full kernel is exactly `P_good`, and

\[
\operatorname{gap}(H_{\rm init}+A)
\ge\frac{a\gamma}{a+B+\gamma}.
\tag{4}
\]

To prove the bound, take a unit vector orthogonal to `P_good` and write
`ψ=x+y`, with `x∈ran(P-P_good)` and `y∈ker P`. Positivity, (3), and
the reverse triangle inequality imply

\[
\langle\psi,(H_{\rm init}+A)\psi\rangle\ge
\gamma\|y\|^2+
\bigl(\sqrt a\|x\|-\sqrt B\|y\|\bigr)_+^2.
\]

Where the bracket is nonnegative, the corresponding two-by-two matrix
has determinant `aγ` and trace `a+B+γ`, giving (4). Where it is negative,
`||y||²≥a/(a+B)`, which gives an even stronger bound. The zero kernel
claim follows by positivity of both summands.

For a clock with `C` equally weighted history vertices, a deterministic
reversible format checker that displays FAIL at least once for every
invalid input supplies `a≥η/C` for FAIL penalty strength `η`. This is
provided the accepted/rejected predicate is diagonal on the initialized
input basis. The inequality for `PAP` then follows by direct averaging
over clock times, including superpositions by linearity.

Crucially, no claim that all pulled-back FAIL projectors commute on
arbitrarily malformed controller states is needed. Such states are
already separated by `H_init`; equation (4) controls their mixing with
the accepted/rejected history sectors.

## 3. Cyclic computation replaces an absorbing halt

A finite unitary cannot send new states into an invariant absorbing
READY subspace. A repeated sweep should instead run a reversible cycle:

1. Compute the classical format predicate and retain its FAIL bit.
2. On valid input, prepare the two-port encoding by a finite unitary
   gate word, possibly approximating the desired encoding.
3. Mark each prepared output record READY and run a nonzero readout
   interval.
4. Clear READY before performing the literal inverse preparation word.
5. Reverse the format and controller computation, returning to the
   starting work configuration and the unchanged logical input.

The machine's direction/phase register distinguishes the forward and
inverse traversals. This describes a reversible cycle, not a fixed
halting state. Its local implementation must still be supplied by the
moving-head interpreter lemma.

The inverse word cancels the actual implemented preparation exactly.
Thus finite gate-approximation error does not accumulate from cycle to
cycle. It changes the prepared isometry `V` to a fixed isometry `V'`.
Each full cycle returns the same logical input, and each active interval
has the same readout matrix element. A terminal clock cut partway through
a cycle is harmless if READY is zero at every intermediate stage whose
readout has not yet been proved correct.

A terminal-only FAIL penalty is insufficient for this protocol, since
the inverse computation eventually clears FAIL. Penalize the dedicated
FAIL symbol whenever it occurs. The correct input never displays it;
every invalid format displays it in each completed checking cycle.
Equation (4) then applies. The actual initialization must pin neutral
ancillas and controller state before the first check; the verifier
cannot certify its own initially untrusted scratch bits.

## 4. An endpoint shuttle keeps the control graph unchanged

The following protocol addresses a boundary whose clock already has
repeated two-site turnaround rewrites but whose bracket is not part
of the computation tape. Attach a work record to the bracket and an
equal-dimensional record to its adjacent interior site. Decorate every
selected boundary turnaround by a SWAP of the entire two records.
The records include their neutral READY bit. This is a species-preserving
unitary and changes neither the neutral clock graph nor its penalties.

Keep this interior boundary site as a shuttle buffer. Perform the
nontrivial encoding in the deeper bulk. Both shuttle records initially
contain fixed dummy states with READY=0, so their automatic interchange
has no effect on the computation. When an output has been prepared,
route its whole record into the shuttle buffer and set READY only after
the last operation that affects its desired readout. Subsequent work
elsewhere commutes with that record's observable. The automatic SWAP
now alternates the prepared READY record with its unready dummy.

Whenever the bracket's READY flag is one, its species operator acts on
the same prepared output record. Whenever the flag is zero, the readout
vanishes. Thus on history states

\[
V_h^\dagger(\Pi_{\rm READY}F_{\rm bracket})V_h
=p_{\rm end}(V')^\dagger F_{\rm output}V',
\qquad p_{\rm end}>0,
\tag{5}
\]

as soon as the prepared record has visited the bracket. Left and right
weights need not agree, and their active time sets need not coincide.
Two extra dummy species records change the required cell length, not
the asymptotic onsite alphabet. A constant number of flag and phase
bits changes that alphabet only by a constant factor.

In a cyclic implementation, clearing READY and recovering each record
must occur at a known shuttle phase before its inverse preparation.
This is a finite scheduling obligation, not an absorbing-state rule.
The full-sweep clock has a regular, known turnaround schedule. An
interpreter must explicitly synchronize its transfer/recovery actions
to that schedule; equation (5) is conditional on that synchronization.

## 5. Polynomial workspace and an exponentially long clock are compatible

The sparse circuit's branch labels, graph enumeration, one-occupation
flip/parity scan, species-word routing, and reversible inverse are
uniform finite algorithms. If their chosen finite-precision circuit
has time and workspace polynomial in `p` and its precision-bit count,
then one may choose a segment length `n(p)=C p^K` with fixed constants
large enough for all work registers. The neutral controller obtains
`p` by a specified decoding of the length, such as an integer-root
encoding with a validity check. It need not have `p` different internal
states.

For all sufficiently large `p`, the `2^n` sweep count exceeds the time
for at least one full computation/readout/uncomputation cycle. Adding
padding to make this true does not increase the onsite alphabet.
Neutral scratch initialization costs local constraints and space; it
does not create extra free witness qubits. These observations establish
resource self-consistency, not the required polynomial-time circuit
bound or its embedding, which must be supplied separately.

Choosing a finite approximation before the thermodynamic limit is
compatible with exact code degeneracy: a unitary gate word gives exact
zero-energy history states, and an exact classical input predicate gives
an exact accepted dimension. Approximation changes compression matrices,
not those two statements. Exact species conservation must be imposed on
every elementary gate, as in the sparse circuit; approximate symmetry
would not suffice for the graph selection rule.

## 6. What constitutes an existential compiler proof

The following is a precise sufficient interface, with no demand for a
numerically optimized alphabet constant.

- A fixed reversible moving-head interpreter is realized by the ordered
  repeated sweep, with full work dimension `D(q)=O(q polylog q)`, including
  heads, neutral scratch, species storage, and READY/FAIL bits.
- Its length-decoding, input checks, preparation, readout, inverse, and
  controller reset are uniform algorithms with a verified workspace
  bound fitting the chosen `n(p)`.
- Its initialized input space has exactly the intended free occupation
  bits; every other controller/scratch degree of freedom is fixed.
- Its FAIL predicate is correct on this space, its quantum gate word
  preserves that bit, and every invalid input is exposed before the
  outer clock ends.
- At every active endpoint time, including during transfers and reverse
  sweeps, the same desired compressed observable is read out; at all
  other times the endpoint flag is zero.

Given those items, the clock theorem and (4) give a positive gapped core
with exactly the desired code, at local dimension `15D(q)` times any
fixed marker/selector overhead. The asymptotic c–d comparison needs only
that the constant is independent of `q`; an explicit minimal numerical
constant is unnecessary.

An actual proof of the local moving-head update and of the complete
uniform program supplies these items. A citation to universal
computation alone does not: it could hide a `q²` tape/head factor,
uninitialized program sectors, or an invalid absorbing readout state.
The preceding lemmas remove those separate spectral and readout
objections once the interpreter contract has been met.

For comparison, the multi-track legal-path and scalar-illegal-state
structure is also visible in [Gottesman–Irani, Sections 5.2 and
5.6–5.8](https://arxiv.org/html/0905.2419). Their initialized witness
track and bracketed segments support the same interface strategy.
The repository's explicit full-sweep clock is preferable for this
project because its local alphabet and full legal graph are already
written down and can be audited directly.

This compiler result concerns a frozen segment core. A microscopic
mobile-cell transport that intertwines the two neighbouring segment
lengths, the local finite-coupling reduction, and the infrared conclusion
are separate requirements. They are not implied by any circuit or
history-space universality theorem.

## 7. The interface is now instantiated

[EXISTENTIAL_CORE.md](EXISTENTIAL_CORE.md) supplies the initialized
finite-core theorem using the exact sparse isometry, indexed algebraic
branch rotations, and the constant-control oracle scheduler. It uses
all-time FAIL penalties and the direct endpoint READY schedule, so the
optional shuttle construction above is not needed. Its core statement
has exact port compressions and exact `2^N` ground-space dimension.
The mobile wrapper and infrared transfer are separate applications of
that theorem, treated by the route composition.
