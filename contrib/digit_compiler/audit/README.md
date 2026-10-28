# Independent audit and an explicit reversible finite-table lift

**Later result in this contribution:** the
[single-physical-head reversible compiler](REVERSIBLE_HEAD.md) closes the
spatial-reversibility gap for a supplied finite source table, using explicit
marked sweeps and a finite history alphabet. It has been composed with the
actual lowered ALU and checked through 7,407,784 forward/inverse physical
microsteps. The multitape lift below remains an independently checked
intermediate construction; it is not the route used for that spatial result.
The root [gate completion](../reversible_gates.py) additionally supplies
stationary-move removal and full-space local permutations. The independent
[gate audit](gate_checks.py) checks every local input of a finite fixture and
five autonomous erasure cycles under a fixed repeated gate word.

## Scope and claim

**[proved]** Every finite deterministic multitape transition table has the
explicit reversible lift in [`logged_lift.py`](logged_lift.py). It adds one
history tape with a **fixed finite alphabet**, not a saved Python integer,
activation record, or variable-sized object. One source transition is four
permutation primitives. The lift has an exact inverse on **all** configurations
of its declared alphabets, including incorrectly initialized history tapes.
Both the ordinary multitape interface and the more economical single-active-
tape interface used by the new arithmetic and memory routines are implemented.

**[checked]** [`checks.py`](checks.py) independently applies the lift to the
actual arithmetic and serial-memory tables. Its results are in
[`results.json`](results.json). This checks a mathematical interface and a
concrete reversible implementation; it does not establish a complete
initialized one-head quantum history Hamiltonian.

The older four multitape primitives inspect symbols at their own, possibly widely
separated tape heads. They are finite multitape operations. **They are not
nearest-neighbor interactions along the physical chain.** The newer
[one-head compiler](REVERSIBLE_HEAD.md) supplies the physical reversible
alternative. Complete-generator compilation, initialized static Hamiltonian
sequencing and illegal-sector checks remain separate obligations.

## Complete permutation construction

Let the deterministic source machine have finite state set `Q`, work alphabets
`Gamma_1,...,Gamma_k`, and `m` defined rules. Give every rule a different fixed
label `t`. A general rule reads

```
(q, a_1,...,a_k) -> (r, b_1,...,b_k; d_1,...,d_k),
```

where every displacement is in `{-1,0,1}`. The new control set is a disjoint
union `Q union {P_t : t is a rule}`. The history alphabet is
`{blank} union {t : t is a rule}`. In particular it has `m+1` symbols independent
of the prime, word width, running time, or tape length.

Each source step is the composition `H E M C`, in the following order.

1. **C, collision:** swap
   `(q, a_1,...,a_k, history blank)` with
   `(P_t, b_1,...,b_k, history t)` for every rule `t`. Fix every other input.
2. **M, movement:** when the control is `P_t`, move work head `i` by `d_i`.
   Otherwise leave all work heads fixed. Preserve control, tapes and history.
3. **E, control release:** when the scanned history symbol is `t`, swap the
   controls `P_t` and `r`, the target of rule `t`. Leave all tape symbols fixed.
   When history is blank, fix the control.
4. **H, history shift:** move the history head right by one cell.

`C` is a product of disjoint transpositions: determinism makes its source
tuples distinct; its targets are distinct because of `P_t`; sources and
targets belong to disjoint control sectors. `E` is one transposition in each
history-symbol sector. `M` is a translation conditional on an unchanged
control. `H` is a translation. Thus **every primitive is a permutation on its
entire space**. No promise on input tape contents is needed for invertibility.

The inverse is `C M^-1 E H^-1`. In the initialized computation the history
head starts at the first unused blank and the control belongs to `Q`.
`C` writes the next rule record and the new work symbols; `M` makes the source
head moves; `E` changes `P_t` into the source target `r`; `H` restores the
next-unused-blank invariant. Induction proves exact agreement with every
defined source step and exact restoration after inverse replay.

This proof uses bi-infinite tapes, or tapes with enough margin that the
promised computation never crosses their ends. A finite physical cell needs
explicit end handling and a capacity proof. The history tape stores **one
symbol per digit-level transition**, so its length is governed by digit
running time, not the earlier word-machine instruction count.

For an asynchronous source state, the active work tape `i(q)` is fixed by
finite control. Each rule then reads and writes only that tape. The same
construction works with all other symbols as spectators; no Cartesian
expansion of the unused alphabets is needed. The `SerialLoggedLift` class
implements precisely these spectator identities. Its list of supported
collision transpositions completely specifies `C`; the default is identity.

## What the checks establish

- A destructive two-tape source fixture erases arbitrary source bits. Its
  complete collision and control-release tables are enumerated: **2,912**
  and **182** inputs respectively. Both are bijective involutions.
- **10,000** randomly selected malformed complete configurations, including
  intermediate controls, overwritten history and shifted heads, are exactly
  restored by forward followed by inverse execution.
- Every supported transposition of the actual ALU and memory tables is
  checked in both directions with each same-symbol spectator filling. The
  independent proof above covers arbitrary spectator combinations; this test
  does not enumerate their exponentially large Cartesian product.
- **60** complete ALU executions, covering all 15 operations, match host
  integer reference answers. Their **9,746** source transitions are executed
  through the lift and reversed to their exact initial tapes and control.
- **14** complete memory executions, including variable-length replacement
  words and empty arrays where allowed, match the original table interpreter.
  All **871** transitions are reversed exactly.

These tests use host integers only to prepare examples, evaluate an
independent expected answer, count transitions, and represent tape positions.
The simulated machine's work values are digit strings. Its runtime operations
are table lookups, local symbol writes, and unit head moves. The Python
dictionary representing a tape is a simulator data structure, not a physical
random-access instruction.

## Interface audit and the polynomial-space issue

The older generator's word-machine is genuinely fixed: 3,115 instruction
addresses and an acyclic call graph. Its journal did **not** itself supply
finite local reversal rules. The construction above closes that specific
gap **for a supplied finite digit-level table**. It must not be applied to the
old word-machine instruction count as though arbitrary integer arithmetic
were a single source transition.

The prior resource notation is

```
A = (N+4q)^2 max(1, ceil(log2 q), ceil(log2 max(u,N))),
B = 128 A,
T <= 18463 B^6,
b = 2 A + ceil(log2(T+1)) + 64.
```

The `2 A` term is necessary. `flag_value` packs `O(q log q)` bits into one
integer; its width cannot be bounded just by `O(log p)`. Source inspection
also shows that:

- `DIV` and `MOD` have nonnegative operands and positive divisors on the
  promised generator inputs;
- `SHL`, `SHR`, `AND` and `XOR` have nonnegative operands;
- negative loop endpoints and differences require signed arithmetic and
  comparisons;
- the **numeric** shift amount is bounded by `2 A` (or a smaller word width).
  It is not enough to bound only the bit width of that shift amount.

The new ALU implements a broader signed interface, which is useful, but its
source-specific complexity analysis can use these stronger promises.
`SHL` and `SHR` iterate over the numeric shift amount, so an unconstrained
`b`-bit shift count would not give polynomial running time in `b`. The
generator's `O(A)` shift bound avoids that problem.

Similarly, linear array scans are polynomial only once both their number of
words and their bit widths are bounded. The existing bound
`(100 T+40000)(2b+4)` is for serialized word-machine data. It does not yet count
the new digit-transition journal, repeated copy/update buffers, virtual-head
markers, or physical single-head traversal. Those must be charged after the
whole bytecode interpreter has been compiled into the fixed table. A bound
of the form `T_digit <= poly(T,b,A)` would make this explicit lift's journal
polynomial; no empirical timing fit is a substitute for that bound.

## Remaining composition obligations

1. Compile the complete interpreter's stack, heap block identifiers, register
   frames, calls, output appends and instruction dispatch into the same finite
   digit-level machine. Independently implemented arithmetic and serial-array
   operations do not by themselves provide that whole-machine table.
2. Apply the [new reversible spatial compiler](REVERSIBLE_HEAD.md) to the
   completed interpreter. The compiler is implemented and proved; its
   composition with the full generator is still outstanding. Directly
   composing the older multitape lift with the deterministic spatial
   lowerer would remain insufficient.
3. Compose copy-output/handoff sequencing and a static initialized quantum
   clock. The new full-space gates already exhibit autonomous computation
   and cleanup on finite erasure fixtures. Their fixed repeated gate word
   is not itself a static initialized Hamiltonian or an output-copy proof.
4. Verify framing and initial blank sectors, and integrate parameter loading,
   the quantum data verifier, output dispatcher and endpoint readout.
5. Count the resulting whole local dimension and establish the fixed-coupling
   infrared claims. Neither task follows from finite arithmetic checks.

## Independent audit of the deterministic spatial lowering

**[checked]** The separate [`single_head.py`](../single_head.py) correctly
lowers an asynchronous table to one physical head on initialized tapes. Its
interleaved encoding uses `k` cells per virtual position, each carrying a
fixed tape tag, one source symbol and one marker bit. Thus it needs only
`2+2k|Gamma|` physical symbols, not the much larger Cartesian product of all
tape tracks. A virtual move is exactly `k` unit physical moves. Each source
step returns to the left anchor, using at most `2S+k+4` microsteps when the
preallocated physical interval has `S` cells including its endpoints.

[`single_head_checks.py`](single_head_checks.py) checks **918** complete
source steps for one, two and three virtual tapes, all active tape choices,
all three move directions, all head-position combinations in a three-position
window, and every scanned symbol. It also exhaustively checks **1,190** local
control/symbol pairs for closure and unit moves. All comparisons and the
claimed time bound pass, using **17,136** physical microsteps. Results are in
[`single_head_results.json`](single_head_results.json).

The initialized-format promise matters. A scan rejects when its selected
marker is missing; it does not validate every tape tag and every marker
multiplicity. Such global format validation remains a separate task.

The same audit contains an explicit **noninjective microstep witness**, even
when the simulated source step merely preserves a symbol and enters its
halting state. A source-update control at position `j` can enter the return
control at `j` without moving; a return control at `j+1` can move left into
the identical configuration. These two predecessor configurations merge.
Therefore composing this spatial lowering with the reversible multitape lift
would not prove a reversible physical controller. The
[new construction](REVERSIBLE_HEAD.md) repairs this using distinct marked
sweep entries and a transition-labelled history alphabet, with the actual
microstep injectivity conditions checked explicitly.

## Reproduction and provenance

```
python3 contrib/digit_compiler/audit/checks.py
python3 contrib/digit_compiler/audit/single_head_checks.py
python3 contrib/digit_compiler/audit/reversible_head_checks.py
python3 contrib/digit_compiler/audit/gate_checks.py
```

The checks use only the Python standard library and run in under a second on
the development machine. Implemented and audited independently of the ALU and
memory authors by a Codex/GPT-6-family worker on 2026-10-28. No human review or
completed counterexample is claimed.
