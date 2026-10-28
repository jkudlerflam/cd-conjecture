# Digit arithmetic and reversible spatial compilation

The sparse-graph candidate now has explicit fixed-alphabet arithmetic and
serial-array transition tables, a one-physical-head reversible compiler,
and full-space local permutation gates for the resulting machine. These
replace several previously abstract computational operations with concrete
rules. The full 3,115-instruction generator still needs its stack, heap,
dispatch and initialization compiled and joined to the quantum dispatcher.
There is no new c-d violation or measured central charge in this contribution.

## Claims and status

**[proved finite constructions; checked implementations]** The arithmetic
machine implements signed integer operations on canonically framed binary
words. The memory machine implements selected-array access and functional
updates on framed binary-word lists. Every transition reads one digit,
writes one digit and moves one tape head by at most one cell. Its finite
control and alphabet do not depend on input width, array size or prime.
The input-format and blank-scratch hypotheses are specified in the component
READMEs; these routines are not whole-input format verifiers.

| Component | Active controls | Rules | Tapes | Symbols per tape |
| --- | ---: | ---: | ---: | ---: |
| Signed arithmetic and comparisons | 4,749 | 28,494 | 8 | 6 |
| GET, functional SET, APPEND, LEN | 97 | 207 | 5 | 6 |
| Preserved-output copy | 8 | 20 | 2 | 6 |

The ALU supports 15 operations, including all 13 arithmetic operations used
by the frozen program. It implements Python floor division, signed modulo,
arithmetic right shift, and infinite-sign-extension bitwise semantics.
Binary long division takes polynomial time in operand width. Numerical
shift counts require their own bound, supplied by the generator's source
promises. [Explicit loop bounds](alu/BOUNDS.md) give an arithmetic-only
contribution below `2^80 A^8` under the existing generator bounds; this does
not include memory, dispatch, routing or history.

**[proved spatial compiler; independently checked]** A tagged, interleaved
encoding converts any supplied asynchronous multitape table into one
physical head, using `2+2k|Gamma|` tape symbols and only unit moves. It uses
98 symbols for the arithmetic source and 62 for the array source. The
[deterministic construction](SPATIAL_LOWERING.md) has a proved finite-space
and polynomial-time overhead, but is itself irreversible.

The subsequent [reversible compiler](audit/REVERSIBLE_HEAD.md) supplies the
missing history routing with that same single physical head. A finite rule
label is carried to a separate history region, stored, and the head returns
to the marked work cell. Every target control fixes its incoming direction;
every target-control/written-symbol pair has at most one predecessor. These
conditions prove partial reversibility even on malformed configurations.
No Python journal or random-access step is needed to run its local inverse.

**[proved local gate completion; checked cycles]** The
[gate construction and proof](LOCAL_GATES.md) are implemented in
[`reversible_gates.py`](reversible_gates.py). Stationary transitions become
right/left excursions. Fresh and updated head sectors make each
collision a disjoint transposition. Head movement is another disjoint
transposition on two adjacent sites. Undefined collisions and zero/two-head
bond sectors have explicit identity action. Thus both gates are full-space
permutations and define exact unitary matrices on their declared finite
alphabets. They can act as the identity on a separate quantum payload.

Repeating exactly the same collision layer and movement sweep also gives
complete forward/reset cycles on 63 initialized destructive-computation
fixtures. The driver never changes to an inverse evolution. The finite
checks see the correct computed answer and then exact restoration of every
work cell, history cell and head state. The gates reflect at undefined
transitions; no absorbing reversible halt is assumed. This is a tested
autonomous fixture, not a theorem that every arbitrary source program
has the needed initialization, readout and reset protocol.

## Evidence from the actual generator

The frozen bytecode is unchanged, with SHA-256
`c9e5dd5329c5200018454406852cdf298428dcba943c02cc13f0e3f002c6ba3f`.

- [Arithmetic replay](generator_results.json) checks the actual arithmetic
  operands used to construct the projective graphs at `p=2,3,5`, additional
  width calculations, and the complete one-edge 906-gate abstract circuit.
  Every distinct operation/operand pair is evaluated by the finite digit
  table and compared with the word-machine answer.
- [Memory replay](memory/vm_results.json) covers all **8,165** memory events
  in the `p=2,3` graph runs and the complete one-edge circuit. The **4,436**
  distinct full-input cases use **75,446,718** forward/reverse/copy digit
  transitions and exactly restore scratch. All 331 actual SET events are
  covered. Graph edge order and the complete circuit match the reference.
- [Physical reversible arithmetic](audit/reversible_head_results.json)
  runs the actual spatially lowered ALU for `0 == 1` and `1 + 1`, using
  **648,266** and **3,055,626** reversible physical microsteps respectively.
  Both complete work tapes agree with the deterministic source. The same
  number of inverse local steps restores all work and history exactly.

The operation-replay tests still use the host VM to select operands, choose
heap blocks and schedule instructions. They validate real generator
interfaces without claiming a compiled whole interpreter. The composed
reversible ALU tests do execute every intervening head-routing step.

## Independent checks and resource scope

The ALU suite checks **15,759** signed cases and **13,919,188** forward digit
steps, with whole-configuration rewinds totaling **11,774,998** steps.
The array suite checks **4,158** cases and **2,714,328** forward/reverse/copy
steps, including variable-width replacement, bounds rejection, counter
carries and 513-bit addresses. Independent audits check local injectivity,
full inverse domains, malformed configurations, and source-step equivalence.
The root and independent spatial suites compare complete configurations
after **5,321** and **918** source steps respectively.

For a one-head source with `Q` controls, `g` tape symbols and `m` rules, the
reversible compiler uses at most `4Q+3m` controls and `2g+m+5` tape symbols.
If it executes `T` source steps on `W` work cells, it needs `W+T+4` total
cells and at most `T(2W+2T+16)` reversible microsteps. The resulting ALU
constants are explicit and very large; no optimization is claimed.

The old word-machine bound is not a digit-history bound. Every digit
transition, serial copy and spatial sweep must be counted before using
these results to choose a cell length. The fixed alphabets preserve the
intended asymptotic dependence on the charged record dimension, but these
component counts are not the final core constant or a violating choice of d.

## Remaining construction

The next finite task is a complete table for the closed interpreter: stack
and register frames, heap allocation and reference lookup, alias-preserving
array replacement, calls/returns, instruction fetch and output appends.
An ALU word and an array word also need a compiled framing adapter. The new
serial SET returns a modified copy of one selected array; it does not itself
update every reference to a heap block.

That complete source must then be initialized from permitted local states,
supplied with the prime and enough tape, and joined to the charged quantum
dispatcher and endpoint READY readout. Illegal sectors, the whole history
kernel and gap, complete local-dimension accounting and fixed-positive-
coupling infrared transfer remain obligations. The separate d=240 fixture
remains a one-mode demonstration and is not changed by these compiler counts.

## Reproduction and provenance

From the repository root, run:

```sh
python3 contrib/digit_compiler/verify.py
```

All checks use the Python standard library. The runner records each check,
timing and source digest in [verification_summary.json](verification_summary.json).
Individual scripts can also be run directly. The root README, STATUS and
REPORT link this continuation of the [autonomous controller](../autonomous_controller/README.md).

Implemented by the coordinating Codex agent and three parallel GPT-6-family
workers on 2026-10-28. The arithmetic, memory and compiler audit were assigned
separately. An independent review found and repaired a missing declared
return-control family in the reversible compiler; the corrected full inverse
domain was checked again. The coordinating agent reran the saved checks.
No human review is claimed.
