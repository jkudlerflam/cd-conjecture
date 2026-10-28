# Explicit finite-state arithmetic on binary tapes

This directory removes the host-integer arithmetic primitive from the generator
interface. It gives one fixed deterministic multi-tape Turing machine for all
15 binary arithmetic/comparison instructions supported by the generator VM
(13 occur in the frozen program; `LE` and `GE` are additional supported entries):
`ADD SUB MUL DIV MOD SHL SHR AND XOR EQ NE LT LE GT GE`.

The implementation is an actual digit transition table, not a Python loop whose
body performs integer arithmetic on the encoded operands. The Python builder
emits the table once. Every execution step reads one symbol, looks up one rule,
writes one symbol and moves one head by -1, 0 or +1. A finite control state
specifies the active tape. Other tape contents and head positions are unchanged.
This asynchronous presentation is an exact factorized standard multi-tape TM;
expanding each rule over the spectator symbols would give the usual synchronous
table without changing its operation.

There are **4,749 nonhalting states, 5 halting states, 8 tapes, 6 tape symbols
and 28,494 transition rules**. These counts do not depend on operand width or on
the projective-plane parameter. Every nonhalting-state/symbol pair has a rule;
unexpected symbols enter `ERROR_MALFORMED`. The machine's arithmetic specification
assumes canonical, finitely framed inputs and blank space beyond their endpoints.
It does not establish halting on arbitrary malformed infinite configurations.

## Encoding and interface

A signed-magnitude word is `P b0 b1 ... bk R` or `M b0 b1 ... bk R`, where the bits
are little endian. `P0R` is the unique zero. Otherwise the last digit is `1`.
The alphabet is `(_, P, M, 0, 1, R)`, with `_` blank. The first two tapes hold the
two operands and their heads start on the sign symbols at position 0. Six scratch
tapes start with `P0R`, heads at position 1. The output is tape 2; input contents
are preserved. Scratch data remain at a successful halt and are recovered by
reversing the recorded computation; ordinary forward arithmetic is not itself
an injective map.

`machine_for(op)` exports `initial_state`, `halting_states`, `alphabet`, `blank`,
`tape_count`, `active_tape` and `transitions`. Rules have the exact form

```
transitions[(state, read_symbol)] = (next_state, write_symbol, unit_head_move)
```

`run_binary(op, left_bits_le, right_bits_le, left_sign=0, right_sign=0)` is a
reference wrapper. It accepts tuples of bits, returns `bits_le` and `sign`, and
contains no encoded-integer decoder or arithmetic implementation. `Runner` gives
access to individual steps and optional history recording. A shared simulator
may instead consume the exported table directly. `transition_table.json` is the
complete serialized table, with each rule a five-entry row; `checks.py` rebuilds
it, records its SHA-256 digest, and checks the parsed rules for exact equality.

## Algorithm and correctness invariants

* **Addition/subtraction:** a left-to-right digit sweep carries one bit of carry
  or borrow in the control. Unequal signs select a magnitude comparison and the
  appropriate subtraction direction. Zeroes are normalized after every result.
* **Magnitude comparison:** sweep from least to most significant digit, replacing
  the remembered comparison whenever the current digits differ; missing digits
  are zero. The last differing digit therefore gives the correct comparison.
  Signed comparisons reverse that order for two negative inputs.
* **Multiplication:** repeatedly read the multiplier's low bit, conditionally add
  the current multiplicand, shift the multiplicand left and multiplier right.
  The invariant is `accumulator + current_multiplicand * remaining_multiplier =
  original_product`.
* **Division/modulo:** binary long division scans a copy of the numerator from
  its most significant digit. Each step shifts the remainder left, appends the
  next digit, shifts the quotient left and conditionally subtracts the positive
  divisor and appends quotient bit 1. Before each iteration `0 <= remainder <
  divisor`; hence after appending one digit at most one subtraction is needed.
  The invariant for the processed prefix is `prefix = quotient * divisor +
  remainder`. Opposite operand signs with nonzero remainder increment the
  quotient magnitude, giving Python floor division. The associated remainder
  has the divisor's sign and satisfies `a = b*(a//b) + a%b`. Division by zero
  enters an explicit error state.
* **Shifts:** a binary counter on scratch tape controls repeated one-digit shifts;
  the counter is decremented by digit transitions. A negative right shift rounds
  toward negative infinity by incrementing the shifted magnitude when a discarded
  digit is 1. This implements Python's arithmetic right shift. Negative shift
  counts enter an explicit error state.
* **AND/XOR:** for each negative input, a finite carry flag converts magnitude
  bits to an infinite two's-complement stream. Once a magnitude endpoint is
  reached its stream sign-extends. One extra sign-extension digit is processed
  beyond the longer input; the result stream is converted back to magnitude with
  another finite carry flag and normalized. The result sign is negative exactly
  when both inputs are negative for AND, or exactly one is negative for XOR.

All loops in these algorithms are edges in the finite transition graph. No
Python algorithm chooses the next iteration, branch, digit position, carry,
borrow or arithmetic result while executing a table.

For maximum operand magnitude width `w`, elementary arithmetic, bitwise and
comparison routines take `O(w)` moves, multiplication and binary division take
`O(w^2)` moves, with constants determined by this table. A shift by numerical
count `b` takes at most `O(b*(w+b+log(b+1)))` moves. These are algorithmic bounds,
not an optimized uniform resource bound for the full generator. In particular,
packed generator flags can have polynomially many bits; they cannot all be
treated as `O(log p)`-bit integers.

[BOUNDS.md](BOUNDS.md) gives explicit conservative constants for these actual
tables and a conditional specialization to the frozen generator. In particular,
its arithmetic contribution is at most `128 T (B+1)(2B+5)` digit transitions
under the stated source-level width and shift bounds; `T` remains the old word
instruction count and is not reused as a digit running time.

## History and the remaining locality distinction

Optional history stores one finite event `(old_state, old_scanned_symbol)` per
step. That event determines the active tape, old symbol, written symbol, move
and next state. Reversal first undoes the head move, checks the written symbol,
restores the old symbol and restores the control. Thus a complete run followed
by its recorded inverse recovers the entire initial configuration, including
all inputs, scratch tapes and head positions. An event can be encoded by one
of 28,494 rule IDs (15 binary digits); its size is independent of operand width.

`checks.py` verifies the forward/inverse operation of **every local transition**
and whole histories for signed arithmetic, including rejected arithmetic inputs.
The append/pop of history is currently an explicit simulator history stack.
This alone is **not** an autonomous reversible nearest-neighbor controller or
a completed qudit Hamiltonian. A physical history tape, its reversible routing,
and the connection to serial memory and quantum program output remain separate
compiler layers. A local permutation lift of the transition alphabet can be
audited independently from these arithmetic correctness checks.

## Verification

Run `python3 contrib/digit_compiler/alu/checks.py` from the repository root.
The check uses Python arithmetic only as an independent expected-value oracle.
It exhausts operand pairs from -16 through 16 (shift counts 0 through 12), then
tests carry/borrow boundaries and seeded random words through 257 magnitude bits.
The reported maximum is 257 because the boundary fixture includes `2**256`.
It checks input preservation, canonical outputs, absence of stray nonblank cells
beyond the output endpoint, and exact complete-configuration history rewinds.
`results.json` contains operation-by-operation counts and the table digest.

Provenance: generated and tested by a GPT-6-family implementation subagent on
2026-10-28. This document does not claim human review or a c-d counterexample.
