# An actual finite local arithmetic layer

This directory replaces eleven aligned-word operations by **4,428 frozen
two-cell rewrite rules** over **900 used cell symbols**, independent of word
width. The rules implement `ADD`, `SUB`, `AND`, `XOR`, `COPY`, `EQ`, `NE`,
`LT`, `LE`, `GT`, and `GE`. A run autonomously computes, copies its answer,
clears its workspace and returns the head to the left endpoint. Both inputs
are preserved. The execution engine only matches adjacent symbol pairs;
there is no integer operation in a transition.

This is a component of the sparse-circuit generator's missing digit-level
compilation. It is **not a complete finite-alphabet implementation of that
generator**, a physical dimension estimate, or a counterexample.

## Interface to the frozen generator

The generator's resource bound supplies a common sufficient signed word
width `b`, including the packed flag's `O(Lw)` digits. These kernels accept
two aligned, little-endian, two's-complement input tracks of any width
`b >= 1`, a blank history track, and initially zero scratch and output
tracks. No control state, rule count or cell alphabet depends on `b`.

`arithmetic.py` exports:

```python
apply(op, a, b, width)            # integer fixture adapter
run = execute(op, a, b, width)  # output, steps, initial_tape, final_tape
run.verify_inverse()           # independent reverse-table verification
```

`ADD` and `SUB` return arithmetic modulo `2**width`, decoded as signed.
They agree with the generator's unbounded integer arithmetic whenever its
operands **and result** fit the chosen signed width. `AND` and `XOR` give
the exact Python signed-integer result when their inputs fit. Comparisons
are signed and return the nonnegative word `0` or `1`, including width one.
At width one a Boolean output of `1` is intentionally interpreted unsigned;
a generator-wide signed convention should use width at least two.
`COPY` copies the first operand to output and ignores the second operand.

The host-side `encode` and `decode` functions bridge test integers to digit
tapes. They are **not** compiled local routines. Likewise, importing these
kernels into a word interpreter does not compile its stack, operand
alignment, call frames, heap addressing, or opcode dispatch. The genuine
local component begins with the encoded, aligned input tape and ends with
the encoded output tape.

## Tape and complete local semantics

A digit symbol `Fabsot` or `Dabsot` has bits `a,b` for the input tracks,
bit `s` for scratch, bit `o` for retained output, and history marker `t`
from the fixed set `_,0,1,L,E,G`. `F` marks the first digit, `D` the other
digits. Boundaries are `L` and `R`. A head symbol is the single cell symbol
`H:operation:phase:state:underlying-symbol`; its components are notation
for a member of a finite alphabet, not separate unbounded fields.

A legal entry tape is

```
H:op:enter:_:L   F a0 b0 0 0 _   D a1 b1 0 0 _   ...   R
```

(The displayed spaces within digit symbols are only for readability.)
The finite table specifies these phases without external switching:

1. `enter` transfers the head from `L` to `F` in the operation's initial
   carry/comparison state.
2. `scan` moves right one cell. It writes the old finite control state to
   the digit's blank history track and writes one computed scratch bit.
   Addition uses carry zero; subtraction computes `a + (~b) + 1`.
   Bitwise operations use one dummy state. Comparison retains the most
   significant unequal-bit result seen so far. When the next cell is `R`,
   the sign-bit ordering is reversed.
3. At `R`, a local rule changes `scan` to `back`. For comparisons it retains
   the final Boolean answer in a finite head-state bit.
4. `back` moves left, verifies the local carry/comparison relation using
   the saved old state, clears history and scratch, and copies the scratch
   bit to output. Comparison instead writes its answer only on `F`.
5. At `L`, the final local rule enters `done`. For comparison, the answer
   bit in `F` determines which rule removes the redundant head-state bit.

Every run takes exactly `2b + 3` local transitions on `b + 2` cells. The
return pass is part of the **forward** finite transition table. There is
no host decision to reverse at the endpoint. The separate inverse-table
option in the test engine is solely a verification facility.

All scratch/history tracks are blank again at completion; input tracks
and delimiters are unchanged. The only retained data are the answer bits
and a fixed `done` control state. Addition's overflow carry is not retained
as output; this is why its public semantics are modular arithmetic.

## Why the transitions are reversible

The old finite state saved at a processed digit uniquely recovers the
predecessor of a scan transition. Conversely, the cleaned digit and old
state in a returning head reconstruct the erased scratch/history and the
next carry/comparison state. The saved verdict bit disambiguates comparison
return rules; the first output digit disambiguates its final exit rule.

`checks.py` exhaustively verifies that every two-cell left-hand side and
right-hand side is unique. It also verifies that every decorated head
symbol chooses exactly one adjacent side, for both the forward and inverse
relation. Since each rule has exactly one head on each side, these finite
checks establish at most one successor and at most one predecessor on
**every one-head tape**, not only on the tested valid trajectories. They
do not establish a valid-initialization energy penalty or exclude unwanted
Hamiltonian components with multiple heads or malformed boundary patterns.

The arithmetic correctness argument is induction over significance: a
scan digit satisfies the binary full-adder equation, bitwise truth table,
or comparison update. The right endpoint handles the sign digit for
comparison. The return sweep copies the computed answer while undoing
each local update. Tests independently compare against integer arithmetic
and check the inverse after all completed runs and after every prefix of
representative runs.

## Artifacts, verification, and remaining work

`build_rules.py` is a finite rule generator. `rules.json` is the actual
executable transition artifact, and `arithmetic.py` loads it directly.
`manifest.json` records its canonical digest and counts: 4,428 rules, 900
used cell symbols, 74 used underlying tape symbols and 90 control states.
The 900 symbols describe only this local arithmetic component, **not**
the final physical dimension after composition with the rest of the machine.
The directed rule relation has overlapping domain and range and includes
stationary endpoint turns. It is not asserted to be a disjoint collection
of transpositions or an involution. Compiling it into the existing
directional sweeping controller's unitary format needs a separate lift;
the 900-symbol count is not the dimension of a completed pair of local
unitaries `F,G`.

`checks.py` checks all operand pairs at widths one through six for all
eleven operations, boundary/random cases through 4,096 digits, exact
whole-run and prefix inverses, and rejection by stalling of dirty scratch,
output and history tracks in a fixture. The latter is a domain check only.

Remaining generator instructions include `MUL`, `DIV`, `MOD`, variable
`SHL`/`SHR`, addressed memory, stack/register transport, calls, journaling
of the whole VM, autonomous output copying and reverse execution of the
whole generator, and the bridge to its dispatcher. One local arithmetic
copy/uncompute sweep does not supply those global mechanisms. Composition
must also prove initialization, unique accepted history, finite range,
coefficient normalization and the required spectral bound.

Reproduce using only the standard Python library:

```sh
python3 contrib/digit_compilation/arithmetic/build_rules.py
python3 contrib/digit_compilation/arithmetic/checks.py
```

Implemented with AI assistance on 2026-10-28. No human review or completed
counterexample is claimed.
