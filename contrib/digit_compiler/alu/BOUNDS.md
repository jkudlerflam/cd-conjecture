# Explicit move bounds for the actual arithmetic tables

These bounds count transitions of `machine.py`, rather than calls to a word-RAM
arithmetic primitive. The constants are deliberately conservative. All counts
assume the canonical finite input framing described in `README.md`. They do not
claim termination or these bounds for arbitrary malformed tapes.

Let `b >= 1` be the larger of the two input magnitude bit widths. For a shift,
let `s >= 0` be the *numerical* shift count, not its bit width. Let `D(op)` count
all transitions from the sign-reading entry state through `HALT`. Then:

\[
\begin{aligned}
D(\mathrm{op})&\le 128(b+1)(2b+5)
&&\text{for every nonshift arithmetic/comparison operation},\\
D(\mathrm{SHL}),D(\mathrm{SHR})&\le64(s+1)(b+s+5).
\end{aligned}
\]

An arithmetic error from division by zero or a negative shift count halts earlier
than the corresponding setup phase. The successful-run bounds already cover
division by zero; a negative shift count takes two entry transitions.

## Primitive loop accounting

Fix an invocation and choose `M` large enough to bound the number of magnitude
digit positions in every input, scratch and temporary word during it, including
unnormalized intermediate results. All heads remain between the header at 0 and
at most one cell past the largest endpoint. Put `H = M+4`. Scans below terminate
on the sign/header or endpoint symbol. Their iteration counts are bounded by
`M+2`, including endpoint processing. Each elementary `action` costs one step.

The following bounds include each macro's own entry and cleanup transitions but
exclude its continuation, which is counted by the caller:

| Actual builder macro | Bound | Accounting |
| --- | ---: | --- |
| `rewind` | `H` | At most one left move per occupied position, then one header transition. |
| `end` | `2H` | A rewind followed by one rightward endpoint search. |
| `clear` | `6H` | Initial rewind; erase sweep; rewind over erased cells; write `0R`; final rewind. |
| `normalize` | `6H` | Seek endpoint; erase it; trim zeros leftward; write endpoint or `0R`; rewind. |
| `zero` | `3H` | Rewind; scan until the first 1 or endpoint; rewind. |
| `copy` | `12H` | Clear destination; rewind source; two transitions per copied bit; endpoint write; rewind both heads. |
| `increment` | `4H` | Rewind; carry scan; at most one endpoint extension; rewind. |
| `decrement` | `8H` | Rewind; borrow scan; normalize. |
| `shift_left` | `9H` | Rewind; one transition per carried digit; endpoint write; normalize. |
| `shift_right` | `12H` | Rewind; first-digit read; three transitions per remaining digit; endpoint replacement; normalize. |
| `compare` | `6H` | Two initial rewinds; two reads per digit column including the endpoints; two final rewinds. |
| `arithmetic` | `20H` | Clear destination; two rewinds; three transitions per digit column; at most four endpoint/carry transitions; normalize; two rewinds. |
| `bitwise` | `20H` | Clear and rewind; three transitions per digit column, including one sign-extension column; endpoint write; normalize and rewind inputs. |
| `sign` | `4H` | Zero test and at most two transitions to replace the sign. |
| `constant` | `7H` | Clear result and write its single output bit. |

For example, the addition/subtraction macro has the bound
`6H + 2H + 3M + 4 + 6H + 2H <= 20H`. This explicitly counts the sequential
single-tape reads and output write; a whole digit column is not treated as one
machine step. Likewise, `shift_right` uses three transitions per copied digit,
because it reads the next digit, moves back to write it, and advances again.

## Composition of the actual control loops

The input sign dispatch costs exactly two transitions.

* Signed addition/subtraction needs at most one comparison, one magnitude
  arithmetic invocation and one final sign normalization, hence at most
  `30H+2`. Bitwise operations and signed comparisons are cheaper.
* Multiplication has setup and finalization at most `37H`, and at most `b`
  multiplier-digit iterations. Each iteration costs at most `57H`: one low-bit
  read, optional `arithmetic` and `copy`, both shifts, and the next zero test.
  Including dispatch, `D(MUL) <= (57b+37)H+2`.
* Long division has setup at most `29H` and exactly as many numerator iterations
  as its canonical magnitude has digits, hence at most `b`. An iteration costs
  at most `65H`: one numerator read, remainder shift and optional increment,
  quotient shift, comparison, optional subtraction/copy/quotient increment,
  and one leftward numerator-head move. The signed DIV/MOD correction and
  finalization cost at most `27H`. Including the endpoint and dispatch,
  `D(DIV),D(MOD) <= (65b+57)H`.
* A shift has setup/finalization/dispatch at most `32H` and exactly `s`
  counter iterations. Each costs at most `28H`: controller entry, decrement,
  one digit shift, optional negative-right-shift correction, and zero test.
  Thus `D(SHL),D(SHR) <= (28s+32)H`.

For nonshift operations choose `M = 2b+1`. This bounds multiplication's shifted
multiplicand and product, the one extra bit in the division remainder, and the
extra sign-extension column of bitwise operations. Therefore `H = 2b+5`, and
all three nonshift bounds are below `128(b+1)H`.

For shifts choose `M = b+s+1`, including the possible temporary extra zero before
normalization. Therefore `H = b+s+5`, and the shift bound is below `64(s+1)H`.
These choices also bound head excursions during clearing and endpoint insertion;
the constants in `H` cover their transient blank cells.

## Conditional specialization to the frozen generator

Let `T` be an upper bound on the number of instructions of the **original
word-RAM generator**, including its arithmetic instructions. Let `A` be the
project's circuit resource parameter. Suppose the separate source/resource audit
establishes that, throughout that generator run,

\[
B=2A+\lceil\log_2(T+1)\rceil+64
\]

bounds every operand's magnitude bit width, and every shift count is nonnegative
and at most `2A`. This paragraph is a **conditional specialization**: the digit
table does not itself prove those source-level operand bounds.

Because `2A+1 <= B+1` and `B+2A+5 <= 2B+5`, the cost of any one arithmetic opcode
is at most

\[
K(A,T)=128(B+1)(2B+5).
\]

There are at most `T` arithmetic opcodes, so their combined digit-TM contribution
is bounded by

\[
\boxed{D_{\rm ALU}\le128T(B+1)(2B+5)
       \le256T(B+3)^2.}
\]

This does **not** identify `T` with the digit running time. It replaces each
word arithmetic instruction by its counted binary-digit computation. It also
does **not** include serial-memory access, VM dispatch/call-stack management,
program emission, a single-head simulation, history-tape routing, or the final
Hamiltonian clock. Those costs must be added or composed with this bound.

Each individual ALU invocation under these assumptions can be laid out on eight
tracks with magnitude width at most `2B+1`, so at most `8(2B+3)` framed tape cells
contain its input/output/scratch words. This is working storage only. A direct
transition journal additionally requires one finite record per executed digit
step; its physical layout and reversible compiler incur separate costs.

The existing [word-resource calculation](../../autonomous_controller/generator/resource_bounds.json)
gives `T <= 18463(128A)^6 < 2^57 A^6`. (Its loop variable named `B` is unrelated
to the operand-width `B` in this note.) For integer `A >= 1`,

\[
\lceil\log_2(T+1)\rceil
\le57+6\lceil\log_2 A\rceil
\le57+6(A-1),
\]

so `B <= 8A+115 <= 123A <= 128A`. Consequently a coarser, explicitly polynomial
corollary is

\[
\begin{aligned}
D_{\rm ALU}
&\le128T(128A+1)(256A+5)\\
&\le4{,}309{,}632\,T A^2
<2^{23}T A^2
<\boxed{2^{80}A^8}.
\end{aligned}
\]

This is still the arithmetic contribution only, conditional on the source-level
resource assumptions just stated. It is not a bound on the finished autonomous
compiler or Hamiltonian history length.

## Numerical checks

`bounds_checks.py` independently compares the inequalities with actual table
runs on signed boundary and seeded random inputs, including widths above 256
bits and shifts through 96. These checks are regression evidence for the stated
accounting; the exhaustive finite transition audit and the loop argument above
are the reason the bound is independent of tested word widths. The check also
evaluates the conditional specialization and verifies its algebraic inequalities
on a grid of resource parameters. Results are in `bounds_results.json`.
