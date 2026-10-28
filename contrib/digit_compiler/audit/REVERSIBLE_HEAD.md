# A finite, reversible, single-physical-head compiler

**[proved]** [`reversible_head.py`](reversible_head.py) compiles an arbitrary
finite deterministic one-head transition table into a reversible one-head
partial transition table. The compiled machine autonomously performs all
work-head updates, history appends and intervening physical sweeps. Its
scanned-symbol alphabet and control set are fixed independently of input
length and computation time. There is no addressed memory, integer operation,
global tape inspection, or externally chosen sweep phase in a runtime step.

**[checked]** The implementation composes with the actual deterministic
spatial lowering of the ALU. It computes `0 == 1` and `1 + 1` by **648,266**
and **3,055,626** reversible physical microsteps respectively. Replaying the
local inverse for exactly the same number of steps restores every work cell,
history cell, marker and head position. This is a genuine composed arithmetic
execution, not only a toy source machine.

This closes the spatial-reversibility issue for a supplied finite source
table. The full word-bytecode interpreter is still to be supplied as such a
table. Initialized quantum history-Hamiltonian legality, autonomous
copy/turnaround control and its full physical dimension are separate tasks.

## Encoding and local rules

For source alphabet `Gamma`, use work symbols `W(a,0)` and `W(a,1)`. Exactly
one work cell is marked `1`, indicating the virtual source head. All work
cells lie to the left of a separator. To its right lie the history records,
one frontier symbol `H`, and blank history capacity. Two outer end symbols
bound the entire physical tape.

A source rule is

```
t = (q,a): (q,a) -> (r,b,d), d in {-1,0,1}.
```

The finite label `t` is stored in control or as one history symbol. It is a
pair drawn from the finite source state/alphabet sets, **not** a coordinate,
word, integer input, or variable-length history object.

Write `plain` for an unmarked work symbol, separator, history record or blank
history symbol. It excludes a marked work cell, the history frontier, and
the two outer ends. These are all the compiled rule families:

| Input | Output | Move |
| --- | --- | --- |
| `F(q), W(a,1)`, source move `d!=0` | `MOVE(t), W(b,0)` | `d` |
| `F(q), W(a,1)`, source move `d=0` | `DEPART(t), W(b,1)` | `0` |
| `MOVE(t), W(c,0)` | `DEPART(t), W(c,1)` | `0` |
| `DEPART(t), W(c,1)` | `SCAN(t), W(c,1)` | `+1` |
| `SCAN(t), plain` | `SCAN(t), plain` | `+1` |
| `SCAN(t), H` | `APPEND(r), history(t)` | `+1` |
| `APPEND(r), blank` | `RETSTART(r), H` | `0` |
| `RETSTART(r), H` | `RETURN(r), H` | `-1` |
| `RETURN(r), plain` | `RETURN(r), plain` | `-1` |
| `RETURN(r), W(c,1)` | `F(r), W(c,1)` | `0` |

`MOVE(t)` exists only for a source rule with nonzero displacement. A source
halting state has no `F` departure. Every unlisted input is undefined, rather
than merged into a common rejecting state. This choice preserves partial
injectivity on arbitrary malformed tapes.

## Proof of reversibility on the full declared configuration space

Two local criteria suffice for a one-head partial TM:

1. Each target control has one fixed incoming move direction.
2. The pair `(target control, written symbol)` identifies at most one source
   `(control, scanned symbol)`.

Given a complete output configuration, the first criterion identifies the
predecessor head position by subtracting that fixed direction. Reading the
written symbol there and using the second criterion reconstructs the unique
predecessor local state. Hence the entire configuration map is injective on
its domain, without initialization assumptions.

Both conditions hold explicitly here. The incoming directions of `F`,
`MOVE(t)`, `DEPART`, `SCAN`, `APPEND`, `RETSTART`, `RETURN` are respectively
`0,d(t),0,+1,+1,0,-1`. The entry to a `SCAN(t)` sweep writes a **marked** work
symbol, whereas the scan loop only writes plain symbols. The entry to a
`RETURN(r)` sweep writes `H`, whereas the return loop only writes plain
symbols. These output sectors are disjoint. At the append, distinct source
rules with the same target `r` write distinct `history(t)` symbols. All
other rules are plainly injective because they preserve the scanned symbol
or keep its old value in `t`.

[`inverse_rule`](reversible_head.py) implements this exact partial inverse.
The inverse runner determines a direction from the current finite control
and reads the single preceding physical cell. It does **not** consult a host
journal, saved head coordinates, source inputs, or snapshots.

## Correctness, capacity and cost

On a valid initialized tape, `F(q)` sits at the unique marked work cell. The
first one or two transitions apply the source rule and move its mark by `d`.
The machine then sweeps right to the frontier, appends `t`, moves the frontier
right by one cell, and returns left to the marked work cell in `F(r)`.
Everything except the source update and the one appended record is preserved.
Induction proves exact simulation for every defined source computation that
stays within its work region and allocated history capacity.

Let the source use a work interval of `W` cells and make `T` source steps.
Allocate `T` blank history cells beyond the initial frontier. The full tape
then has

```
S = W + T + 4
```

cells, including both ends, the separator and initial frontier. Each source
step takes at most `2S+8` physical microsteps. Therefore

```
T_reversible <= T (2W + 2T + 16).
```

Storage is `W+T+4`; reverse replay takes the same number of microsteps. Thus
this is a polynomial-overhead, fully spatial reversible simulation. If the
source is the earlier deterministic multitape lowering, substitute its
explicit `T <= T_digit(2W+k+4)` bound. The remaining whole-generator resource
question is to prove and instantiate `T_digit`; the physical reversible
compiler introduces no new superpolynomial dependence.

For a source with `Q` controls, alphabet size `g`, and `m` defined rules, a
safe count is

```
reversible controls <= 4Q + 3m,
reversible tape alphabet <= 2g + m + 5.
```

For the actual lowered ALU, the code uses `Q=47,589`, `g=98` and
`m<=4,663,134`, giving at most **14,179,758 controls** and **4,663,335 tape
symbols**. These large constants are explicit, independent of the prime,
and not optimized. They are **not** the final qudit dimension. The source
rule labels and the factorized local evaluator specify this finite table
exactly without exporting its enormous Cartesian product. Each evaluator
call depends only on its finite control and scanned symbol.

## Checks and the remaining boundaries

[`reversible_head_checks.py`](reversible_head_checks.py) checks a complete
40-control, 20-symbol fixture table: **265 defined rules**, all **800** inverse
control/symbol pairs, fixed incoming directions, and unique incoming pairs.
It additionally round-trips **10,000** malformed complete configurations
through individual defined local transitions. A destructive source that
erases either `0` or `1` produces the same work output while preserving its
distinguishing history: every binary string through length six is run to
completion and then restored, **127 cases and 14,463 steps each way**.

The two actual lowered ALU runs add **3,703,892 forward and 3,703,892 inverse
physical microsteps**, with exact arithmetic answers and complete cleanup.
Results are saved in [`reversible_head_results.json`](reversible_head_results.json).

The forward sweeps are autonomous, but the test harness chooses when to
invoke the inverse runner in the tests above. The subsequently implemented
[`reversible_gates.py`](../reversible_gates.py) supplies an explicit
stationary-move conversion and full-space local collision/movement
permutations, described below. Its fixed repeated gate word also exhibits
autonomous computation/erasure cycles on finite fixtures; a physical
copy-output/handoff controller remains to be composed. The initialization
adapter lays out the tagged source
tape and blank history capacity; it is not an initialization algorithm or
an illegal-sector spectral proof. The initial configuration may have a
predecessor outside the initialized computation, so the result is not a
claim that the unpenalized history graph begins there.

## Full-space gates and independent autonomous-cycle audit

**[proved]** The root gate module replaces every stationary transition into
control `r` by a move right into a fresh bounce control `I_r`, then a move
left into `r` preserving the neighboring symbol. The existing fixed incoming
direction property prevents these from colliding with nonstationary incoming
transitions. One spare neighbor suffices. All resulting moves are `+1` or
`-1`.

Introduce fresh and updated head phases. For each normalized partial rule,
the collision gate swaps its fresh source `(q,a)` and updated target `(r,b)`;
all unmatched inputs are fixed. These are disjoint transpositions by the
proved local injectivity. A separate two-site movement gate swaps an updated
head and its fresh neighbor position in its fixed incoming direction.
Zero-head and two-head sectors are fixed. Hence both gates are permutations
on their **complete local spaces**, including malformed inputs.

**[checked]** [`gate_checks.py`](gate_checks.py) independently checks the
normalized finite table, **440** source-rule/bounce-neighbor contexts, all
**407** collision inputs, and all **165,649** two-site movement inputs.
The latter includes **121** zero-head, **8,712** one-head and **156,816**
two-head inputs. All maps are exact involutions and movement preserves both
physical symbols and head count.

Independently of the root cycle tests, five erasure inputs return exactly to
their initial complete configuration under repeated application of the same
gate word: collisions at every site, then bond movements from left to right.
No inverse is invoked during this evolution. The inputs `""`, `"0"`, `"1"`,
`"010"`, `"10110"` have cycles of **38, 70, 70, 158, 278** fixed sweeps.
Each orbit first visits the correct erased output. These are finite-fixture
checks; the initialized state can have predecessors, so the orbit may
include an additional excursion beyond the chosen initial configuration.
An undefined partial transition is a reflecting point of the permutation
dynamics, **not** an absorbing reversible halt.

The externally repeated, state-independent gate word is explicit. Making
that gate word a static initialized quantum history Hamiltonian, retaining
a copied output, and composing the complete generator remain separate
requirements. The finite-cycle result alone proves none of those claims.

Reproduce with:

```
python3 contrib/digit_compiler/audit/reversible_head_checks.py
python3 contrib/digit_compiler/audit/gate_checks.py
```

**Provenance.** Implemented by a Codex/GPT-6-family audit worker and reviewed
independently by the arithmetic worker on 2026-10-28. The independent review
caught and repaired a declared-control mismatch: the return states now
exist for every source control, matching the complete inverse and the
`4Q+3m` bound. No human review or completed counterexample is claimed.
