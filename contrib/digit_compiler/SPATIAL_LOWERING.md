# One moving head for the digit machines

## Deterministic simulation theorem

**[proved; independently checked on finite examples]** Let a deterministic
asynchronous machine have a fixed finite alphabet `Gamma`, `k` work tapes,
and a finite control set. In each state `q` it reads and writes just the
current cell of tape `i(q)` and moves that head by `-1`, `0`, or `1`.
[`single_head.py`](single_head.py) compiles it into a deterministic machine
with one physical head and alphabet size

```
2 + 2 k |Gamma|.
```

This theorem assumes a correctly initialized tagged tape, one virtual-head
marker per work tape, and enough preallocated space for the source execution.
It is a theorem about classical deterministic computation. This lowering by
itself is not reversible and cannot be used directly as a quantum gate.

## Tape and local transitions

Reserve a virtual interval `lo,...,hi`, containing `L` positions. Its physical
encoding has `S=kL+2` cells. The first and last cells have symbols `LEFT` and
`RIGHT`. Between them, position `p` is the block

```
(0, a_0(p), h_0(p)), ..., (k-1, a_(k-1)(p), h_(k-1)(p)).
```

Here `a_i(p)` is a single source symbol and `h_i(p)` a marker bit. For each
`i`, exactly one marker is set. The tape tag is a constant-size neutral symbol;
it does not depend on the prime, source operand width or cell length.

A source-step boundary has the physical head at `LEFT`, with its source
state in finite control. The local rule evaluator receives only that control
and the currently scanned physical symbol. It has no reference to the tape,
the head position, another cell, or a decoded source operand.

1. `SCAN(q)` walks right until it encounters the marked cell tagged `i(q)`.
   It then consults the fixed source rule `(q,a)->(r,b,d)`.
2. If `d=0`, write `b`, preserve the marker, and enter `RETURN(r)`.
   Otherwise write `b`, clear the old marker, and walk exactly `k` physical
   cells in direction `d`. The remaining distance, at most `k-1`, is a finite
   control label. Set the marker at the destination, which must be an
   unmarked cell of the same tape tag. Then enter `RETURN(r)`.
3. Walk left to `LEFT`, preserving all cells. Start the next `SCAN(r)`, or
   enter the declared halt state if `r` is a source halt state.

Hitting an endpoint during a virtual move, failing to find the selected
marker, or encountering an undefined source rule rejects. A separate
initialization theorem is needed to exclude other malformed tapes; this
simulator does not certify every tag and marker multiplicity globally.

## Correctness and cost

Assume the encoding invariant at a source-step boundary. The first sweep
finds exactly the active virtual head. A shift of `k` cells changes the
virtual position by exactly one and preserves the tape tag. Only the active
source symbol and its marker move; all other source symbols and virtual
heads are unchanged. The final sweep restores the physical head to `LEFT`
and retains precisely the next source state. Induction proves exact
simulation at every source-step boundary, including all scratch contents
and head positions.

One macrostep uses at most `2S+k+4` physical microsteps. The two long scans
have total length at most `2S`; updating, marking and changing phase cost at
most four additional transitions, with the bounded virtual move costing at
most `k`. This deliberately loose bound covers zero moves and halting.

If the initial data and virtual heads occupy an interval of length `L0`,
then a source run of at most `D` digit steps fits inside an interval of length
`L0+2D`: no virtual head can move farther than `D` cells in either direction.
The spatial simulation therefore takes at most

```
D [2(k(L0+2D)+2) + k + 4]
```

microsteps on `k(L0+2D)+2` cells. This is a bound for a supplied digit-level
source program. It does not replace `D` by the old word-machine instruction
count; the compiled interpreter and all of its storage operations must first
be counted at the digit level.

The implementation explicitly enumerates its finite controls and gives a
total local evaluator on that control set and alphabet. Its present counts
are:

| Source component | Physical symbols | Physical controls | Complete local table entries |
| --- | ---: | ---: | ---: |
| All arithmetic operations | 98 | 47,589 | 4,663,722 |
| Selected-array operations | 62 | 602 | 37,324 |

The local tables are specified by the evaluator, not saved as millions of
JSON rows. The source arithmetic and memory tables are separately saved in
full. These are intermediate classical counts, not the controller constant
or the microscopic dimension of a counterexample.

## Reversibility boundary

The [independent audit](audit/README.md) exhibits two distinct predecessors
that merge in a return sweep. Thus reversibility of the simulated source
would not make this microstep table reversible. The audit's separate
history-recording lift is a full-space permutation on multitape
configurations, but its heads can be far apart. One must preserve
reversibility in the physical routing itself before composing these facts
into a local quantum controller. The subsequent
[reversible one-head compiler](audit/REVERSIBLE_HEAD.md) now supplies that
repair, and [local permutation gates](reversible_gates.py) supply its gate
embedding. The noninjectivity warning still applies to using this
deterministic intermediate directly.

## Checks and provenance

The root [checks](checks.py) compare every source-step boundary with an
independent table interpreter for all 15 arithmetic entries, selected-array
GET/SET/APPEND/LEN, and random finite tables. They compare complete tape
contents, controls and head locations, and test rejection when capacity is
insufficient. The separate [audit](audit/single_head_checks.py) enumerates
918 source-step fixtures spanning one to three tapes and all virtual move
directions, and checks closure of 1,190 local control-symbol pairs.

Implemented by the coordinating Codex agent and reviewed independently by a
parallel Codex worker on 2026-10-28. No human review is claimed.
