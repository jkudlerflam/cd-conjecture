# Independent review of the one-head reversible history compiler

Reviewer: the arithmetic implementation subagent, separate from the author of
`audit/reversible_head.py`, 2026-10-28. This is a code/mathematical review, not
human review and not a proof of the full initialized Hamiltonian construction.

## Local criterion checked

For a deterministic one-head transition

`(q, a) -> (r, b, d)`,

the following two properties suffice for a partial global inverse on arbitrary
configurations of the declared finite control/tape alphabets:

1. The incoming displacement `d` depends only on the target control `r`.
2. The map `(q,a) -> (r,b)` is injective on defined rules.

Indeed, from a target configuration with head at `x` and control `r`, the unique
possible predecessor head is `x-d(r)`. Reading that cell gives `b`; property 2
then determines at most one predecessor `(q,a)`. This argument does not assume
one work marker, a well-formed history, or a correctly initialized tape. Those
promises are needed to simulate the desired computation, not for the inverse
of a defined transition.

The review traced each control phase in `reversible_head.py`:

| Target phase | Incoming displacement | Why the written symbol has at most one preimage |
| --- | ---: | --- |
| `F(r)` | 0 | Only `RETURN(r)` on a marked work cell enters it. |
| `MOVE(t)` | Source displacement in event `t` | Event `t` identifies the old source control and read symbol; the write is fixed. |
| `DEPART(t)` | 0 | Nonzero-displacement events enter from `MOVE(t)` and retain the scanned work symbol; zero-displacement events enter only from their unique `F` rule. These event classes are disjoint. |
| `SCAN(t)` | +1 | Entry writes a **marked** work cell; scanning self-loops write only unmarked work, ordinary history, separator or blank symbols. |
| `APPEND(r)` | +1 | The written history record contains the complete source transition `t`; its target is checked to equal `r`. |
| `RETSTART(r)` | 0 | Only `APPEND(r)` on blank enters it, writing the frontier. |
| `RETURN(r)` | -1 | Entry writes the frontier; scanning self-loops exclude the frontier and marked work cells. |

The distinction between the marked entry symbols and the allowed sweep symbols
is essential. Allowing `SCAN` to pass marked work cells, or `RETURN` to pass the
frontier, would merge an entry transition with a sweep transition and break
injectivity. The reviewed implementation excludes both.

No collision or motion rule consults a host history stack, a remote scanned
symbol or a runtime integer value. The finite event is carried in the control
on the outward trip, deposited in a physical history cell, and recovered from
that same cell by the inverse trip. Original irreversible source rules can
therefore share a target: their distinct event records keep their inverse
histories separate.

## Finite-control bookkeeping finding

The first reviewed version of `materialize()` created `F(q)` for every source
state, but created `APPEND(q)`, `RETSTART(q)` and `RETURN(q)` only for states
occurring as transition targets. Nevertheless `inverse_rule(F(q), marked_work)`
returned `RETURN(q)` even when `q` was a non-target initial state. That inverse
control was outside the declared materialized control set.

An independent two-symbol erasing-source fixture exposed two such extra inverse
entries. All defined forward entries of that fixture were already injective and
had the correct incoming displacement. The issue was the claim of a complete
partial inverse on the declared alphabet, not agreement on initialized runs.

The author repaired this by including all three history/return phases for
**every** source state. The existing bound `4|Q|+3m` already allows this choice. An equally
valid but more cumbersome repair would guard the forward and inverse routines
against these undeclared controls. The independent check should compare the
domain of `inverse_rule` against the image of the materialized forward table,
including image points that have no predecessor, rather than checking only
forward-then-inverse on defined forward entries.

After the repair, an independent exhaustive check compared the inverse domain
with the forward image on two fixtures. They had respectively 154 and 476
control/symbol pairs and 44 and 156 defined rules. Together the fixtures cover
irreversible erasure, all three motion directions, and an initial state absent
from the source transition targets. Every pair passed, including pairs with no
predecessor. Each fixture's finite collision map was also completed by explicitly
matching unused domain and range pairs, and the resulting complete map was
checked to be a permutation. These fixture checks support, rather than replace,
the general phase analysis above.

## Scope of the conclusion

The phase analysis establishes the reversible partial-TM criterion for the
repaired finite-control bookkeeping. The source-step history append
uses one record cell, and the outward/return sweeps are linear in the allocated
work-plus-history interval. Its length is set by the **digit-level** source step
count, not by the earlier word-RAM instruction count.

A finite injective partial collision map can be completed to a permutation by
matching unused domain and range pairs; the fixed-direction head motion is also
invertible on an unbounded tape. This observation is not itself the remaining
Hamiltonian construction. Finite endpoints, a concrete choice of completion,
autonomous turnaround/readout, preparation of the work/history format, and
exclusion or energy penalties for unwanted sectors still require composition
with the core construction. In particular, a halting partial transition cannot
simply be declared an absorbing fixed point of a reversible dynamics that has
already entered it from a different state.
