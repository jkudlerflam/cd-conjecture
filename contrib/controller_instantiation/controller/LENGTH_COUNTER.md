# Explicit blank-to-length reversible controller

Status: **implemented finite transition table, with a direct correctness
proof and finite exhaustive checks**. This is a second controller primitive,
separate from the preloaded-program dispatcher. It makes a concrete step
toward the latter's missing initialization contract.

The machine starts with static marks `LEFT BULK^(n-2) RIGHT`, both dynamic
bits zero at every site, and one head at LEFT in state `F:ADVANCE`. For
every `n>=4`, it reaches a deterministic midpoint with:

- the binary integer `n-1` on the BULK sites, little-endian from the left;
- a single frontier bit at RIGHT;
- its head back at LEFT.

It then exactly reverses the computation, returns both dynamic tracks to
zero, and repeats. There is no absorbing halt. The midpoint can be used as
a reversible subroutine handoff; the current standalone table immediately
uncomputes, so retaining the answer for a following computation needs an
explicit composition or reversible copy operation.

## 1. Actual alphabet, states, and collision table

The neutral alphabet is

\[
\Gamma=\{\mathrm{LEFT,BULK,RIGHT}\}\times\{0,1\}_{\rm counter}
\times\{0,1\}_{\rm frontier},\qquad B=12.
\]

There are eight forward control labels:

| Label | Next movement |
|---|---|
| ADVANCE | right |
| BACK | left |
| CARRY | right |
| RETURN | left |
| SCAN | right |
| FINAL_BACK | left |
| FINAL_CARRY | right |
| FINAL_RETURN | left |

Each has a reverse label with the opposite movement direction. Thus
`K=16`, and the direct-sum fresh/updated head implementation has

\[
\boxed{(1+2K)B=396}
\]

states. If this controller overlays an otherwise inert quantum record
of dimension `D(q)`, the product is `396 D(q)`. This counts only this
primitive; it is not the dimension of a composed initialized core.

[length_counter_table.json](length_counter_table.json) contains all 192
entries of its local collision permutation. There are 54 prescribed
entries, with the others supplied by an explicit permutation completion
within each static-mark sector. Consequently static marks are conserved
on every local input, including malformed sectors.

## 2. Forward motion and locally injective increment

The head first moves one site right in ADVANCE, sets the new frontier, and
enters BACK. It returns to LEFT, enters CARRY, increments the binary
counter, then returns to LEFT. SCAN walks right until it finds the
frontier, clears it, and enters ADVANCE for the next new site.

The binary increment uses these literal local rules, preserving the
frontier bit in each case:

\[
(\mathrm{CARRY},1)\mapsto(\mathrm{CARRY},0),\quad
(\mathrm{CARRY},0)\mapsto(\mathrm{RETURN},1),\quad
(\mathrm{RETURN},0)\mapsto(\mathrm{RETURN},0).
\]

The two ranges entering RETURN have different counter bits, so they do
not merge. RETURN crosses exactly the trailing ones that CARRY changed
to zeros. At LEFT it enters SCAN. The rule has no dependence on the
numerical counter value and no alphabet that grows with its width.

Other possible local merges are excluded explicitly:

- ADVANCE enters BACK with frontier `1`, while BACK's scan loop has
  frontier `0`.
- RETURN enters SCAN at LEFT, while SCAN loops only on BULK with frontier
  `0`.
- SCAN enters ADVANCE at BULK, distinct from the initial-state bounce at
  LEFT.
- FINAL phases have separate labels, so the terminal increment is not
  merged with an ordinary iteration.

On reaching RIGHT, ADVANCE sets its frontier and enters FINAL_BACK. The
machine returns to LEFT and performs one final increment in separately
tagged CARRY/RETURN states. The final arrival at LEFT is the midpoint.

## 3. Correctness and exact time

When the frontier first reaches physical site `j`, counting LEFT as site
zero, exactly `j` increments have occurred by the end of that iteration.
This follows by induction from the initially zero counter. Frontier
motion never revisits a processed site in the forward computation.
Therefore the final counter is `n-1`.

There are `n-2` binary storage sites. For every `n>=4`,

\[
n-1\le 2^{n-2}-1,
\]

so no increment overflows. The counter occupies only `O(log n)` nonzero
bits, plus the constant-track frontier mark. No travel history is stored.

If `v_2(j)` is the exponent of two dividing `j`, increment number `j`
travels right through `1+v_2(j)` bits and returns the same distance. For
`1<=j<=n-2`, a complete iteration takes

\[
1+2j+2(1+v_2(j))
\]

moves. The final iteration at RIGHT omits the SCAN and takes
`n+2+2v_2(n-1)` moves. Thus the exact forward time is

\[
\boxed{T(n)=n^2+3n-4-2\operatorname{popcount}(n-1).}
\]

The full forward/inverse cycle takes exactly `2T(n)` moves. The primitive
therefore uses quadratic time and a fixed alphabet, with no dependence
on an externally supplied program or precomputed value of `log n`.

## 4. Actual reverse execution, including endpoint bounces

For each nonterminal forward collision

\[
(q,a)\mapsto(q',b)
\]

the table also includes

\[
(\bar q',b)\mapsto(\bar q,a),
\]

where `bar q` moves in the opposite direction to `q`. After each reverse
move, this collision undoes the preceding forward collision. The final
forward RETURN at LEFT enters its reversed RETURN label directly. This
is exactly the inverse-ready state needed before retracing the last
leftward move.

At the other end, reversed ADVANCE arrives at LEFT and enters the original
forward ADVANCE label. These two endpoint transitions occupy ranges not
used by the other prescribed rules. They are not one-way bounces whose
targets have been incorrectly left fixed. The generator checks every
local domain and range for collisions before completing the permutation.

As with the dispatcher, the collision permutation gives an on-site
fresh/updated involution `F`; direction-tagged movement gives the fixed
two-site involution `G`. The quantum oracle is identity for every movement
here. An ordered `G` sweep and on-site `F` sweep execute one machine step.
The resulting local gates are defined on the entire Hilbert space.

## 5. Reproduction and scope

```sh
python3 contrib/controller_instantiation/controller/length_counter.py --output contrib/controller_instantiation/controller/length_counter_table.json
python3 contrib/controller_instantiation/controller/length_counter_checks.py
```

[length_counter_results.json](length_counter_results.json) records all 192
local pairs and inverses; complete blank-to-parameter-to-blank cycles for
every `n=4,...,128`, totaling 1461216 elementary steps; and 8149 arbitrary
nonoverflowing binary increments with different frontier positions.
Every tested cycle agrees with the exact time formula. At `n=128`, the
counter reaches 127 after 16750 moves and is blank again after 33500.

This establishes an explicit finite means of obtaining a size parameter
from a blank neutral tape. It does **not** yet extract the graph prime
from a chosen length convention, generate the sparse-encoding circuit,
load or verify the dispatcher's program, initialize all quantum input
sectors, transport the ports, or implement READY readout. Those are
separate remaining finite-controller contracts. The initial static marks
and unique head are also assumptions of this primitive, to be enforced by
the initialized-core Hamiltonian.

The table, implementation, and proof were developed with GPT-6-family
assistance and are supplied for independent checking.
