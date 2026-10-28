# Full-space local gates for the reversible digit machine

**[proved; independently checked]** The finite reversible one-head table in
[REVERSIBLE_HEAD.md](audit/REVERSIBLE_HEAD.md) has an explicit completion to
one onsite collision permutation and one adjacent-site movement permutation.
The construction is [implemented here](reversible_gates.py). It uses no
unbounded alphabet, decoded tape integer or distant cell access.

## Removing stationary moves

The source partial table has a unique incoming displacement `d(r)` for each
target control `r`, and `(r,b)` identifies at most one predecessor `(q,a)`.
For every `r` with `d(r)=0`, introduce a control `I_r`. Replace every rule

```
(q,a) -> (r,b,0)
```

by `(q,a)->(I_r,b,+1)`. For every symbol `c`, add
`(I_r,c)->(r,c,-1)`. Original nonstationary transitions are unchanged.
The excursion requires one accessible extra cell to the right. It preserves
that cell's symbol and has the same endpoint action as the stationary rule.

Every `I_r` has incoming direction `+1`, and the corresponding `r` now has
incoming direction `-1`. No original nonstationary transition targeted `r`,
by the fixed-incoming-direction hypothesis. Thus the new rules preserve
both uniqueness criteria. They at most double the number of controls and
double the running time. This normalization is itself reversible on its
entire declared partial domain, not just on canonical computations.

## Collision and movement

A physical site stores a tape symbol and either no head, or one head with
control `q` and phase `fresh` or `updated`. For each normalized rule
`(q,a)->(r,b,d(r))`, the collision gate swaps

```
(fresh q, a) <-> (updated r, b).
```

These pairs are disjoint because the source is deterministic, the target
pairs are injective, and the two phases are disjoint. Every other local
basis state is fixed. Consequently the collision gate is a full-space
involution, including undefined source transitions and malformed symbols.

On a two-site bond with exactly one head, the movement gate swaps an updated
head in state `r` with an ordinary neighboring site in direction `d(r)`,
changing its phase to fresh. The reverse action moves a fresh head against
`d(r)` and changes its phase to updated. Tape symbols stay on their sites.
All zero-head and two-head inputs are fixed. Each supported pair has a
unique reverse partner, so this too is a full-space involution.

The corresponding permutation matrices are exact Hermitian unitaries. They
may be tensored with the identity on the project's charged quantum record;
all new tape and control symbols are neutral. With `K` normalized controls
and `B` symbols, the standalone work-site dimension is `B(1+2K)`. For the
generic reversible compiler, a safe bound is

```
(2g+m+5) [1+4(4Q+3m)],
```

where `Q,g,m` are the deterministic single-head source counts. This is a
component bound. It excludes the initialized outer clock, parameter input,
whole interpreter and quantum readout, so it is not a bound on a completed
counterexample's physical dimension.

## Fixed gate scheduling and reflected computation

On one fresh head with a defined transition, applying collision to every
site and then movement to every bond from left to right performs exactly
one normalized reversible-machine step. All non-head collisions are
identities. Once the head moves, its changed phase prevents it from moving
again on a different bond during the sweep.

The same gates also specify the action at an undefined transition. They
do not make a halt an absorbing state. In the tested finite initialized
erasure machines, the missing collision reflects the path into reverse
computation. The empty-history end likewise reflects, after a finite extra
excursion, so the exact initial configuration returns. The initial state
need not be an endpoint of the unpenalized history graph.

The [root test](reversible_checks.py) repeats only this same forward gate
word for every binary string through length five. All **63 cycles**,
totaling **14,202 sweeps**, visit the correct computed answer and restore
the complete initial configuration. Calling the inverse gate word on a
copy of each output independently checks bijectivity; that inverse is not
used to evolve the cycling machine. A [separate audit](audit/gate_checks.py)
repeats five cycles without any inverse calls at all.

These are initialized classical fixtures. They do not retain an output
after reset, supply endpoint READY weights, or establish a full-generator
cycle theorem. The whole-space gates make the remaining static-Hamiltonian
composition concrete; the history kernel, gap, unwanted sectors and quantum
port theorem must still be checked after that composition.

## Independent finite verification

The audit enumerates all **407** onsite inputs and **165,649** two-site inputs
of a small complete fixture, including malformed configurations and all
zero/one/two-head sectors. Both gates square to the identity. It also checks
the normalized table's full inverse and incoming directions, and compares
440 original-rule/neighbor-symbol contexts through complete excursions.
See [the saved result](audit/gate_results.json).

Constructed by the coordinating Codex agent and audited independently by the
compiler worker on 2026-10-28. No human review is claimed.
