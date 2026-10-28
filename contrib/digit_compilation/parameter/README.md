# Compiled finite parameter reader

Status: **explicit finite local table and whole-space gates, with checked
parameter output and autonomous reverse cleanup**. This instantiates the
finite walk left open in
[`PARAMETER_ORACLE.md`](../../autonomous_controller/audit/PARAMETER_ORACLE.md).
It is one controller primitive, not a compiled circuit generator or a
counterexample to the conjecture.

## Interface and counted resources

For `n >= m >= 2`, begin with static marks `LEFT BULK^(n-2) RIGHT`, zero
digit and end tracks, one fresh `F:LAUNCH` head at RIGHT, and the existing
indexed helper `H_(0,0)` at RIGHT. Every other Sigma tensor factor may be
arbitrary, including loaded `T` records and malformed `H` records. The
full sparse-code family has `m >= 10`, so the `m >= 2` condition imposes no
additional restriction on that family.

After exactly **`2m` controller sweeps**, the head is back at RIGHT in
`R:BACK`, the helper is again `H_(0,0)`, and all other Sigma factors are
back at their original sites. The digit track contains `p`, little-endian
from RIGHT, with zero padding to width `m`. The end track has its only
one at distance `m-1` left of RIGHT. After another **`2m` sweeps**, the
entire initialized configuration is restored. No host-side handoff,
counter comparison, stored trajectory, or program tape drives this cycle.

There are eight directional control labels and twelve neutral symbols:

\[
K=8,\qquad B=3\cdot2\cdot2=12.
\]

The standard no-head/fresh/updated implementation used by the existing
controller therefore has the explicitly counted local dimension

\[
\boxed{(1+2K)B D(q)=204 D(q)},\qquad
D(q)=q+3+2(q+1)m.
\]

This standalone count includes its three static marks, both output tracks,
and the complete head overlay. It does not include a sweeping clock,
graph/circuit generator, format verifier, readout, or boundary enforcement.
Shared tracks may reduce a composed implementation's count, but no such
reduction is used here.

The actual neutral table is saved in
[`transition_table.json`](transition_table.json): **96 entries**, of which
52 are prescribed. Its values and control labels are independent of
`p`, `q`, `m`, and `n`.

## Literal local collision and movement gates

Write `W` for the existing six-gate word

```
PARAM, COPY_DIGIT, PARAM, END, COPY_END, END.
```

These are the complete-alphabet matrices in the original parameter oracle.
For a valid helper `H_(j,0)`, `W` writes `p_j` and `[j=m-1]` into the
physical digit/end tracks and restores the helper bit to zero. For a
malformed `H_(j,c)`, it XORs `c XOR p_j` and `c XOR [j=m-1]` instead;
for non-H records it is identity. Thus it is an explicit permutation on
every local state. Its inverse is the reverse six-gate word.

The complete on-site collision is the fixed composition

\[
P=W_R^{-1}\,C\,W_F.
\]

`W_F` applies `W` only under input controls `F:START` at RIGHT or `F:OUT`
away from RIGHT. `W_R^{-1}` applies `W^{-1}` only under output controls
`R:START` at RIGHT or `R:OUT` away from RIGHT. Each control condition is
unchanged by the corresponding gate, so both controlled gates are
permutations on the full alphabet. There is no growing operation hidden
in this factorization: each active controlled word has six given gates.

`C` preserves the static mark, digit and end. Its prescribed control
transitions are below; all digit values are included.

| Mark | End | Input | Output |
|---|---:|---|---|
| LEFT or BULK | either | F:LAUNCH | F:START |
| LEFT or BULK | either | R:START | R:LAUNCH |
| LEFT or BULK | 0 | F:OUT | F:OUT |
| LEFT or BULK | 1 | F:OUT | F:BACK |
| LEFT or BULK | 0 | F:BACK | F:BACK |
| LEFT or BULK | 0 | R:OUT | R:OUT |
| LEFT or BULK | 1 | R:BACK | R:OUT |
| LEFT or BULK | 0 | R:BACK | R:BACK |
| RIGHT | 0 | F:START | F:OUT |
| RIGHT | 1 | F:START | F:BACK |
| RIGHT | 0 | R:OUT | R:START |
| RIGHT | 1 | R:BACK | R:START |
| RIGHT | 0 | F:BACK | R:BACK |
| RIGHT | 0 | R:LAUNCH | F:LAUNCH |

These pairs have disjoint domains and ranges. The remaining entries are
completed separately within each `(mark,digit,end)` sector: first fix the
unused common domain/range, then pair the remaining indices in increasing
order. This is an actual finite bijection, including every malformed
neutral sector; the saved table makes its completion directly inspectable.

The directional controls and two-site movement words are:

| Control | Direction | Quantum movement word |
|---|---:|---|
| F:LAUNCH | left | identity |
| F:START | right | identity |
| F:OUT | left | SHIFT departure, Sigma-only SWAP |
| F:BACK | right | Sigma-only SWAP, SHIFT_DAG arrival |
| R:LAUNCH | right | identity |
| R:START | left | identity |
| R:OUT | right | Sigma-only SWAP, SHIFT_DAG arrival |
| R:BACK | left | SHIFT departure, Sigma-only SWAP |

`SHIFT` is the original cyclic index shift on both H and T records. It
does not change the species carried by T. Sigma-only SWAP leaves the
physical digit/end tracks and static marks in place.

The full on-site involution `F` maps an updated head through `P` to a
fresh head, and a fresh head through `P^{-1}` to an updated head. No-head
states are fixed. The full two-site involution `G` moves a correctly
oriented fresh head with its movement word and moves its updated partner
back with the exact inverse word. All zero-head, two-head and other
orientation sectors are fixed. Consequently `F^2=G^2=I` on the full
Hilbert spaces. `F` preserves each one-site species charge; `G` preserves
each two-site total species charge, including malformed H/T sectors.

One ordered left-to-right G sweep followed by a right-to-left F sweep
executes one head move and collision. The updated-head tag prevents a
second move within the same sweep. This is the existing controller's
fresh/updated mechanism, with explicit new finite control entries.
Each sweep applies `n-1` G gates and `n` F gates, so the full cycle has
`4m(2n-1)` scheduled local gate applications, including identity actions.

## Why output and cleanup are independent of quantum input

The first two identity moves go one site left, then return to RIGHT and
write digit zero. They supply a reversible entrance/exit without a
stationary head or a one-way halt. For each subsequent outward move,
`SHIFT` acts on the tracked helper tensor factor and the adjacent swap
moves that factor left. Inductively it is exactly `H_(j,0)` at the site
distance `j` from RIGHT. The local collision writes the prescribed bit
and delimiter. It enters BACK exactly when `j=m-1`.

Thus every stopping decision depends on the deterministically generated
trusted end bit. An arbitrary H-like record elsewhere is never selected
by its shape: the helper is followed through prescribed adjacent swaps.
In particular no untrusted quantum record controls the head schedule.

On the rightward leg, the movement word is the inverse of the matching
outward word. Each displaced Sigma factor returns to its original site;
the helper index decreases to zero. The digit/end tracks remain at their
physical sites. The collision at RIGHT changes to `R:BACK`, giving the
midpoint after `2 + (m-1) + (m-1) = 2m` sweeps.

The reverse labels invert these steps. They first carry the helper back
to the delimiter, erase that site's two outputs with `W^{-1}`, then erase
each preceding site's outputs while returning to RIGHT. The two identity
launch moves are reversed, and the final RIGHT collision returns to
`F:LAUNCH`. The complete cycle is `4m` sweeps. This proof is a tensor-factor
identity on the entire promised-helper input subspace, so it also applies
to superpositions entangled with arbitrary external systems.

The forward quantum word, excluding the two identity launch moves, is
exactly the prior `10m-4` primitive-gate extraction word. The classical
finite head implementation is what is added here. The claimed initialized
trajectory requires the stated pins; no input-independent trajectory is
claimed for uninitialized neutral tracks or a missing helper.

## Reproduction and remaining interface obligations

```
python3 contrib/digit_compilation/parameter/walk.py
python3 contrib/digit_compilation/parameter/checks.py
```

The checks verify the saved 96-entry table, exhaust full on-site F columns
for the actual `p=2` and `p=3` family alphabets, exhaust every Sigma pair
and one-head selector in a small complete G alphabet, test the zero/two
head sectors, exhaust all untrusted three-site payloads under the helper
pin, and run larger adversarial cases. Actual micro-sweeps are compared
with the direct head step. Output, exact payload restoration, identical
schedules across payloads and full reverse cleanup are all checked.

The midpoint is a reversible subroutine interface, not a halt: its next
step begins cleanup. A generator must explicitly compose with that state
to consume the bits, then return for uncomputation. Enforcing the unique
RIGHT head, blank tracks and helper pin in a complete initialized core is
also a composition obligation. The bounded table above removes the
parameter-reader compilation gap only; graph/circuit generation and the
finite-coupling infrared argument remain separate.

Developed and checked with GPT-6-family assistance. No independent human
review is claimed.
