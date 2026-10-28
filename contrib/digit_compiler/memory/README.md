# Serial, fixed-alphabet array memory

This component replaces addressed lookup with an explicit digit-level tape
computation. The same **97 active control states and 207 transition rules** work
for every address width, array length, and word width. The tape alphabet is
`_ ^ 0 1 # $`, with `_` blank. The control state chooses one of five tapes; each
rule reads/writes only that tape's current cell and moves its head by `-1`, `0`,
or `1`. There is no host array indexing, integer parsing, or integer arithmetic
inside the transition interpreter.

The four entries implement `GET`, variable-width functional `SET`, `APPEND`,
and `LEN`. This is an array-storage primitive, **not yet a complete VM heap**.
The caller must supply the selected array on its input tape. Heap allocation,
heap-reference resolution, alias-preserving replacement of a modified block,
and a spatial single-head implementation remain to be composed.

## Input and output contract

The five data tapes are:

| Tape | Initial content | Purpose |
| --- | --- | --- |
| `M` | `^word#word#...$` | Immutable array, with nonempty binary words |
| `A` | `^address$` | Unsigned little-endian address bits |
| `V` | `^word$` | Replacement/appended word |
| `C` | entirely blank | Counter work tape |
| `R` | entirely blank | Result work tape |

All heads begin at coordinate zero. `^$` is the empty array. A binary word is
nonempty; its bits are opaque to memory. For example `^10#001#$` contains the
two words `(1,0)` and `(0,0,1)`. Array indices start at zero. The address may
have redundant high zero bits. Negative addresses have no representation.

`GET` returns `^word$`. `SET` returns the entire modified array on `R`, so the
replacement word may have a different width. `APPEND` returns a new array with
`V` at the end. `LEN` returns the canonical little-endian binary array length
inside `^...$`, including `^0$` for an empty array. `APPEND` and `LEN` do not
inspect the address. `GET` and `LEN` do not inspect the value. All successful
base executions return the five heads to coordinate zero.

An out-of-range address reaches a missing rule at the array end and rejects;
it cannot silently wrap. The rules also reject malformed symbols when those
symbols are encountered. **This is not a whole-input validator:** for example,
`GET(0)` need not inspect later words. Proper initial tape framing is a caller
invariant, to be supplied by the generator's initialization verifier.

Signed arithmetic values can be stored as an opaque word containing a sign bit
followed by little-endian magnitude bits. Connecting that representation to
the ALU's `P/M` sign delimiter requires a finite symbol-conversion adapter;
the memory routine never interprets or changes the sign convention.

## Why the address lookup is actually serial

`GET` and `SET` first copy the address one bit at a time from `A` to `C`.
Suppose the input address is `a`. At the beginning of lookup iteration `j`:

1. the memory head scans the first bit of word `j`, or the final `$`;
2. the counter contains the fixed-width binary expansion of `a-j`;
3. for `SET`, the output already contains exactly the first `j` input words.

The zero test traverses every counter digit and returns to its `^` delimiter.
If the counter is nonzero, the memory head traverses one complete word to the
next `#`, and a borrow-carry loop subtracts one from `C`. The borrow loop changes
each low `0` to `1` until it changes the first `1` to `0`; a reverse sweep then
returns to `^`. It is called only for a nonzero counter, so it cannot underflow.
These rules preserve the invariant and advance `j` by exactly one.

When the counter is zero, `GET` copies that word. `SET` traverses that word,
copies the replacement from `V`, then serially copies the remaining input
array. Encountering `$` instead of the selected/skipped word rejects. No word
number is ever converted to a Python index or physical head coordinate.

`LEN` similarly scans each input word, incrementing a binary counter at each
`#`. Its carry loop replaces low `1` digits by `0`; if it reaches the terminal
`$`, it writes a new `1` and moves the terminal right one cell. The final counter
is copied to the result tape. `APPEND` is just serial input copying followed by
serial value copying.

Consequently a valid array with `n` words, total serialized length `B`, address
width `w`, and replacement width `v` takes `O(B + (n+1)w + v)` base steps for
`GET`/`SET`, `O(B+v)` for `APPEND`, and `O(B+n log(n+1))` for the stated simple
`LEN` implementation. The bound depends on the length of the array, not on the
possibly enormous integer encoded by an out-of-range address. These are bounds
on multitape transitions; routing between spatially separated tape heads is
not included.

## Exact rollback and clean output

The implementation has a separate history tape. Each forward transition writes
one symbol identifying its rule, then performs that rule. There are **208
history symbols including blank**, independent of all inputs. A symbol contains
only one of 207 fixed transition IDs; it does not contain an address, word,
Python object, or unbounded integer. Reading the preceding history symbol
determines the old control state, active tape, old/new symbol, and head move.
Reversing that head move and restoring the old symbol recovers the preceding
configuration exactly. Rejection does not destroy the current state, so a
rejected computation can also be rolled back.

After a successful computation, the **8-state, 20-rule** machine in
`copy_table.json` copies the framed result to a fresh preserved-output tape `P`,
then returns the `R` and `P` heads to their initial delimiters. Reversing the
base computation erases `R`, `C`, and the history and restores all input heads.
Only the preserved output remains. `clean_run` checks this complete initial
configuration equality, not just the logical output.

The interpreter's one recorded base step is a reversible simulation macro;
`machine.py` does not claim it is already one nearest-neighbor physical bond
transition. The independent [compiler audit](../audit/) supplies a local
collision-permutation lift of these finite rules. A spatial layout and routing
compiler, whole-Hamiltonian legality constraints, and autonomous composition
of compute/copy/reverse still need to be supplied. The Python driver chooses
these phases at present.

## Files and reproduction

- `machine.py` generates the rules and executes only current-cell transitions.
- `transition_table.json` is the complete base table, with machine-readable
  tape indices, entry/halting states, alphabet, and transition metadata.
- `copy_table.json` is the complete preserved-output copy table.
- `checks.py` compares against ordinary reference array semantics, checks
  exact reversal of every prefix of selected computations, and regenerates
  the stored tables for equality.
- `results.json` records the latest successful check results.

Run from the repository root:

```sh
python3 contrib/digit_compiler/memory/machine.py
python3 contrib/digit_compiler/memory/checks.py
```

The checked suite contains **4,158 semantic cases**, **2,714,328 base forward,
reverse, and copy transitions**, **620 exact consecutive-prefix rollback
checks**, and **10 malformed-input rejections**. It includes variable-width
replacements, empty arrays, binary-counter carries through length 256,
out-of-range addresses, redundant high address zeros, and address widths up
to **513 bits**. An independent generic interpreter should also replay the
exported tables; this supplies a cross-check against implementation-specific
interpreter errors when the component is integrated.

These files were produced by a GPT-6-family worker for the c-d conjecture
project on 2026-10-28. Automated checks and agent review do not constitute
human review or a counterexample proof.
