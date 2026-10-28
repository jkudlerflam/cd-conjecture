# An explicit finite neutral dispatcher

Status: **implemented and checked, conditional on a preloaded valid program**.
This replaces the unspecified interpreter transition table for one concrete
instruction-tape interface. It does **not** supply the classical algorithm
that computes this tape from the cell length, its microscopic initialization
and acceptance test, or the full counterexample Hamiltonian.

A separate [explicit blank-to-length controller](LENGTH_COUNTER.md) now
computes `n-1` from zeroed dynamic tracks, then exactly uncomputes it. Its
192-entry collision table, 16 control labels, and 396-state head alphabet
are implemented and checked. This is progress toward initialization; it
has not yet been composed with the circuit generator and dispatcher.

The executable [dispatcher.py](dispatcher.py) generates all local classical
collision entries. The checked-in [transition_table.json](transition_table.json)
is an actual finite table, not pseudocode. Its gates on the quantum factor are
named oracle matrices, whose definitions belong to the charge-preserving
gate implementation. No quantum value enters a classical control label.

## 1. Exactly what the interpreter accepts

An addressed circuit is a sequence `(opcode, left_address)` of two-site
gates. All addresses refer to a fixed data interval, with zero-based sites
`0,...,w-1`; `left_address < w-1`. One-site gates act on the left tensor
factor and the identity on its right neighbor.

`lower_gate_stream` turns this into a word over a fixed alphabet:

- `MOVE_R`, `MOVE_L`: move a neutral data-cursor mark one site;
- one symbol for each oracle: execute that two-site matrix at the cursor.

For example, `[(H,0),(CNOT,0),(CV,2)]` becomes
`[H,CNOT,MOVE_R,MOVE_R,CV]`. The lowerer is ordinary host-side code.
It is **not** a claim that the physical Hamiltonian has already generated or
verified this word.

If the addressed circuit has `g` gates and width `w`, this word has length
at most `gw`. Thus it preserves a polynomial-in-`p` space bound for the
sparse-encoding circuit family. It must not be used to expand the separate
physical-length-`n` transport walk into `O(n)` instructions: such an expansion
would create the space circularity already identified in ORACLE_SCHEDULER.

The physical tape has the following neutral program track:

```
START instruction_1 ... instruction_m END DATA ... DATA
```

Every site also has two neutral bits `P,D`, for the current program and
data cursor. `D` initially marks the first DATA site; `P` is initially zero
everywhere. The single neutral interpreter head starts at START, pointing
right. The quantum register is a separate factor on **every** site; it may
contain arbitrary superpositions. Program-region quantum records are not
touched by this dispatcher.

Each oracle step returns the head to its marked instruction, clears `P`,
and advances to the next instruction. At END the direction of program
execution reverses. The same word is traversed backward, with inverse
cursor moves and adjoint quantum gates. At START the initial head state,
cursor bits, program, and quantum input are all restored. There is no
absorbing halt. The END visit is a deterministic computation midpoint.

## 2. Explicit finite counts

Let `r` be the fixed number of oracle families. The program-symbol alphabet
consists of `r` oracles, two cursor moves, and START/END/DATA. Therefore

\[
B=4(r+5)
\]

counts **all** neutral tape symbols, including malformed combinations of
the two cursor bits. It does not depend on the program length, tape width,
species count, or matrix entries of the oracle families.

For each execution sign `+,-` the finite control contains:

- one ADV state;
- GO and RET for each of the `r+2` instructions;
- EXEC and BACK for each of the `r` oracle instructions;
- SHIFT for each of the two cursor moves.

Thus

\[
K=2[1+2(r+2)+2r+2]=8r+14.
\]

Every control label has **one fixed direction**. Fresh/updated head phases
therefore give exactly

\[
D_{\rm dispatcher}=(1+2K) B D(q)
\]

single-site states, where `D(q)` is the quantum-record dimension. This is a
small improvement on the generic `(1+4K)BD(q)` bound that allocates both
directions for every opcode/control label. The direction is still present
in the label and is essential to the ordered-sweep proof.

The checked-in menu has 17 oracle names:

```
SWAP X H Z CNOT CV CV_DAG TYPE SHIFT SHIFT_DAG QUERY LOAD
EXCHANGE ROT ROT_DAG CH INIT
```

This gives

\[
\boxed{B=88,\quad K=150,\quad D_{\rm dispatcher}=26488D(q).}
\]

The collision permutation has `BK=13200` entries, of which 1860 are
prescribed semantic rules. The others are an explicit bijective completion.
These numbers do **not** include the outer clock, initializer, verifier,
endpoint transport/readout, or any length-to-program computation. In
particular, 26488 is not the constant `C` of a completed counterexample.

## 3. The collision table and local unitaries

A control label determines the next movement direction and, only in an
EXEC state, its quantum oracle or adjoint. After that move, a collision acts
on the arrived-at neutral symbol and control label. The semantic rules are:

| Arriving control and neutral symbol | Collision output |
|---|---|
| ADV at instruction `k`, `P=D=0` | Set `P=1`, enter GO carrying `k` and sign |
| GO at any `P=D=0` symbol | Continue right |
| GO for oracle `k` at `DATA,D=1` | Enter EXEC; next right move applies `U_k` or `U_k†` |
| EXEC at unmarked DATA | Enter BACK, with identity left move |
| BACK at marked DATA | Enter RET |
| GO for cursor move at `DATA,D=1` | Clear `D`, enter the appropriate left/right SHIFT |
| SHIFT at unmarked DATA | Set `D=1`, enter RET |
| RET at any `P=D=0` symbol | Continue left |
| RET at matching instruction `k,P=1` | Clear `P`, erase opcode using `k`, enter ADV |
| ADV+ at END | Enter ADV− |
| ADV− at START | Enter ADV+ |

This is a partial **injection on local pairs**, not only a reversible
trajectory from one input. Its potentially overlapping ranges are separated
by actual local bits:

- Entry into GO has `P=1`, whereas its scan-loop range has `P=0`.
- Entry into RET has `D=1`, whereas its scan-loop range has `D=0`.
- Distinct instructions erased into ADV remain distinguished by their
  immutable program symbol.
- The endpoint-bounce ranges have START/END program symbols, disjoint from
  instruction-to-ADV ranges.
- Every other phase has a distinct control label.

The generator rejects any duplicate domain or range. It then completes
unmatched pairs to a permutation **separately for each read-only program
symbol**. First it fixes pairs unused in both domain and range; then it
matches the remaining domain/range lists in increasing index order. Thus
the program track is conserved on all inputs, including malformed ones.
There is no assumption that unused states can simply be fixed when that
would conflict with a prescribed image.

Call this finite permutation `C`. With ordinary quantum symbol `a`,
ordinary neutral symbol `gamma`, and a control label `k`, define the head
states `H^0_(k,gamma,a)` and `H^1_(k,gamma,a)`. The on-site unitary is

\[
F\lvert H^1_{k,\gamma,a}\rangle
=\lvert H^0_{C(k,\gamma),a}\rangle,
\qquad
F\lvert H^0_{k,\gamma,a}\rangle
=\lvert H^1_{C^{-1}(k,\gamma),a}\rangle.
\]

It fixes ordinary sites. Therefore `F^2=I` on the entire local Hilbert
space. It never changes a quantum record or species charge.

For a right-directed control `k`, the two-site unitary `G` has the block

\[
\lvert H^0_{k,\gamma,a},\delta,b\rangle
\mapsto
\sum_{a',b'} (U_k)_{a'b',ab}
\lvert\gamma,a',H^1_{k,\delta,b'}\rangle,
\]

with its adjoint on the reverse block. The left-directed block is the
corresponding reversed-head orientation. Other blocks are identity. The
controls' fixed directions, together with the fresh/updated phase, make all
blocks disjoint. `G` is an involution on the whole two-site space. If every
oracle conserves species charges, so do `G` and `F` on all inputs.

One ordered left-to-right `G` sweep followed by the on-site `F` sweep
executes exactly one of the machine steps above. The partial collision
rules are sufficient for every valid preloaded program; the arbitrary
completion is never used along those trajectories.

## 4. Cleanup and quantum-input independence

The dispatcher writes only the neutral cursor bits. Its control and travel
depend on immutable program symbols and those cursor bits. The data record
under the interpreter head is never copied to finite classical control.
Consequently, for every time `t`, one has exactly

\[
\mathcal U^t J_{c_0}|\psi\rangle=J_{c_t}V_t|\psi\rangle,
\]

where `c_t`, the running time, and the END midpoint are independent of the
entire quantum input, including superposed species and branch registers.
At the midpoint `V_t` is the requested circuit. At the end of the cycle it
is the identity, and all neutral tracks and the head equal their initial
values. This cleanup is by literal adjoint execution, not measurement or
discarding a history register.

Each instruction takes at most `2n+4` elementary moves on an `n`-site tape.
The full cycle is bounded conservatively by `2(m+1)(2n+6)` moves. There are
no stored travel histories. The classical program costs `m` sites and the
two marker bits already counted per site.

## 5. Verification and remaining obligations

Run from the repository root:

```sh
python3 contrib/controller_instantiation/controller/dispatcher.py --output contrib/controller_instantiation/controller/transition_table.json
python3 contrib/controller_instantiation/controller/checks.py
python3 contrib/controller_instantiation/controller/integration_check.py
```

[results.json](results.json) records:

- all 13200 collision inputs and inverses, including read-only program
  conservation outside the valid sector;
- all 4200 eligible one-head ordered sweeps for widths 2 through 8 and all
  150 control labels;
- 2271 exhaustive short legal forward/inverse programs;
- 100 random addressed scripts using the full 17-family menu;
- 20 coherent four-qubit input trials, with H, CNOT, and a non-self-adjoint
  controlled square-root-X gate; direct forward-circuit comparison and
  complete cleanup error are below `3e-16`;
- equality of the **entire classical trajectory** for those quantum inputs.

The independent [integration_check.py](integration_check.py) connects the
actual sparse-isometry local gate generator to this instruction alphabet.
For its one-edge fixture it reconstructs all 36140 gate addresses from the
lowered program, with no mismatch. The resulting concrete program has 84763
instructions, needs 85053 physical sites including its 288-site DATA
interval, and has a conservative full-cycle bound of 28838747136 moves.
[integration_results.json](integration_results.json) records its SHA-256
fingerprint and exact opcode counts. This long complete trajectory was
**not** simulated; the gate-stream/operand interface was checked, while
the circuit's coherent test and the dispatcher's small-program tests are
separate checks. This fixture is not a violating example.

The remaining microscopic obligations are explicit:

1. Generate or verify the required program from `p`/physical length with a
   fixed reversible neutral machine, and include its alphabet and states
   in the final dimension. The host-side lowerer does not do this. Indeed,
   the table here conserves its program symbol on **every** input, so it
   cannot generate a nonblank program from an initially blank program
   track. A separate generator/controller extension or an energetic
   verification of a preloaded program is mathematically required.
2. Force the initialized accepted sector, including the exact program and
   cursor pattern, with the clock/verifier construction. Allowing arbitrary
   program words introduces unwanted sectors and is not acceptable.
3. Attach compact endpoint-transport routines. The DATA interval here is
   right of the program region, so its leftmost quantum record is not yet
   at the physical left boundary.
4. Implement READY/readout and verify every elementary outer-clock vertex,
   including intermediate scan positions, not merely completed interpreter
   steps. The deterministic END checkpoint alone is not the required
   endpoint-readout projector.
5. Compose the initialized-core spectral gap and the finite-coupling
   infrared argument. This interpreter neither proves these claims nor
   produces a numerical violating `d,p,s`.

This module was developed with GPT-6-family assistance. Its complete table
and reproducible checks are supplied so that the local reversibility and
resource claims can be audited without invoking a universal-computation
existence theorem. An independent agent reviewed the collision-range
separation and emphasized the unavoidable preloaded-program restriction.
