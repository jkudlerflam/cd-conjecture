# A frozen program for the entire sparse-circuit generator

**Result.** There is now an executable, fixed **3,115-instruction word-machine
program** that takes the prime `p`, constructs its ordered projective-plane
incidence graph, constructs the two-port encoding circuit, lowers every gate
to nearest-neighbor operations, and emits exactly the symbols accepted by the
existing dispatcher. It has 34 functions, 27 primitive instruction families,
at most 38 integer registers per activation, and an acyclic call graph of
maximum depth nine. These counts are independent of `p`.

This is a concrete intermediate compiler target. **It is not yet a finite
alphabet, nearest-neighbor reversible Turing machine**, an initialized history
Hamiltonian, or a dimension estimate for the final counterexample. In
particular arithmetic on arbitrary-length words and addressed array accesses
have not been expanded into tape-head moves and digit operations. No final
controller factor is inferred from the 3,115 instruction addresses.

## What executes

- [`program.py`](program.py) is source in a small, closed language: integer
  arithmetic, indexed arrays of integers, conditionals, loops and calls to
  its 34 declared functions. It is read as syntax; its functions are never
  executed by the host Python interpreter.
- [`machine.py`](machine.py) compiles that syntax into
  [`bytecode.json`](bytecode.json). The frozen JSON is the actual executable
  artifact. The compiler rejects unrecognized expressions, imports, and
  external function calls.
- The same file contains the complete instruction interpreter. Runtime
  values are integers; array-valued source variables hold integer heap
  addresses. There are no graph, matrix, circuit, sorting, or synthesis
  instructions hidden in the interpreter.
- [`program_manifest.json`](program_manifest.json) records the instruction
  count, fixed instruction families and SHA-256 of the executable.

The instruction set is `PUSH, GETREG, SETREG, NEW, GET, SET, APPEND, LEN,
JMP, JZ, CALL, RET, POP, EMIT`, together with the integer operations
`ADD, SUB, MUL, DIV, MOD, SHL, SHR, AND, XOR, EQ, NE, LT, GT`.
`NEW` constructs a literal array of at most seven integer entries; `APPEND`
adds one entry. `EMIT` appends one pair of integers to a separate output
track. Conditional jumps inspect only classical integers, never the quantum
payload. The prime is an initialized classical input. Primality checking is
not part of the program: the model family is parameterized by a chosen prime.

The source faithfully reproduces the allocation and instruction order of
[`sparse_circuit.py`](../../controller_instantiation/circuit/sparse_circuit.py):
explicit Boolean comparison networks, branch preparation, charged-record
queries, sorting controls, Fock signs, inverse decoding, cleanup and both
port stages. The output mode selects an abstract circuit, local gate stream,
or literal dispatcher command stream. It never enumerates logical Fock states.

Abstract opcode integers 1 through 9 mean `X, CX, CCX, H, CH, R, QUERY,
CSWAPD, P`. Local opcode integers 1 through 17 have exactly the existing
`LOCAL_OPCODES` order. Dispatcher symbols 18 and 19 mean `MOVE_R` and
`MOVE_L`. Rotations still refer to the existing exact, `p`-dependent local
matrix families. The classical generator does not calculate their real
matrix entries or introduce approximate rotations.

## Reversal that has actually been implemented

The interpreter can append an explicit history record for each instruction:
the old instruction pointer, consumed stack words, number of produced words,
and overwritten register/array entries or returned activation record.
`undo()` is its literal inverse. Array allocation rolls back the next address;
array appends, output appends, calls, returns and all overwritten values are
restored. A returned activation record contains at most 38 registers.

This is a **journal-reversible word computation**. The journal is variable
length and consumes polynomial space. Finite local transition tables for
writing, reading and erasing it are not supplied here. The test driver also
chooses when to run forward and when to invoke the inverse. Autonomous
copy-output/reverse/handoff control remains to be composed. Consequently a
successful rollback test does not establish that illegal initial tapes have
positive energy, or that the resulting history graph has a unique accepted
component.

## Exact comparisons

[`checks.py`](checks.py) compares against the earlier implementation, which
is imported only by the tests and never by the generator. All seven checks
passed in [`results.json`](results.json):

| Check | Exact result |
| --- | --- |
| `p=2` incidence graph | All 21 edges in the same canonical order; full rollback |
| One-edge abstract circuit | All 906 operations and wire allocations match; full rollback |
| `p=2` abstract circuit | All 287,798 operations match; 11,027,317 VM steps |
| `p=3` abstract circuit | All 1,588,624 operations match; 60,109,113 VM steps |
| One-edge local stream | All 36,140 nearest-neighbor gates and addresses match |
| One-edge dispatcher stream | All 84,763 command symbols match |
| Output-producing prefix | After 500 emitted symbols, rollback restores the exact initial state, including empty output and heap |

The `p=2` and `p=3` **local** streams are too large for these tests and were
not expanded. The tests compare the complete abstract streams at those
primes, and the complete local and dispatcher streams only in the one-edge
fixture. Two matching implementations do not by themselves prove the
mathematical encoding theorem; earlier circuit tests check its semantics.

## Source-derived polynomial resource bound

Use the earlier notation

\[
L=N+4q,\qquad
w=\max\{1,\lceil\log_2q\rceil,\lceil\log_2\max(u,N)\rceil\},
\qquad A=L^2w.
\]

Every loop of this generator has at most `B=128 A` iterations on its promised
input. Here is the bound's basis, rather than an extrapolation from tests:

1. The abstract gate list has at most `128 A` entries, and the neutral wire
   count is at most `32 A`, by the existing
   [source-loop resource argument](../../controller_instantiation/circuit/RESOURCE_BOUND.md).
   The new generator has the same allocation and append structure.
2. All other construction lists are graph vertices/edges, words, flags,
   predicates, copies of gate blocks, or routing plans; their lengths are
   bounded by the gate-list bound, `L^2`, or `Lw`, hence by `128 A`.
3. Local lowering uses data width `D<=33 A`. Its gathering scans and cursor
   moves have at most `D` iterations per invocation. No loop iterates over
   the integer value of a packed flag: that value is only shifted and queried.
4. The prime-dependent graph loops have at most `q` iterations at each level.
   The integer-width loop has at most `w` iterations. Constant bounds, prefix
   counts and rotation indices also fit under `128 A`.

[`resource_analysis.py`](resource_analysis.py) recursively counts bytecode
cost from the restricted source, charging **both** conditional branches and
allowing **every** loop `B` iterations. The call graph is acyclic. For the
whole `main` function, including literal dispatcher output, it gives

\[
T\le193+4010B+7469B^2+5021B^3+1306B^4+424B^5+40B^6
\le18463B^6.
\]

The coefficients are saved in [`resource_bounds.json`](resource_bounds.json).
This deliberately loose bound counts word instructions, not digit operations
or local head transitions. It bounds all supported output modes.

At most seven heap words are allocated by any instruction, at most two words
are emitted, and a history record can be serialized in at most 64 words.
Heap block headers and the fixed stack/frame overhead fit within a total
`100 T + 40000` words. This retains all allocated arrays: no garbage collector
or space recycling is assumed.

The packed neutral flag integer can occupy `O(Lw)` bits. Thus the word size
must **not** be assumed to be only `O(log p)`. A sufficient signed word width
for values, addresses and bookkeeping is

\[
b=2A+\lceil\log_2(T+1)\rceil+64.
\]

The flag's largest bit position is below `2qw<=2A`; other numeric loop/data
indices are bounded by `128A`, heap addresses by `T`, and program addresses
are fixed constants. A delimiter-encoded serialization of the data, output
and journal therefore fits within `(100 T + 40000)(2b+4)` bits. These are
bounds on **serialized word-machine data**, not a claimed final physical
cell length: tape-level instruction routines, their workspace, verification,
clock sectors and endpoint readout still need to be specified and counted.

Since `A<=300 p^7` for projective-plane inputs, the bound is polynomial in
the numerical prime. In particular it supplies a concrete polynomial-space
input for subsequent digit-level compilation without any growing program
counter alphabet or an exponential Fock-space enumeration.

## Reproduction

```
python3 contrib/autonomous_controller/generator/machine.py
python3 contrib/autonomous_controller/generator/resource_analysis.py
python3 contrib/autonomous_controller/generator/checks.py
```

All three scripts use only the standard Python library. A test sink can hash
and compare output as it is emitted; the journalled model instead retains
the output track so that its inverse is defined. This testing optimization
is not a claim of zero output-tape space in the physical construction.

**Provenance.** Implemented with AI assistance by Codex on 2026-10-28.
No human review or completed counterexample is claimed.
