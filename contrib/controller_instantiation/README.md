# From the microscopic construction argument to executable local rules

**Status: implemented finite components and a repaired initialization
contract; the complete autonomous high-channel Hamiltonian is not yet
instantiated.** No new c-d violation or measured central charge is claimed.

This work answers Pace's request to make concrete progress on the missing
microscopic Hamiltonian. It starts from commit `5514132`. Author: OpenAI
Codex, GPT-6 family, with three parallel circuit/controller/audit workers,
2026-10-28 UTC. The workers exchanged interfaces and independently checked
the new local dispatcher and physical hopping projection. No human review
is claimed. The distinction between exact finite identities, floating-point
checks, and unproved thermodynamic conclusions is retained throughout.

## 1. The circuit is now an executable gate sequence

The [circuit generator](circuit/README.md) explicitly implements the sparse
two-port encoding. There is no remaining unspecified Boolean function
inside its quantum transducer: the graph lookup, canonical charged words,
Fock signs, routing controls, and inverse decoder are expanded into
NOT/CNOT/Toffoli networks and the previously specified charged operations.
The entire sequence lowers to 17 fixed two-site opcode families.

For the first projective plane, p=2, the exact generated counts are:

| Quantity | Value |
| --- | ---: |
| Species | 14 |
| Logical fermions | 21 |
| Circuit record dimension | 317 |
| Circuit sites | 81,342 |
| Abstract gates | 287,798 |
| Nearest-neighbor gates after localization | 1,213,761,392 |

The billion-gate count comes from the actual routing algorithm; that stream
was counted without being materialized or simulated. Its large workspace
and time do not themselves increase the circuit's local alphabet, but they
make the need for a compact autonomous controller concrete.

The one-edge fixture is small enough for complete coherent execution of
its **36,140 local gates on 288 sites** using a sparse state representation.
Both logical input columns produce the desired two-port isometry and return
all auxiliary sites to blank, with amplitude discrepancy below 1.2e-16.
Small graphs also reproduce both port matrices, including the many-mode
right parity dressing. Exact integer tests cover 2,592 p=2 branch cases.
The p=2 instance is a construction check, not a violating parameter choice.

The [source-level resource bound](circuit/RESOURCE_BOUND.md) also replaces
the circuit/dispatcher layer's unspecified polynomial constants. With
L=N+4q, w the maximum species/branch bit width, and A=L^2w, the actual
generator uses at most 32A neutral records, 128A abstract gates and
70,000A^2 local gates. The preloaded dispatcher fits within 2,310,035A^3
sites. For projective planes this is below 2^50 p^21. This is an explicit
capacity bound for that layer only; it excludes autonomous generation,
verification, endpoint transport and the full core. Checks include exact
generated counts for p=3, without expanding its 32.9-billion-gate local stream.

## 2. A finite dispatcher table replaces one existence claim

The [controller directory](controller/README.md) contains an actual
13,200-entry local collision permutation. It reads an instruction tape,
moves a neutral data cursor, calls the specified quantum gates, and then
executes their adjoints in reverse. Its trajectory never depends on a
quantum value. The complete neutral configuration returns to its start.

For the 17 preparation opcodes its exact resource count is

\[
B=88,\qquad K=150,\qquad
D_{\rm dispatcher}=26{,}488\,D_{\rm circuit}.
\]

Every local transition has a unique predecessor, including all malformed
classical states after an explicit permutation completion. The head holds
only its current quantum record; no hidden second species register is
present. All 13,200 collisions and 4,200 one-head sweeps were checked, as
were 2,271 short program cycles, 100 full-menu addressed scripts, and
20 coherent quantum trajectories with non-self-adjoint gates.

The actual one-edge gate sequence lowers to 84,763 dispatcher instructions
on 85,053 sites. Reconstructing every operand address from that program
matches the gate stream exactly. The conservative bound of 28,838,747,136
steps for its full forward/inverse cycle was calculated, not simulated.

**This count is for a preloaded dispatcher, not the complete core.** Its
program symbols are conserved on every input, so the table cannot generate
its own nonblank program from blank initialization. That is a specific
missing operation, not an issue that can be solved by choosing a smaller
perturbation. The controller's TYPE opcode is the bit/helper conversion;
the verifier's TYPE-query families are additional operations and are not
silently included in this 17-family count.

### A blank-to-length controller is also explicit

The separate [length counter](controller/LENGTH_COUNTER.md) starts with
zero dynamic bits and the static endpoint marks. Its 192-entry collision
permutation computes the physical length minus one, then autonomously
uncomputes and repeats. It uses 12 neutral symbols and 16 direction-fixed
controls, for a standalone head-overlay factor of 396. Its exact forward
time is n^2+3n-4-2 popcount(n-1), and no travel history is stored.

All 125 lengths from four through 128 pass complete cycles, totaling
1,461,216 elementary steps; 8,149 additional binary increments pass.
This supplies the first parameter-reading primitive from blank. It is
not yet linked to prime decoding, program generation, or the dispatcher:
the standalone table immediately uncomputes its result at the midpoint,
so a reversible handoff or copy must be incorporated in the composition.

## 3. A real input-format ambiguity was found and repaired

The [adversarial audit](audit/README.md) shows that TYPE and species queries
alone cannot distinguish an ordinary neutral zero bit from a one bit.
Allowing a supposedly fixed zero-bit input would therefore introduce
unwanted accepted states.

The accepted untrusted input is now specified to contain only ordinary
charged data/reference records and the **unique blank**. A separately
pinned boundary helper and separately pinned verifier scratch remain
trusted. After verification, an explicit involution converts designated
blanks to working zero bits; inverse preparation restores them. No growing
alphabet factor is introduced. The acceptance space is exactly 2^N for
this format. The original construction notes have been amended to say so.

Exact permutation tests cover 32,532 local query columns and 33,504 helper
round trips, including malformed loaded/unloaded symbols. A separate
all-sector initialized history fixture enumerates all 22,500 basis states
of a d=150 two-site core and all 14,800 propagation edges. Exactly the two
intended histories survive its initialization penalties; illegal clocks
and every head-number sector are included. It is a complete finite test
of that mechanism, not the desired high-channel core.

## 4. The outer microscopic bond can now be generated and inspected

The [local assembler](local/README.md) evaluates the complete physical
two-site bond, given actual core and port matrices. It includes every
selector, wrapper, filter, helix, compensated hopping, marker-motion and
chemical-potential term, on the entire local Hilbert space.

For an explicit one-fermion positive core of length three, it produces a
**d=240** microscopic model. The entire 57,600 by 57,600 two-site bond is
saved in [coordinate form](local/fixture_bond.npz), with 122,192 entries.
Every column was evaluated: the matrix is Hermitian to below 7e-15 and
conserves U(1) charge and marker number exactly in its sparsity pattern.

An independent test contracts the **actual physical hopping entries**
through explicit endpoint reduced operators at all four combinations
of cell lengths four and five. It recovers the intended hopping matrix
with maximum relative discrepancy below 5.8e-13. Negative controls show
that dropping X# or omitting the compensators gives the wrong result.
This checks the local formula's implementation and normalization. It is
one logical mode per cell and does not violate c-d. Its chosen finite
coupling has not been certified to lie in a thermodynamic conformal phase.

## 5. The next finite construction that must be completed

The decisive remaining deliverable is a **single initialized autonomous
controller**, with one complete alphabet count, that performs:

1. Determine the physical core length and decode p.
2. Generate the exact instruction tape and data layout from that parameter,
   using fixed finite rules and workspace bounded by the selected cell.
3. Check every untrusted record using the repaired input format.
4. Execute the circuit, transport the outputs to the physical endpoints,
   and expose READY only during verified clean intervals.
5. Reverse all work and reset, with correctness checked at every elementary
   history-clock vertex, not only completed program steps.

The present dispatcher does not cover items 1–3 or endpoint readout; the
separate length primitive supplies part of item 1 before composition. The
host-side Python circuit generator cannot be substituted for a physical
controller executing item 2. Pinning a different arbitrary instruction
at every site would break the required uniform-bond construction.

After this controller is counted, the actual microscopic dimension can be
computed and compared with N+1 for a definite p. The fixed-positive-coupling
infrared argument remains a further requirement, even after every local
matrix is specified. Current finite checks do not resolve that limit.

## Reproduction

Each linked component README gives its exact commands and generated result
files. Python's standard library suffices for circuit, controller and graph
checks; the full matrix assembly uses NumPy. Large gate streams need not
be written to disk. The checked-in tables, manifests and fixture matrix
are all small enough for direct repository review.

Run all new finite checks with

```
OPENBLAS_NUM_THREADS=1 python3 contrib/controller_instantiation/verify_all.py
```

The coordinator reran all nine checks successfully. The combined outcome
is recorded in [verification_summary.json](verification_summary.json).
