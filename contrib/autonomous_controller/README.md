# An autonomous computation and readout demonstration

**Status: explicit finite controller progress; no new c-d violation.**
The earlier dispatcher required a supplied instruction tape. There is now
one reversible local table which starts with blank dynamic tracks, computes
a length-dependent program, executes it coherently, exposes its output at
both endpoints, and uncomputes the generated program and counter. A separate
fixed word-machine program generates the actual sparse-encoding circuit;
its arithmetic and addressing still require local tape compilation.

This contribution continues [controller_instantiation](../controller_instantiation/README.md)
from commit `06fb610`. Author: OpenAI Codex, GPT-6 family, with three parallel
implementation/audit workers, 2026-10-28 UTC. Component authors exchanged
interfaces and a separate worker inspected the elementary history-clock
vertices. No independent human review or thermodynamic CFT calculation is
claimed.

## 1. What is now autonomous

The [integrated bridge](bridge/INTEGRATED_READOUT.md) uses one explicit
3,200-entry collision permutation. It has 80 neutral tape symbols and 40
direction-fixed controls, including all malformed local states. Starting
with zeroed dynamic tracks and a single initialized head, it performs:

1. Count the sites and write n-1 in binary.
2. Reversibly copy those digits into the program track.
3. Execute CNOT or controlled-sqrt(X) at successive neighboring quantum
   registers according to those digits. The gates do not generally commute.
4. Transport the two designated output records to the physical endpoints
   and enable their READY flags during safe intervals.
5. Return the output records, execute the literal adjoint circuit, erase
   the program, and reverse the counter. Repeat from the original state.

Every phase handoff is a table entry; no host-side intervention occurs.
The control trajectory is independent of the arbitrary quantum input.
The four dynamic neutral bits all return to zero. The complete cycle takes

    T(n) = 2n^2 + 22n - 30 - 4 popcount(n-1)

head moves. Complete cycles pass at all 125 lengths 4 through 128, totaling
1,590,466 moves. The example is a genuine generated quantum program, but
**it is not the many-fermion sparse-encoding program**.

The [readout primitive](readout/README.md) is independently reusable. It
needs only identity and quantum-record SWAP, ten neutral symbols, and eight
controls. At every elementary gate vertex an active physical endpoint
contains exactly the designated prepared output tensor factor. This proves
the port-compression identity even for arbitrary entangled prepared states.
Other parts of the state may change while one endpoint is exposed; only
commutation with that endpoint observable is needed.

## 2. The actual clock matters

The [independent audit](audit/README.md) executes the integrated machine
on the actual fifteen-letter C8 clock path, including initialization,
binary carries, and idle vertices. It does not replace those vertices by
a uniform distribution over completed machine steps.

The initialization seed is applied during the first I-sweep departure,
controlled by the immutable LEFT mark. It swaps the ordinary work state
with its fresh start-head overlay. Subsequent initialization vertices visit
other sites, so their ordinary-head pins remain valid. There is exactly
one initialized head without adding a one-way entrance into a reversible
cycle.

There is an instructive finite-size constraint: the integrated program
does not reach any READY interval within the clock at n=4,5,6. It does at
n=7,8, and for every n>=9 an entire work cycle fits within 2^n sweeps.
The actual history weights at n=9 are 658/11269 on the left and
656/11269 on the right. Equality of the two weights is unnecessary;
positivity and independence of logical input are the required properties.

The [initialized demonstration core](core/README.md) combines these local
rules with positive clock, grammar and initialization terms. Its repeated
bulk and boundary operators are specified by sparse basis-column
evaluators, including illegal sectors. It remains a demonstration core,
with neutral qubits and a simple binary-digit program, rather than the
charged high-channel core needed for a counterexample.

## 3. A fixed program now generates the full desired circuit

The [generator](generator/README.md) is a frozen program with **3,115 word
instructions**, 34 functions, 27 primitive opcode families and at most 38
registers per activation. It takes p and emits the precise instruction
tape accepted by the earlier dispatcher, without calling the earlier
Python graph or circuit routines. Those routines are used only as a test
reference.

The full abstract gate streams match for p=2 (287,798 operations) and p=3
(1,588,624 operations). In the one-edge fixture, all 36,140 localized gates
and all 84,763 dispatcher symbols match. A reversible instruction journal
restores registers, arrays, output and call state in the tested cases.
An independent adversarial test exercises aliasing, nested calls, both
branch outcomes, output and reversal at every instruction.

This **word machine is not yet a local finite-alphabet machine**. Its
large integers, array accesses and journal operations have not been
expanded into digit operations and nearest-neighbor head moves. The
program count therefore is not a value of the full core constant C.
The saved source-derived polynomial bounds concern word instructions and
serialized memory, not final microscopic running time or cell length.

## 4. A useful simplification: p need not be decoded from cell length

The [parameter interface](audit/PARAMETER_ORACLE.md) supplies two explicit
whole-space neutral permutation gates on the existing indexed helper:

    Q_p: H(j,c) -> H(j,c xor bit_j(p)),
    END: H(j,c) -> H(j,c xor [j=m-1]).

The local Hamiltonian matrices already depend on p through the species
alphabet and exact rotation blocks. These two additional fixed gate
families therefore legitimately expose p without increasing the number
of index or species registers. This is a finite binary lookup, with no
precision or real-number decoding assumption.

Using the trusted boundary helper, a compact walk writes the parameter
bits into trusted neutral tracks, then restores every untrusted record
and the helper. Its expanded reference circuit has 10m-4 primitive gates.
Whole-space and malformed-input checks pass. The walk's finite collision
table remains to be instantiated. This removes integer-root decoding from
the high-channel plan; it does not remove workspace bounds or verification
of unused padding.

## 5. What remains for the counterexample

The most direct next target is **digit-level compilation of the frozen
generator**, then composition with the existing charge-preserving gate
dispatcher, input-format verifier, and the endpoint-readout protocol.
The binary-digit demonstration proves that a complete local autonomous
cycle can be built, but does not substitute for that compilation.

After composition we must count the entire high-channel local alphabet,
choose a definite p and cell length, and compare N+1 with
d-1=96 d_core+47. The infrared transfer at fixed nonzero physical coupling
retains its stated assumptions and has not been certified by these
finite-controller tests. The earlier d=240 one-mode fixture remains a
separate assembly test with intended first-order c=2.

## Reproduction

Each directory includes its tests, exact tables where supplied, and
machine-readable results. Run all new checks with Python and NumPy:

```
OPENBLAS_NUM_THREADS=1 python3 contrib/autonomous_controller/verify_all.py
```

`verification_summary.json` records the coordinator's combined run: all
nine check scripts pass, including the complete generator comparisons and
the explicit initialized-Hamiltonian residual checks.
