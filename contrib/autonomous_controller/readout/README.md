# An explicit reversible endpoint-readout cycle

**Status: finite local rules implemented and checked.** This is a separate
transport/readout primitive, conditional on a prepared quantum state and
verified output-location marks. It does not initialize those marks,
prepare the sparse high-channel encoding, or complete its Hamiltonian.

Author: OpenAI Codex, GPT-6 family, 2026-10-28 UTC. The coordinating agent
wrote the table and checks; a parallel audit is recorded in the enclosing
contribution. No human review is claimed.

## Contract and explicit count

There are n>=4 sites and two output slots 0<i<j<n-1. Immutable neutral
marks identify LEFT, RIGHT, OUT_L at i, OUT_R at j, and BULK elsewhere.
The one neutral READY bit is initially zero everywhere. The head starts
at OUT_L with control LEFT_OUT. The quantum state is arbitrary, including
entanglement between the outputs and all other records.

`transport.py` supplies a complete 80-entry collision permutation, with
20 prescribed entries. It has ten neutral symbols (five marks times two
READY values) and eight direction-fixed control states. The ordinary plus
fresh/updated head space has dimension

    (1+2*8)*10*D = 170 D,

where D is the single quantum-record dimension. This is the standalone
primitive's factor, not an additional mandatory factor to multiply into
every other module or a count of the finished core.

The complete table is saved as `transition_table.json`. Its completion
preserves the static marks even on malformed inputs. The two quantum
opcode families are identity and whole-record SWAP. No second charged
record is held in finite control. Neutral metadata stays at its physical
site when quantum records move.

## The cycle

| Phase | Head travel | Quantum action | Arrival transition |
| --- | --- | --- | --- |
| LEFT_OUT | i to 0 | successive SWAPs | Set left READY |
| FIND_RIGHT | 0 to j | identity | Enter RIGHT_OUT |
| RIGHT_OUT | j to n-1 | successive SWAPs | Set right READY |
| READ_BACK | n-1 to 0 | identity | Clear left READY |
| LEFT_RETURN | 0 to i | inverse left SWAPs | Enter FIND_END |
| FIND_END | i to n-1 | identity | Clear right READY |
| RIGHT_RETURN | n-1 to j | inverse right SWAPs | Enter FIND_START |
| FIND_START | j to i | identity | Return to LEFT_OUT |

All phase changes occur in the collision table. No host-side mutation,
external halt signal, or instruction string is used. The eight distances
sum to exactly 4(n-1), independent of i and j. The two SWAP walks have
disjoint supports; their inverses restore every quantum record. The
READY bits and all control state also return exactly to their start.

## Local reversibility and full-space gates

Every phase has a distinct control label. A phase's scan loop encounters
only its listed interior marks; its entry image has a different mark or
READY value. For example, FIND_RIGHT enters on LEFT with READY=1, whereas
its scan images are BULK or OUT_L with READY=0. LEFT_OUT's entry image
is OUT_L, whereas its loop image is BULK. These separations prove that
the 20 prescribed local updates have disjoint domains and ranges.

The finite unused domains/ranges are bijectively paired within each static
mark sector. Fresh/updated head phases turn the permutation C into an
onsite involution F, using C on updated heads and C^-1 on fresh heads.
The two-site involution G moves one head in its fixed direction while
applying its quantum opcode; the reverse block is its adjoint. Zero-head
and two-head sectors are fixed. Thus F and G are defined on the entire
local Hilbert space, including malformed metadata and head sectors.

An outward ordered G sweep makes exactly one head move. The return F
sweep supplies exactly one collision. The original sweep-interface proof
therefore applies without adding species-dependent control states.

## Endpoint compression at every elementary clock vertex

Define the physical left port to be a quantum endpoint operator O,
multiplied by READY and by the absence-of-head projector. At any elementary
G/F vertex t where this flag is active, the record at site 0 is exactly the
original record at i. Consequently, for the prepared isometry V,

    V_t^* [READY_0 P_no-head O_0] V_t
        = a_L(t) V^* O_i V,       a_L(t) in {0,1}.

The same identity holds at the right endpoint with j. If a head is at
the endpoint, the port vanishes. Setting or clearing READY occurs only
while that head is present, so intermediate collision/sweep vertices
cannot expose a partially transported record.

The right flag remains active while the LEFT_RETURN phase undoes the
left transport. This is safe: that phase acts on the complement of the
right endpoint and therefore commutes with its observable. It is not
necessary that the entire prepared state remain frozen. The two port
identities hold even when the endpoint records are entangled.

For a history with orthogonal clock vertices and weights w_t>=0, summing
gives exactly p_L V^*O_i V and p_R V^*O_j V, where
p_L=sum_t w_t a_L(t) and similarly on the right. The schedule and weights
are independent of the logical quantum input. One complete cycle has
positive exposure at each end. Additional work-identity clock vertices
cannot invalidate this identity. The uniform sample weights in results
are **not** advertised as the full binary clock's visibility fractions.

## Checks and limits

Run:

```
python3 contrib/autonomous_controller/readout/transport.py
python3 contrib/autonomous_controller/readout/checks.py
```

The checks verify all 340 one-site and 115,600 two-site states for a
two-state quantum record, including every head sector and READY value.
They follow every elementary sweep vertex for all 364 output geometries
on lengths 4 through 15: 402,948 vertices. Exact record permutations
verify endpoint safety and whole-cycle cleanup. A separate arbitrary
complex four-dimensional code in six quantum bits verifies coherent
operator compressions; errors are below 2e-16. Removing READY produces
large discrepancies, so the projection test is sensitive to flag misuse.

Not supplied here: verified initialization of OUT_L/OUT_R, a local
handoff from the high-channel preparation to this routine, its return
into inverse preparation, or an all-sector spectral proof of the fully
composed initialized core. A cycle cannot acquire an arbitrary one-way
lead-in under a finite permutation; integration must splice actual local
transitions with unique predecessors. The adjacent bridge contribution
investigates that composition for a smaller autonomous example.
