# Adversarial audit of the initialized controller interface

Author: Codex parallel compiler-interface audit agent, 2026-10-28.
The scope is the exact finite compiler and input-format claims in
`contrib/exact_route` and `contrib/torus_bottleneck/EXISTENTIAL_CORE.md`.
This is not a proof of the infrared limit or a completed microscopic
counterexample.

## Concrete gap and a constant-alphabet repair

The original format-verifier description combines three TYPE-bit
queries with indexed species queries. Its TYPE encoding is

| Record | TYPE value | Species-query output |
| --- | --- | --- |
| `D_s` | 0 | bits of `s` |
| blank | 1 | zero |
| ordinary bit 0 or 1 | 2 | zero |
| `H_(j,c)` | 3 | zero |
| `T_(j,c,s)` | 4 | zero |

Consequently those observables cannot distinguish ordinary bit 0 from
ordinary bit 1, or two different H records. A verifier built exclusively
from these queried functions cannot require an arbitrary specified
neutral record. In particular, allowing an ordinary zero-bit slot in
the stated initial accepted format without an additional value query
would also accept its one-bit version. The resulting history space
would have unwanted degeneracy. This is an actual distinguishability
failure, rather than a matter of numerical precision or controller size.

A sufficient explicit repair is to choose the accepted *untrusted*
`Σ_q` input slots as follows:

1. Each logical occupation slot permits exactly its two ordinary
   charged symbols `D_s,D_t`.
2. Each charged reference slot permits exactly its specified `D_s`.
3. Every other untrusted slot permits the unique blank symbol.
4. The single right-boundary `H_(0,0)` helper is separately fixed by
   the existing initialization pin and is excluded from the untrusted
   scan. Trusted verifier scratch remains in the independently pinned
   neutral Q tracks.

For this format, every accepted fixed symbol is an individual
equivalence class of the query signature. Each occupation slot's
accepted equivalence classes contain precisely its two desired charged
symbols. The accepted tensor-product space therefore has dimension
`2^N`; every other Σ record at every other slot fails. This proof
includes malformed H and T states, which have distinct TYPE codes.

After acceptance, a designated neutral scratch slot is converted by
the explicit one-site permutation

\[
U_{\mathrm{blank},0}
=|\mathrm{bit}\,0\rangle\langle\mathrm{blank}|
+|\mathrm{blank}\rangle\langle\mathrm{bit}\,0|
+\sum_{x\notin\{\mathrm{blank},\mathrm{bit}\,0\}}
|x\rangle\langle x|.
\]

It is an involution on the entire alphabet and preserves all species
charges. The existing bit-to-`H_(0,c)` conversion and SHIFT produce
other required neutral helper inputs if needed. These operations are
part of the preparation word and are reversed by its literal inverse.
No extra indexed register or opcode family growing with q is required.
They are executed also on rejected inputs, where their arbitrary-input
unitarity suffices because FAIL is already retained.

The earlier theorem statements should explicitly adopt this format,
or supply additional neutral-label queries. The supplied finite query
checks establish this repair; they do not retroactively identify an
unspecified format with this one.

## What was checked

`format_query_checks.py` constructs the complete alphabet, including all
loaded and unloaded indices and neutral values. It verifies:

- the neutral-bit alias above, and the absence of aliases among the
  repaired accepted references;
- TYPE and species queries as charge-preserving involutions on **all**
  two-site basis states, including H/H and H/T malformed pairs;
- SHIFT and blank-to-zero conversion as full-space involutions or
  inverse permutation pairs;
- the actual helper transport, indexed query, copy, unquery, and inverse
  transport on all inputs of the tested short tapes, including arbitrary
  malformed records and either initial trusted-copy-bit value;
- exact accepted dimension four for a two-occupation example with a
  charged reference and a blank padding slot.

The exact integer/permutation checks cover 32,532 local-query columns
and 33,504 complete helper round trips. There are no floating-point
tolerances in these checks. Results are in `format_query_results.json`.

Run:

```text
python3 contrib/controller_instantiation/audit/format_query_checks.py
```

## Remaining controller requirements

The fresh/updated G-block construction is a valid local unitary
interface. Its disjoint domains and ranges prevent a second move in
one outward sweep. This is stronger than a simulation on a selected
one-head input: every pair containing two interpreter heads is fixed,
so the full matrix remains unitary.

That result does not provide a controller collision table. Any actual
F collision must be a bijection from updated to fresh neutral control
and local-metadata states. For example, setting a right-moving head to
left-moving at an endpoint while also fixing a left-moving head there
would merge two inputs and is invalid. A swap of those two local cases,
or a larger disjoint permutation cycle, repairs the issue. Completion
must preserve the prescribed valid computation and not merely fill
unassigned matrix entries arbitrarily.

READY correctness must be checked at every elementary outer-clock
vertex, including the intermediate return sweep that applies F at
one site at a time. The existing head-absence projector is a sound way
to suppress transitions while a flag is set or cleared at an endpoint.
It remains necessary that all quantum preparation has already finished
before either endpoint is exposed and that all later exposed motion
acts identically on its data record.

The initialization and positive FAIL-compression gap argument has no
identified defect in this audit once its stated hypotheses hold.
However, that conclusion still depends on a correct actual finite
controller, a complete check on all untrusted positions, a data-independent
cycle that fits in the outer clock, and the corrected accepted input
format above. Existing tests of the sweep interface do not establish
those missing transition-table properties.

## An all-sector initialized toy core

`tiny_history_checks.py` supplies a fully specified two-site example of
the history initialization mechanism. Each site has 15 clock states and
10 work states: two ordinary charged records, plus fresh/updated and
left/right head overlays of each record. The head executes the explicit
two-record SWAP opcode. Endpoint collisions flip direction by a full
permutation, and the right seed swaps an ordinary record with its
specified starting-head overlay.

The script enumerates all 22,500 physical basis states and all 14,800
decorated propagation edges, including every illegal clock word and
every zero-, one-, or two-head work sector. Each component's exact
Laplacian zero mode is constant; a positive diagonal pin removes that
zero mode whenever it is nonzero at any component vertex. This gives
an exact integer computation of the kernel, without selecting legal
histories in advance or performing floating-point diagonalization.

There are 100 legal-clock history components before initialization.
The initial I-sweep pins require no head anywhere and fix the first
charged input to zero, leaving the second input free. Precisely two
components survive, each with the required 30 clock vertices. Every
illegal clock sector has a positive defect penalty. Thus this explicit
toy core has an exactly two-dimensional zero-energy space on its entire
Hilbert space.

This verifies the interface and initialization in a finite example. It
does not instantiate the high-channel preparation program, its format
verifier, or its READY protocol. Results are in
`tiny_history_results.json`; run with ordinary Python 3:

```text
python3 contrib/controller_instantiation/audit/tiny_history_checks.py
```

## Audit of the explicit preloaded dispatcher

The concurrently supplied `../controller/dispatcher.py` is an actual
finite rule table for one useful part of the compiler. Its prescribed
collision entries avoid the common reversible merges:

- A fetched instruction enters GO with its program mark set; GO's
  scanning self-loops have that mark clear.
- A marked DATA cursor distinguishes entry into execution or cursor
  shifting from the corresponding scanning states.
- The retained program mark and immutable opcode identify the RET
  exit, allowing the finite opcode memory to be cleared injectively.
- Endpoint bounces belong to the partial injection, whose completion
  fills unused domain and range elements within each program-symbol
  sector. They are not one-way maps merged with fixed output states.

Its control labels each have a fixed direction. Therefore the count
`(1+2K)B` includes both fresh and updated phases without an additional
independent direction bit; it is compatible with the earlier
`(1+4K_h)D` convention. For its 17-oracle menu, `B=88`, `K=150`, and
the actual dispatcher factor is 26,488.

There is also a precise remaining obstruction: this table preserves
each site's read-only `program` symbol on **every** input, including
malformed inputs. Consequently it cannot generate a nonblank
instruction tape from the original blank-initialized classical track.
A separate program-generating control phase, or a different initialized
sector with a complete verifier, is required. The measured 26,488
factor is thus a dispatcher-only constant. It must not be substituted
for the complete core's constant until those extra controls, tracks,
flags, and activation rules have been counted.

Likewise, its forward/inverse cycle alone does not implement the
endpoint READY protocol. A clean readout interval has to be inserted
with injective control updates and checked at elementary clock times.

## Final independent component review

The later length-counter table and local microscopic-bond assembler
received a read-only review against their source formulas and recorded
checks. No additional mathematical defect was identified in that review.

For `../controller/length_counter.py`, the reverse control label has
the opposite direction to its corresponding forward label. On arriving
at the previous tape site, its prescribed collision reverses the earlier
forward collision. The terminal `FINAL_RETURN` bounce and initial
`ADVANCE` bounce occupy otherwise unused domains and ranges. They thus
close a reversible orbit without a local merge. The binary increment
distinguishes its RETURN-entry bit one from its RETURN scanning zeros;
frontier bits survive the increment. The capacity inequality
`n-1 <= 2^(n-2)-1` applies for every promised `n>=4`. The standalone
factor 396 correctly counts that primitive and is not a complete core
count. Its initial left-end head and static marks still have to be
enforced and connected to the outer history initialization.

For `../local/bond.py`, the site encode/decode order agrees with the
stated selector, wrapper, A-helix and B-helix tensor factors. The wrapper
contains its oriented boundary terms, and the selector gates them on
internal bonds. Both X# and R# activate the charged hopping. Its right
and left compensators act on the respective A and B helix factors; the
reverse hopping term has the required adjoints. Marker motion uses the
forward vacancy-conversion phase and its conjugate on the reverse
transition. The one-mode fixture's right creation operator agrees with
the one-mode `f^dagger P` convention. The chemical potential formula
uses the intended free-fermion energy `(k_F-tan(k_F))/pi`.

The complete matrix column test establishes Hermiticity and charge
conservation of the instantiated local matrix. The separate endpoint
contraction checks its hopping normalization through actual matrix
entries. These facts do not establish a periodic-chain CFT or the
global selected ground sector. The source and overall README preserve
that distinction and clearly separate the circuit, dispatcher, length
counter, and fixture dimensions from the still-unavailable complete
high-channel core dimension.
