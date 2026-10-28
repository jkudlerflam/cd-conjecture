# Sparse compiler audit: separate computability from polynomial efficiency

Status: **[proved resource observation; incomplete Hamiltonian audit]**.
AI-assisted independent audit by a GPT-6-family agent, 2026-10-28.

The parameter-query gates already supply the selected prime `p` through
finite matrices on the existing helper alphabet. Consequently the
existential dimension comparison does **not** require a polynomial bound
on the runtime or storage of the fixed classical circuit generator.
For each fixed `p`, any total finite program has finite runtime and
workspace; the selector can choose a sufficiently large finite segment
length. The charged work alphabet remains
`D(q)=q+3+2(q+1)m=O(q log q)`, multiplied by a constant neutral-machine
factor. Enlarging a segment changes neither local dimension nor this
asymptotic comparison.

This observation removes a resource-efficiency prerequisite; it does not
supply a missing machine table. The following distinctions matter:

1. A fixed universal/interpreter program may be preloaded by a finite
   boundary-recognizing prefix automaton. The prefix length is a constant
   independent of `p`, so this has constant alphabet overhead. A
   **prime-dependent gate stream** cannot be supplied the same way without
   counting its growing prefix alphabet.
2. A finite one-head reversible classical machine acts on a neutral tape
   tensored with **one** charged record. Its head overlays that record;
   nothing in universal classical simulation requires a second charged
   species register or a `q^2` single-site factor.
3. Long scans of unused padding and endpoint transport must remain compact.
   Storing a separate history entry for each step of an arbitrary
   `poly(n)`-time scan can exceed the `n`-site capacity. The already
   specified reversible marked walks and binary violation counter address
   this issue; simply choosing larger `n` does not cure a circular `n^2`
   workspace requirement.
4. A sufficient timing condition is `T_core(p)+C(p)n^a < 2^n` with fixed
   `a` and finite `C(p)`. Once the `p`-dependent core uses finite workspace,
   some sufficiently large `n` satisfies both capacity and this condition.
   No computable polynomial upper bound in `p` is needed.
5. Universal computation alone does not prove the initialized history
   kernel, accepted-input dimension, clean two-endpoint READY readout,
   exact cyclic reset or mobile infrared limit. Those interfaces remain
   necessary.

The [existing scheduler](../../exact_route/ORACLE_SCHEDULER.md) already
argues the constant-factor existential compilation. The digit compiler
adds inspectable finite rules; its unfinished whole interpreter is an
implementation gap and cannot, by itself, refute the existential argument.
Conversely its component state counts are not a certified numerical value
for the complete constant `C`. This audit does not claim a concrete
violating prime or microscopic dimension.

Primary-source check: Morita, Shirasaki and Gono's publisher abstract states
that any irreversible machine can be converted to a one-tape reversible
machine, first with 32 tape symbols and then two. The abstract does not
supply the complete instrumented quantum-oracle table or this project's
constant; no such count is inferred here.
[Original publication](https://globals.ieice.org/en_transactions/transactions/10.1587/e72-e_3_223/_p).

The near-term search has moved to the much smaller explicit d=23 frozen
cell candidate and a one-vacancy mobile lift. That route can avoid this
large compilation task entirely if its low-energy stability can be shown.
