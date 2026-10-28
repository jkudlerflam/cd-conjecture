# Independent audit of the autonomous-controller continuation

This folder contributes two concrete results: a parameter-input shortcut,
and independent checks of the endpoint-readout/clock and journal-machine
interfaces. It does not certify the full high-channel controller or its
infrared limit.

## Parameter input

[PARAMETER_ORACLE.md](PARAMETER_ORACLE.md) constructs two local 0/1
permutation matrices exposing the chosen prime through the existing
indexed neutral helper. Since the Hamiltonian family already depends on
that prime, this is permitted without a growing control alphabet. The
extraction word restores every other quantum record, including malformed
ones, and can be inverted to clear its output. It removes the need to
recover the prime by taking an integer root of the physical cell length.
The compact extraction algorithm still needs its own finite-table lowering.

## Readout under the actual outer clock

The transport controller preserves the output's original tensor factor
at each active endpoint. That is enough for clean readout on arbitrary
entangled inputs. It is **not** necessary to freeze the entire prepared
state: undoing left-port transport can occur while the right port remains
active, since that operation does not act on the physical right endpoint.

More formally, if `P_t` is the payload permutation and endpoint `e` is active
only when it contains original output factor `o`, then

\[
P_t^\dagger A_e P_t=A_o.
\]

For any input isometry `V` and any local operator `A`, each active vertex
therefore has the identical compression `V^dagger A_o V`. The flag and
head-absence projector make every other vertex contribute zero. This is
an operator statement for an arbitrary input code, not only a basis-word
or product-state test.

The root tests enumerate all elementary G/F vertices. The independent
`sweep_audit.py` also runs the literal C8 clock path: initial I sweep,
outward G updates, return F updates, all binary carry/turn/return idles,
and the truncated final controller cycle. Those idle vertices affect
readout weights and are included in the recorded exact integer counts.
The readout input is supplied with its head at the designated interior
output mark; these tests do **not** claim to prepare that state using the
outer clock activation.

## Initialization orientation and cycle entry

The existing C8 activation seeds at RIGHT, whereas the standalone blank
length/bridge machine starts at LEFT. They cannot be joined by merely
adding a one-way transition into a closed machine cycle: a finite
permutation has no incoming edge to an invariant cycle from outside it.
The root construction supplies a concrete orientation fix: decorate
`IZ -> LI` with a seed involution on the departed site's work register,
controlled on its static LEFT mark. Only the first initialization edge
sees LEFT. It swaps ordinary with the chosen fresh start-head overlay;
all other head types are fixed. The right activation becomes identity.

The first I vertex pins LEFT ordinary **before** that seed; every later
I vertex pins a different site, untouched by the seed. Consequently the
pulled-back initialization projectors still require exactly the original
all-ordinary head input, while the computation begins with one fresh head
at LEFT. The existing right-end helper pin is unchanged. This is a full
work-space involution preserving static marks and species, so it can
decorate the existing clock rule without changing clock legality.

The audit exhausts the head-type initialization basis through six sites,
and then follows the seeded **integrated** blank-to-program-to-readout
machine through every actual C8 clock vertex for lengths four through ten.
Every completed sweep agrees with the independently stepped macro machine,
including oracle order and adjoints. Counter idles are included in the
recorded readout weights. A complete work cycle is guaranteed to fit for
all `n >= 9`: its time is bounded by `2n^2+22n-30`, which is below `2^n`
at nine and remains so by induction. Smaller lengths can truncate before
the first readout and should not be assumed to have positive port weights.

Similarly, preparation, execution, readout and cleanup are reversible
subroutines only when their handoffs are actual injective local rules.
A host program calling separately closed cycles is not that composition.
The bridge work addresses this by splicing its phase boundaries directly
into one partial collision permutation.

## The instantiated finite core and its gap

The independent source audit of [core.py](../core/core.py) and
[local_gates.py](../bridge/local_gates.py) found no mismatch between the
full-space gates and the checked initialized trajectory. The full work
dimension is 12,960 and the fifteen-letter clock gives `d_core=194400`.
The sparse column functions define every matrix entry, including illegal
clock/work sectors. This is a demonstration core, with arbitrary quantum
qubits as its input; it is not the sparse high-channel code.

For `n >= 4`, clock and static grammar select a unique clock path and the
unique static word `LEFT OUT_L BULK^(n-4) OUT_R RIGHT`. Every propagation
decoration is a whole-space involution. The initial pins then leave exactly
the `n` quantum qubits free, giving an exact `2^n`-dimensional zero space.

An explicit positive gap follows without measuring a giant matrix. Let

\[
C=(2n+4)2^n+n-4,\qquad \gamma=1-\cos(\pi/C).
\]

On the legal clock/static sector, propagation is unitarily equivalent to
the half-weight path Laplacian, whose gap is `gamma`. Let `P` be its history
kernel and let `A` be the sum of I-gated initialization pins. There is at
most one I clock letter at a legal vertex, so `||A|| <= 1`. The seed acts
only at LEFT after its pin, so all pulled-back pins are precisely the
original commuting input projectors `P_bad,i`. Therefore

\[
PAP=\frac1C\sum_iP_{{\rm bad},i}
\ge\frac1C(P-P_{\rm good}),\qquad AP_{\rm good}=0.
\]

Applying the proved positive-pin inequality in
[COMPILER_INTERFACE, section 2](../../torus_bottleneck/COMPILER_INTERFACE.md)
with `a=1/C`, `B=1` gives

\[
\boxed{\Delta_n\ge\frac{\gamma}{1+C+C\gamma}>0.}
\]

Clock or static defects cost at least one and are preserved by propagation,
so the displayed (smaller-than-one) bound also holds for the full core.
The audit rechecked the inequality's proof: the quadratic-form lower bound
is `gamma ||y||^2 + (sqrt(a)||x||-sqrt(B)||y||)_+^2`; its positive branch
has determinant `a gamma` and trace `a+B+gamma`, while its other branch is
stronger. No unproved commutativity on arbitrary malformed work inputs is
needed.

The script also checks the complete 12,960-state work initialization pin,
and samples Hamiltonian columns in both directions of every clock rule,
including every oracle control with fresh/updated heads and complex CV
entries. These column samples supplement the operator proof; they do not
enumerate the roughly 38 billion two-site basis columns.

## Fixed bytecode is still a word machine

The generator source has no callback into the existing encoding builder,
and its interpreter uses integer arithmetic, array addresses and a fixed
bytecode. The journal inverse is correct on reached semantic states,
including frame returns and heap aliasing. Independent adversarial programs
check every intermediate reversal, not just their final checksum.

Its instruction count, opcode count and registers per frame are **not**
a local-state count. The unbounded integer operands, addressed heap, return
frames and saved predecessor/overwritten-value history require tape bits
and a local reversible interpreter. The host-directed halt/reverse command
also requires an explicit controller handoff. These are remaining lowering
tasks; the current word-machine implementation does not claim otherwise.

## Reproduction and provenance

Run the two scripts:

```
python3 contrib/autonomous_controller/audit/parameter_checks.py
OPENBLAS_NUM_THREADS=1 python3 contrib/autonomous_controller/audit/sweep_audit.py
```

The latter uses NumPy and reads the existing literal C8 clock rules without
importing its unrelated spectral utilities. Result files state the checked
finite ranges. This audit was performed by a
GPT-6-family agent separate from the new readout, bridge and generator
implementers; no independent human review is claimed.
