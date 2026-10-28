# Local digit compilation and a quantitative infrared certificate

**Result: concrete progress on the sparse first-order route, with no new
counterexample claim.** This continuation supplies the previously missing
finite parameter-reader walk, explicit local arithmetic rules for eleven
operations, integration with the actual sparse-circuit generator, and a
proved Gaussian matrix bound that makes one infrared requirement explicit.

The conjecture concerns a one-site translation-invariant, nearest-neighbor
Hermitian chain with microscopic local dimension d and total unitary
Virasoro central charge c. The candidate retains the predicted
N+1=(p+1)(p²+p+1)+1 modes and aims for d=O(p² log p). None of the new
component counts is the local dimension of a completed violating chain.

## What is established

| Contribution | Established result | Remaining limitation |
| --- | --- | --- |
| [Parameter reader](parameter/README.md) | A literal 96-entry collision table; eight controls; complete local gates; parameter output after 2m sweeps and exact reset after 4m | Its midpoint must be joined to the generator and initialization |
| [Digit arithmetic](arithmetic/README.md) | 4,428 fixed adjacent-cell rules implement eleven operations at arbitrary width, with autonomous scratch cleanup | Aligned operands are assumed; the existing sweeping-controller unitary lift remains |
| [Generator integration](integration/README.md) | Actual local-rule execution replaces eight opcodes in the frozen generator and reproduces graph/circuit outputs | Memory, dispatch, operand transport and other arithmetic remain at word level |
| [Gaussian certificate](audit/GAUSSIAN_CERTIFICATE.md) | A full relative stiffness perturbation below 5/13 keeps every reference-dimension-three-or-higher vertex irrelevant | Requires actual filling/momentum margins and a controlled microscopic low-energy reduction |

**Finite construction and checks.** The parameter reader costs exactly
204 D(q) states as a standalone head/record system, where
D(q)=q+3+2(q+1)m. Its count includes both output tracks and all head phases;
it excludes the outer clock and remaining controller. Whole-space
involution, charge conservation, malformed payload restoration and exact
primitive gate-word checks are provided.

The arithmetic rules operate on bits, not unbounded words. A b-digit run
takes 2b+3 nearest-neighbor transitions and restores its inputs and blank
scratch. Its 900 used cell symbols are not a completed physical F/G gate
alphabet. Exhaustive tests cover 60,060 small signed cases, with additional
checks through 4,096 digits and exact inverses of intermediate states.

The short integration run executes 7,616 local arithmetic calls, producing
both p=2 and p=3 incidence graphs and all 906 gates of the one-edge circuit.
The separately recorded larger run executes 1,175,426 local arithmetic
calls at width 128 and exactly reproduces all 287,798 gates of the full
21-mode p=2 circuit. Its 304,435,334 local digit transitions are actually
executed, not extrapolated from operation counts.
The integration documentation explicitly separates implemented digit rules
from host-side word operations and operand encoding.

**Proved matrix bound; conditional physics application.** For the full
Gaussian energy matrix G, including its density-current block, set

    delta = ||G_ref^(-1/2) (G-G_ref) G_ref^(-1/2)||.

The scaling-dimension form obeys

    Delta_G >= sqrt((1-delta)/(1+delta)) Delta_ref.

Thus delta<5/13 suffices to keep every surviving vertex with Delta_ref>=3
above dimension two. The derivation covers unequal, noncommuting velocities
and all vertex charges simultaneously. The bound is sharp as a matrix
inequality, not a necessary phase boundary. It has been incorporated into
the existing [infrared audit](../research_continuation/IR_AUDIT.md).

## What the next complete implementation must supply

The remaining finite tasks are multiplication/division/modulo and variable
shifts; addressed memory and stack/register transport; a sweep-compatible
arithmetic lift; and reversible handoffs joining parameter input, the real
generator, format verification, charged execution and endpoint readout.
Their combined alphabet must then be counted and a definite violating
prime chosen. Standalone component factors cannot be multiplied or added
without specifying their shared tracks and controls.

The major remaining physics assumption is still local massive-sector
elimination at fixed nonzero coupling. The existing defect-energy bound
does not itself prove the all-order connected-excursion estimate with
moving cell projectors. This continuation found no new resonance
obstruction, but does not settle that assumption or measure c.

## Reproduction and provenance

```
OPENBLAS_NUM_THREADS=1 python3 contrib/digit_compilation/verify_all.py
python3 contrib/digit_compilation/integration/checks.py --full
```

The first command needs NumPy only for the Gaussian audit; the controller,
arithmetic and integration checks use the standard library. The optional
full integration takes several minutes. Each subdirectory retains its
executable data and result report. `verification_summary.json` records the
coordinator's reruns.

Produced by Codex (GPT-6 family) and three parallel research workers on
2026-10-28, continuing the clean checkout at `588aa43`. The coordinator
independently reviewed the matrix proof and local rules and reran the new
checks; a separate worker reviewed the parameter-reader implementation.
No human review, completed high-channel core, final violating d, or
thermodynamic central-charge measurement is claimed.
