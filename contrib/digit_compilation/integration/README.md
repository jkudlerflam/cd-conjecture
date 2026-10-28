# Executing the actual generator with local digit kernels

**Status: exact finite integration checks; a hybrid implementation.** Eight
integer opcodes of the frozen sparse-circuit generator now execute through
the fixed adjacent-cell rules in [arithmetic](../arithmetic/README.md).
The result returned by those rules drives every subsequent branch and
memory operation. It is not a comparison against an unused arithmetic
subroutine.

The other VM operations, operand placement, and the remaining integer
operations still execute at word level. This test does not compile memory,
program dispatch, or autonomous transfers between subroutines. No final
local dimension follows from its success.

`checks.py` subclasses the existing word interpreter without changing its
source or the frozen 3,115-instruction program. For ADD, SUB, AND, XOR,
EQ, NE, LT and GT, it encodes the two input words on aligned binary tracks,
executes only the frozen local rewrite table through compute/copy/cleanup,
and decodes the returned result. The old interpreter runs independently in
lockstep. Instruction pointer and stack are compared at every instruction;
the complete semantic state and instruction journal are compared at the
end. Separate comparisons use the original graph and circuit implementation.

The short run covers all 21 edges at p=2, all 52 edges at p=3, and all
906 abstract gates of the one-edge encoding. It performs 7,616 actual local
arithmetic calls at signed width 16. All three runs also reverse their
word-level journals and recover the exact initial state. Each arithmetic
call already clears its own digit-level scratch by autonomous table rules;
the surrounding VM journal is still an abstract memory structure.

The optional `--full` run also executes the entire p=2, 21-mode abstract
encoding with signed width 128. The saved passing run performs 1,175,426
local arithmetic calls (304,435,334 local digit transitions) and matches
all 287,798 abstract output gates to the existing circuit. It took about
137 seconds on the development host. This larger case retains no word journal and does
not claim a full-machine rollback. Its report is saved separately so the
fast test cannot overwrite it.

ADD and SUB in the local kernel wrap at the chosen width. Matching the
unbounded VM at every step verifies that no such overflow occurs in these
tested workloads. The selected widths are test parameters, not a proof of
the word-size bound for all primes.

## Reproduction

Only the Python standard library is required:

```
python3 contrib/digit_compilation/integration/checks.py
python3 contrib/digit_compilation/integration/checks.py --full
```

The short run takes less than a second on the development host. The full
run takes several minutes. Results are saved as `results.json` and
`full_results.json`, respectively. In both reports the unresolved host
arithmetic is counted explicitly. A test passing does not turn those
operations, stack/array access, or dispatch into local transitions.

Developed and independently integrated by the coordinating Codex agent
(GPT-6 family), 2026-10-28. No human review or central-charge measurement
is claimed.
