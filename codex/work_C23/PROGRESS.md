# C23 progress

2026-10-27. Read the project STATE, route-A background, repository STATUS, C8 progress and the relevant C7/C8/C13/C18/C19/C20 results. C21 and C22 progress files were absent on the initial read. No subagents or network used. Only this folder is writable for this task.

## Completed source audit

- C19 supplies an exact finite matrix cell and second-order target, not a spatial compiler. Its reflection dimension is 16(n+2)2^n; its port source weights, rather than its eigenvalues alone, are the required interface.
- C8 wrapper and helix equalizer prove first-order identities conditional on a positive frustration-free initialized core. Neither proves equality of route-A second-order resolvents at the two lengths. A non-frustration-free C19 simulator would also need an extension of the wrapper theorem.
- C18 gives a 15-letter clock and two idle ports, with an initialization gap that can be smaller than the clock gap. The exact positive static-Bell interface has scalar source resolvents and cannot realize C19.
- C20 gives M(n)=1+8n+32 binom(n,2)+96 binom(n,3)+128 binom(n,4), normalized detuning 1/(4M), and a conditional local-vertex/gluing test. C13 supplies the physical-coordinate dictionary at the reference projection and qualifies finite-ring uniqueness and the accumulation point.
- Read the collaborator continuation, exact/sweeping interfaces and graded repair. The sparse first-order family has growing alphabet and does not close fixed-d route A. Its first-order order-of-limits hierarchy cannot be copied to a second-order payload.

## Next

Derive explicit finite-cluster bounds, a conditional hierarchy separating structural errors from tunable higher orders, and n=2,3 numbers. Check the resonance list and Jacobi visibility using a small local Python calculation. Write RESULT.md with an interface ledger, honest unknown resource functions, ranked gaps and concrete calculations. No actual microscopic bond or full accumulation Hamiltonian is currently specified.

## Finite checks completed

`check_assembly.py` now verifies the Jacobi Gram identities with exact rational arithmetic for n=2,3, enumerates all 49 and 217 C20 vertex pairs, and checks explicit rational momentum scaffolds by integer arithmetic. It gives rho^2=1/172800 and 1/10536960. The target has displaced centers and order-one velocities in both examples.

Derived a scalar-rescaling obstruction: if the two lengths have source gaps d0,d1, their scalar resolvent factors 1/(d_p+d_q) cannot all be made equal by endpoint rescaling unless d0=d1. This obstructs using first-order helix equalization as a second-order proof, not route A itself.

Derived a finite-bond cubic bound 27 lambda^3/(4 Delta^2) for C19 and the absence of mixed-bond cubic contributions on the infinite nearest-neighbor line. Higher-order interactions need not be Gaussian, but the C20 detuning conditions control their weak-coupling RG fate conditionally.

The arithmetic hierarchy check uses explicitly assumed completed-core resources. Its very small logical scales and very large clock coefficients are feasibility arithmetic only, not parameters of a constructed Hamiltonian. Actual clock strength, wrapper gap, accuracy and cell length are interdependent. A scale proof cannot treat them as freely independent until the compiler is supplied.

Read newly appearing C21 and C22 progress. C21 is testing boundary moments and local spectral ladders; C22 is testing a looped excitation bus. Neither progress report yet contains a completed fixed-alphabet target. The integrator will not duplicate these calculations.

## Completed handoff, 2026-10-28

`RESULT.md` is complete. It contains the parts/alphabet/bond ledger and checked-versus-assumed interfaces; C19's matrix construction and logical signs; all available selector, clock, wrapper, visibility and logical-scale formulas; the C20 n^-4 margin and a simultaneous velocity lower bound; explicit n=2,3 targets; cubic, clock, mixed-process and leakage estimates; the second-order equalization obstruction; conditional scale choices and their resource-dependence limitation; and the accumulation/ring qualifications. Remaining gaps are ranked with concrete calculations.

Verdict: no general incompatibility killing route A was found, but the present pieces are not a counterexample. The spatial spectral core and its compatibility with the mobile wrapper and second-order interface are the principal missing finite ingredients. The full normalized accumulation bond cannot yet be identified. Only the bare selector's z=2 limit is established.

Validation: `python3 check_assembly.py` passed exact Jacobi rational checks for n=2,3, complete low-shell integer resonance checks (49 and 217 conjugate pairs), the scalar equalization determinant check, and all conditional hierarchy inequalities. `checks.json` records inputs and outputs. Markdown file links were checked to exist, numerical table values were compared to the record, and prose contains no em dashes. These are algebra/resource checks, not ED or CFT measurements.

Final concurrent-source read: C21 now reports a product-subspace boundary-moment commutator obstruction, delayed Jacobi filters and a failed unitary-decorated pole ladder; no RESULT file was yet present. C22 still has an unvalidated looped-bus candidate. RESULT identifies these as concurrent reports rather than silently treating them as a completed compiler.

Only work_C23 was modified. No network, external messages, subagents, cluster computation, changes to project STATE.md, or edits to another worker's files. Resume by reading RESULT Sections 6.3, 7 and 11; do not rederive the inherited filters or clock.
