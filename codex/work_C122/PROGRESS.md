# C122 progress

2026-10-29.

[Work record] Read the project memory, every update of specs/_bound_background.md, REPORT sections 9.68 through 9.81, and the C118/C120 derivations, results, checks and source. Settled constants are inputs. All writes and Python execution remain in work_C122. A local STATE.md will carry the handoff because writes to the project memory are forbidden by this task.

[In progress] Audit the half-time operator comparison, explicit temperature dependence of imaginary-time locality, and the proposed low-energy restriction. Implement full conditional-min-entropy SDPs using NumPy and SciPy with primal and dual feasibility checks on small rings.

[Open] No new universal beta/L estimate is established. The C120 zero-temperature identity already places the full doubling-energy slope in log[mu(Y)/mu(X)].

[Proved] Full operator order is too strong. Every Y<=cX satisfies log c>=beta Tr(sigma V). On XX, the Gibbs variational principle yields log c>=(2beta/L)csc(pi/L)-2log 2, diverging linearly in L at beta=tau L. This is an obstruction to the order route, not the scalar target.

[Proved] Derived an explicit volume-independent expansional constant by counting support-growth events on the folded rung circle. Its logarithm is exponential in beta M. The locality tail has factorial spatial decay at fixed beta, but the width it guarantees grows exponentially in beta M. Derived a local-congruence approximation bound with the required factor q in its error.

[Numerical] Completed 70 full finite-temperature SDPs and 14 ground-projector SDPs in 91.06 seconds. The widest numerical primal-dual bracket for the log ratio is 3.11e-6. Beta-zero and pure-state analytic values lie in the brackets. Widening the seam changes the bounded ground offset while leaving the ground-energy slope. No low-energy truncation or selected-sector Gibbs approximation was made.

[Proved] Quantified the gap-only obstruction to a uniform low-energy restriction. The spectral-tail bound contains d^(2L)exp(-beta delta), and even exact ground projection leaves the Casimir energy difference. Stated the remaining dual-witness transport lemma in PROOF Section 7.

[Work record] Wrote RESULT.md, PROOF.md, CHECKS.md and local STATE.md. Primary literature was checked with the browsing tool. Independent diagonal-SDP, full-order, recurrence and width checks pass. No cluster connection, large local computation, external write or independent agent review was performed.

[Work record] Final validation passed Python syntax, prose restrictions, saved-result counts and every saved primal/dual feasibility check. All ground-state analytic values remain inside the reported brackets. The universal scalar target is left explicitly open.
