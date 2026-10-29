# C103 progress

2026-10-29.

[Proved work record] Read the supplied background including the round-one update, project memory, REPORT.md Sections 9.68 to 9.71, and the C98 to C101 results. Read the free-fermion theorem and the C101 normalization and stiffness derivations. All writes are confined to work_C103.

[Open] The general inequality A <= pi M log_2 d is not established. Standard infrared CFT matching alone has not supplied a microscopic bound on its coefficient.

[Proved analytic progress] The raw Casimir estimator is not a nonincreasing function of length, even for a nearest-neighbor quasi-free chain with d=8 obtained by blocking a range-three hopping chain. An explicit dispersion and exact finite-ring energy give approach to A from below.

[Proved analytic progress] The normalized thermal energy-current variance is not a nonincreasing function toward the infrared, already for the nearest-neighbor XX chain. Its first low-temperature correction is negative. This does not contradict a c-theorem for a differently constructed function.

[Assumed; proved conditional implication] The useful positive mechanism remains a complete conserved-current description with a lower bound on the stiffness metric and a microscopic diamagnetic budget. C101 provides actual examples. Neutral sectors and the sharp microscopic budget remain unaccounted for in the general setting.

[Proved completed work] RESULT.md and PROOF.md now give the derivations, all additional assumptions, and the precise scope of the failed monotonicity proposals. A positive-current-moment improvement of the optical bound is included. At Heisenberg, the uncorrected diamagnetic budget exceeds the sharp budget, exposing a concrete gap in the strongest positive route.

[Numerical] python3 checks.py passes. It verifies the blocked free-ring energies and dispersion, the XX thermal-current integral, the C101 formula evaluations, and the optical-moment and twist-curvature identities on small XXZ rings. CHECKS.md records the results. No large computation or network shell command was used.

[Open final status] The conjecture is neither proved nor falsified. The scale-free conditional argument requires current completeness and an independent sharp microscopic optical budget. No unknown crossover scale is inserted into the final conjectured inequality. Local STATE.md records the handoff, since the parent project memory lies outside the permitted write folder.
