# C111 progress

2026-10-29

[Proved, work record] Read the project memory, specs/_bound_background.md including both updates, REPORT.md Sections 9.68 through 9.77, contrib/fable_bound_modular/README.md, specs/_smooth_lemma.md, and C105's STATE.md, RESULT.md, PROOF.md and numerical implementation. All edits in this task remain in work_C111.

[Proved, obstruction] The weighted Hamiltonian is not invariant under h -> h + o tensor I - I tensor o. The uniform periodic Hamiltonian and M are invariant. This already prevents a bound for every representative without additional control.

[Proved, route audit] C105's positive-interval obstruction covers arbitrary interval lengths and overlaps, not only equal-length windows. Changing lengths alone does not evade it. A deformed trial vector gives an upper ground-energy bound, hence the wrong direction for the requested deficit bound.

[Proved] The stronger counterexample is complete in PROOF.md. A four-state spectator has a unique gapped product ground state, optimized M = 1, and a uniquely minimizing bond representative up to scalar shifts. Its all-A trial sector has the exact response D = sum [sqrt(delta^2 u^2 + lambda^2 v^2) - delta u]. With lambda = 1/10 fixed, D/Q tends to 1/(200 delta). Stacking with XX gives d = 8, an optimizing bond with M = 3/2, A = 1, and liminf Def/(M Q) >= 1/(300 delta). Uniqueness of the minimizing gauge is asserted only for the four-state spectator.

[Proved] The quadratic modular profile also violates the proposed bound. The exact spectator trial deficit is (lambda^2/delta) log N + O_delta(1), while Q = 2 log N + O(1). This cannot be hidden in a bounded endpoint remainder. The uniform infrared CFT is unchanged. The proposed general matching of the logarithmic deficit to A alone is therefore invalid for arbitrary bond representatives, even ones attaining M.

[Numerical] checks.py passed the bond-spectrum, three-site swap, full spectator energy and gap checks, with largest Hilbert-space dimension 729 across the entire suite. Scalar evaluations verify the predicted smooth response and approach to the logarithmic coefficient. XX, Heisenberg, XXZ(1/2), SU(3) and four seeded random bonds were checked on small rings without sector restriction. The final diagnostic profile includes a second harmonic because a pure first harmonic gives a vanishing response in several small critical rings.

[Proved, infrastructure record] The one Typhon connection attempt returned "Connection timed out during banner exchange". No job was submitted. The large Gaussian XX calculation is prepared in large_gaussian.py and cluster_job.sbatch and remains unrun. No heavy computation was moved locally or onto a login node.

[Proved, completion record] RESULT.md summarizes the counterexample, PROOF.md gives all derivations, NUMERICS.md records checks and their limitations, and STATE.md is the handoff memory. The universal smooth-modulation lemma is disproved as written. The separate velocity-weighted bound remains open. Parent files were not modified.
