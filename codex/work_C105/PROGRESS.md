# C105 progress

2026-10-29

[Assumed] The target is a crossover-free upper bound on A = sum_a c_a v_a, with the infrared matching specified in the background. All work is confined to work_C105. The parent project STATE.md is read-only for this task, so this folder will carry the session memory.

[Proved, bookkeeping] Read the project memory, specs/_bound_background.md with both updates, contrib/fable_bound_modular/README.md in full, and REPORT.md sections 9.68 through 9.75. No delegation was requested or used.

[Heuristic] Separate the raw lattice bond profile from its smooth stress-energy component. The source README already reports alternating XX contributions larger than x^(-2). Also separate a trial-state upper bound on ground energy from the lower bound needed here.

[Proved, infrastructure observation] Both Typhon SSH attempts failed before authentication with "Connection timed out during banner exchange". No job was submitted. No heavy numerical calculation has run locally or on a login node.

[Heuristic] Active analytic directions are positive interval decompositions with general overlaps, telescoping local positive certificates, and the response formula for a local deformation. Numerical scripts will distinguish ramp, quadratic, and periodic smooth weights, and will report raw and filtered local energies separately.

[Proved] Found an exact infinite-half-line spectral representation for the XX deficit, using the Meixner-Pollaczek recurrence and the exact half-line BW relation. Verified the needed spectral relation against Eisler, arXiv:2410.16433, and the polynomial formulas against DLMF. The source links and convention conversion are in PROOF.md.

[Numerical] Evaluated the half-line profile through bond 8193 with two quadrature orders. Runtime was 0.12 seconds. The raw scaled deficit grows strongly on odd bonds while the filtered result approaches 1/(24pi). No large finite matrix was used. The predicted 1/(x log x) staggered asymptotic is explicitly tagged heuristic because the uniform remainder estimate is incomplete.

[Proved] Derived an XX lower-energy certificate from congruence of the uniform hopping operator with its absolute value. It produces a logarithmic allowance for the quadratic weight. Positive interval sums and fixed-range local-marginal SDP certificates have separate rigorous obstructions. Neither obstruction excludes nonlocal certificates.

[Numerical] Ran the four requested models at small sizes, tested ramp and quadratic weights, optimized over all positive interval decompositions, and checked the XX certificate. The proposed centered-energy interacting extension fails for Heisenberg. All completed numerical suites were subsecond. Periodic smooth weights were not needed for the final checks and were not implemented.

[Proved] Replaced reliance on the numerical Heisenberg failure by an exact four-site singlet witness. The weighted singlet block has eigenvalue -sqrt(7)/4, which violates the proposed extension. The proof and constants are in PROOF.md Section 3.

[Assumed, unresolved] The general interacting theorem is not established. The adiabatic route reaches an inverse-energy susceptibility with no available crossover-free bound. No gap estimate or failed candidate is presented as a proof.

[Proved, deliverables] Wrote RESULT.md, PROOF.md, NUMERICS.md and the local STATE.md, along with executable Python checks and an unsubmitted SLURM file. All edits are inside work_C105.
