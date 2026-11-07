# C125 progress

2026-11-06.

[Assumed] Scope and rules are those in specs/_bound_background.md and the C125 request. All new files and Python execution stay in work_C125. No cluster access or network shell command is used.

[Proved, inherited] Read project STATE.md, the full bound background and reconnection note, REPORT sections 9.80 through 9.86, and C118 RESULT.md and PROOF.md. C118's XX bound and the previous route obstructions are inputs, not new results here.

[Heuristic, under audit] Reflection positivity must be distinguished from entrywise positive weights and from positivity of a conventional one-site quantum transfer matrix. A two-site transfer construction is being checked explicitly. The proposed linear Matsubara spectrum bound may fail because lattice commutators have tails outside a Lieb-Robinson cone.

[Heuristic, under construction] A positive crossed spectrum with many low-frequency modes and a tunable ultraviolet crossover appears capable of satisfying the stated microscopic entropy ceiling while leaving the infrared Casimir coefficient unbounded. This will be checked against the full torus scaling form.

[Proved] Completed the positive-transfer audit and an eigenvalue compactness construction for exact even-ring traces L>=4. The Renyi identity and unscaled monotonicity hold in the positive class. A finite Racah recurrence proves the required positive coefficients for every SU(2) symmetric integrable representation. The on-site RP criterion fails for negative XXZ anisotropy and N>2 same-representation SU(N) models.

[Proved] Derived the exact XX crossed product. Locally detected high-frequency transfer gaps grow logarithmically, refuting the proposed global spectrum condition. Their real momentum tilt diverges for every nonzero tilt. Derived the corrected low-frequency and logarithmic high-frequency Matsubara bounds from the C109 Lieb-Robinson estimate.

[Proved] Completed an ordinary finite-dimensional positive spectral countermodel with a full compact-CFT infrared limit. It obeys stronger forms of all hypothetical inputs, has R<=beta M/2, and has A=n v unbounded as n varies. Theta is concave below an arbitrarily long crossover. Explained the missing rotation relation, the separate null edges for unequal velocities, and sufficient extra ground-state scale comparisons.

[Numerical] Completed small ED, finite-Trotter spectra, large free-fermion rings and exact crossed-product tests, all requested Bethe-data signs, countermodel covariance checks, and exact rational checks of the recoupling formulas. Recorded actual failures of R monotonicity and Theta convexity at Delta=-0.9, and negative SU(3) D(4). CHECKS.md records the initial large-beta quadrature issue and its correction, as well as the removed exhausted enumeration cases.

[Work record] Wrote RESULT.md, PROOF.md, CHECKS.md and a local STATE.md. Prepared run_typhon.slurm without submitting it. No outside write, network shell command, Typhon connection, or independent agent/referee review was used. Final limitations are stated in RESULT.md and STATE.md.
