# C112 progress

2026-10-29

[Proved, bookkeeping] Read the project STATE.md, the complete bound and smooth-lemma specifications, REPORT.md sections 9.68 through 9.77, the modular README, and C105's derivations and numerical record. All new files stay in work_C112. The project STATE.md is outside the writable task scope, so session memory will be written here.

[Proved, infrastructure observation] The first Typhon connection timed out during banner exchange. A second short connection is being checked. No heavy calculation has run locally or on a login node.

[Proved within CFT] The stress-tensor coefficient is 1/(24 pi). The finite-ring response relative to its uniform ground energy also contains a Casimir subtraction. A single fundamental cosine has zero response in the CFT vacuum. The full derivation is being written in PROOF.md.

[Proved] The stated deficit depends on discrete-divergence additions to the bond, whereas A and M do not. An independent-site gapped model gives an exact unbounded ratio in this gauge class. This does not violate a bound on A. Numerical tests will distinguish this obstruction from tests in a fixed representative.

[Proved, infrastructure observation] The second Typhon attempt also timed out during banner exchange. No cluster connection or submission succeeded. Prepared a Python SLURM array script and larger guarded grids instead of moving the work locally.

[Proved within CFT] Completed two derivations of kappa = 1/(24pi), one from the circle diffeomorphism and one from Virasoro second-order perturbation theory. The two sectors of each CFT and different velocities add to A. The Gaussian and cosine coefficients agree. Verified the stress-tensor transformation convention against Moosavi's primary paper using the browsing tool.

[Proved] Derived the exact XX thermodynamic particle-hole susceptibility and its q^2/(48pi) coefficient for the cosine energy per site. The ring Casimir term must be retained at fixed Fourier harmonic. Fundamental n = 1 modulation gives zero CFT adaptation deficit.

[Proved] Completed the independent-site counterexample for every d >= 2, with M = delta and exact ratio limit t^2/(2delta^2), while A = 0. Added a separate exact second-order counterexample at the gauged XXZ ferromagnetic endpoint. A spin-coherent variational bound establishes failure at fixed amplitude and a fixed bond as L increases.

[Numerical] Completed the small local suites in 3.43 seconds and an independent scalar check in less than one second. CFT comparisons, random full-space bond checks, numerical M brackets, exact onsite spectrum checks and the endpoint perturbation calculation are recorded in the JSON outputs. All assertions passed. Tiny negative subtraction residuals are retained. These tests are not the unrun moderate-size survey.

[Proved, deliverables] Wrote RESULT.md, PROOF.md, NUMERICS.md, local STATE.md, executable checks and the unsubmitted SLURM script. The literal lemma is disproved, the specified-minimizer variant remains unresolved, and the limitation from inaccessible Typhon is stated explicitly. No files outside work_C112 were modified.
