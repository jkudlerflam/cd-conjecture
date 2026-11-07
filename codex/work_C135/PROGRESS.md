# C135 progress

[Proved, inherited conventions] Work uses the unhalved trace norm and gauge-optimized half-spread M. All edits and computations remain in work_C135. The parent STATE.md was read but cannot be edited under this task's explicit restriction.

[Numerical, preparation] The installed environment contains NumPy 2.4.6, SciPy 1.17.1, CVXPY 1.9.2, Clarabel and SCS. No installation is needed.

[Proved, audit] C132's correction applies: the conditional-entropy continuity argument does not disprove a uniform trace-norm route. F1 to F4 require no ring saturation assumption.

[Heuristic, under verification] Sine-window smoothing of a finite LTI fermionic covariance may give a uniform quadratic-chain LTI rate. An explicit weighted XX dual certificate is being derived independently of optimal-dual numerics.

[Assumed workflow] Run small falsification computations locally. Put large dense SDPs, large ring scans and extensive multistart searches in resumable SLURM-ready scripts. Prepared jobs will be distinguished from executed computations.

[Proved] F1 to F4 and support completion are written in PROOF.md sections 1 and 2. Ordinary seam-crossing windows of a twisted ring need not be LTI. A twisted-translation invariant ground state in a fixed seam-free linear ordering is LTI.

[Proved] Ring saturation is false for the real reflection-symmetric XX bond with field b=1/2. E_LTI(3)=-5/12, while the best untwisted ring at length >=3 has density -3/8. An exact rational Toeplitz covariance at n=4 gives -7553/20000, below -3/8 even when all on-site symmetry twists are allowed. Integer principal minors prove covariance feasibility. A uniform bound on rings L>=8 closes the infinite tail, so this is not a finite-ring numerical comparison.

[Proved] A finite-covariance sine taper proves the LTI rate for every quadratic qubit bond, with K=pi²/2. For XX an explicit telescoping X made of weighted bond translates has K=pi. Its positive decomposition uses extended fermionic modes. This builds on the weighted XX inequality already in C105, with a new LTI application.

[Numerical] The analytic XX certificate passed one-particle checks for n=2 through 32 and full spin telescoping checks through n=8. Small SDP benchmarking reaches n=8 quickly after magnetization reduction. The full local scan is running with n<=8 for symmetry-reduced models, n<=6 for generic bonds, and rings through L=10 or 12. Larger sizes are reserved for SLURM.

[Numerical, completed] The scan now has 62 bonds, including 30 random bonds, 738 ring energies, 345 LTI values at n<=8, 48 additional free covariance values at n=9,10, 372 hull-distance upper bounds and 161 distance lower bounds. Generic real bonds reach n=6. Generic complex bonds reach n=4, plus one n=5 entry; their first n=6 solve was stopped as too expensive locally. All larger requested entries are marked NA rather than inferred. The free spectra extend to L=20 without exponential diagonalization.

[Numerical, completed] A real covariance formulation replaced the poorly conditioned early complex covariance SDP. Full-spin data remain archived in the result records. Independent covariance and spin computations agree. The largest problematic real random n=6 solve was refined, and independent direct gauge-dual SDPs tightened loose scale bounds. Eleven final LTI entries and 28 distance lower entries remain optimal_inaccurate and are explicitly exploratory.

[Numerical, completed] Four-start eight-angle searches ran with 100 evaluations per start for n=2 to n=4 and 300 for each ring-upper objective at L=6,8,10. All exhausted their budgets. The largest-n all-angle search remains in the unrun SLURM continuation. The local maxima do not establish asymptotic drift or an unbounded family.

[Proved and numerical, final audit] PROOF.md, RESULT.md, NUMERICS.md, TABLES.md and STATE.md are written. Independent spin checks pass the exact counterexamples, the covariance smoothing identities and the explicit XX certificate. Table validation finds no reversed brackets and no corrected LTI lower bound above a compatible computed ring energy. Python syntax checks pass. Four SLURM-ready continuation scripts are present; none was submitted. No file outside this folder was edited.
