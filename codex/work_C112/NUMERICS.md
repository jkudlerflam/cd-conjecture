# C112 numerical record

[Numerical] All computation is Python with NumPy and SciPy. The completed local suite is lightweight. The main suite took 3.43 seconds and the independent scalar check took less than one second. The largest interacting sector had 2520 basis vectors, the largest unrestricted generic-bond Hilbert space had 256, and the largest free one-particle matrix had 128 rows. No heavy local calculation was performed.

## Reproduction

[Numerical] From this folder run

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 smooth_checks.py --suite all --output small_results.json
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 scalar_checks.py

The second command reads small_results.json and writes scalar_results.json. additional_results.json is an earlier standalone output of the additional suite, also included in the final small_results.json. It is retained as a work record.

[Numerical] The main output contains 24 XX cases, 162 XXZ cases, 36 generic-bond cases, 15 additional SU(d) and gauged-ferromagnet cases, and the exact independent-site counterexample checks. The 36 generic cases come from nine bonds, two representatives and two profiles. The seeded draw is reproducible, using seed 1122026. All matrices defining those representatives are stored in the output.

[Definition] For every finite-ring calculation the deficit is

    D_L = E_0(H)*mean(w)-E_0(H_w),
    G_discrete = sum_i (w_(i+1)-w_i)^2/w_i.

The ratio is D_L/(M G_discrete). This removes the uniform finite-ring Casimir energy. The CFT prediction retains the harmonic-mean-length correction. G_discrete is used for the requested lattice ratio, while a spectral derivative supplies the continuum functional used in the CFT comparison. They are not silently interchanged.

## Solvers, sectors, and verification

[Numerical] General bonds use sparse full-space Hamiltonians with no symmetry restriction. Sparse eigensolvers request the two lowest eigenvalues with a deterministic nonsymmetric starting vector. The stored residual is ||H psi-E psi||. Matrices of dimension at most 64 use dense diagonalization.

[Numerical limitation] Standard XXZ checks use the balanced magnetization sector at even L. For Delta <= -1 the exact polarized ground energy is also included, and the bondwise lower bound proves it is the full ground energy. For the other standard XXZ and SU(d) models, the quoted balanced-sector results are numerical sector checks, not a new proof that every modulated ground state stays in that sector. Generic bonds and the small independent-site spectrum checks are unrestricted. The small-amplitude gauged ferromagnet's balanced-sector selection is established analytically in PROOF.md Section 5.2.

[Numerical] The free XX calculation uses the half-filled antiperiodic fermion vacuum at lengths divisible by four. Its uniform spin-ring energy is -csc(pi/L). The thermodynamic particle-hole integral is evaluated independently in scalar_checks.py, verifying the closed form before taking the long-wavelength limit.

[Numerical] M is known analytically for the standard spin, swap and onsite examples. For generic bonds, minimization over all traceless Hermitian one-site divergence generators supplies an upper bound. A dual matrix with equal partial traces, zero trace and unit trace norm supplies a lower bound. Both are stored; the largest width was 6.10e-5. The brackets have floating-point accuracy and are not certified interval-arithmetic results. The random table reports ratio intervals where that width is visible.

[Numerical] Assertions checked nonnegative adaptation deficit to numerical tolerance, primal/dual ordering for M, the exact onsite formula against full-space ED, the independent XX susceptibility integral, the ferromagnetic trial-energy bound, and invariance of uniform ground energies under the numerical gauge change. All passed. The largest recorded many-body residual was 4.40e-14. The SU(4) four-site deficit is -8.88e-16 from subtraction; this is recorded rather than clipped to zero.

[Proved, guard check] Invoking smooth_checks.py with --cluster outside a SLURM allocation exits before calculation with

    Heavy calculations require a SLURM allocation.

This expected rejection was tested. No should_not_exist.json file was created.

## Prepared but unrun Typhon work

[Proved, infrastructure observation] Two calls using ssh -o LogLevel=ERROR with a 12-second connection timeout and batch authentication failed before login with

    Connection timed out during banner exchange

Neither cluster discovery nor submission succeeded. No remote directory, job ID, or result is claimed. The task explicitly allowed Typhon SSH; no other network command was used. A browsing tool was used for the primary CFT source.

[Assumed, proposed execution] cluster_job.sbatch is a Python SLURM array script. When submitted from this folder on Typhon, its four tasks run the XX, XXZ, random, and additional suites. It requests one CPU, 16 GiB and eight hours per task. No partition was guessed because discovery failed. The resource request has not been benchmarked and may require adjustment to the available CPU partition.

[Assumed, planned grid] The guarded larger suite contains free XX matrices through L = 1024; standard XXZ rings at L = 16,20,24; 24 generic bonds at each d, with d = 2 lengths 12,16, d = 3 lengths 7,9, and d = 4 lengths 6,7; SU(3) rings at L = 12,15 and SU(4) rings at L = 8,12; and the gauged ferromagnetic endpoint at L = 16,20,24. The generic-bond cases remain unrestricted. They may be costly, and none of these sizes is represented as completed data.

[Assumed, follow-up] Even a completed version of this grid would only test an envelope over sampled bonds. Establishing continuum behavior requires increasing L relative to each bond's crossover scale. Near Delta = -1, the standard small-ring response strongly undershoots the CFT value, illustrating this issue. The exact counterexamples do not depend on such an extrapolation.
