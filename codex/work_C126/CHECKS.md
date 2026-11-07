# C126 checks and reproduction

[Numerical] Run these commands from this directory.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 algebra_checks.py

The first writes checks_local.json and TABLES.md. The second writes algebra_checks.json. Both completed successfully. The first uses only small sparse ED sectors and two-site contractions. The largest local sector is 184756 states. The local ED suite took about two seconds in this environment.

[Numerical] checks.py explicitly asserts the fractional-contraction identity, every reported variational inequality, the Heisenberg b=4a^2/3 relation, the 3/32 lower bound, and eigenpair residuals below 1e-8. The actual maximum residual is 1.06e-11. It uses 17 deterministic initial profiles for the fractional angles and records the final gradient. The best-found value is not a global optimization certificate.

[Numerical] algebra_checks.py independently checks the following quantities.

| Check | Maximum error or outcome |
| :--- | :--- |
| Actual on-site seam conjugation, complex oriented bonds | 8.89e-16 |
| Balanced-sector energies versus small full-space ED | 1.34e-15 |
| Heisenberg two-site correlation tensor | 3.34e-16 |
| Heisenberg fractional coefficient | 3.89e-16 |
| Fractional functional versus an eight-spin trial vector | 1.78e-15 |
| Odd-prime cyclic logarithm on Ising cat marginals | 1.67e-16 |
| XX odd and signed-even sewing, through L=30 | Below 1e-10, asserted |
| Paired-XY signed-cycle trial versus an eight-spin physical energy | 7.51e-15 |
| Ordinary Ising spin-flip defect | 2J |
| Unrestricted Ising replica-permutation defect | 0 |

[Proved or assumed, distinguished] The fractional formula and its Heisenberg lower bound are analytical statements in PROOF.md. The ED ground-sector choices at larger sizes use the standard balanced-sector property of these antiferromagnets. The full-space checks verify that choice only at the stated small sizes. Lanczos residuals check eigenpair accuracy, not a proof that every other conserved-charge sector was compared at the largest sizes.

[Pending] run_typhon.slurm runs checks.py --cluster and the algebra checks inside SLURM. The cluster flag is rejected outside a SLURM allocation. The extended suite includes SU(3) at L=9 and 18, using a matrix-free Hamiltonian for the 17,153,136-state 18-site sector. The file was prepared but not submitted, and that result is not included in any table.

[Work limits] No network command, SSH connection, external write, or independent referee run was performed. Primary literature was consulted only with the browsing tool. A project-level memory update is forbidden by this task, so STATE.md here supplies the handoff.
