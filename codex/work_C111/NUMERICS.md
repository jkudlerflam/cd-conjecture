# C111 numerical record

[Proved, reproduction] Run the small checks from work_C111 with

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py

The only output data file is checks_results.json in this folder. The script uses NumPy, SciPy and Python's standard library. Random-bond generation has fixed seed 111. Every many-body calculation uses the full Hilbert space without a charge-sector restriction. The maximum dimension is 729, and the complete final run took 0.102790 seconds.

[Numerical] All assertions passed. The largest eigenvector residual was 3.079887e-12. This is a residual check, not a certified lower spectral bound. The counterexample itself does not rely on numerical eigenvalues.

## Algebra checks

[Numerical] The spectator's 16 by 16 bond matrix has minimum and maximum eigenvalues -1 and 1 for delta = .05, .01, .001. The independent stacked bond has minimum and maximum eigenvalues -1.5 and 1.5. Gauge optimality is proved by explicit marginal witnesses in PROOF.md, not inferred from these spectra.

[Numerical] Full spectator rings at L = 3 and 4 give E0 = (-.75-delta)L and first excitation gap 2 delta, within 2e-11. Weighted-ring energies obey the all-A variational comparison in every tested case. The three-B-site matrix P12+P23 has minimum eigenvalue -1.

## Smooth response

[Numerical] For w_i = 1+.3 cos(2 pi i/L), the spectator trial deficit was evaluated by the rationalized expression to avoid cancellation. At L = 16384, the results for D/Q are

    delta        measured D/Q          limiting lambda^2/(2 delta)
    .05          .09999999867          .1
    .01          .49999986552          .5
    .001         4.9998668475          5
    .0001        49.8676439003         50

Here lambda = .1. The slower convergence at smaller delta is visible and is consistent with the expansion parameter lambda |v|/(delta u). The proof takes the smooth limit at each fixed delta and requires no common matching scale.

[Numerical] For w_j = j(N-j)/N, the last logarithmic slope, from N = 16384 to 65536, is

    delta        measured slope         limiting lambda^2/delta
    .05          .1999637266            .2
    .01          .9991213508            1
    .001         9.9140662555           10
    .0001        92.4499890692          100

These are one-dimensional scalar sums of the exact formula, not large Gaussian or many-body computations. At the smallest gap the accessible slope has not reached its asymptotic value. Equation (6.2) of PROOF.md establishes the limit analytically.

## Standard and random bonds

[Numerical] The final profile is

    w_i = 1 + .25 cos(2 pi i/L) + .15 cos(4 pi i/L).

The second harmonic avoids a vanishing response that occurs for a pure first harmonic in several of these small critical rings. The reference deficit is E0(H_1) sum w_i/L - E0(H_w), so it removes the constant-weight finite-size offset. The following are the ratios of that deficit to M_raw Q at the largest tested size.

    model            L       ratio
    XX               8       .01663505
    Heisenberg       8       .02750763
    XXZ(.5)          8       .02246963
    SU(3) swap       6       .05056613
    random_0         8       .15205041
    random_1         8       .01690693
    random_2         8       .14733860
    random_3         8       .09388785

[Numerical, limitation] For the random bonds, M_raw is the centered raw half-spread, normalized to one. No gauge minimization was performed, and no conformal infrared theory is asserted for them. These tests verify the implementation and measure finite-size response. They do not support a uniform-in-bond constant or test every infrared sector. The exact family in PROOF.md is what disproves such a constant.

## Cluster work

[Proved, infrastructure record] The connection command was

    ssh -o LogLevel=ERROR -o BatchMode=yes -o ConnectTimeout=8 typhon 'hostname'

It exited with code 255 and output

    Connection timed out during banner exchange

[Proved, unrun work] large_gaussian.py prepares quadratic open XX calculations at N = 512 through 16384 using a tridiagonal one-particle matrix. It also reports the exact spectator trial deficit and a lower bound on the total stacked deficit. cluster_job.sbatch requests one compute task and one hour. A valid Typhon partition must be selected on submission. No job identifier exists because the connection failed before submission. The large script requires SLURM_JOB_ID and SLURMD_NODENAME and rejects hostnames containing "login". The local code was syntax-checked without executing that calculation.
