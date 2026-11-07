# C124 checks

2026-10-30. All statements below are numerical checks unless tagged otherwise.

## 1. Reproduction and implementation

[Numerical] Run these Python commands from work_C124 with one BLAS thread. They write their output files only beside the scripts.

```sh
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python3 extra_checks.py
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python3 interval_checks.py
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 python3 transfer_checks.py
```

[Numerical] The environment has NumPy 2.4.6 and SciPy 1.17.1. The python alias is absent; python3 works. The main survey took 11.4 seconds and the interval survey took 25.0 seconds. All computations were small local checks. No cluster calculation was attempted.

[Proved provenance] xx_kernels.py and ed_kernels.py are local copies of C123's finite-grid and symmetry-block routines, with their top-level survey removed. Their formulas and conventions are inherited rather than rederived. The new scripts independently test those kernels against full matrices and use them for the requested model checks. No script imports a module from another work folder or writes outside this one.

## 2. Commuting quantum separators

[Numerical] checks.py tests a physical d=4 site consisting of two virtual qubits. The bond acts only on R_i tensor L_(i+1). The real stoquastic example is

    k=-0.4 X tensor X-0.23 Z tensor Z-0.11 Z tensor I.

The complex example adds 0.19 Y tensor X+0.07 X tensor Z. Both have nontrivial quantum virtual-bond eigenstates. The adjacent-bond commutator is checked directly. Full physical rings at L=2 and 2L=4 are compared to the independently assembled matrix B of PROOF (5), at beta=0.2,2,6 and N=1,2.

[Numerical] All 12 cases satisfy q_s=B_ss and q_D(s,t)=B_st B_ts, with maximum relative residual 3.50e-14. All obey G<=0. The measured G range is [-1.847361858,-1.085823211]. Their one-sector partition ratio is R=0 to numerical precision.

[Numerical] extra_checks.py tests two central sectors with local dimensions 2 and 1, so d=3. The outgoing virtual bond matrices depend on the next central sector and include a complex off-diagonal term. This separately tests the sums over intermediate sectors in (5). Nine cases use beta=0.2,2,6 and N=1,2,3, on rings 2 and 4. The separator residual is at most 1.95e-14. The partition traces agree with the independently formed two-state transfer matrix, and R<=log 2. The G range is [-1.098043029,-0.298156138].

[Numerical] A three-dimensional Fourier rotation mixes the central sectors, so the computational measurement basis is no longer the adapted one. All nine rotated cases satisfy the proved arbitrary-basis bound G<=log 2+2N log 3. Their G range is [-1.098556614,-0.523851951]. These checks do not assert that G<=0 holds in every basis.

## 3. Diagonal and critical-spectator boundary factors

[Numerical] The diagonal counterexample is checked at L=4,16,64,256,1024,8192, a=0.25,1,4 and N=1,2,4,8, giving 72 cases. At a=1,N=1,

| L | C_L/L | G |
| :--- | ---: | ---: |
| 4 | 0.789744104 | -0.635309779 |
| 16 | 1.222119405 | -0.689258651 |
| 64 | 1.344461881 | -0.692903109 |
| 256 | 1.375789573 | -0.693131922 |
| 1024 | 1.383665289 | -0.693146227 |
| 8192 | 1.385965623 | -0.693147166 |

[Proved comparison] The respective limits are 2 log 2 and -log 2. These values follow from PROOF (13), independently of the numerical samples.

[Numerical] transfer_checks.py independently enumerates all diagonal matrix units needed for ||u_S||^2 at L=2,3,4 and verifies that only the all-zero configuration carries the top transfer eigenvalue. It reproduces (12) to floating-point precision. Thus the calculation does not rely on a cutoff that would inadvertently merge the slow excited levels with the ground eigenspace.

[Numerical] The d=4 critical XX plus on-site-qubit example is evaluated in 45 cases, at L=4,16,128,1024,8192, a=0.25,1,4 and N=1,2,4. At L=8192,a=1,N=1,

    spectator boundary factor = 11353.830381229,
    full boundary factor >= 11354.516086497,
    full G = -0.007441897.

[Proved interpretation] The full boundary lower bound uses C_XX>=G_XX at N=1, not an assumption about the XX top eigenspace. The divergent factor is compatible with bounded G because the numerator and denominator thermal contributions cancel.

## 4. XX at large L

[Numerical] Sixteen independent full-Hilbert-space calibrations use L=4,6, beta=0.4,4, N=1,2 and either one or two measured sites. The determinant and dense recorded histories agree to a maximum relative error 8.92e-15. Eighteen overlapping charge-block ED and determinant values of G agree to a maximum absolute error 3.34e-13.

[Numerical] The main free-fermion survey has 46 cases. It uses L=4,16,128,1024,8192, a=0.25,1,4 and N=1,2,4, with one additional N=6 case at L=8192,a=1. Both physical parity sectors and the bordered zero-mode determinants are retained. The inherited C123 bound

    G<=2 log 2+[2 log 2+(10pi/3)a]N

passes every case. At L=8192,a=1,

| N | G/N |
| :--- | ---: |
| 1 | 0.685705292 |
| 2 | 0.839919367 |
| 4 | 0.876720683 |
| 6 | 0.877622285 |

[Numerical] These positive values also reject an extension of the commuting conclusion G<=0 to the noncommuting stoquastic XX chain. They reproduce C123's scale dependence. No new large-L extrapolation is required for the commuting theorem or the boundary-factor counterexample.

## 5. Full-spectrum interacting ED

[Numerical] Ninety main-survey cases use XX, Heisenberg, and XXZ at Delta=-0.9,0.5 on L=4,6 and SU(3) on L=3. For every model and size, a=0.25,1,4 and N=1,2,4 are tested, plus a=1,N=8. Doubled rings have lengths up to 12 for qubits and 6 for SU(3). The largest retained charge block has dimension 924. Every eigenvalue of every sector is included, with color-permutation multiplicities treated explicitly. No ground-projector approximation is made.

[Numerical] At a=1,N=8,

| Model | L | G/N |
| :--- | ---: | ---: |
| XX | 6 | 0.940026731 |
| XXZ, Delta=-0.9 | 6 | -0.303410460 |
| XXZ, Delta=0.5 | 6 | 1.446196265 |
| Heisenberg | 6 | 1.889413953 |
| SU(3) | 3 | 2.979479717 |

[Numerical] Every main-survey row satisfies the inherited pinching bounds and (1). Summing the explicitly enumerated numerator histories agrees with the independently propagated nonselective channel to 3.15e-14. G<=0 fails somewhere in each model family, including XXZ at Delta=-0.9 when a=4. These are failures of an extension outside the commuting hypotheses, not failures of the proved theorem.

[Proved basis observation] On even rings an alternating on-site Z rotation makes the spin-chain off-diagonal matrix elements nonpositive and commutes with the seam projectors. Thus the XX, Heisenberg and XXZ checks also test their stoquastic representatives without changing any recorded weight. No stoquastic claim is made for the SU(3) antiferromagnet.

## 6. Full thermal interval increments

[Numerical] interval_checks.py adds 156 cases, some overlapping the main survey. Free fermions use L=128,8192, a=0.25,1,4 and all N=1,...,6. ED uses the five model families at their larger main-survey circumference and all N=1,...,8, at the same three a values. Tau is held fixed while N changes. At a=1, the observed increments Gamma_n=G_(n+1)-G_n range as follows.

| Model | L | Smallest increment | Largest increment |
| :--- | ---: | ---: | ---: |
| XX, free fermions | 8192 | 0.878263322 | 0.994133441 |
| XX, ED | 6 | 0.942389248 | 1.063371133 |
| XXZ, Delta=-0.9 | 6 | -0.428252952 | -0.261438979 |
| XXZ, Delta=0.5 | 6 | 1.450112863 | 1.557021416 |
| Heisenberg | 6 | 1.900881848 | 1.935315385 |
| SU(3) | 3 | 3.025650986 | 3.378926417 |

[Numerical] transfer_checks.py independently contracts the full temporal maps on XX, Heisenberg, both XXZ examples at L=4, and SU(3) at L=3, for N=1,...,4 and a=1. It does not truncate to the top eigenvector. Fifteen overlapping G values agree with the symmetry-block ED to 1.87e-13. Its directly evaluated numerator moments agree with explicit histories to 2.90e-13. This checks the definitions entering PROOF (15)-(16).

[Open scope] These finite ranges do not prove a uniform c_1(d) in PROOF (17). No explicit unproved universal constant is inferred from the sample maxima. In particular, a top-eigenvalue comparison alone has not been substituted for the finite-N quotient.

## 7. Limitations and files

[Numerical scope] The commuting and direct-sum full-matrix checks use small rings, while their all-size claims are analytic. Interacting ED is too small to establish uniform reconnection at arbitrary L. SU(3) L=6 versus 12 was not computed, and no heavy local computation or cluster submission was performed. All SU(3) vacuum tests here use lengths divisible by three.

[Work record] Raw scalar outputs are results.json, extra_results.json, interval_results.json and transfer_results.json. The run transcripts are run.txt, extra_run.txt, interval_run.txt and transfer_run.txt. The independent parallel attempt was not accessed.
