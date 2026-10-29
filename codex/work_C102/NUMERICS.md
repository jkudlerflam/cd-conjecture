# C102 numerical checks

2026-10-29

[Numerical] Reproduce with `python3 checks.py`, then `python3 summarize.py`, from this folder. The first script writes numerics.json and asserts the consistency checks. The second writes this report. Python, NumPy, and SciPy are the only computational tools. No network commands or cluster jobs were used.

[Numerical] The sparse ring matrix is assembled from the two-site bond in a base-d product basis. XXZ uses total Sz=0 and its known zero-field ground-state sector. Spin-1 SU(2) multiplets all contain Sz=0, so the spin-1 sector calculation includes a global ground state. Generic rotations and gauges use the full space. Independent small full-space comparisons check both sector implementations.

[Numerical] Recorded 136 eigenpairs, with largest matrix dimension 12870. The largest eigenpair residual was 2.286e-11. These are modest local checks. No large-volume computation was performed.

## Estimators and Bethe input

[Proved] For lengths L<K, the bulk-canceling two-size estimator is

    A_pair=(6/pi)[E(K)/K-E(L)/L]/[L^(-2)-K^(-2)].

[Numerical] The second estimator fits E(L)/L=e+a/L^2+b/L^4 and sets A_fit=-6a/pi. Where the exact bulk energy is known, the third estimator is A_L=(6/pi)L[L e-E(L)]. The polynomial fit is a diagnostic and is not a controlled treatment of logarithmic or noninteger-power corrections.

[Assumed, Bethe input] XXZ uses A=pi sin(gamma)/(2 gamma) and delta=cos(gamma). Its bulk energy is evaluated independently by quadrature of

    e(delta)=delta/4-sin(gamma) integral_0^infinity
              sinh((pi-gamma)x)/[sinh(pi x) cosh(gamma x)] dx.

[Numerical] The integrand is implemented with decaying exponentials to avoid overflow. Bethe testing here means this bulk-energy quadrature and the closed velocity formulas, not a finite-size Bethe-root solver. TB and ULS use the inherited thermodynamic data specified in PROOF.md.

## XXZ

[Numerical] Every row uses L=6,8,10,12,14,16, in spin units J=1. All optimized M values are 1/2.

| delta | Bethe A | Pair 14,16 | Polynomial fit | Bulk subtraction L=16 |
| :-- | --: | --: | --: | --: |
| -0.75 | 0.429534903 | 0.426823284 | 0.432156977 | 0.428425857 |
| -0.5 | 0.649519053 | 0.653213647 | 0.649924977 | 0.651125519 |
| 0.0 | 1.000000000 | 1.010446443 | 0.999286231 | 1.004516223 |
| 0.5 | 1.299038106 | 1.316577608 | 1.299055390 | 1.306750746 |
| 0.75 | 1.437574203 | 1.461828945 | 1.441813049 | 1.449139682 |

[Numerical] The decisive delta=-1/2,0,1/2 midpoint defect from bulk subtraction is

| L | Midpoint defect |
| --: | --: |
| 6 | 0.026377321020 |
| 8 | 0.025771916362 |
| 10 | 0.025618211270 |
| 12 | 0.025577307072 |
| 14 | 0.025571698106 |
| 16 | 0.025578090905 |

[Proved conditional] The Bethe limit is 1-9 sqrt(3)/16=0.0257214207425. All sampled signs agree. The finite-size defect is not monotone in L. The bulk Jensen defect is 0.00644964023537; at L=16 the ring Jensen defect is 0.102357202698. Their difference after the extensive subtraction has the positive sign in Eq. (1) of PROOF.md.

## XYZ

[Numerical] The convention is h=sum_a J_a S_a tensor S_a. Full spaces at L=6,8,10,12 are used. The first two bonds are rotated XX chains. The last two are massive anisotropic free XY chains. The middle bond is a rotated easy-axis XXZ chain, with a massive thermodynamic sector.

| (Jx,Jy,Jz) | Thermodynamic A input | Pair 10,12 | Polynomial fit |
| :-- | --: | --: | --: |
| (1, 1, 0) | 1.000000 | 1.019774652 | 0.999086654 |
| (1, 0, 1) | 1.000000 | 1.019774652 | 0.999086654 |
| (1, 0.5, 0.5) | 0.000000 | 0.869411678 | 0.808980800 |
| (1.5, 0.5, 0) | 0.000000 | 0.037839767 | -0.073335779 |
| (0.5, 1.5, 0) | 0.000000 | 0.037839767 | -0.073335779 |

[Numerical] The gapped (1,1/2,1/2) model still gives an apparent A near 0.87 at these lengths. It is not evidence for a nonzero asymptotic Casimir coefficient. The two anisotropic XY endpoints have A=0, while their midpoint is (1,1,0), with A=1. This additional convexity failure uses the gapped extension of A, unlike the principal XXZ example.

## Spin-1 bilinear-biquadratic

[Numerical] Here h_b=X1+b X1^2, X1=S_i.S_{i+1}. For b<1 use L=4,6,8,10. For b>=1 use L=3,6,9, compatible with the period-three ULS vacuum. The last column is a fit diagnostic, not an asymptotic determination.

| b | Thermodynamic A input | Last pair | Polynomial fit |
| :-- | --: | --: | --: |
| -1.0 | 9.424777961 | 10.020640000 | 9.528693527 |
| -0.5 | 0.000000000 | 6.747825436 | 6.426105459 |
| 0.0 | 0.000000000 | 2.616596983 | 2.435141309 |
| 0.3333333333333333 | 0.000000000 | -0.000000000 | 0.000000000 |
| 0.5 | 0.000000000 | 0.453941792 | 4.390686106 |
| 1.0 | 4.188790205 | 4.522578012 | 4.196024607 |
| 1.5 | unresolved here | 4.874365688 | 4.534527216 |
| 2.0 | unresolved here | 5.056935426 | 4.726589741 |

[Numerical] For b=1,1.5,2, the midpoint defect is 0.084608969 from the last pair and 0.073220042 from the fit. Both support the XXZ direction, but three small lengths do not establish the limiting sign. The b=1/2 fit is particularly unstable, and the nonzero Haldane-phase values are finite-size artifacts. AKLT at b=1/3 supplies an exact zero check.

## Seeded critical bonds

[Numerical] Four independent triples were generated from XXZ with random J and anisotropy, a common on-site SU(2) rotation for each triple, and independently chosen scalar and Hermitian coboundary gauges. Each midpoint is the literal average of the two full matrices. Each entire segment is known to be critical. All computations use full spaces at L=6,8,10. Seed=102 fixes the run. This is a constrained critical-family sample, not a claim that unconstrained random Hermitian bonds are critical.

| Triple | J endpoints | delta endpoints | Exact midpoint defect | ED defect L=10 |
| --: | :-- | :-- | --: | --: |
| 0 | 0.905834, 0.622725 | -0.409396, 0.669652 | 0.019217354 | 0.018388488 |
| 1 | 0.778992, 0.844660 | -0.715465, 0.230206 | 0.029174203 | 0.030711084 |
| 2 | 1.007015, 1.328810 | -0.643066, 0.549840 | 0.045550884 | 0.045695189 |
| 3 | 1.124991, 0.935466 | -0.605911, 0.678090 | 0.044695036 | 0.044209763 |

[Proved] Gauge-optimized M for each transformed bond is inherited from its XXZ source, J/2. The JSON field raw_half_width is the unoptimized width and is not substituted for M. [Numerical] Rotations, scalar shifts, and coboundaries preserve the checked ring energies to 1.199e-14 after accounting for the scalar.

[Numerical] Two additional randomly rotated and gauged critical transverse-field Ising bonds use h=-J XX-J(ZI+IZ)/2, with Pauli matrices. Their exact A=J comes from c=1/2 and v=2J. These are an independent free critical class, with on-site field terms retained in the bond.

| Bond | Exact A | Pair 8,10 | Bulk subtraction L=10 |
| :-- | --: | --: | --: |
| random_Ising_0 | 0.689998876 | 0.695114628 | 0.691990312 |
| random_Ising_1 | 1.316636459 | 1.326398192 | 1.320436460 |

## Opposite sign and extreme-bond checks

[Proved] Hopping amplitudes t=exp(+/-i pi/3)/2 have A=1 each, whereas their midpoint has A=1/2. The nonzero straight segment stays critical. L=6,12 remove the endpoint phase twists. [Numerical] All six energies obey E=-A csc(pi/L) within 6.217e-15.

[Proved] The exposed extreme bonds XX and YY with Pauli matrices have E=-L and A=0. Their midpoint has M=1 and A=2. [Numerical] Rings L=4,6,8 reproduce the classical endpoint energies within 3.553e-15.

## Validation

[Numerical] Script validation status: passed. All error checks require <1e-9; the eigenpair residual check requires <1e-8. The analytic proof, not the grid, establishes the curvature sign.

| Check | Value |
| :-- | --: |
| XX_bulk_quadrature_error | 0 |
| XX_ring_formula_max_error | 2.22044604925e-15 |
| AKLT_exact_energy_max_error | 6.2172489379e-15 |
| XXZ_sector_full_error_L6 | 3.10862446895e-15 |
| spin1_sector_full_error_L4 | 3.5527136788e-15 |
| largest_sampled_XXZ_f_second_derivative | -0.0444460548596 |
| free_hopping_ring_formula_max_error | 6.2172489379e-15 |
| classical_ring_energy_max_error | 3.5527136788e-15 |
| max_eigenpair_residual | 2.28606602205e-11 |
| max_unitary_gauge_scalar_identity_error | 1.1990408666e-14 |

[Scope] No broad classification of random critical bonds, finite-volume Bethe-root solution, large spin-1 extrapolation, or proof for arbitrary projector bonds was completed or claimed. None is needed for the exact XXZ obstruction. Only files in work_C102 were modified.
