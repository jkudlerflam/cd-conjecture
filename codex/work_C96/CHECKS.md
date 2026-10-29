# C96 checks

[Numerical] These are finite checks of proved formulas and bounds. No uniform tagged-kernel estimate or full d=13 trace is inferred.

| N | v | Bath dimension | Kernel discrepancy | First moment after c_dens subtraction |
| :-- | :-- | :-- | :-- | :-- |
| 3 | 2 | 6 | 8.327e-17 | 0.153084597654 |
| 3 | 3 | 10 | 1.388e-17 | 0.108734414363 |
| 4 | 4 | 35 | 9.021e-17 | 0.109951240625 |
| 5 | 5 | 126 | 7.182e-16 | 0.110084393219 |
| 6 | 6 | 462 | 2.501e-15 | 0.110257956948 |
| 7 | 7 | 1716 | 2.09e-15 | 0.110213935979 |
| 8 | 8 | 6435 | 2.907e-15 | 0.110223764747 |

[Numerical] Complete determinant basis, eigenvalue, and form-factor checks are included through N=4. The largest full-basis error is 6.788e-15. At every checked size the contact spectral weight agrees with the projector identity.

| J | T | Resummed connected coefficient | Static covariance upper bound | Operator upper bound |
| :-- | :-- | :-- | :-- | :-- |
| 0.03 | 0.2 | 0.124917048217 | 0.124929759775 | 0.150100142118 |
| 0.03 | 0.03 | 0.122187147846 | 0.12228778971 | 0.146265232116 |
| 0.03 | 0.003 | 0.07532472274 | 0.0757668011067 | 0.0780051389195 |
| 0.03 | 0.0003 | 0.073411598272 | 0.0738622133487 | 0.0752105435816 |

[Numerical] The material sources include a complex Pauli matrix and do not commute. Matrix positivity and the imaginary-time covariance inequality pass at five times for each row. The independent thermal determinant formula (33a) agrees with direct bath correlations at N=v=3 and T=0.5,0.1, with maximum error 1.419e-15.

| N | First moment of periodized 1/(1+r squared) | Weaker norm, b=3 |
| :-- | :-- | :-- |
| 32 | 9.609263146 | 39.03798417 |
| 64 | 10.99315647 | 39.12445179 |
| 128 | 12.37884897 | 39.14612054 |
| 256 | 13.76499277 | 39.15154099 |
| 512 | 15.15124948 | 39.1528963 |
| 1024 | 16.53753443 | 39.15323515 |

[Scope] This table concerns the explicit test kernel in PROOF equation (29), not the tagged contact kernel. Its first moment diverges logarithmically while the weaker norm stays bounded. A deterministic scalar-cell comparison also satisfies (19).

| T / J squared, t=1 | Exact R / J squared | Quadratic R2 / J squared |
| :-- | :-- | :-- |
| 10 | -0.125778065754 | -0.125 |
| 1 | -0.132788725896 | -0.125 |
| 0.1 | -0.188568479789 | -0.125 |
| 0.01 | -0.243062278506 | -0.125 |
| 0.001 | -0.249300603131 | -0.125 |

[Numerical] The fourth-order coefficient at T=0.2 is -0.0312507093417; the finite-difference value at J=0.02 is -0.0312491455384. The exact residual obeys its uniform bound on every tested row.

[Check correction] The first run stopped at the provisional assertion weaker_norm<36, with value 39.0379841724 at N=32. That test threshold omitted the absolute row sum. The analytic bound is 36+pi*coth(pi)<40, as proved by telescoping the second difference. The corrected check passes. No contact-kernel theorem depended on the provisional threshold.

[Reproduction] From work_C96 run python3 -B tagged_kernel.py and python3 -B checks.py. Maximum bath dimension is 6435. No parent script is executed, no large interacting calculation is run, and all outputs stay in this folder.
