# C95 checks

[Numerical] Small quadratures and closed-form checks, with t=1. These do not test the full d=13 thermal trace.

Half-filling density-mode coefficient: 0.07857579749228.

| T | s | Material density | Density-mode coefficient |
| :-- | :-- | :-- | :-- |
| 0.1 | 0.5 | 0.4531311445 | 0.05950989653 |
| 0.1 | 0.636619772 | 0.499159968 | 0.07817453349 |
| 0.1 | 0.75 | 0.534298529 | 0.09349008419 |
| 0.03 | 0.5 | 0.4541737461 | 0.06001669039 |
| 0.03 | 0.636619772 | 0.4999249471 | 0.07853959854 |
| 0.03 | 0.75 | 0.5348850323 | 0.09371134544 |
| 0.01 | 0.5 | 0.4542642343 | 0.06006069047 |
| 0.01 | 0.636619772 | 0.499991666 | 0.07857177476 |
| 0.01 | 0.75 | 0.5349364359 | 0.09373158876 |
| 0.003 | 0.5 | 0.4542745159 | 0.06006569036 |
| 0.003 | 0.636619772 | 0.49999925 | 0.07857543544 |
| 0.003 | 0.75 | 0.5349422813 | 0.09373389842 |

[Numerical] At v=N, the adjacent-sector contact jump is O(1/N). The exact formulas are evaluated directly.

| N | Crossing tilt | N times contact jump |
| :-- | :-- | :-- |
| 15 | 0.610116148317 | 0.334810697525 |
| 30 | 0.623453021115 | 0.343030521344 |
| 60 | 0.630056102101 | 0.347165637797 |
| 120 | 0.633342671909 | 0.349240687161 |
| 240 | 0.634982382048 | 0.350280233271 |
| 480 | 0.635801364228 | 0.350800530268 |
| 960 | 0.636210639684 | 0.351060812076 |

[Numerical] Exact two-particle-hole form-factor sum for the uniform contact. The N=4,6 and (N,v)=(9,4) values agree with the signed sums of C94's independently computed gap-resolvent rows.

| N | Vacancies | Signed contact-kernel sum |
| :-- | :-- | :-- |
| 3 | 3 | 0.02391975308642 |
| 4 | 4 | 0.024536904631338 |
| 6 | 6 | 0.024972535151528 |
| 9 | 4 | 0.0082386441046654 |
| 15 | 15 | 0.025263029011575 |
| 30 | 30 | 0.025304385619553 |

[Numerical] Gapped two-level bath diagnostic, J=0.01. R2 is the Taylor coefficient at J=0 multiplied by J squared. The remainder is not uniformly smaller than R2 as T/J decreases.

| T/J | R / J squared | (R-R2) / abs(R2) |
| :-- | :-- | :-- |
| 10 | -0.0609393962137 | 0.0249696565621 |
| 1 | -0.0472106499578 | 0.244629600675 |
| 0.1 | -0.000841816281982 | 0.986530939488 |
| 0.01 | -2.56808351164e-23 | 1 |

[Numerical] All checks passed. The signed contact-sum bound, density-mode bounds, and two-level identities are proved separately in PROOF.md. No spatial first-moment bound or all-order d=13 estimate is inferred from these checks.
