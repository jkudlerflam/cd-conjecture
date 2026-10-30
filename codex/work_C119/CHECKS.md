# C119 checks

2026-10-29.

[Numerical] Run inside this folder:

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py

The output is checks_small.json. The completed run took less than one second.

[Numerical] The C106 spin Hamiltonian was built directly from its 16 by 16 bond and diagonalized on lengths 3 and 4, with largest Hilbert dimension 256. Full traces at beta=0, 0.7 and 5 agree with the two-grid parity formula to 4.9738e-14. This checks the common parity projection and the closing bond.

[Numerical] The XX correction in PROOF.md equation (5) was checked at y=0.5,1,2,4 and L=32,64,128,256. After subtracting the predicted L^-2 term, the remainder scales as L^-4 until roundoff from bulk subtraction becomes visible. For example at y=1 its predicted coefficient is -0.1949786389, and L^4 times the remainder at L=32,64,128 is -3.7293,-3.6339,-3.6112. The largest length is not used to infer an additional asymptotic coefficient.

[Numerical] The C106 curvature correction was checked at a=1.2 and epsilon=2sin(pi/8)/sqrt(a^2-1), so all tested L=64,128,256,512 align the Fermi points. At tau=1 and 2 the predicted L^-2 coefficients are -5.3998627643 and 30.1021344810. The measured coefficients at L=512 agree within 0.01 times max(1, absolute coefficient), the explicit check tolerance. The residual follows L^-4 over the resolved range. These parameters test the correction formula, while the separate small-epsilon scan tests the nonuniform crossover.

[Numerical] Slow XX contributions computed from the modular channel are

| y | r_XX(y) | leading 4 exp[-pi/(2y)] |
| :--- | ---: | ---: |
| 0.25 | 0.007448881679 | 0.007469770927 |
| 0.125 | 0.00001394929646 | 0.00001394936942 |
| 0.0625 | 4.864622684e-11 | 4.864622684e-11 |

At beta=L/v, the value is 0.6326341534 for v=1,0.1,0.001. This checks the distinction between fixed beta/L and the inverse-velocity scale.

[Numerical] The two-spin-XX stack at tau=1 has ratios 0.3657304331, 0.5577695327, 0.6169021483 and 0.6286622745 for epsilon=0.5,0.1,0.02,0.005 and velocities 1/(1+epsilon),epsilon/(1+epsilon). These agree with the C116 finite-size data within their finite-size corrections.

[Numerical] C106 crossover checks use L=32,128,512,2048, epsilon=0,1e-5,0.001,0.01,0.1, a=1.2 and beta=L/M. At epsilon=0 the value tends to 0.66326794, not the independent-spin-stack value 2r_XX(1)=1.26526831. This is the shared-parity effect. At unaligned finite epsilon, negative values occur and are allowed. No positive lower bound was assumed.

[Numerical] For the massive a=0 model with t=1, m=0.1 and tau=0.5, R at L=32,128,512 is 0.568344052,0.0151449200, approximately 5.8e-11. The last value is at the scale of numerical bulk-subtraction error and is not used to resolve an exponential prefactor. Smaller masses retain order-one crossover contributions.

[Numerical] At fixed y=1, the compact-boson ratios for K=1,10,100,1000 are 0.6326341534,1.8519234821,3.0032160286,4.1545085751, consistent with the derived logarithmic growth. For the XXZ sequence gamma=pi-eta, the theta limits at tau=1,4,16 are 0.4441688024,1.2654086932,1.9586593040. At eta=0.001 the corresponding CFT evaluations are 0.4443974952,1.2655681144,1.9588184010. These are CFT evaluations, not interacting finite-ring diagonalizations.

[Numerical] The candidate pointwise scaling ceiling r_XX(y)<=pi y/4 was checked at 200 logarithmically spaced y from 0.02 to 100. The largest excess was 7.11e-15, numerical roundoff. This is not an analytic proof of that ceiling.

[Numerical] A deterministic 360-point C106 scan uses L=4,8,16,32,64, eight epsilon values from 1e-5 to 0.5, a=1.001,1.2,2, and beta=sL/M for s=0.1,1,4. All points obey R<=3log 2+(10pi/3)s. The largest left-minus-right difference is -3.126638926. The rigorous proof, not this finite scan, excludes divergent parameter sequences.

[Work record] The extended scan in run_typhon.slurm uses no large matrices but samples many more sizes and parameters. It is guarded by SLURM_JOB_ID, was not run locally, and was not submitted. No Typhon connection was attempted.
