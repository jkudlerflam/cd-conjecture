# C109 progress

[Proved, record] Read the project STATE.md, the full supplied _bound_background.md, REPORT.md sections 9.68 through 9.77, and C104 RESULT.md before starting the analysis. C108 RESULT.md was not present at the initial read.

[Assumed] Work uses the nearest-neighbor, finite-dimensional, unitary infrared setting of the background. All writes and Python scripts stay in work_C109. No heavy numerics are planned.

[Heuristic, investigation] Test discrete symmetry seams, replica permutation seams, vacuum stress and energy-current correlators, and possible charged extensions. The immediate questions are the missing factor of ring length in a seam norm estimate and the velocity powers in the stress correlators.

[Proved, scope] The project-level STATE.md is outside the authorized write directory. A local STATE.md will record the outcome instead.

[Proved, current route] The equal-space stress contributions are B2/(4 pi^2 x^4) for the energy density and B4/(4 pi^2 x^4) for the energy current, where Bp=sum c_a v_a^p. This agrees with C100. The current coefficient is not B2. A local norm controls the correlation value, not its asymptotic coefficient.

[Proved, discrete route] A single on-site symmetry seam has norm difference at most 2M. Its CFT energy shift is 2 pi sum v_a x_a(g)/L. The estimate therefore retains L. A cyclic seam on n replicas is unitarily a ring of length nL, so its energy difference is E0(nL)-nE0(L). The n=2 route is C98 ring gluing.

[Proved, new obstruction] The d=2 critical anisotropic XY family has M=sqrt(1+gamma^2), c=1/2, and v=2 gamma. Thus v/M tends to zero inside the Ising critical line. A velocity floor cannot repair the vacuum-current route.

[Proved, conditional on dynamical infrared matching] A nearest-neighbor Lieb-Robinson estimate gives v_a<=4 e M. Hence one unitary minimal-model factor obeys cv<4 e M, one Ising factor cv<=2 e M, and one standard Z_N parafermion factor cv<8 e M. Sector multiplicity remains uncontrolled.

[Proved, completed] Wrote PROOF.md and RESULT.md with conditional single-factor bounds, the precise seam-budget gap, replica equivalence to gluing, the covariance-coefficient obstruction, the d=2 Ising velocity-floor counterexample, and the requirements for a charged lift. Added local STATE.md because the parent memory is outside the write scope.

[Numerical] All checks in checks.py passed. Maximum matrix dimension was 64. Replica spectra agree to 5.8e-15, Ising twist shifts to 9.4e-15, and independent XY ring energies to 9.8e-15. No heavy numerics were performed.

[Proved, final limitation] Arbitrary neutral-sector multiplicity and general coset contributions remain unbounded by the available arguments. No general A_neutral<=C M f(d) theorem or counterexample is claimed.

[Proved, final audit] Checked the source recursion for the 4 e M Lieb-Robinson ceiling, the XY dispersion normalization, the parafermion order/disorder dimensions, and the continuum current-extension formula. All required files are nonempty and pass the prose-style scan. The general bound remains open for the reasons stated in RESULT.md.
