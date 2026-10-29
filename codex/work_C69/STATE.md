# C69 local state

2026-10-28. Finished the investigation requested for work_C69. All generated files remain in this folder. The parent project memory cannot be updated within the allowed write scope.

[Proved conditional] Equal-spatial-width conformal boundary-state quenches have half-line entropy rate pi A/(12 ell), with A=sum c_a v_a. Ancillary SIE gives |Delta S|<=18 M log d Delta t for the optimized bond half-spread M, including arbitrary local coboundaries. Controlled lower matching at ell_0 gives A<=[216/(pi alpha)]ell_0 M log d. RESULT Section 3 also gives a finite-time version with an explicit error budget.

[Proved restricted obstruction] Scale-covariant entropy rates carry an inverse length. Multiplying by that length loses it on the SIE side too. This is a failure of this inference, not a general no-observable theorem and not a counterexample chain.

[Proved crossover theorem] If s(T_*)>=alpha pi T_* sum(c_a/v_a)/3, then c<=3 xi_H log d/(pi alpha)<=3 xi_max log d/(pi alpha), xi_H=[c/sum(c_a/v_a)]/T_*, xi_max=v_max/T_*. In the prompt's literal pi/6 convention, C'=6/(pi alpha) and xi_cross=u_max/T_*. Keep the central-charge and velocity conventions consistent. Thermal wavelength is not automatically an arbitrary operator's correlation length.

[Proved spatial theorem] Lower CFT matching for S(2n)-S(n) at n=xi_cross gives c<=3 xi_cross log d/(alpha log 2), using only subadditivity and finite on-site dimension.

[Numerical] checks.py is a small standard-library formula check. It passed and generated CHECKS.md. No new many-body computation, network request, or external source verification was performed.

[Assumed/open] Universal A<=C M log d and universal exclusion of all observables remain unresolved. A sufficient missing lower bound is E(h)>=kappa A for the single-bond ancillary entangling capacity. No uniform preparation length is implied by the existence of an infrared CFT.

File map: RESULT.md is the full derivation; PROGRESS.md is the completed work record; checks.py reproduces CHECKS.md.
