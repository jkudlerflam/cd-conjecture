# C69 progress and handoff

2026-10-28. Investigation complete. RESULT.md contains the derivations and explicit assumptions. All writes are confined to work_C69. No network access, delegated agents, large computation, cluster job, commit, or push was used. The project-level STATE.md was read first and was not edited because it is outside the writable folder. Local STATE.md records this handoff.

## Results

[Proved conditional on the stated CFT preparation] A product of regularized conformal boundary states with common spatial extrapolation length ell has half-line entropy increment

    Delta S(t) = sum_a (c_a/6) log cosh(pi v_a t/(2 ell)),

and late-time rate pi A/(12 ell). The strip-map derivation is in RESULT Section 3. Equal temporal preparation lengths instead give a rate proportional to c. Ordinary common-temperature reservoir currents lose the velocity weighting.

[Proved from the published ancillary SIE theorem] The entropy increment across one lattice cut is bounded by 18 M log d times elapsed time. The optimized half-spread M is legitimate because a bond coboundary adds only terms local to the two sides of the cut. If the lattice entropy growth matches the CFT rate from below with factor alpha at ell_0, then A <= [216/(pi alpha)] ell_0 M log d. A finite-time version needs matching only over T_obs = q ell_0/v_min, with alpha = 1-2 log 2/(pi q)-delta for relative increment error delta. No infinite plateau or derivative convergence is needed in that version.

[Proved restricted obstruction] For a fixed-shape scale-covariant entropy functional, dilating all lengths and times by lambda divides its physical entropy rate by lambda. Multiplication by the preparation length also multiplies the SIE upper bound. An explicit family of formal entropy profiles obeys the rate, entropy-density, and speed inequalities for arbitrarily large A by increasing the preparation length. It is not a microscopic counterexample and does not exclude every possible observable.

[Proved crossover theorem] Under the one-temperature lower estimate s(T_*) >= alpha (pi T_*/3) sum c_a/v_a, finite-dimensional entropy capacity implies c <= 3 xi_H log d/(pi alpha), where xi_H = [c/(sum c_a/v_a)]/T_*. The explicit upper length xi_max = v_max/T_* gives the same bound with xi_H replaced by xi_max. The literal pi/6 convention in the prompt gives C' = 6/(pi alpha) with xi_cross = u_max/T_*. The conversion between chiral and nonchiral counting, or u=v/2, is stated in RESULT Section 7.

[Proved spatial variant] If S(2n)-S(n) >= alpha (c/3) log 2 at n = xi_cross, then c <= 3 xi_cross log d/(alpha log 2). A one-site increment assumption improves this to c <= 3 log d/[alpha log(1+1/xi_cross)]. No thermal-to-spatial crossover identification is assumed.

[Proved lattice estimate; assumed continuum visibility] A conservative overlapping-bond LR count gives v_LR = 6 e M. If each sector has a nonzero algebraic real-time signal in suitable local observables, v_a <= v_LR. This supplies c <= [18e/(pi alpha)] (M/T_*) log d, but by itself only A <= 6 e c M. The sharper LR constant quoted in older notes is not required.

## Verification

[Numerical, elementary formulas] Ran python3 checks.py successfully. CHECKS.md records the output. The maximum finite-difference error in the strip entropy derivative was 4.30618e-10. The finite-time fraction at q=4 is 0.889682199924. Scaling covariance, the harmonic-mean velocity comparison, the pi/6 conversion, the LR series majorant, and the spatial logarithm inequality passed. No new microscopic spectra or CFT identifications were computed.

[Proved review correction] The current norm is stated for a chosen centered representative g, as ||i[g_12,g_23]|| <= 2||g||^2. It is not incorrectly assigned to an arbitrary original h using the optimized M. Current gauge dependence remains separate from the gauge-invariant SIE argument.

## Open question

[Assumed/open] Neither a universal bound A <= C M log d nor a universal no-observable theorem was established. A concrete sufficient missing theorem is E(h) >= kappa A with dimension-independent kappa, where E(h) is ancillary single-bond entangling capacity. SIE would then give C = 18/kappa. Alternatively, uniform controlled boundary-state preparation at a bounded spatial length would give C = 216 ell_0/(pi alpha). The infrared CFT assumption supplies neither premise. The sharp pi/log 2 coefficient remains conjectural.
