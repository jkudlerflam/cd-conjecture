# C100 progress

2026-10-29. Completed the requested thermal and transport investigation.
Read the supplied bound background and C100 specification, project memory,
and C69's established entropy and propagation-speed conclusions. All
writes remained in work_C100.

[Proved CFT calculations.] Derived thermal free energy and entropy,
current variance, current midpoint correlators, local current spectra,
the thermal Drude coefficient with explicit normalization, the
energy-density structure factor, its positive first moment, and its
inverse-frequency static susceptibility. A = sum c_a v_a occurs in
Var(J)/L = pi A T^3/3, D_th = pi A T/3, and the q^2 coefficient
A/(12pi) of the static density susceptibility.

[Proved lattice result.] For a gauge-optimized bond of norm M,
0 <= <[h_-q,[H,h_q]]>/(2L) <= 8 M^3 q^2 at all T >= 0.
The derivation uses stationarity to produce two powers of q. Positive
low-energy spectral matching at T_* with lower factor alpha gives
A <= 48 M^3/(pi alpha T_*^2), with an additional positive q_*^2
sum c_a v_a^3 term retained in the stronger version.

[Proved implications of other matching assumptions.] Entropy capacity
gives A <= 3 v_max^2 log d/(pi alpha_s T_s). C69's conservative
speed bound yields 108 e^2 M^2 log d/(pi alpha_s T_s), under its
real-time identification assumption. A direct ring bound
|E0(L)-L e0| <= 2M yields A <= 12 M L_*/(pi alpha_E) when a
finite-size Casimir lower bound is available at L_*.

[Scope.] No bound without a crossover scale was proved. The static
susceptibility is a remaining possible route, but local norm and positive
spectral moments do not control its inverse-frequency moment. A family
of analytic thermodynamic functions demonstrates the insufficiency of
entropy capacity, convexity, and leading high-temperature data alone.
It is not a family of constructed chains. Modular invariance is applied
only within the CFT, not promoted to a microscopic duality.

[Verification.] checks.py passed all final assertions on four seeded
six-site random chains and continuum quadratures. Maximum spectral
identity residual 1.776e-15. An initial Fourier expectation omitted the
v^2 converting density to current; this was corrected and the full
script rerun successfully. No infrared inference is made from these
small calculations. Literature checks confirmed the Drude normalization
and thermal sum-rule context.

Deliverables are RESULT.md, PROOF.md, checks.py, CHECKS.md, and STATE.md.
No changes outside the assigned folder, no cluster jobs, and no commits
or pushes. The parent project STATE.md was left untouched under the
folder restriction; local STATE.md carries the handoff.
