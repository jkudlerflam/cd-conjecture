# C100 result. Thermal and transport routes

2026-10-29. The CFT energy-current variance and thermal Drude coefficient
do measure A = sum_a c_a v_a. None of the comparisons established here
removes the crossover scale. The proposed bound A <= pi M log_2 d
remains unproved.

The main new lattice result is a positive energy-density spectral-moment
bound. Its comparison with the CFT gives an explicit conditional bound
without requiring an exact lattice Drude peak.

## Main theorem

[Proved, no CFT assumption.] Choose a centered bond representative g
with ||g|| = M, set H = sum_i g_i and
h_q = sum_i exp(-iqi) g_i, and take a Gibbs state or a ground state.
Then, for nonzero momenta |q| <= pi,

    0 <= m1(q,T) := <[h_-q,[H,h_q]]>/(2L)
                 <= 8 M^3 q^2.                       (1)

The moment is positive because it is a positive-frequency absorption
integral with weight omega(1-exp(-omega/T)). The upper bound follows
from stationarity, the nearest-neighbor commutator support, and
||[g_j,[g_i,g_(i+1)]]|| <= 4 M^3. PROOF.md gives the full derivation.

[Proved within the assumed CFT.]

    m1_CFT(q,T) = pi A T^2 q^2/6
                 + q^4 sum_a c_a v_a^3/(24pi).       (2)

[Assumed matching; proved conclusion.] If a positive low-energy part
of the lattice moment at (q_*,T_*) is at least alpha times (2), with
0 < alpha <= 1, then

    A + q_*^2 sum_a c_a v_a^3/(4pi^2 T_*^2)
        <= 48 M^3/(pi alpha T_*^2).                  (3)

The simpler bound A <= 48 M^3/(pi alpha T_*^2) follows. The matching
condition is on integrated spectral weight and allows peak broadening.
It is stronger than specifying the zero-temperature CFT. Taking the
infrared limit T_* -> 0 makes this bound weaker, not stronger.

## Route comparison

Put B = sum_a c_a/v_a. All CFT entries below refer to the continuum
stress tensor. Microscopic matching at a specified temperature or size
is an additional hypothesis.

| Quantity or route | CFT value | Lattice control | Outcome |
| :--- | :--- | :--- | :--- |
| Total energy-current variance per length | pi A T^3/3 | Local norm ||j_i|| <= 2M^2 only gives Var(J)/L <= 4LM^4 | Spatial summation requires correlation control; no uniform bound follows |
| Local current midpoint correlator | pi^2 c T^4/4 | At most 4M^4 | Measures c, with a matching temperature |
| Local zero-T current spectrum | c omega^3/(24pi^2) | Total connected weight at most 4M^4 | Measures c, with a matching frequency cutoff |
| Thermal Drude coefficient, Re kappa = pi D_th delta + regular | D_th = pi A T/3 | Transport sum rule gives D_th <= 16M^3/T under the stated transport limits | Needs both matching and a nonzero lattice Drude component |
| Positive density first moment | Equation (2) | Equation (1) | Gives (3), without assuming a Drude delta function |
| Static density susceptibility | pi B T^2/3 + A q^2/(12pi) | Its inverse-frequency moment is not bounded by local norm or positive moments | A new susceptibility estimate would be needed |
| Thermal free energy and entropy | f=e0-pi B T^2/6; s=pi B T/3 | s <= log d; high-T variance per site <= 3M^2 | Bounds B at a specified matching temperature |
| Modular relation and ring energy | e0 L-E0(L) = pi A/(6L)+o(1/L) | Absolute ring correction at most 2M | Needs a specified matching size L_* |

The CFT Drude convention agrees with [Gawedzki, Langmann, and
Moosavi](https://link.springer.com/article/10.1007/s10955-018-2025-x).
The conductivity sum-rule framework is also discussed by
[Shastry](https://arxiv.org/abs/cond-mat/0508711); the nearest-neighbor
bound used here is derived directly in PROOF.md.

## The bound containing log d

[Assumed matching; proved implication, using C69's entropy result.] If

    s(T_s) >= alpha_s pi B T_s/3,

then

    A <= 3 v_max^2 log d/(pi alpha_s T_s).            (4)

This follows from s <= log d and A <= v_max^2 B. Under the local
real-time CFT identification in C69, its conservative propagation bound
v_max <= 6 e M gives

    A <= 108 e^2 M^2 log d/(pi alpha_s T_s).          (5)

When both hypotheses hold, take the minimum of (3) and (4), or (5)
if only the microscopic speed estimate is available. Neither estimate
uniformly dominates the other. No optimality of the numerical constants
is claimed.

[Assumed finite-size matching; proved implication.] If at L_*,

    e0 L_*-E0(L_*) >= alpha_E pi A/(6L_*),

the ring-energy bound gives

    A <= 12 M L_*/(pi alpha_E).                      (6)

## Why the onset scale survives

[Proved limitation of the listed arguments.] An infrared asymptotic
formula supplies some low-temperature or large-size matching regime.
It does not control its endpoint in lattice units. Local current norms
bound amplitudes but not spatially integrated susceptibilities or
inverse-frequency moments.

[Proved thermodynamic insufficiency.] The analytic family
f_tau(T) = -(log d) sqrt(T^2+tau^2) obeys entropy capacity, positive
specific heat, pressure convexity, and a bounded energy interval, but
its low-T quadratic coefficient grows as 1/tau. PROOF.md also gives
a two-scale version with fixed leading high-T variance. These are
counterexamples to an inference from those thermodynamic constraints,
not constructed spin chains and not counterexamples to the conjecture.

[Proved CFT distinction.] Modular invariance relates the continuum
Casimir and thermal coefficients. With multiple velocities it acts
sector by sector on different moduli. It does not identify the
microscopic high-temperature chain with a large continuum ring.

[Open.] A bound on the q^2 coefficient of the static energy-density
susceptibility could escape these particular failures. That coefficient
contains A/(12pi) in the CFT. Its control by M log d, without a
frequency cutoff or a gap, remains unproved. This investigation is
not a theorem excluding all thermal or dynamical proofs.

## Verification and files

[Numerical checks.] Small six-site random-chain calculations verify
continuity, the double-commutator spectral identity, positivity, and
the bounds. CFT midpoint integrals and Fourier transforms agree by
quadrature. All final checks pass, with identity residuals below
2e-15. These are normalization checks, not infrared evidence.

- PROOF.md contains the derivations, gauge qualification, transport
  limit qualification, asymmetric-speed formulas, and analytic checks.
- checks.py reproduces CHECKS.md with NumPy and SciPy.
- PROGRESS.md and STATE.md record the handoff.

All created files are in work_C100. No files outside this folder were
modified, including the project-level STATE.md. No cluster job, commit,
push, or large numerical calculation was performed.
