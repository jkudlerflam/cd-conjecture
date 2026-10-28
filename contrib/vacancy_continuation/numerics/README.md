# Full microscopic vacancy/XX control with two criticality diagnostics

## Claim and status

[numerical] A literal three-state, one-site translation-invariant,
nearest-neighbor Hermitian chain supports the expected charge and payload
modes at half occupancy. Full microscopic exact diagonalization at fixed
J/t=0.2, 0.1 and 0.05, on rings L=8,12,16, gives half-chain entropy slopes
consistent with total Virasoro c=2. Independent flux and number-sector
responses calibrate the two velocities, whose sum agrees with the ground
energy Casimir coefficient. The largest ground-state sector has 900,900
states; the largest response sector has 2,018,016 states.

This is a control of the vacancy mechanism, not a counterexample:
d=3 and c=2 saturates d-1. It is not a simulation of the 24-state bond,
its tiny scales, its cell selection or its high-channel stability.
Finite-size deviations above 2 below are not evidence for a violation.
No thermodynamic theorem is claimed.

## Literal model

[explicit finite definition] Each site has basis {0, down, up}, where 0
is vacant. With t=1 and J>0, the only nonzero two-site bond entries are

    h = -sum_(a=up,down) (|0,a><a,0| + |a,0><0,a|)
        -J (|up,down><down,up| + |down,up><up,down|).

The payload is the ferromagnetic XX chain in occupied-site coordinates.
This sign makes the full microscopic bond stoquastic and selects its
zero-momentum finite-ring ground branch. It avoids a separate finite-ring
momentum frustration that arises for the antiferromagnetic sign at some
occupied numbers. Both signs have the intended critical XX payload, but
no claim that their finite-ring gluing is identical is used here.

Total occupied number N and magnetization M=(N_up-N_down)/2 are exact
U(1) charges. The principal runs use N=L/2 and M=0. Matrix construction
uses the complete physical basis at these charges. No charge-sea,
squeezed-space or logical-code projection is made during diagonalization.
The partial trace is also taken over contiguous *physical* sites.

[checked finite matrix] `xx_vacancy.py` independently rebuilds the L=8
Hamiltonian from the literal 9-by-9 two-site matrix and finds zero entry
error. Charge sectors are built exactly, sparse Hermiticity is checked,
and saved eigenvector residuals are below 1e-7. The ground-state solver
starts with a deterministic random vector, with no momentum restriction.

## Entropy

[read, CFT diagnostic] We use the periodic single-interval formula

    S_L(ell) = (c/3) log[(L/pi) sin(pi ell/L)] + constant + corrections.

See [Calabrese and Cardy, equation (1)](https://arxiv.org/html/hep-th/0405152v3).
The coefficient adds for independent compact bosons with unequal speeds.

[numerical] At ell=L/2, fitting S versus log L removes the small-ell
shape corrections. The first fit uses all three rings; the second uses
only L=12 and 16. These are finite-size effective slopes, not statistical
confidence intervals.

| J | c from L=8,12,16 | c from L=12,16 | Center-shape c at L=8 | L=12 | L=16 |
|---|---:|---:|---:|---:|---:|
| 0.20 | 2.025754 | 2.011260 | 2.505935 | 2.241363 | 2.106726 |
| 0.10 | 2.033771 | 2.016165 | 2.527493 | 2.253410 | 2.112954 |
| 0.05 | 2.038487 | 2.019092 | 2.539793 | 2.260632 | 2.116735 |

The center-shape estimator is
`3[S(L/2)-S(L/2-2)]/[-log cos(2pi/L)]`. Its sizable drift illustrates why
one ring or a short-interval slope would be misleading. All individual
entropies and fits with different short-distance exclusions are retained.

[checked independent calibration] `calibrate.py` uses the same physical
partial-trace code on the fully polarized vacancy chain, whose charge
sea is free. For L=8,12,16 every entropy agrees with the independently
constructed fermion correlation-matrix answer to less than 9e-14.
Energies agree with -2/sin(pi/L) to less than 8e-15. The finite-size
half-chain estimate is c=0.979755, approaching the known one-mode value;
it is not artificially normalized to exactly 1.

## Independent response and Casimir check

[derived Gaussian diagnostic] Normalize each conserved charge Q so a
Gaussian mode has zero-mode and flux energies

    E(Q+delta Q)-E(Q) = pi v (delta Q)^2/(2 K L),
    E(Phi)-E(0) = v K Phi^2/(2 pi L).

At zero magnetization, spin reversal prohibits the mixed quadratic
charge-spin term. The two diagonal responses therefore determine
separate v and K. For charge we use N+/-2, keeping M=0; for spin we use
M+1, keeping N fixed, and spin reversal supplies M-1. Thus

    (v/K)_charge = L[E(N+2)+E(N-2)-2E(N)]/(4pi),
    (v/K)_spin   = 2L[E(M+1)-E(M)]/pi.

The twist is applied to *every microscopic transport term*. A swap of
local states a,b has phase exp[i Phi (q_a-q_b)/L]. Charge uses
q=(0,1,1), and spin uses q=(0,-1/2,+1/2). In particular, the spin twist
acts on vacancy hopping as well as on the occupied XX exchange.
We compute `vK=2pi L[E(Phi)-E(0)]/Phi^2` at Phi=.04,.02 and cancel the
leading Phi^2 error by Richardson extrapolation. The largest change on
halving Phi at L=16 is below 8e-7.

[numerical] The following velocities come from full microscopic responses
at L=16. Separately, ground energies at L=8,12,16 are fitted to
`E(L)=e_inf L - pi V/(6L) + B/L^3`. Two c=1 modes predict V=v_c+v_s.
The three-parameter fit interpolates three data points, so its agreement
with independently measured responses is a consistency check, not an
estimated precision or proof. Omitting B and using the two larger rings
gives V=2.17928,2.10399,2.06506 respectively, showing the finite-size
systematic effect.

| J | Response v_charge | Response v_spin | Response sum | Casimir V | V/(v_charge+v_spin) |
|---|---:|---:|---:|---:|---:|
| .20 | 1.918739 | .234759 | 2.153498 | 2.149341 | .998070 |
| .10 | 1.959688 | .118790 | 2.078478 | 2.076089 | .998851 |
| .05 | 1.979940 | .059748 | 2.039688 | 2.038236 | .999288 |

The spin speed decreases linearly with J while the charge speed remains
near 2. Both compressibilities and stiffnesses are positive at all tested
points. The energy coefficient is velocity weighted; it is not converted
into c by setting both speeds to the same value.

## Analytic weak-coupling calibration

[proved first-order finite-ring input] In the zero-twist charge branch,

    E_charge(L) = -2/sin(pi/L),
    alpha_L = 1/2 - 2/[L^2 sin^2(pi/L)],
    E_XX(N) = -2/sin(pi/N),          N=L/2,
    E_leading = E_charge + J alpha_L E_XX.

The numerical remainder `(E0-E_leading)/(J^2 L)` ranges from -0.0174
to -0.0189 across the recorded points. It tends smoothly with decreasing
J, supporting ordinary second-order charge relaxation in this control.
No infinite-volume estimate is inferred from that small-size test.

[derived first-order continuum prediction] Write
`alpha=1/2-2/pi^2`. The physical payload velocity is
`v_spin=4J alpha+O(J^2)`. The bulk energy at density rho is

    e(rho) = -(2/pi) sin(pi rho)
             -(2J/pi)[rho^2-sin^2(pi rho)/pi^2]+O(J^2).

Its compressibility gives `(v/K)_charge=2-8J/pi^2+O(J^2)` at half
occupancy. The charge twist gives `(vK)_charge=2+O(J^2)`, hence

    v_charge = 2-4J/pi^2+O(J^2),
    v_charge+v_spin = 2+J(2-12/pi^2)+O(J^2).

The charge velocity correction is necessary. Equivalently, expanding
`alpha_L=alpha-2/(3L^2)+...` supplies it in the finite-ring Casimir
coefficient. Using unrenormalized charge speed 2 would miscalibrate the
O(J) energy fit. The first-order sums 2.156829,2.078415,2.039207 agree
with the response trend.

## Reproduction and provenance

Run from the repository root with Python 3, NumPy and SciPy:

```bash
OPENBLAS_NUM_THREADS=1 python3 contrib/vacancy_continuation/numerics/xx_vacancy.py --L 8
OPENBLAS_NUM_THREADS=1 python3 contrib/vacancy_continuation/numerics/xx_vacancy.py --L 12
OPENBLAS_NUM_THREADS=1 python3 contrib/vacancy_continuation/numerics/xx_vacancy.py --L 16
OPENBLAS_NUM_THREADS=1 python3 contrib/vacancy_continuation/numerics/responses.py --L 12
OPENBLAS_NUM_THREADS=1 python3 contrib/vacancy_continuation/numerics/responses.py --L 16
OPENBLAS_NUM_THREADS=1 python3 contrib/vacancy_continuation/numerics/calibrate.py
OPENBLAS_NUM_THREADS=1 python3 contrib/vacancy_continuation/numerics/analyze.py
```

The default J list is .2,.1,.05. All outputs are written beside the
scripts. The recorded L=16 entropy run took 6 seconds and its response
run 68 seconds on the current execution host; all calculations ran on the current execution host. Runtime varies by machine.
An optional figure can be regenerated with
`python3 contrib/vacancy_continuation/numerics/plot.py` after installing
Matplotlib. [View the figure](vacancy_control.png).

No large state vectors are saved. `summary.json` contains the combined
numbers; per-size files retain all energies, entropies and residuals.

Author: GPT-6-family Codex worker, 2026-10-28. Independent algebra/code
review by the sibling all-n audit worker checked the basis, physical
partial trace, first-order formula and charge-velocity correction.
A separate full 3^8 tensor embedding and direct Schmidt decomposition at
L=8, J=.1 reproduced all seven bipartitions within 4.45e-16; see
[the independent audit](../all_n_audit/audit_numerics.md).
The primary CFT source was opened. No human review is claimed.
