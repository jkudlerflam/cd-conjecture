# Independent review of the three-state vacancy control

**[checked code and analytical normalization; numerical verification]**
The microscopic control and its reported interpretation are consistent.
It tests one charge boson plus one squeezed XX payload boson at d=3,
which saturates c=d-1. Its finite-size fits do not establish c>2 and must
not be presented as a counterexample. The saved summary correctly makes
that distinction.

## Hamiltonian, Hilbert space and entropy

**[checked]** The ternary code uses 0 for vacancy, 1 for down and 2 for up.
The basis has exactly L!/[N_up! N_down! (L-N_up-N_down)!] states. The
vacancy mask selects exactly one empty endpoint; the exchange mask
selects exactly the (up,down) pairs. Each legal directed transition
appears once per physical bond. The Hamiltonian is therefore

    h = -(|0,up><up,0|+|0,down><down,0|+h.c.)
        -J (|up,down><down,up|+h.c.).

The minus sign for XX exchange is deliberate and consistent with every
energy formula. On the even squeezed rings considered, it is the usual
gapless XX control. All terms preserve both species numbers. The
matrices are real, symmetric and irreducible in the tested sectors, with
nonpositive off-diagonal entries for J>0, so the finite-sector ground
state is unique and positive. The sign convention avoids an accidental
ground-momentum ambiguity in the numerical control.

**[checked]** The entropy routine correctly blocks the Schmidt matrix by
subsystem up/down counts, with complementary subsystem numbers fixed by
the global sector. Every tensor coefficient is represented once. Its
unit-trace check is appropriate, and discarding singular values below
1e-15 is harmless at these dimensions. The real-valued Schmidt matrix
is correct for the untwisted real ground states used here; it would need
a complex dtype if reused for twisted ground states. The flux routine
does not presently call it.

The separate `calibrate.py` comparison with the full free-fermion
correlation matrix is a meaningful entropy calibration. It tests the
physical tensor factorization and boundary convention independently of
the interacting central-charge fit.

## Weak-J energy and velocity normalization

**[proved]** At half occupancy, with N=L/2 and the tested zero-momentum
squeezed XX ground state,

    E_charge(L)=-2 csc(pi/L),
    alpha_L=1/2-2/[L^2 sin^2(pi/L)],
    E_XX(N)=-2 csc(pi/N).

The first derivative of the energy is alpha_L E_XX(N). With
alpha=1/2-2/pi^2, its expansion is

    E(L) = -(2/pi+J alpha/pi)L
           +[-pi/3+J(2/(3pi)-2pi alpha/3)]/L
           +O(L^-3)+O(J^2) at fixed L.

In particular, the 1/L term includes the finite-size correction to
alpha_L. It would be incorrect to retain the bare charge velocity 2
while comparing the O(J) Casimir coefficient with the payload velocity.

The bulk energy density at general occupied density rho is, to first
order in J at zero magnetization,

    e(rho)=-(2/pi) sin(pi rho)
           -(2J/pi)[rho^2-sin^2(pi rho)/pi^2]+O(J^2).

At rho=1/2 this gives charge v/K=2-8J/pi^2+O(J^2). A uniform charge
twist does not change the first-order contact expectation on the fixed
adiabatic branch, so charge vK=2+O(J^2). Consequently

    v_charge=2-4J/pi^2+O(J^2),
    v_spin=4J alpha+O(J^2),
    v_charge+v_spin=2+J(2-12/pi^2)+O(J^2).

The physical spin velocity includes the coordinate conversion from
squeezed length N to physical length L. The worker's final analysis
uses this corrected sum.

**[checked]** For a U(1) flux phi, the finite-size ground-energy curvature
gives vK=2pi L [E(phi)-E(0)]/phi^2. The charge generator is occupied
number; the spin generator is (N_up-N_down)/2. The response code uses
those exact transported charges. Its number-sector formulas are

    (v/K)_spin = 2L [E(Sz=1)-E(Sz=0)]/pi,
    (v/K)_charge = L [E(N+2)+E(N-2)-2E(N)]/(4pi).

These factors are correct: the spin change is one and the balanced
charge change is two. Spin interchange at zero magnetization excludes
mixed charge-spin quadratic responses. The two flux values and the
quadratic extrapolation test the finite-flux error. Finite ring and
finite number-step errors remain and are correctly visible in the data.

## What the two diagnostics support

**[numerical; L=8,12,16]** Half-chain size scaling gives apparent c from
2.0258 to 2.0385 as J ranges from .2 to .05. Restricting the slope to
L=12,16 gives 2.0113 to 2.0191. Single-ring chord fits have appreciably
larger upward bias, even when their fit residual is small. This is
consistent with finite-size drift toward two, and does not justify
treating the excess above two as a violation.

**[numerical consistency diagnostic]** Independent stiffness and
compressibility calculations resolve both a fast charge velocity and
a slower positive spin velocity. The energy Casimir coefficient agrees
with their sum at the reported sizes. This is useful independent
support for two modes. It is not a precision determination of their
central charges: the fit with bulk, 1/L and 1/L^3 terms uses exactly
three sizes and thus has no residual degrees of freedom. The supplied
two-size alternative and analytic weak-J coefficient are essential
cross-checks; the small agreement error alone is not a statistical
uncertainty estimate.

The calculation is a full microscopic three-state chain calculation in
fixed symmetry sectors. It validates an accessible version of the
vacancy mechanism, not the 24-state bond or its 23-mode payload.

## Independent numerical spot-check

`audit_numerics.py` checks the basis dimension and both exact directed
transition counts at L=8, J=.1, then embeds its ground vector in all 3^8
tensor-product states. It compares complete Schmidt decompositions for
every cut against the worker's symmetry-block entropy calculation and
the saved values. The maximum entropy discrepancy across all seven cuts
is 4.45e-16; the recomputed ground energy differs from the saved value
by 2.67e-15. `numerics_audit.json` records the numerical results.
This optional spot-check needs numpy and scipy, unlike this folder's
standard-library-only C39 certificate checks.

    OPENBLAS_NUM_THREADS=1 python3 contrib/vacancy_continuation/all_n_audit/audit_numerics.py

Review by the independent GPT-6-family all-n audit worker, 2026-10-28.
No human review is claimed. The numerical worker independently derived
the same first-order charge-velocity correction before this review was
completed. No numerical-worker source or data were edited.
