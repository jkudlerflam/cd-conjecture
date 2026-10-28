# Independent review of the density-window continuation

**Verdict [proved finite statements; conditional continuum application]:**
The finite-volume estimate and reference momentum bound in
[`DENSITY_WINDOW.md`](../DENSITY_WINDOW.md) are correct for L>=2. They
remove the need to tune the occupied-particle chemical potential to exact
half filling, subject to the same stated microscopic matching and
phase-stability inputs. I found no counterexample. This review does not
promote the density estimate to a proof of criticality.

## Exact hopping energy and twists

At fixed N, an internal cyclic-translation eigenstate supplies a twist
of the hard-core hopping problem. Its Jordan-Wigner fermions have the
additional parity boundary phase, zero for odd N and pi for even N.
The N lowest one-particle energies form a consecutive momentum sea,
allowing a choice at tied endpoints. The geometric sum of their
exponentials has magnitude sin(pi N/L)/sin(pi/L). Its real part cannot
exceed that magnitude. Thus every internal twist has energy at least

    -2t sin(pi N/L)/sin(pi/L).

A completely polarized internal word has cyclic eigenvalue one and
realizes the symmetric sea with that energy, for both parities. This
argument also gives zero for N=0,L and works at odd L. No assumption on
the internal interaction or q enters the purely kinetic minimum.

## Constant in the finite-volume density bound

The operator inequality |h_lift|<=M n_x n_(x+1) gives the claimed
sector energy brackets. With N0=floor(L/2), comparison to N0 yields

    A_L[1-cos(pi delta)]
      <= (3/2)JM + A_L[1-sin(pi N0/L)].

Use A_L>=2t/pi, 1-cos(pi delta)>=4delta^2, and
1-sin(pi N0/L)<=pi^2/(8L^2). Dividing by 4A_L gives precisely

    delta^2 <= 3pi JM/(16t) + pi^2/(32L^2).

This verifies the odd-ring correction as well as the interaction
coefficient. At JM/t<=1/1000 the limiting deviation is at most
0.024270324..., strictly smaller than 1/40. The statement concerns each
ground-state number sector, so degeneracy or superposition of sectors
causes no problem. Equation-of-state smoothness does not follow.

## Momentum margin and finite-coupling scope

For a reference vertex of dimension below three, its unwrapped physical
momentum has magnitude less than 2pi rho sqrt(3). Its distance to a
nonzero reciprocal vector is therefore larger than
2pi(1-21sqrt(3)/40). Its distance to zero is at least
19 delta_int/(40b) in every nonforward payload channel. Pure marker
harmonics obey the stronger 19(2pi)/(40b) bound. Gaussian zero-harmonic
currents are retained; derivative descendants cannot reduce dimension.
This establishes the stated minimum of the two bounds.

A nonzero reciprocal-vector resonance contributes m^2/rho^2 to the
reference dimension, exceeding 1600/441>3 in the window. Thus the
half-filled reference is not essential to the locking exclusion.

A minor wording qualification was sent to the author: K_gap=rho^-2
is the *free reference* value at that density. The finite-J stiffness
is dressed and must be controlled by the stated relative Gaussian
matrix budget. No exact interacting equality is needed. Also state
L>=2 explicitly, since a one-site periodic two-site bond needs a
separate convention.

## Independent diagnostics and provenance

`density_review_checks.py` independently sorts the twisted one-particle
spectra in 10,710 cases for L=2,...,64 and every N. The untwisted
hard-core sector attains the formula to 3.3e-14 absolute error. A
separate sample of nonuniform Fermi points tests 89,532 harmonic instances
at b=1,2,5,24 and densities including the closed-window endpoints. All
bounds pass. These are diagnostics accompanying the analytic review,
not proofs of the continuum dictionary or a measured central charge.

Run `python3 contrib/vacancy_continuation/local_audit/density_review_checks.py`.
Only the Python standard library is required. Independent GPT-6-family
Codex worker, 2026-10-28; no human review is claimed.
