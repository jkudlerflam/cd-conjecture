# Mobile length selection: tail control and the remaining marker response

Date: 2026-11-17. Author and AI assistance: OpenAI Codex, GPT-6 family.
No independent human or agent review is claimed.

## Claim and status

**[Proved, awaiting independent review.]** In the exact legal vacuum
sector of the literal d=22 nearest-neighbor model, with lambda=f=0,
the weak marker-changing comparison has the stated first-order margin.
Two new results keep the terminal tail rather than using a bulk norm:

1. A nonperturbative bound reduces each length word to explicit scalar
   vacancy Hamiltonians. The added potentials scale as q^(2k+2), with
   the actual marker contacts retained. It uses the existing allowance
   q<=1/4, J/t<=1/80, t/2<=s<=3t/4, without a D or length factor.
2. For **every** n>=1 and every finite common vacancy sector, the clock
   part of the second-order correction favors n+2 cells of n bits over
   n+1 cells of n+1 bits. No charge-gap bound is used for this sign.

**[Unresolved.]** These results do not prove mobile length selection
uniformly in volume. The scalar response to changing the marker pattern
does not cancel. It opposes selection in the checked free baths, although
the favorable clock term is larger. An all-order, uniform comparison
of the explicitly reduced scalar problems remains necessary.

Fixed mu_sel=-(n+2)E_n survives all bounded checks. No retuning is found
necessary in them. Its all-n thermodynamic validity is neither proved
nor refuted. No smaller J/t or added penalty is prescribed to cover
an uncontrolled remainder.

## Exact scope

All forbidden-arrow flags and charged-cell flags are fixed to zero.
The ordered marker pattern is fixed separately in each sector. Payload
bits may be fixed to zero because every payload string is isospectral
at f=0. Clocks and legal vacancies remain fully quantum and entangled;
marker positions move in physical space. Both terminal and marker
onsite energies retain their two contact halves.

The nonperturbative theorem uses the positive length tilt and its
minimizing canonical blocks. A positive tilted difference transfers
to a **target-minimizing** block; transfer to arbitrary densities
is not assumed. Some larger-n checks below are canonical diagnostics
outside that ensemble and are labeled accordingly.

In this restricted real vacuum sector, positivity improves the previous
global cyclic-boundary allowance to **zero**: the physical cyclic and
untwisted rooted ground energies agree exactly. The usual single
96t/N comparison allowance is still the available general statement
once complex charged/port processes are restored.

## What the weak comparison tests

Both words occupy N=(n+1)(n+2) material sites, at the same vacancy
number V. Their frozen difference is

    (n+1)(mu_sel+E_(n+1))-(n+2)(mu_sel+E_n)
       = (n+1)E_(n+1).

Its first-order mobile value is J alpha_(N,V)(n+1)E_(n+1).
This margin is much smaller than the balanced (n-1,n+1) defect cost
and includes the otherwise missed change in marker count.

The common FF clock u_k proportional to (1,q,...,q^k) makes an exact
reference with the entire marker-dependent vacancy bath retained.
The previous clock-pin channel supplies additive coercivity in that
reference. A positive two-dimensional matrix estimate then gives the
scalar sandwich in [PROOF.md §3](PROOF.md). Compressing the target
onto its **frozen** clock eigenvectors gives a sharper upper bound.
Comparing these two scalar bounds is a concrete sufficient test.

The second-order calculation separates scalar marker response from
clock excitation. The clock form factors obey

    F_(k+1)(theta) <= (512/81)q^2 F_k(theta)
                         <= (32/81)F_k(theta).

This is a pointwise analytic inequality, valid for all lengths and
momenta. The response matrix is positive, so the longer clocks obtain
less energy lowering and the intended word gains a favorable margin.
The proof does not assume the scalar response has the same sign.

## Bounded numerical checks

**[Numerical corroboration, not interval certificates.]** At q=1/4,
J/t=1/80, the exact finite matrices give:

| Target n | N,V | First-order margin | Observed margin | Ensemble status |
| --- | --- | ---: | ---: | --- |
| 1 | 6,6 | 1.57208384e-6 t | 1.57571292e-6 t | Both minimize s/t=0.6 |
| 2 | 12,2 | 4.19735303e-7 t | 4.20052172e-7 t | Canonical diagnostic only |
| 2 | 12,4 | 3.42813705e-7 t | 3.43301808e-7 t | Canonical diagnostic only |

At n=1,N=V=6, the exact second-order coefficients are

    marker/scalar: -3.25703e-6/t,
    clock:        +2.65078e-5/t.

The nonperturbative scalar lower bound for the competitor minus the
frozen-clock scalar upper bound for the target is about **4.08e-8 t**.
It is positive but much more conservative than the observed margin.
The scalar minimizers in that check are also V=6. This bounded check
does not establish positivity for repeated words or infinite volume.

The n=1 crossing in mu moves from -0.01001928594 at first order to
about -0.01001828537 at J/t=1/80. The chosen value -0.01044901688
therefore remains on the intended side. The finite-volume crossing
formula and its limitations are stated in PROOF §7.

The scans keep positive weak-pair margins through J/t=1 for n=1 and
through 0.2 for n=2,V=2. These larger values are diagnostics, outside
the common theorem allowance. No thermodynamic threshold is inferred.

## Reproducing and checking

From the repository root, with NumPy, SciPy and mpmath installed:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B contrib/mobile_length_selection_20261117/checks.py
```

The script forces one BLAS thread and writes only [checks.json](checks.json).
The recorded run takes about five seconds. The largest full clock/gap
matrix has dimension 110,565; the largest saved eigenpair residual is
below 7.3e-12. Checks include literal source matrices for both words at
n=1,2, the analytic form-factor bounds, separately resolved second-order
coefficients, finite margins and crossings, scalar comparison bounds,
and direct physical 22-state versus rooted ground energies.

An initial comparison used a 1e-15 tolerance between an existing
double-precision E_2 and the new high-precision value; their propagated
arithmetic difference was 1.61e-15. The cross-implementation tolerance
was changed to 1e-14. No formula, sign, or reported margin was changed.
The previous contributions and parent bond are imported unchanged.

## Simplicity and relation to the project

[SIMPLICITY.md](SIMPLICITY.md) distinguishes independent controls,
literal processes and microscopic coefficient ratios. This contribution
adds no Hamiltonian terms. It removes no active term from the proposed
full high-c mechanism either. The tiny physical selector coefficient
and the full architecture's other hierarchies remain.

The result advances the relative-response question left by the earlier
[clock control](../mobile_vacuum_clock_20261117/README.md) and
[charged-rail pin](../small_d_mobile_matching_20261116/README.md).
Arbitrary wrong-length words, charged sectors, nonzero ports/readout,
and illegal sectors still require separate comparisons. No unitary
compact CFT or total nonchiral Virasoro central charge is established;
the proposed c=n+1 remains conditional.
