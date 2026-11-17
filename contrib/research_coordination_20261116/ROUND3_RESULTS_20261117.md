# Third-round results: simple qubit phases and mobile length selection

Date: 2026-11-17. OpenAI Codex, GPT-6 family.
Status: summary of the three completed chats, not an independent proof
review or a new numerical calculation. Source contributions and their
qualification of each claim are linked below.

This records the outcomes of the
[third-round assignments](LATEST_REVIEW_20261117.md). All three chats are
idle at this checkpoint. The four-term calculation was stopped at the
user's battery-saving request with resumable checkpoints retained; this
summary does not restart it or assign further computation.

## Four-term qubit model: encouraging diagnostics, unresolved c

The [V=4 bulk study](../v4_bulk_diagnosis_20261117/README.md) kept the
fixed Pauli Hamiltonian

```text
H = sum_j [(9/16) X_j Y_(j+1) + (7/16) Y_j X_(j+1)
           + Z_j Z_(j+1)] - (7/4) sum_j Z_j.
```

It has four displayed monomials and coefficient spread four. The earlier
qualification about the invariant coupling hierarchy still applies.

**[Numerical.]** The main new result is a comparison of L=48,96,144
using connected density/bond correlations and calibrated entropy fits.
The central density modulation amplitude falls from about 0.154--0.158
at L=96 to 0.128 at L=144. Normalized connected structure factors and
long-distance correlations also fall. This is compatible with decaying
boundary-induced oscillations on the tested sizes; it does not exclude
a smaller nonzero plateau or a crossover. An independent L=96 start
reproduces the previous state, with overlap magnitude 0.999999999875.

The entropy slope c_ent is 2.8184 between 48 and 96 sites. Between
96 and 144 sites it is 2.4542 at the last completed accuracy stage and
2.5669 using four saved sweeps of the interrupted higher-accuracy stage.
Window choices at the latter stage give roughly 2.34--2.65. Exact free
controls at the same larger sizes return 0.995 and 1.954 for known c=1
and c=2, validating normalization but not supplying a universal correction.

The largest completed stage fails the bond-dimension entropy gate, and
the interrupted stage has not passed the strict convergence requirements.
Thus **2.57 is a provisional finite-size entropy estimate, not an
established thermodynamic Virasoro central charge**. Positive-velocity
mode evidence, a further converged size sequence and a robust parameter
interval remain missing. No new microscopic counterexample is established.

## Two-term model: evidence favors c=1, with ordering at stronger coupling

The [two-term study](../two_term_qubit_phase_20261117/README.md) analyzed

```text
H(g) = sum_j [2 X_j Y_(j+1) + g Z_j Z_(j+1)].
```

**[Numerical.]** At g=1/2 and g=1, calculations through L=20 show
approximately inverse-length neutral gaps and decreasing density/dimer
fluctuations, favoring c=1. An exact basis map also identifies previously
saved L=128 data at g=1, consistent with c=1 but with appreciable fit
dependence. The equal-coefficient point g=2 trends toward density order;
g=4 is strongly ordered. No c>1 evidence was found.

**[Proved/derived, as qualified in the source.]** A local conserved-charge
obstruction excludes the standard regular homogeneous difference-form
Yang--Baxter route, not every possible form of integrability. A separate
obstruction excludes faithful term-by-term free-Majorana mappings.
An exact duality plus an established gap-stability theorem gives a
gapped bulk at sufficiently large |g|, with a non-explicit threshold.
It does not certify any of the moderate numerical couplings.

Near g=0 the slow Majorana hopping appears already at first order,
with v_B^(1)=8|g|/pi. This line is not the previously tuned cancellation
curve. That observation is not an all-order criticality theorem.
A global c<=1 bound for this line and its precise phase boundary remain
open. Its few terms and moderate coefficients do not currently improve
the evidence for a counterexample.

## Fixed d=22: sharper relative estimates, but selection remains open

The [mobile length study](../mobile_length_selection_20261117/README.md)
tested the weak competitor identified in the assignment: n+2 cells of
n bits versus n+1 cells of n+1 bits, with the same total material length.
It works in the exact legal vacuum sector at lambda=f=0, retaining the
terminal, marker contact halves, quantum clocks and vacancy dynamics.

**[Proved, awaiting independent review.]** The first-order difference is
J alpha_(N,V)(n+1)E_(n+1). New nonperturbative bounds reduce the comparison
to scalar vacancy Hamiltonians while preserving the exponentially small
clock-tail scale. For every n and finite common vacancy sector, the
clock contribution at second order favors the intended cell length.
In this restricted real sector, physical cyclic and untwisted rooted
ground energies agree exactly, removing the earlier boundary allowance.

The remaining scalar response to the changed marker pattern does not
cancel and opposes selection in the tested baths. The favorable clock
term is larger in those tests, but a uniform all-order comparison is
still missing. The positive-tilt theorem and transfer to target-minimizing
canonical blocks do not imply validity at arbitrary prescribed density.

**[Numerical.]** At n=1, q=1/4, J/t=1/80 and N=V=6, the observed
positive selection margin is about 1.58e-6 t. The conservative scalar
comparison yields about 4.08e-8 t, also positive. Limited n=2 canonical
tests likewise retain the intended ordering. The existing selector
chemical potential needs no retuning in these tests. These finite
calculations are not a thermodynamic certificate.

No new terms or stronger coupling hierarchy were introduced. The full
model's existing small coefficients and many processes remain. Arbitrary
wrong-length words, charged/port/readout effects, illegal sectors and
the final CFT still require control. Unbounded c at fixed d is unproved.

## Publication and overall status

The underlying results were already committed and pushed through
4a9a12f6 (qubit bulk study), a700ebd1 (two-term study), and 608b00e3
(mobile length study). This coordination update only records their
outcomes. Its local documentation links and patch were checked; the
source numerical campaigns were not rerun.

The four-term qubit model remains the more promising simple c>1
candidate, while the two-term model currently favors ordinary c=1 or
order. The mobile construction has stronger control of one weak
selection mechanism but remains conditional. No universal B(d) or
verified fixed-d unbounded-c construction was obtained this round.
