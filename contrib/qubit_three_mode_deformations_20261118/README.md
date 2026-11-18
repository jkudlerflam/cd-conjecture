# Full three-mode deformation map and a crystal-preemption certificate

Date: 2026-11-18. Focused round, assignment 2.

## Claim

For analytic microscopic Hermitian nearest-neighbour qubit deformations
near the m=3 root-of-unity parent, preserving its complete order-zero,
order-two and order-four kernels and leaving the last kernel scalar
through order five, the complete first retained operator is

```
K6/(128/3) = P[sum |0,-1,1,a><0,-1,1,a|_(j,j+3)
              +(3/16)G8+zeta T1+zeta* T1†
              +eta B+eta* B†-mu N]P.                       (1)
```

P forbids separations 1,2,4; G8 has the conditioning specified in
[NORMAL_SPACE.md](NORMAL_SPACE.md). The seven real controls are complex
a,zeta,eta and real mu. An exact physical normal chart shows that five
of them exhaust the normal complement to the full rank-one parent
manifold; a is the complex displacement tangent to that manifold but
transverse to the resonant curve. There is no independent repulsion
direction hidden among the other two-site bond entries at this order.
All particle-number-changing terms are retained.

For the complete number-conserving leading slice a=eta=0, define
`gamma=|Re(zeta)|+|Im(zeta)|`. On every compatible periodic L>=12:

```
0<gamma<=1/10, mu>=(7/10)gamma
    => ground space = the three period-three crystals;

mu>(7/10)gamma
    => gap_L >= mu-(7/10)gamma.                            (2)
```

This is a new all-length result for the leading interacting operator,
proved by summing positive twelve-site blocks with exact rational
certificates. The highest of the three bare one-particle well minima
is above `(71/100)gamma`. Therefore raising mu to fill that third well
already lies strictly inside the gapped crystal region. Three equal
global one-particle minima are possible only at zeta=0, the known
ferromagnetic threshold. See [SATURATION.md](SATURATION.md).

No microscopic unitary compact Virasoro c>=3 is established, and no
universal c<=2 bound is claimed. The theorem does not exclude dressed
three-mode mechanisms below the bare third threshold, general nonzero
a or eta, larger hopping, other lower constraints, or higher-order
effects. Transferring it to a fixed microscopic epsilon still needs a
uniform remainder estimate.

## Status

- **[Proved]** Exact five-real-dimensional normal complement and
  exhaustion of the seven-real-parameter leading family under the
  stated lower-kernel hypotheses. The local sixth coefficient itself
  uses the prior full-sector calculation, whose scope is retained.
- **[Proved, exact computer-assisted arithmetic]** The saturation and
  gap inequality (2) for the displayed leading model, and the bare
  third-well comparison. Fifteen rational PSD checks have largest
  real matrix dimension 74; this is not finite-size extrapolation.
- **[Numerical verification]** Six complete microscopic elimination
  comparisons at L=10,12,13, with arbitrary complex higher pair jets
  and all normal controls, agree within 2.89e-15. Independent block
  embeddings on L=15,18 agree within 1.78e-15.
- **[Numerical, bounded escape test]** An explicit real six-Pauli-term
  family with a=2,zeta=1/100,eta=0,mu=1/2 includes both a positive D3
  repulsion and its charge-changing companion. On L=15,18, it has
  number fluctuations but more than 99.99997% single-residue support
  even after conditioning on at least two particles. This does not
  certify a mixed three-liquid state or a thermodynamic phase.

## How it was established and checked

- Author and AI assistance: OpenAI Codex; assigned configuration
  GPT-6 Astra, Max reasoning. No independent human or agent review is
  claimed, and no subagent was used.
- The local coordinate proof factors an arbitrary physical bond into
  an exact parent with the same facilitation/pair coordinates plus
  chemical, hopping and source terms. It closes the physical-control
  exhaustion gap in the earlier alpha/tau studies.
- The all-length positivity proof uses a separately written open-block
  occupation construction and exact Fraction LDL arithmetic. Global
  embedding multiplicities are proved and independently numerically
  checked against the prior periodic occupation implementation.
- The microscopic Schur comparisons reuse the earlier elimination
  routine; this dependence is explicit and source hashes are recorded.
- The earlier exact-parent ground-space context agrees with the opened
  [Bravyi--Gosset paper](https://arxiv.org/html/1503.04035). No external
  theorem is used to infer a CFT or replace the new PSD certificate.

## Reproducing

From the repository root, with NumPy, SciPy and threadpoolctl:

```bash
python3 -B contrib/qubit_three_mode_deformations_20261118/check.py
python3 -B contrib/qubit_three_mode_deformations_20261118/escape_check.py
```

The scripts write `checks.json` and `escape_results.json` in this folder.
Observed runtimes were about **1.6 seconds** and **0.14 seconds**, with
one worker. The exact arithmetic uses the Python standard library,
not a floating-point PSD tolerance. `requirements.txt` lists dependencies.
No DMRG, background production or large scan was run.

The optional `block_certificate.py --length 12 --gamma .1 --exact`
prints the local thresholds and certificates. Its floating-point
threshold estimates are diagnostic only; the proof uses the rational
vertices, not those estimates. The construction works on complete local
spaces, including blocked-hop diagonals.

## Simplicity and next step

[ESCAPE.md](ESCAPE.md) gives the literal six Pauli monomials, the four
independent real family controls, exact rational points, and the full
bounded escape test. At epsilon=1/5 its coefficient spread is 4.61045,
but the perturbative phase identification is uncontrolled there. At
epsilon=1/20 the spread is 142.334, with a much smaller emergent scale;
neither point is advertised as a proven phase.

The decisive next calculation is a mixed-residue interface energy and
motion calculation in (1) with a!=0, including every number-changing
process, followed by the matched relative stiffness and allowed locking
amplitudes. The finite escape test makes the energetic preference for
one residue the immediate issue. Merely adding a repulsive diagonal
entry, counting three wells, or imposing N conservation would miss it.

## Relation to the repository

This addresses assignment 2 of the
[focused round](../research_coordination_20261116/FOCUSED_ROUND_20261118.md).
It extends the complete-control interpretation of the
[cancellation classification](../qubit_cancellation_classification_20261116/README.md)
and the [resonant-liquid detunings](../qubit_resonant_liquids_20261116/README.md),
while keeping the [exact resonance hierarchy](../qubit_resonance_hierarchy_20261116/README.md)
distinct from interacting conformal liquids. Shared README/STATUS/REPORT
and coordination pages were not edited. Proposed summary text is in
[COORDINATOR_SUMMARY.md](COORDINATOR_SUMMARY.md).
