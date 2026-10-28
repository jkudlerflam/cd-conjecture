# Frozen-cell continuation: thermal reduction and spectral obstruction

## Claim and status

**[proved; independently checked by the coordinating AI agent]** The frozen
23-letter candidate has exponentially small defect corrections to its exact
homogeneous selected-cell free energy:

    0 <= f_cell(T)/b - f_full(T) <= C T exp(-delta/T).

The statement includes every wrong cell length and every illegal-arrow
sector, at fixed nonzero microscopic coupling. It uses the already proved
local sector bounds, extends their boundary-energy argument to partition
functions, and sums an explicit positive overcount of all defect patterns.
See [the complete proof](DEFECT_PRESSURE.md).

**[proved]** Under the same hypotheses, the ground-space degeneracy grows
at least linearly with ring length on incompatible lengths `b` not dividing
`L`. Thus the full all-ring finite-size spectrum has a genuine obstruction
to a single fixed compact-CFT partition function, despite the thermal
reduction. See [the precise obstruction and proof](ALL_RING_OBSTRUCTION.md).

Neither theorem establishes the interacting payload CFT. No counterexample
to the c-d conjecture is declared.

## Application to the candidate

C33 certifies a 23-band reference at local dimension `d=23`, selected cell
length `b=24`. The inherited selector proof establishes the local hypotheses
for every `0<|lambda|<=2^-20` and `K=10`. Conservative exact rational
witnesses in this folder give

    delta >= 8.83523286277916e-8  (decimal display),
    C = 1200,
    T <= 1.104404107847395e-8    (sufficient displayed threshold).

The JSON file records the exact rational bounds underlying these displays.
These temperature scales are only defect bounds. They are far above the
extremely small payload bandwidth at the certified parameter point and do
not assert that the payload CFT is visible at those temperatures.

**[conditional implication]** If the exact selected-cell chain satisfies

    f_cell(T)=u-(pi/6) sum_a(c_a/v_a^cell) T^2+o(T^2),

then the full microscopic chain has the same coefficient per physical
length, with `v_a^physical=b v_a^cell`. Independently resolving the modes or
proving a scaling-limit theorem is still needed to identify the individual
central charges. One scalar thermal coefficient alone does not count them.

**[research judgment]** This makes the frozen d=23 route materially shorter
for a bulk interpretation of the conjecture: only the selected finite-cell
microscopic-to-infrared step remains, rather than the additional gapless
vacancy matching. The target is `c=23>d-1=22`. The strict all-ring spectral
interpretation remains obstructed and would still motivate the mobile
constructions.

The original paper was opened directly: [Latorre and Sierra, arXiv:2403.17242v2,
Section 2](https://arxiv.org/html/2403.17242v2#S0.SS2). It does not impose ground-state
uniqueness or explicitly forbid spontaneous translation breaking. It also
does not resolve the two infrared interpretations above. We do not add an
unstated assumption to its conjecture or claim that a thermal coefficient
alone fulfills its spectral hypothesis.

## Reproducing the finite checks

Run from the repository root:

    python3 contrib/decisive_attack_20261028/frozen/check_pressure.py

Standard library only, under a second on the development machine. It writes
[pressure_witnesses.json](pressure_witnesses.json), evaluates the selector
using rational alternating-series intervals, verifies all positive constants,
and tests the analytic convex-envelope estimate at 211 representative
lengths including its crossover and long tails. The proof covers every length;
the finite checks are diagnostics and positivity witnesses.

The elementary interval formulas come from the existing
[`check_global_sectors.py`](../../counterexample_campaign/route_a/check_global_sectors.py).
No central charge is measured by this script.

## Provenance and relation to other work

Author: Codex GPT-6-family subagent, 2026-10-28. The coordinating Codex agent
independently checked the finite-temperature inequalities, the combinatorial
sector overcount, and the all-length exponential tail against the original
microscopic-sector argument. No human review is claimed.

The existing [scope note](../../direct_attack_20261028/SCOPE.md) established
that fixed ground-state degeneracy does not alter a logarithmic entropy
coefficient. This continuation adds the full thermal defect sum and proves
that the incompatible-ring degeneracy is in fact unbounded. It does not
modify the mobile-chain infrared work being pursued by other agents.
