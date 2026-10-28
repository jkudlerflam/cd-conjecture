# Direct proof and counterexample campaign, 2026-10-28

## Claim and status

**The conjecture remains open.** This continuation compares the two current
mobile constructions and attacks their missing microscopic-to-infrared steps.
It does not reinstate the withdrawn claim that the d=24 counterexample is
verified.

Starting snapshot: `c4b223b`. Four parallel Codex agents worked in one shared
checkout: the coordinating agent, a vacancy-matching worker, an infrared
worker, and an independent alternative-route worker. Results are AI-assisted,
not human-refereed. Each subfolder records the scope and checks of its claims.

## Current route ranking

This ranking is **[research judgment]**, not a proof or numerical confidence
estimate.

1. **For refuting c <= d-1 at a concrete small dimension: the d=24 vacancy
   construction.** Its n=23 finite band certificate and frozen-sector
   comparison already exist. The main extra task is the finite-J,
   common-bulk-subtracted comparison uniformly over competing material
   patterns, followed by the low-energy matching needed to retain all 24
   modes. Repeating finite certificates will not resolve that task.
2. **For the strongest current analytic route and unbounded c at fixed d:
   C39/C42/C44/C45, d0=4464.** The finite spatial construction and coupling
   schedule now exist. Its main missing input is now the retarded-action
   multichannel infrared theorem and its physical-ensemble application.
   The new abstract thermal matching theorem below handles an infinite
   excited cell space under explicit hypotheses; their complete microscopic
   instantiation is still required. C46 supplies useful partial theorems,
   but explicitly does not prove all of S. The large
   alphabet is a disadvantage for the first numerical violation, not a
   mathematical objection to a fixed-n existence argument.
3. **For a direct general proof:** no current route provides the needed
   bound independent of the conformal crossover scale. Restricted-class
   theorems remain useful, but do not establish the full conjecture.

The precise finite-J task should use a one-sided sector inequality or a
length-weighted remainder. An absolute error proportional only to the number
of arbitrarily long wrong cells was disproved in the preceding
[matching attack](../counterexample_audit_20261028/matching_attack/README.md).

## New work

- [Vacancy matching](vacancy/README.md): an explicit local change of basis
  gives a volume-independent ultraviolet remainder with one spatial moment.
  The initial finite-depth proof has cutoff power 18; a
  [translation-covariant single exponential](vacancy/SINGLE_EXPONENTIAL.md)
  improves that power to 10. At s=J/t and cutoff delta=s^a, 0<a<1/10, its
  physical remainder is o(J).
  This addresses a specific technical gap in the earlier matching route.
  It does not prove the remaining low-frequency energy comparison.
- [Scope clarification](SCOPE.md): the original conjecture does not state
  ground-state uniqueness. A fixed number of ordered sectors does not
  multiply the entropy central charge. The frozen construction still has
  an unresolved all-ring spectrum issue and the payload infrared assumption.
- [Hidden free-duality exclusion](alternatives/RESULT.md): six local Pauli
  anticommutation obstructions show that a uniform nearest-neighbor qubit
  Hamiltonian admitting a faithful term-by-term map to single Majorana
  bilinears already rotates to the standard quadratic or commuting Ising
  class. An exact check covers all 4,096 Pauli supports. This excludes a
  particular shortcut to c>1, not general interacting qubit chains.
- [Finite-cell massive matching](infrared/FINITE_CELL_MATCHING.md): a finite
  spectral embedding and a published fermionic determinant expansion give
  an exact small logical action correction, linear in the leakage, after a
  locally block-diagonal reference. The bounds include spatial/time moments
  and do not divide by the slow logical bandwidth. They require the stated
  strengthened coefficient norms. The independent
  [audit](alternatives/IR_AUDIT.md) found and repaired three bookkeeping
  issues. This advances fixed finite-cell matching; the subsequent logical
  RG/CFT limit remains open.
- [Arbitrary-cell thermal matching](infrared/INFINITE_CELL_PROGRESS.md):
  an operator-valued auxiliary representation and convergent time-block
  polymer expansion remove the finite-cell dimension dependence. With an
  explicit gap, heat-trace bound and small exponential interaction norm,
  the full thermal action correction is linear in leakage, with time/space
  moments, scalar density and source control. Its kernels approach the
  periodized zero-temperature kernels exponentially in inverse temperature.
  A separate [audit](alternatives/INFINITE_CELL_AUDIT.md) checked the proof,
  including the thermal trace and continuous-time limits. C46 supplies a
  candidate heat bound on its specified rooted space. Physical sectors,
  fixed length and the final multichannel RG/CFT still require proof.
- [Length-cutoff stabilization](alternatives/LENGTH_CUTOFF.md): projected
  perturbative coefficients agree exactly below order
  `4 C(K-b)`, where `C(r)=sum_(j=1)^r ceil(j/4)` and lengths through K
  are retained. The quadratic growth comes from the finite supply of holes
  in nearby cells. This is a coefficient identity, not a convergent
  cutoff-removal theorem. The exact combinatorial checks pass in 870 finite
  cases and 24 explicit paths.
- [RG initialization](vacancy/RG_EXTENSION.md) and
  [degree reparameterization](vacancy/DEGREE_REPARAMETERIZATION.md):
  bandwidth normalization changes every pinned kernel norm by the same
  factor 1/W, regardless of field degree. Finite microscopic degree and
  analyticity recover the powers of the small coupling required by a
  published RG tree estimate. The full marginal-flow and counterterm
  induction is the remaining task; power counting alone does not close it.

Further worker results are added only after their statements and scope have
been reviewed. The active goal is to solve the conjecture, not to stop at this
partial reduction.

## Checking and reproduction

The coordinator independently checked the vacancy proof's filter signs,
causal-cone gate count, interaction-norm multiplicity, and cutoff exponents,
and reran the finite-matrix and exact Pauli-support checks. A separate agent
reviewed both the UV proof and the finite-cell massive matching argument.
The proof uses the standard finite-range Lieb-Robinson estimate for bounded
interactions and no many-body gap. Its detailed assumptions and limitations
are in [UV_NORMAL_FORM.md](vacancy/UV_NORMAL_FORM.md).

Research provenance: Codex GPT-6-family coordinating and parallel agents,
2026-10-28. No human review or numerical measurement of a violating central
charge is claimed.
