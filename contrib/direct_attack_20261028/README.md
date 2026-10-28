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
2. **For unbounded c at fixed d: C39/C42/C44/C45, d0=4464.** The finite
   spatial construction and coupling schedule now exist. Its main missing
   input is the stated Wilsonian decoupling and infrared stability statement
   S, with the physical ensemble and all cell lengths retained. C46 supplies
   useful partial theorems, but explicitly does not prove S. The large
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

- [Vacancy matching](vacancy/README.md): an explicit finite-depth local
  change of basis gives a volume-independent ultraviolet remainder with one
  spatial moment. At s=J/t and cutoff delta=s^a, 0<a<1/18, it is o(J).
  This addresses a specific technical gap in the earlier matching route.
  It does not prove the remaining low-frequency energy comparison.
- [Scope clarification](SCOPE.md): the original conjecture does not state
  ground-state uniqueness. A fixed number of ordered sectors does not
  multiply the entropy central charge. The frozen construction still has
  an unresolved all-ring spectrum issue and the payload infrared assumption.

Further worker results are added only after their statements and scope have
been reviewed. The active goal is to solve the conjecture, not to stop at this
partial reduction.

## Checking and reproduction

The coordinator independently checked the vacancy proof's filter signs,
causal-cone gate count, interaction-norm multiplicity, and cutoff exponent.
The proof uses the standard finite-range Lieb-Robinson estimate for bounded
interactions and no many-body gap. Its detailed assumptions and limitations
are in [UV_NORMAL_FORM.md](vacancy/UV_NORMAL_FORM.md).

Research provenance: Codex GPT-6-family coordinating and parallel agents,
2026-10-28. No human review or numerical measurement of a violating central
charge is claimed.
