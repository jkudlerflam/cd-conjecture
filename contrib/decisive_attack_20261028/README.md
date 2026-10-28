# A shorter frozen-chain route and concrete infrared checks

## Status and current direction

**The c-d conjecture remains unresolved.** This campaign develops a shorter
route to a bulk counterexample and checks what the single-mode demonstrator
can actually establish. It does not identify a measured violating central
charge or assert a new conformal-limit theorem.

**[research judgment]** For the usual infinite-volume bulk interpretation,
the frozen d=23,n=23 construction now deserves priority alongside the mobile
d=24 construction. Its certified 46 Fermi points already exist. The theorem
below controls all frozen defects in the full low-temperature pressure, so
one can focus on the homogeneous finite-cell chain's infrared theory.
This avoids the extra mobile length fluid for this particular route.
The distinction from an all-ring finite-size conformal claim must remain:
the pressure theorem does not eliminate forced defects or their exact
translation degeneracy on incompatible ring lengths.

The decisive remaining task is a controlled source-dependent infrared
analysis of the small, high-degree retarded logical interaction produced
by massive matching. A full-rank chiral current algebra could give the
lower bound c>=23, using C46's existing Sugawara argument. A thermal proof
must instead retain the vacuum terms and identify each mode's contribution.
Neither conclusion follows just from bounded action coefficients.

## Results

- **[proved; independently checked under the imported sector hypotheses]**
  [Frozen defect pressure](frozen/DEFECT_PRESSURE.md): the exact full
  microscopic free energy satisfies
  `0 <= f_cell(T)/b - f_full(T) <= C T exp(-delta/T)`.
  Wrong cells of every length, illegal-arrow sectors and arbitrary physical
  ring lengths are included. The proof uses trace almost-additivity and
  explicit positive defect activities. Rational witnesses at b=24 give
  C=1200 and a conservative delta greater than 8.8e-8 for the stated
  selector interval. The temperature threshold is extremely small; this
  is an existence estimate. No payload infrared asymptotics are assumed
  in the pressure inequality itself.
- **[proved algebra; matching conditional on the existing theorem]**
  [A weak-XX demonstrator](demonstrator/README.md): the ten-state model has
  exact second-order anisotropy
  `|Delta_XXZ|=3|Jp/Jz| kappa^2/[4(4D^2-kappa^2)]`.
  Choosing a small fixed nonzero Jp/Jz puts it arbitrarily near XX.
  The cell gap, parity, 45-mode auxiliary embedding and strong coefficient
  norm requirements are verified explicitly. The intercell coupling is
  chosen afterwards. Twelve exact rational resolvent checks pass.
- **[proved; independently checked]**
  [Branch-resolving local sources](sources/BRANCH_SOURCES.md): trigonometric
  interpolation constructs a finite-range observable for each distinct
  Fermi crossing. An analytic source-matching lemma bounds every fixed
  finite source derivative by the leakage norm, under the existing
  finite-cell theorem. This supplies initial sources, not their infrared
  Ward identities. Exact checks cover 2920 interpolation entries through
  46 branches; an independent Fourier-space check also passes.
- **[proved; coordinator checked]**
  [Exact-encoding obstructions](exact/PROOF.md): shallow circuit dressing
  retains the separated-boundary obstruction to an exactly invariant
  hopping code. An exponentially local exact dressing must have tails
  comparable to the hopping scale. The candidate's tiny hopping scale
  remains compatible with this inequality; it does not refute the candidate.

The separate [two-flavour selector audit](demonstrator/N2_SELECTOR_AUDIT.md)
finds an unsupported global-sector inference, a missing infinite-length
endpoint in the tail bound, and saved zero-coupling eigenvalues inconsistent
with explicit zero modes. It supplies a selected-sector variational upper
bound and keeps those findings separate from the reported numerical c data.

## What would complete the proposed route

For fixed n=23, preserve the certified band margins while choosing the
conversion strength small enough to control all coefficient norms; choose
the final port coupling afterwards. Apply finite-cell massive matching
with physical sources. Then prove the retarded-kernel marginal-flow and
source comparison, and obtain a sufficient conformal/current or thermal
identification. The full frozen-sector pressure theorem transfers any
proved selected-cell thermal coefficient to the complete microscopic chain.

No limit in n or cell size is needed for a single d=23 violation. Constants
may depend arbitrarily badly on this fixed finite cell. Couplings must be
fixed and nonzero before the thermodynamic limit. An argument demanding
conformal finite-size spectra on every ring still needs separate treatment
of incompatible lengths; the mobile construction is relevant to that goal.

## Reproduction and provenance

From the repository root, run:

```sh
python3 -B contrib/decisive_attack_20261028/frozen/check_pressure.py
python3 -B contrib/decisive_attack_20261028/demonstrator/verify_tuning.py
python3 -B contrib/decisive_attack_20261028/sources/check_interpolation.py
python3 -B contrib/decisive_attack_20261028/exact/checks.py
python3 -B contrib/decisive_attack_20261028/exact/source_checks.py
```

The first three use the Python standard library; the last two need numpy.
The interpolation check uses exact fractions and can take longer than the
other checks. Outputs are written only beside their scripts. Finite checks
supplement the written proofs and do not constitute an infrared simulation.

Four GPT-6-family Codex agents worked in parallel, with the coordinator
reviewing the defect-pressure and exact-encoding arguments and rerunning
the certificates. A separate agent audited the branch-source theorem.
Other active chats independently own `direct_attack_20261028` and
`thermal_ensemble_20261028`; this campaign preserves their work and links
their results where used. No human review is claimed.
