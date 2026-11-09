# Next calculation after the microscopic L=128 downturn

AI-only recommendation, OpenAI Codex (GPT-6 family), 2026-11-09. This is a
prospective diagnostic protocol, not a new result. No simulation was launched.

**Do one fixed-Hamiltonian convergence study of the unrounded C181153
microscopic chain at L=128, in both parity sectors.** First resolve local
solver/sweep error at χ=128, then continue the same states to χ=256 and
384. This directly distinguishes optimizer bias from a persistent entropy
downturn. A new rounded bond or larger V changes the physical question
before the present numerical uncertainty is resolved.

The present negative result is substantive: the selected 64→128 slope is
0.7865; all four parity pairings give 0.738–0.931; the L=128 central-half
profile fits are 0.749 and 1.076. Restoring a slope of two needs about
0.1402 extra central entropy. Meanwhile, the even branch still changes its
mean entropy by 0.00405 per final sweep, and the sector energy difference
is comparable to visible energy drift. See the full [review](../d2_tenpy_20261109/REVIEW.md).

## Freeze the physics, improve only the numerical state

Use `contrib/d2_tenpy_20261109/run.py`, `--name c181153`, with the existing
uniform open-chain fields and exact parity conservation:

```
Jx=0.5002549063202599  Jy=-0.4997450936797401
Jz=10  Dz=0.006244799657946725  Dy=0  hz=-19.97589125
```

1. Optimize both sectors with `--strict-local`, first at χ=128. The existing
   flag disables TeNPy's adaptive relaxation of the Lanczos `P_tol`, retains
   `P_tol=1e-13`, and uses a reorthogonalized Krylov space up to 160. This
   targets an actual difference from the original runs. Save states with
   `--save-state`; continue with `--resume-state`, which checks the exact
   Hamiltonian, L and parity. The old JSON profiles alone cannot resume an
   MPS. If no matching trusted checkpoint exists, a fresh strict run is
   necessary. Use `--checkpoint-every 10` together with `--save-state`
   for long runs: atomic sweep-boundary snapshots survive interruption,
   and resumption canonicalizes the loaded state before rebuilding its
   environments. Do not relabel a new Hamiltonian as a continuation.
2. From the strict χ=128 states, continue to χ=256 and then 384. Keep the
   field fixed. Use a bounded sweep budget, for example 80 additional
   sweeps per stage; exhausting it means “unconverged,” not “passed.”
   Record energy, central entropy, density and symmetry defect after each
   sweep. The wrapper now retains `eng.update_stats`, per-update Lanczos
   counts/cap flags, and a directly measured untruncated local residual
   from one extra Hamiltonian application, with per-iteration summaries.
   `P_tol` is a dimensionless Krylov convergence criterion, not that
   absolute residual. The residual-pass status explicitly provides no
   global ground-state certificate.
3. Strictly recheck the L=64 baseline before the final size comparison,
   and use a second initialization for the leading L=128 sector at the
   final χ. For this wrapper use a different `--init-density`, for example
   0.22 versus the default 0.12: changing `--seed` alone does **not** change
   its deterministic parity-sector product state. Preserve both sector
   records. If the parity ordering is smaller
   than observed numerical changes, quote the envelope over both sectors.

These stages form one convergence experiment, not a parameter scan.
Checkpointing avoids repeating earlier work. No infinite-chain result or
additional L is needed to answer this first numerical question.

## Prospective decision rules

Before interpreting a phase, require a numerical-quality pass: directly
measured final local residuals at most 1e-8 in physical Pauli energy units,
equal to **8e-7 of |t2|=0.0125** here, and no unresolved Krylov-cap failures;
final four total-energy changes
below 1e-8; central-entropy range below 5e-4; canonical norm error below
1e-8; and central entropy changing by at most 0.005 between the final two
χ values and between independent starts. The entropy tolerance corresponds
to about 0.043 in a single-size contribution to the two-size c proxy.
Check that the reflection defect falls below 1e-3, or explicitly diagnose
near-degenerate symmetry-breaking variational states before accepting the
comparison. These are practical acceptance tolerances, not rigorous
ground-state error bounds. A failed quality test leaves the result unresolved.

Retain the predetermined three-central-cut slope. Use the two relative
windows of half-width L/8 and L/4 and the L=128 central-half profile fit as
preselected sensitivity checks; retain the 16-site density homogeneity test.

- **Numerical suppression identified:** tightening the solve at unchanged
  χ or increasing χ lowers the variational energy and materially restores
  the missing entropy. Report which change caused the recovery. The
  original negative size trend remains a result of the original settings.
- **Two-mode finite-size reading recovered:** after the quality pass, the
  three-cut slope, both relative-window slopes and central-half profile
  coefficient all exceed 1.6, with density variation below 25%. This
  restores support for a c≈2 candidate; it still does not prove an infrared
  counterexample.
- **Earlier c≈2 reading fails the improved test:** after the quality pass,
  the three-cut slope and both relative-window slopes remain at or below
  1.3, with no continuing entropy recovery as χ increases. This supports a
  physical crossover explanation and warrants deprioritizing this field.
  It is not a theorem that c≤1 or that the entire parity family is closed.
- **Anything between these outcomes, or disagreement between estimators:**
  inconclusive. Preserve the data and limits rather than selecting a
  favorable window or silently moving the field.

## Why the alternatives come later

| Alternative | Usefulness | Why it is not the next decisive calculation |
|---|---|---|
| Rounded parity bond with the independent Davidson solver | Checks implementation independence and whether decimal simplification preserves finite-size behavior | Changes couplings and has its own unresolved local-residual/cap issues; a positive 32→64 result cannot explain the unrounded L=128 downturn |
| Higher V at matched effective parameters | Tests whether finite-V number-changing vertices destroy the effective-model liquid | Changes the Hamiltonian, requires density/valley matching, and reduces the retained energy scale, making solver conditioning harder before it has been controlled |
| Additional transverse-DM scans | Could find a different phase window | Existing calibrated finite-chain points read near one; it does not resolve the current parity-family discrepancy |

If a well-converged fixed-H test retains the downturn, a larger-V comparison
becomes informative. Match the effective hopping, valley offset and actual
bulk density, re-extract at the physical field, and use the validated solver.
Changing V while merely copying a bare field or a far-vacuum fit would repeat
the chemical-potential mismatch already identified in C177/C186.
