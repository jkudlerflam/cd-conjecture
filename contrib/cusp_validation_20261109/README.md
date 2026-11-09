# Independent cusp reduction check

Status: exact local perturbative coefficient derived and independently checked, 2026-11-09. AI-assisted by GPT-6 Codex. No counterexample or thermodynamic conclusion is claimed.

This contribution independently checks the low-energy reduction of the nearest-neighbour qubit candidate in C173. It does not estimate the central charge.

## Result

The first omitted term creates or annihilates a pair at distance three, with coefficient `2 Delta^3/(3V^2)`. It survives at the candidate with imaginary hopping. At `V/Delta=40` its coefficient is `1/30` of the leading hopping scale, so the full local remainder begins at relative order `Delta/V`, not `(Delta/V)^2`. Direct low-energy atomic pairing remains off resonance by momentum; this finding does not invalidate the second-order bound-pair result. It identifies the additional elimination and stability argument needed to transfer that result to the microscopic chain.

The [derivation and scope correction](PROOF.md) prove the coefficient by local paths. The independent [matrix checks](check_reduction.py) pass on 45 full-spin perturbation matrices through L=12. The [saved data](results.json) also compare short-ring spectra, with their vacuum and finite-size limitations stated explicitly.

## Reproduction

```bash
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 python3 contrib/cusp_validation_20261109/check_reduction.py --max-L 12 --spectrum-L 12 --output contrib/cusp_validation_20261109/results.json
```

Requires numpy and scipy. The recorded run took about 10 seconds on the local machine with at most two numerical-library threads. No repository modules are imported.

## Registered checks, before calculation

1. The coefficient at physical order `Delta^3/V^2`, with `w = w0 Delta^2/V`, should consist of a projected distance-three pair creation term with coefficient `2/3` and a projected distance-two pair term with coefficient `-2 Re(w0)`. At the candidate `w0 = i/2`, the latter cancels.
2. Direct finite-basis perturbation matrices must agree with this formula to `1e-12` for open chains through at least L=10, and for both imaginary and generic complex hopping.
3. Full-chain and second-/third-order effective spectra will be compared at the actual candidate chemical potential `mu=(mu_eff+6)/(2V)`, `mu_eff=-2.2`. The old C173 spectrum check used a vacuum-side chemical potential with the opposite sign of the leading shift, so it does not itself validate finite density.
4. Spectral errors will be reported as measured; no finite-size spectral match will be treated as a thermodynamic CFT result.

A preliminary open-chain calculation showed strong edge filling at the actual chemical potential, because an open chain with uniform fermion chemical potential has different boundaries from the uniform-field spin bond. The spectral comparison therefore uses periodic spin rings. Open-chain matrices are retained only as local coefficient checks. This switch was made before inspecting periodic spectrum output.

## Microscopic repair audit

A proposed transverse-DM perturbation can generate the desired distance-two repulsion using a genuine nearest-neighbour qubit bond. Its complete reduction also generates a density-assisted particle creation term at the same perturbative order. The [full derivation and tuned model](FACILITATED_REPAIR.md) show why a free adjustment of the density interaction is insufficient. Six independent full-spin matrix checks pass exactly; use `check_facilitated.py` to reproduce `facilitated_results.json`.

The [all-particle-number ring diagnostic](repair_ed.py) also verifies its matrix against the direct microscopic reduction, then checks 48 finite systems through L=20. At `mu_eff=-2.27`, the repulsion-only model stays in the vacuum on the L=20 ring, but the complete deformation already has negative ground energy at `lambda/Delta=0.2`. The omitted term therefore materially changes the onset. [Saved results](repair_ed_results.json) report energies, densities and particle-number variances; no central charge is extracted.

A further scan uses the coordinator's small-r valley-retuning prescription and retains all blockaded particle-number sectors through L=24. The finite-level pair-versus-atoms combination becomes positive near `lambda/Delta=0.15`, but the positive margin falls sharply with L. [Details and limitations](FACILITATED_REPAIR.md#6-small-r-valley-retuning-promising-finite-spectra-unresolved-limit) and [corrected-tuning data](repair_ed_retuned_L24_results.json) are recorded. This is a promising concrete deformation to investigate, not evidence sufficient for infinite-volume unbinding or `c=2`.
