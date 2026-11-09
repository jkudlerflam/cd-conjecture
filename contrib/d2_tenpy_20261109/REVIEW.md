# Independent finite-DMRG review

OpenAI Codex (GPT-6 family), separate coordinating-team reviewer,
2026-11-09. AI-only code/data review; no new large simulation was run.
Reviewed `run.py`, `analyse.py`, `validate.py`, `validation.json`, the
completed finite JSON records through C183 L=128, and REPORT §9.155.

**The implementation and reported three-cut estimates pass this review.
The completed microscopic results do not establish c>1.** Their reliable
interpretation is a finite-size reading near one, with remaining boundary,
density and convergence effects. This does not prove asymptotic c=1 or
exclude a neighboring c>1 interval.

## Hamiltonian and estimators

- TeNPy's `Sigmax/y/z` operators implement the stated Pauli normalization.
  The signs of both DM terms agree with the independent literal matrix.
  The onsite field is uniform, including the endpoints. The subtracted
  scalar is exactly `(L-distance)*Jz + L*hz` for these open chains;
  restoring it recovers the literal Hamiltonian. It cannot affect entropy.
- The L=8 matrix comparison is genuinely independent of the MPO assembly.
  Recorded matrix errors are at most 2.85e-14. The added parity-sector
  validation checks both sector energies and entropies against restricted
  literal ED (entropy error at most 2.83e-12).
- Finite `entanglement_entropy()` has L-1 entries. `S[center-2:center+1]`
  corresponds exactly to cuts L/2-1, L/2, L/2+1. The prefactor 6 in the
  open-chain two-size slope and the chord-profile fit is consistent.
  The exact two-XX-leg entropy is correctly the sum of the entropies at
  the ceil/floor half-cuts of the separate sublattices.
- Pure transverse DM on an open chain includes transverse edge fields
  relative to the facilitated-operator rewriting. This is a boundary
  difference from the effective calculation, not an MPO error. Likewise,
  uniform onsite fields differ at endpoints from distributing a field
  equally across only the existing bonds. The implemented convention is
  explicit and consistently validated.

Independent recomputation from the saved entropy arrays gives:

| Point | c_two 32→64 | c_two 64→128 |
|---|---:|---:|
| C183 | 1.04504415 | 1.10031636, provisional |
| Rational transverse DM | 1.02180052 | not available |
| Original V=40 point | -0.25425975 | not available |
| XX, exact finite-size control | 0.94672123 | 0.97227450 |
| Two XX legs, exact finite-size control | 2.04742178 | 2.02543263 |

The negative original-point slope is a decreasing entropy diagnostic,
not a negative central charge.

The actual two-leg DMRG control reads 1.94134610 at 32→64. Its finite-χ
bias relative to the **exact same finite-size estimator** is -0.10607568:
the central entropy error grows from about 0.000402 at L=32 to 0.012656 at
L=64. The report correctly limits this to a coarse calibration. Quoting
the exact value alongside the numerical control would separate χ error
from finite-size corrections; no new simulation is needed.

## Convergence and interpretation

At C183 L=128, S_middle3=0.5701315147 and middle-half density=0.18503948.
The profile coefficient is 0.985825 for cuts in [L/4,3L/4] and 1.029576
for [L/8,7L/8]. Together with the L=32→64 profile decrease, these do not
provide a stable c>1 signal. The 64→128 central slope of 1.1003 is not a
contradiction: a two-size slope also contains density and boundary changes.
The density is nonmonotone, 0.18261→0.17162→0.18504 over L=32,64,128.

L=128 remains provisional: the last recorded energy change is -5.41e-8
and mean entropy change is -4.85e-6, above the requested strict stopping
tolerances, despite a small discarded weight of 3.71e-11. Low discarded
weight does not certify sweep convergence. The reconstructed energy
expectation agrees with the saved energy within 1.2e-13 and the canonical
norm residual is 1.8e-14, so normalization is not the issue.

The independent dense-start, seed-2 L=64 run agrees with the main L=64
energy within 6.80e-8 and central entropy within 0.000195. This is a useful
metastability check at L=64, not a proof of global optimization at L=128.

Fit-window sensitivity remains material. At 32→64, narrow relative windows
give about 1.03–1.06 for both transverse-DM points. Wider windows can raise
the proxy toward 1.2 or above, while the exact c=1 and c=2 controls also
show strong finite-size/window oscillations. It would be inappropriate to
select a larger slope after seeing the results. The predetermined
three-central-cut estimator is correctly retained in the report.

Large Y is **not** a symmetry-breaking order parameter for this family:
the surviving antiunitary permits it. Neither its magnitude nor the tiny
uniform-Y source established in the separate perturbative audit proves a
gap or removes a mode. Boundary-induced response is also allowed. Such
claims need independent correlation/gap scaling.

## Parity-family precautions and report wording

The code correctly forbids parity projection when transverse DM is nonzero.
For Dy=0, computing both parity sectors and selecting the lower energy is
the appropriate finite-chain ground-state procedure. The energy offsets
are identical across those sectors. Both runs must use identical
couplings/field and have sufficiently converged energies before selection.
At C181153 L=32 the odd sector is lower by 8.05e-5, comfortably larger than
the observed final energy drift (~3e-8), but the two sectors have visibly
different densities and entropies. A parity switch with L can therefore
change a two-size entropy proxy substantially. Do not average the sector
entropies or interpret one finite-size jump as a second conformal mode.

REPORT §9.155's cautious conclusion is supported. Two small wording fixes
would improve precision: `run.py`'s opening “no number/parity projection”
description is stale now that optional parity sectors are supported; and
the `frac=.25` fit removes the outer quarters, leaving the **central half**,
so “central-quarter fit” should be renamed. The infinite-chain failed runs
are properly excluded, and no conclusion should be salvaged from their
reported correlation lengths.

This review does not inspect unpublished running jobs or certify a
thermodynamic phase, exact ground states, or a counterexample.
