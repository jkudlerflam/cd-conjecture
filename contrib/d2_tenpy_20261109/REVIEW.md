# Independent finite-DMRG review

OpenAI Codex (GPT-6 family), separate coordinating-team reviewer,
2026-11-09. AI-only code/data review; no new large simulation was run.
Reviewed `run.py`, `analyse.py`, `validate.py`, `validation.json`, the
final saved finite JSON records through C183 L=128 and both C181153 parity
sectors at L=32,64,128, and REPORT §9.155. “Final saved” means the requested
run returned; several runs stopped at the sweep limit without converging.

**The larger L=128 microscopic C181153 check does not confirm the earlier
c≈2 signal: every parity-sector pairing gives a 64→128 entropy slope below
one.** This negative size trend must not be hidden behind the convergence
caveat. The runs are incompletely optimized, so it also does not establish
asymptotic c=1 or a gap. The transverse-DM points separately continue to
read near one. The implementation and reported three-cut estimators pass
this review; no robust microscopic c>1 conclusion survives the larger-size
test currently available.

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
| C181153, lower reported energy at each L | 2.14523037, provisional | 0.78647328, unconverged |
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
The new C181153 results are evaluated in detail below. Do not average the
sector entropies or interpret one finite-size jump as a second conformal
mode.

REPORT §9.155's cautious conclusion is supported. Two small wording fixes
would improve precision: `run.py`'s opening “no number/parity projection”
description is stale now that optional parity sectors are supported; and
the `frac=.25` fit removes the outer quarters, leaving the **central half**,
so “central-quarter fit” should be renamed. The infinite-chain failed runs
are properly excluded, and no conclusion should be salvaged from their
reported correlation lengths.

This review does not inspect unpublished running jobs or certify a
thermodynamic phase, exact ground states, or a counterexample.

## Earlier C181153 32→64 signal and its sector/density checks

All six saved records use the identical literal nearest-neighbour
qubit Hamiltonian, with no blockade or fixed-particle-number projection:
Jx=0.5002549063202599, Jy=-0.4997450936797401, Jz=10,
Dz=0.006244799657946725, Dy=0, hz=-19.97589125. Parity conservation is exact.

| L, sector | Shifted energy | Three-cut entropy | Middle-half density | Final discarded weight |
|---|---:|---:|---:|---:|
| 32, even | -0.777184222709 | 1.3121480504 | 0.21609504 | 6.20e-11 |
| 32, odd | -0.777264754544 | 1.3532321408 | 0.25494289 | 1.27e-10 |
| 64, even | -1.584442503872 | 1.6010588714 | 0.21825608 | 1.23e-9 |
| 64, odd | -1.584393384639 | 1.5844057687 | 0.20162920 | 9.27e-10 |

Selecting the lower energy therefore selects odd at L=32 and even at L=64.
The energy margins are 8.05e-5 and 4.91e-5. These exceed the observed last
four-sweep energy spans (about 3.4e-8 at L=32 and 9.3e-7/7.9e-7 at L=64),
so the selection is stable relative to the visible drift. There is no
residual-based global ground-state error certificate, however.

The selected two-size slope is 2.14523037. A useful sector sensitivity test
is to compute every pairing, without claiming the excited-sector curves
are ground-state central charges:

| Sector at 32 → sector at 64 | Three-cut c proxy |
|---|---:|
| odd → even, lower energies | 2.14523037 |
| odd → odd | 2.00107828 |
| even → even | 2.50086125 |
| even → odd | 2.35670916 |

Every pairing remains well above one. Thus the parity switch by itself
does **not** manufacture the positive slope. At L=64, both sectors also
give a central-half profile coefficient near 2.05 (2.05039 and 2.05516);
all four recorded fit intervals yield 2.05–2.40. For the selected lower
energies, relative central windows of half-width 0.10, 0.125, 0.20 and
0.25 times L give two-size proxies 2.1849, 2.2880, 2.2718 and 2.3254.
This is stronger finite-size evidence than a single large profile fit,
and is distinct from the near-one transverse-DM results.

Two qualifications remain substantial. First, the selected bulk density
falls from 0.25494 to 0.21826, and the local profiles have sizable Friedel
oscillations. The even-to-even comparison has nearly matched mean density
(0.21610→0.21826) and still gives a large slope, which argues against mean
density drift as its sole cause, but does not remove boundary/shell effects.
Second, both L=64 sector runs exhaust the requested 40-sweep budget (the
engine records 41 iterations), rather than satisfying strict convergence. Their
last energy changes are -4.86e-7 and -7.75e-7; their last mean-entropy
changes are +9.48e-6 and -9.30e-5. Small discarded weights, canonical norm
errors near 1e-14, and energy-expectation agreement near 2e-14 establish
good internal normalization, not completed optimization. Further sweeps,
larger χ and an independent start remain useful checks.

The signed even-minus-odd difference changes sign. Its absolute value
times L is 0.00257702 at L=32 and 0.00314363 at L=64 in the literal Pauli
energy units used by these JSON files. Here the leading effective scale
is |t2|=(Jx-Jy)^2/(2V)=1/80=0.0125, so the dimensionless values
L|E_even-E_odd|/|t2| are 0.20616150 and 0.25149047. That is compatible with
a small finite-size splitting but is not a thermodynamic gap test from
two sizes. It is also **not** the fixed-N curvature
E(N+1)+E(N-1)-2E(N) measured in C181's number-conserving effective model.
Do not equate the two diagnostics.

The density homogeneity check was independently recomputed as well.
C181's actual implementation convolves by a 16-site box and takes the
central half of the **shortened convolution array**. This is narrower
than the physical central half of the original chain at these small L.
The following table uses precisely that implementation:

| L, sector | Total occupation | Smoothed minimum | Smoothed maximum | Smoothed mean | (max-min)/mean |
|---|---:|---:|---:|---:|---:|
| 32, even | 8.01013785 | 0.19988589 | 0.21690158 | 0.21120616 | 8.0564% |
| 32, odd | 9.00497405 | 0.23603224 | 0.25494289 | 0.24862074 | 7.6062% |
| 64, even | 16.01522889 | 0.21767702 | 0.22577849 | 0.22170471 | 3.6542% |
| 64, odd | 15.02329348 | 0.20100978 | 0.20852687 | 0.20472875 | 3.6717% |

Both L=64 sectors pass C181's 25% finite-size homogeneity threshold by a
wide margin. As a stricter boundary sensitivity check, selecting smoothed
windows by their actual center coordinate in the original physical
central half gives:

| L, sector | Minimum | Maximum | Mean | Relative range |
|---|---:|---:|---:|---:|
| 32, even | 0.19988589 | 0.26474797 | 0.22521211 | 28.8004% |
| 32, odd | 0.23603224 | 0.29855801 | 0.26128167 | 23.9304% |
| 64, even | 0.21767702 | 0.23503437 | 0.22292165 | 7.7863% |
| 64, odd | 0.20100978 | 0.21780018 | 0.20594406 | 8.1529% |

The L=32 even-sector result highlights why the precise window definition
matters, while both L=64 sectors pass even the wider physical window.
The L=64 outer quarters have mean densities about 0.28222 (even) and
0.26785 (odd), above their central-half means 0.21826 and 0.20163. Thus
these microscopic open chains show **edge enhancement**, not the boundary
depletion seen in the number-conserving effective-model data. The latter
boundary picture must not be transferred without matching the microscopic
open-chain reduction. The computed interior profiles do not diagnose a
central dense droplet and pass the registered screen, but this is a finite
homogeneity check, not a proof against thermodynamic phase separation.

This positive 32→64 reading justified the larger-size check. The following
L=128 results supersede it as the assessment of whether that signal persists.

## L=128 check: c≈2 is not sustained, with unresolved optimization

Both parity-sector jobs have now returned their final files:

| Sector | Shifted energy | Three-cut entropy | Middle-half density | Total occupation | Final discarded weight |
|---|---:|---:|---:|---:|---:|
| even | -3.199287936003 | 1.6863152143 | 0.22006166 | 29.88745228 | 2.65e-8 |
| odd | -3.199296931209 | 1.6919158280 | 0.21392232 | 29.20462129 | 3.78e-9 |

The odd sector has the lower last reported variational energy, but its
8.9952e-6 advantage is comparable to the even-sector last-sweep change
(+5.7700e-6). Sector selection is therefore not certified. Importantly,
**changing the selection does not rescue the entropy-growth result**:

| Sector at 64 → sector at 128 | Three-cut c proxy |
|---|---:|
| even → odd, lower last energies | 0.78647328 |
| even → even | 0.73799342 |
| odd → odd | 0.93062538 |
| odd → even | 0.88214551 |

For the nominal lower-energy selection, the relative central windows used
above give 0.7716, 0.7719, 0.8328 and 0.9105. At L=128, central-half profile
fits fall to 0.74875 (even) and 1.07585 (odd), from approximately 2.05 in
both L=64 sectors. Other recorded intervals span 0.75–1.23 and 0.91–1.37,
respectively. Thus the negative trend is present in multiple estimators,
not merely in one parity switch or one chosen cut. The exact finite-size
controls at 64→128 read 0.97227 and 2.02543.

Mean density drift is now modest: the nominal ground-sector comparison
changes 0.21826→0.21392, and even→even changes 0.21826→0.22006, yet both
slopes fall sharply. The 16-site-smoothed C181 homogeneity test still passes:

| Sector | Smoothed minimum | Maximum | Mean | Relative range |
|---|---:|---:|---:|---:|
| even | 0.21518091 | 0.22312232 | 0.21901457 | 3.6260% |
| odd | 0.20919129 | 0.21720855 | 0.21291663 | 3.7654% |

The wider physical-central-half version gives 3.6252% and 3.7644%. These
data do not diagnose a droplet as the cause of the entropy change.

The physical interpretation remains open for concrete numerical reasons:

- **Possible physical crossover:** entropy slopes fall below one across
  all sector pairings, profiles move toward one, density is relatively
  stable and the interior stays smooth. These are negative indications
  for the earlier two-mode interpretation, consistent with crossover
  toward fewer gapless modes or an eventual gap.
- **Incomplete/local optimization:** both L=128 jobs exhaust the same
  40-sweep budget. The even run's final mean-entropy change is +0.004055;
  its last eight mean entropies rise from 1.4638 to 1.4860. The odd run's
  final change is +7.57e-5, with +0.00205 accumulated over the last eight
  records. The saved entropy profiles also have reflection defects
  max|S(l)-S(L-l)| of 0.00995 and 0.02636, versus below 0.001 at L=64.
  The finite Hamiltonian has inversion-conjugation symmetry, so these
  defects are an optimization/near-degeneracy diagnostic for a putative
  unique ground state. They are not evidence of a new bulk phase.
- **Finite χ remains untested:** both runs use χ=128. Their relatively
  small discarded weights do not rule out entanglement suppression or
  trapping in a low-entropy variational state. No fixed-L χ sequence or
  independent optimized L=128 start is supplied here. Additional sweeps
  and larger χ are needed to distinguish that possibility from the
  physical crossover indicated by the current data.

Recovering c_two=2 for the nominal 64→128 pair would require raising the
L=128 central entropy by about **0.1402**, from 1.6919 to 1.8321. This is
a substantial change, not something justified merely by citing a small
stopping-tolerance violation. Conversely, the current histories provide
no rigorous bound excluding such a correction after better optimization.

The latest saved files have canonical norm errors around 1e-14 and energy
expectation differences below 4e-14. Earlier intermediate norm warnings
should not be attributed to these final files: the remaining issue is
state optimization, not their present normalization. The nominal L=128
parity splitting gives L|gap|=0.00115139 in physical energy units, or
0.09211091 after division by |t2|=0.0125, but the near-degenerate energy
ordering and unfinished sweeps make it unsuitable for a gap claim.

The updated conclusion is **an unconfirmed initial c≈2 reading, with a
negative larger-size check and a clear convergence task outstanding**.
Neither a counterexample nor a microscopic c=1 theorem follows.
