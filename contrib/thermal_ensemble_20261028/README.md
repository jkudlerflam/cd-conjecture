# Physical ensembles and the remaining thermal infrared step

## Claim and status

**The c-d conjecture remains open.** This campaign proves ensemble and
thermal comparison results needed by the C44/C45 mobile construction.
None establishes its microscopic Virasoro central charge.

[research judgment] For the original periodic-spectrum question, the d=24
vacancy and d0=4464 wrapper remain the main mobile candidates. The former
has a smaller explicit target and an unresolved dressed-defect comparison;
the latter now has a clearer analytic route through all-length massive
matching and physical pressure transfer. For an infinite-volume bulk
interpretation, the parallel [frozen d=23 route](../decisive_attack_20261028/README.md)
also deserves priority after its new full-pressure defect bound. Forced
defects on incompatible rings still qualify that interpretation. A direct
general proof lacks a bound independent of the conformal crossover length.
This ranking is research judgment, not a probability estimate.

## Results with proofs

- [Pinned physical marker](ensemble/PINNED_MARKER.md): a rank-q local pin
  changes a finite-d nearest-neighbor partition function by at most a
  constant factor times `exp(O(beta))`, uniformly in ring length. Exact
  all-word cell coordinates and a bounded root seam transfer the full
  interacting physical pressure to rooted cells summed over marker counts.
  Length-dependent bonds and arbitrarily long cells are retained. For
  the literal C44 alphabet, `d/q=4464/1116=4`.
- [Length fugacity](exact/LENGTH_FUGACITY.md): the physical pressure is the
  interior zero of the rooted pressure as a function of length penalty.
  This applies to interacting cell chains with a trace-class neighborhood
  and bounded cut strength. Orthogonal prefix-length subspaces prove
  existence of the selected fixed-physical-length renewal pressure.
  Applying C46's thermal comparison after the length tilt gives an
  exponentially small full-reference versus selected-reference pressure
  difference. The partial trace need not preserve fixed physical length.
- [Tilted auxiliary embedding](TILTED_EMBEDDING.md): length is preserved
  exactly in the auxiliary representation. Linear length costs give a
  uniform tilted gap, heat trace, fixed-order derivatives, and an
  exponential bound on invalid embedding states. The unbounded length
  operator stays in the semigroup. The separate
  [complex-tilt extension](../direct_attack_20261028/infrared/TILTED_MATCHING.md)
  treats derivatives of the massive influence action.
- [Scalar thermal criterion](thermal/THERMAL_CRITERION.md): a precise
  finite-temperature vacuum-scale comparison implies `o(T^2)` even when
  convergence to the infrared fixed point is logarithmic. Those scalar
  estimates are not yet derived for the encoded interacting chain.
- [Exact retarded quadratic control](thermal/QUADRATIC_MEMORY.md): Gaussian
  elimination preserves the thermal coefficient with its complete scalar
  normalization. The frequency residue changes dressed velocities;
  discarding it gives an incorrect coefficient in a checked example.

Every theorem states its hypotheses and proof. An all-word coordinate
representation does not by itself verify the stronger positive onsite
and small intercell decomposition on every malformed sector.

## Remaining work

The separate [massive theorem](../direct_attack_20261028/infrared/INFINITE_CELL_PROGRESS.md)
controls thermal periodization and scalar terms before integrating logical
fields. The ensemble results here then apply once the microscopic
all-sector decomposition and the pressure-zero location are established.
The latter has an explicit low-temperature bracket including all scalar
shifts and the retained bandwidth.

The decisive analytic task is the multichannel retarded logical RG,
including normalized marginal flows and either positive-level current
algebras or scalar vacuum estimates relative to the dressed Luttinger
reference. A small fixed pressure error is not `o(T^2)`. A thermal
coefficient is velocity weighted and does not alone identify the
unweighted central charge. Work continues on these obligations and the
full microscopic sector application.

## Checks and provenance

Four parallel OpenAI GPT-6-family Codex agents worked in this folder,
coordinated with two other research chats and the GitHub notebook. The
coordinator rederived the pin, pressure-zero, orthogonal-prefix,
tilted-trace and scalar estimates and reran all checks. A separate worker
and the direct-attack coordinator reviewed the tilted embedding. No human
referee or formal proof assistant was used.

Python with numpy is sufficient. From the repository:

    python3 -B contrib/thermal_ensemble_20261028/ensemble/checks.py
    python3 -B contrib/thermal_ensemble_20261028/exact/checks.py
    python3 -B contrib/thermal_ensemble_20261028/check_tilted_embedding.py
    python3 -B contrib/thermal_ensemble_20261028/thermal/checks.py

Saved controls include 13 generic complex physical/rooted matrix pairs
and 39 thermal comparisons, 18 tilted reference inequalities, four renewal
controls, 18 invalid-embedding comparisons, and 30 heat-kernel derivative
bounds. All pass. These check algebra and the stated comparisons; they
do not measure a central charge for either candidate.
