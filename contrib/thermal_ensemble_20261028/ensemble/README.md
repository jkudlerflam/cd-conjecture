# Interacting physical-to-rooted pressure transfer

## Claim and status

**[proved]** A nearest-neighbor physical chain with a conserved marker
number has the same thermodynamic pressure as its complete all-length
rooted cell representation, summed over marker counts. This holds with
length-dependent interacting bonds. The finite-volume partition functions
differ by factors bounded independently of physical length and cell count.

The proof in [PINNED_MARKER.md](PINNED_MARKER.md) pins the full marker
subspace at one physical site, removes/restores two bounded physical bonds,
and adds the bounded rooted seam. For C44 the marker rank is 1116, so the
finite prefactor is `d/q=4`. Arbitrarily long cells and illegal words are
retained. The common-bond factorization used in C46 is unnecessary.

This is a coordinate and ensemble result, not a microscopic CFT theorem.
A low-temperature theorem proved only in C45's well-formed invariant
sector still needs a treatment of the other sectors. The new
[positive completion](ILLEGAL_SECTORS.md) supplies an explicit repair:
add positive endpoint pins on illegal selector pairs, leaving every
legal sector unchanged, and absorb the positive multi-token readout
into the onsite operator. The resulting full all-string cell space
has the needed gap and a linear length cost for heat-trace control.
The local dimension remains 4464. The
physical theorem cannot be applied directly to a constrained projected
code by replacing its local Hilbert space with a smaller dimension.

## Checks and provenance

The proof was derived by a GPT-6-family Codex agent; the parent agent
independently proposed the pinned-marker method. A separate agent
independently checked the pin/seam proof and reran its checks. No human
review is claimed. [checks.py](checks.py) constructs generic complex Hermitian
marker-conserving physical bonds and independently constructs rooted
operators by reanchoring root-marker moves.

The saved [checks.json](checks.json) records 13 full-matrix cases at local
dimensions 2, 3, 4; marker ranks 1, 2; and physical lengths 3 through 7.
All basis-counting, Hermiticity, seam-support, marker-conservation and
39 finite-temperature partition comparisons pass.

Reproduce with `python3 contrib/thermal_ensemble_20261028/ensemble/checks.py`
from the repository root, using a Python with numpy. It takes under one
second in the tested environment. No scipy or repository model imports
are needed. These checks validate the coordinate identities and bounds,
not the central charge of C44.

`python3 contrib/thermal_ensemble_20261028/ensemble/illegal_checks.py`
uses only the standard library and checks 3,280 selector cells, 2,080
marker reanchorings, 176 core pair rewrites, 12 endpoint source types,
and the positive-completion support. Its output is
[illegal_checks.json](illegal_checks.json).

## Relation to the program

The theorem removes cyclic gluing as an additional pressure-transfer
assumption when the **full** rooted microscopic ensemble is controlled.
Covariant local Schrieffer-Wolff transformations preserve its traces at
each fixed cell count and total length, so the physical comparison may be
made before that transformation. The remaining massive matching and
infrared analysis can then be performed in rooted coordinates.

The companion length-fugacity calculation identifies the physical pressure
with an interior zero of the rooted cell pressure, where its trace-class
hypotheses hold. A selected-code zero can serve as a comparison quantity
without separately identifying it as the pressure of a constrained
physical projected Hamiltonian.
