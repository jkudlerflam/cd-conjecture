# Counterexample campaign, 2026-10-28

Status: active research, no established counterexample claimed.

This campaign uses four simultaneous Codex workers: one coordinator and
three independent research workers. It starts from collaborator commit
`eff9cfc`, including C33's new finite-g, n=23 Fermi-point certificate.
The original conjecture concerns the total Virasoro central charge of
unitary conformal infrared sectors, not an effective central charge.

## Priority

The most concrete existing candidate is the explicit 23-letter route-A
bond in `../explicit_candidate/`. C33 now certifies 46 simple crossings
and positive quadratic and quartic nonresonance margins at n=23. This is
a finite-matrix statement, not yet a microscopic CFT determination.
Frozen cells create a separate ring-sector problem. The existing mobile
wrapper has d=2736 and therefore needs substantially more modes.

We are testing a cheaper mobile completion. Add one vacancy state to the
23-state physical alphabet, move it by a species-independent adjacent
swap, and apply a weak copy of the original bond only to pairs of
occupied sites. Ordered occupied labels form a squeezed chain. At zero
internal coupling the charge positions are free spinless particles;
first-order averaging gives a scalar multiple of the original bond on
the squeezed labels. If finite-coupling spin-charge separation and the
n=23 payload infrared inference can both be established, this would give
a particularly concrete d=24, c=24 counterexample.

That conditional statement is a research target, not a result. A
finite-ring charge gap must not be used as a uniform thermodynamic gap.

## Independent work

- `route_a/`: global cell-length and illegal-sector bounds for the
  explicit 23-state bond.
- `infrared/`: finite-coupling infrared stability, including the new
  vacancy completion and its allowed mixed vertices.
- `sparse_core/`: audit of the alternative species-graph controller.
- `doping_check/`: independent finite-system squeezed-space checks.
- Coordinator: source verification, candidate formulation, independent
  checks, synthesis, and GitHub synchronization.

All workers are AI assistants in the GPT-6 family. No human review is
claimed. Results must distinguish proofs, finite numerical checks,
standard continuum inputs, and open hypotheses. Scripts and exact
parameters are retained alongside their reports.

## First checkpoint

[proved finite identities; numerical cross-checks] The vacancy-lift bond is
implemented in `doping/bond24.py` and a moderate-scale 576-by-576 fixture
is saved. All three exact U(1) checks and Hermiticity pass. Independent
squeezed-space checks cover q=2,3 rings with full dimensions 24 to 270,
including the ring twist and complete first-order charge-ground projection.

[proved] The original chain's global sector proof in `route_a/README.md`
now excludes all wrong-length and illegal sectors on compatible rings
for a fixed nonzero interval of port couplings. This does not supply the
finite-coupling infrared conclusion of the vacancy construction.
