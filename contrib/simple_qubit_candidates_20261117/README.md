# Simpler qubit candidates: calibrated full-chain diagnostics

Date: 2026-11-17. OpenAI Codex, GPT-6 family. **Campaign in progress;
no new microscopic counterexample or thermodynamic central charge is
established.** Completed raw runs are under `outputs/`; `analysis.json`
contains convergence gates and all parity/size comparisons. A run marked
`done` means its bounded schedule finished, not that every gate passed.

## Claim and model class

We investigate Hermitian, one-site translation-invariant,
nearest-neighbour chains of physical qubits in the four-coupling family

```
H = sum_j [(XX-YY)/2 + kappa/(4V)(XY-YX) + (V/4)ZZ]_(j,j+1)
    + (-V/2+a/V) sum_j Z_j.
```

Coefficients are fixed independently of L and the full onsite field acts
on every open-chain endpoint. The target is a unitary compact conformal
bulk with total Virasoro c>1. No c_eff reinterpretation is used.
[MODEL.md](MODEL.md) gives exact source conventions and infrared limits.

The older C190/C195 production point remains the evidence benchmark,
with V=40,kappa=1,hz=-19.976047250613814. Its previously audited
256-to-512 main-window reading is 1.982737, but its last entropy
increment 0.010375 remains above the registered 0.01 threshold.
The [earlier registration](../qubit_convergence_20261116/REGISTRATION_20261116.md)
and conclusions are unchanged. The epsilon=10^-3 sixth-order construction
is retained as a [mechanism study](../qubit_resonant_liquids_20261116/README.md),
not the practical candidate for this simplicity-first search.

## Simplicity comparison

There are **four independent real coupling parameters, three after
overall scale fixing**. In the original common Pauli basis there are
six nonzero monomials: XX,YY,XY,YX,ZZ,Z. The same common onsite rotation
U=exp(-i pi Z/8) at all sites gives **four monomials**:

```
(1/2+D)XY+(1/2-D)YX+(V/4)ZZ+hz Z.
```

The total field is counted once, and constants and telescoping terms
are removed. No unequal coefficients are grouped to reduce the count.
Both literal coefficient spreads are reported with Delta=1:

| Point | Original-basis spread | Common-rotation spread | Generated leading scale | Evidence role |
|---|---:|---:|---:|---|
| V=40, production field | 3196.17 | 40.46 | 1/80 | Strongest existing size/convergence evidence; still provisional |
| V=8,kappa=1,a=24/25 | 124.16 | 8.2773 | 1/16 | Preferred simpler follow-up under current diagnostics |
| V=12,kappa=1,a=24/25 | 284.16 | 12.3548 | 1/24 | Intermediate-hierarchy comparison |
| V=4,kappa=1,a=1 | 28 | 4 | 1/8 | Older entropy anomaly with strong density modulation; exploratory |

The rotation shortens the expression but encodes the small hopping in
the difference of the two cross couplings. Invariant ratios
(Delta/2)/D=2V and Jz/D=V^2 remain. A generated small kinetic scale or
mode velocity is distinct from the microscopic coefficient spread.
The [sibling term analysis](../qubit_minimal_terms_20261117/README.md)
independently confirms the count and leaves some simpler faces open;
it does not establish that all four structures are universally necessary.

The explicit rational V=8 follow-up bond, in the common rotated basis,
is

```
H = sum_j [(17/32) X_j Y_(j+1)+(15/32) Y_j X_(j+1)+2 Z_j Z_(j+1)]
    -(97/25) sum_j Z_j.
```

It is exactly the original family with pairing scale one, D=1/32,
V=8 and a=24/25. Its finite-size evidence must not be read as proof
that this exact rational point has a c>1 thermodynamic limit.

## Numerical status

The initial a=1 pilots at V=8,12 were substantially denser than the
older candidate. Sequential rational-field bracketing selected the
common a=24/25,kappa=1. At L=32,64 both candidate ground branches pass
the residual, sweep-stationarity, discarded-weight and chi-step gates.
The main-window 32-to-64 slopes are about 2.4; the exact two-leg control
itself reads 2.389 there. Narrow and wide windows show large finite-size
effects, so those slopes alone do not identify a CFT.

The L=128 extension through chi=192 is complete: the V=8 and V=12
main-window slopes are 2.059737 and 2.109011. Both fail the 0.01 chi
gate, with entropy steps of 0.0417--0.0514; neither is a converged CFT
measurement. The preferred V=8 point is continuing to chi=256. The
one-/two-leg controls are compared with exact correlation-matrix
entropies, whose corresponding slopes are 1.004026 and 1.858274. The
L=128 two-leg DMRG control still has energy drift 1.22e-6 and an entropy
error -0.003144; its failed gate is retained and the exact value is the
calibration reference. The
two-leg control has range two and is used only for calibration, not as
an admissible nearest-neighbour candidate. The old
V=8 infinite-chain reading near 1.4 is retained as unresolved contrary
evidence, with its different fitted coefficients and older-solver
limitations stated in the source audit.

`analysis.json` records densities, density modulation, valley weights,
connected separation-two pair correlations, parity gaps, every parity
pairing of size estimators, and the extrapolation flag on any secondary
density-corrected estimate. Balanced valleys do not rule out pairing.
No Luttinger exponent is fitted from a few short distances. Symmetry-
allowed two-pair charge locking remains a physical concern on this line.
The raw valley convention assigns both boundary momentum bins to one
half; `valley_half_boundary` splits those bins equally and removes the
spurious roughly 0.16% imbalance. Neither convention diagnoses charge
locking. The stricter registered stopping rule and the four local
diagnostic gates are reported separately; `schedule_limit_hit` remains
visible even when the looser local gates pass.

## Verification, reproduction and remaining tests

- **[Checked numerically]** Independent literal Pauli assembly at L=8
  agrees with the real computational-gauge MPO; both parities agree with
  ED in energy, every entropy cut, density and a pair observable. Full
  endpoint fields are explicitly checked. This validates conventions
  and local numerics, not an infrared theory.
- **[Numerical]** The unchanged C190 parity solver supplies all candidate
  and control states. Source files and old-data hashes are recorded in
  `source_audit.json`. Its shared solver is not an independent DMRG
  implementation. No human scientific review is claimed.
- **[Derived and numerically checked]** `simplicity.py` verifies the
  common rotation and reports both counts/spreads. The term-removal
  results are attributed separately to the sibling contribution.

Install `requirements.txt`, then from the repository root:

```bash
python3 -B contrib/simple_qubit_candidates_20261117/audit_sources.py
python3 -B contrib/simple_qubit_candidates_20261117/validate.py
python3 -B contrib/simple_qubit_candidates_20261117/simplicity.py
python3 -B contrib/simple_qubit_candidates_20261117/run.py <config.json>
python3 -B contrib/simple_qubit_candidates_20261117/analyze.py
```

`campaign.py` runs an explicit list of configs with at most two processes,
each with one BLAS thread. Warm continuations require their earlier
checkpoint, regenerated by its config; checkpoints stay ignored in
`runtime/`. `REGISTRATION.md` gives the initial diagnostic rules and
`DECISIONS.md` records subsequent evidence-based choices. Per-sweep
runtime, local residual and truncation data remain in every output.

To establish a simpler c>1 phase still requires converged larger sizes
and chi, an infrared diagnostic independent of a small-ring entropy
slope, control of pairing/charge locking, and a parameter neighbourhood
whose phase is verified rather than inferred from raw density. The
absent old production environment blocks its specific very-large-chi
test; it has not been used as a reason to omit these smaller calculations.
