# An actual microscopic bond matrix and its reusable assembler

**Status:** implemented finite local operator, with exhaustive local checks
and finite-cell projection checks. This directory does **not** instantiate
the proposed high-channel controller or establish a c-d violation.

`bond.py` implements the entire repeated two-site interaction in
[MICROSCOPIC_HAMILTONIAN.md](../../research_continuation/MICROSCOPIC_HAMILTONIAN.md).
The input `CoreData` contains actual core-bond/boundary/port column functions,
charge representations, and finite cell length. It contains no unspecified
controller hidden inside an implementation stub. Supplying the intended
high-channel core remains an obligation of the larger project.

Every microscopic site is ordered as `(selector, wrapper, helix_A, helix_B)`.
The wrapper is a direct sum of three vacancies and three occupied tags,
each holding one core record and a filter qubit. `encode` and `decode` fix
the basis ordering. `MicroscopicBond.column(j)` returns **every** nonzero
entry in column `j` of the two-site bond, including actions on illegal
selector/tag sectors. The `base`, `hopping`, `move`, and `field` components
can also be queried separately. This representation does not require a
dense matrix, and the core operators themselves may be sparse callables.

The implementation includes the wrapper's oriented endpoint terms, both
selector cut patterns X# and R#, the two independent helix pins, compensated
charged ports, the complex vacancy-conversion phase, and the charge-offset
correction to the marker chemical potential. All are two-site operators
with coefficients independent of the eventual periodic ring length.

## Fully instantiated fixture

To exercise this assembler without assuming the missing controller, use
the explicit qubit core of length three:

\[
h_c=|s\rangle\langle s|+|11\rangle\langle11|,
\qquad |s\rangle=(|01\rangle-|10\rangle)/\sqrt2,
\qquad a_c=b_c=0.
\]

Its complete zero-energy code is the vacuum and the normalized uniform
one-particle state. The core gap is 1/2. Each endpoint raising matrix
compresses to the logical creation matrix divided by sqrt(3), including
the required one-mode right parity convention. The physical U(1) number
acts as logical number with zero reference offset.

The wrapper selects physical lengths four and five. With U=1,
vacancy hopping 1/8, helix pitch 1/4, and s=10^-7, the full site has
**d=240**. The file [`fixture_bond.npz`](fixture_bond.npz) contains the
entire 57,600 by 57,600 Hermitian bond in compressed coordinate form:

```python
import numpy as np
data = np.load("contrib/controller_instantiation/local/fixture_bond.npz")
# h[data["row"][j], data["column"][j]] = data["value"][j]
# data["shape"] gives the two-site matrix shape.
```

There are 122,192 stored entries. This is an actual finite microscopic
matrix with a completely specified basis, not an effective Hamiltonian.
It tests the assembly mechanism with **one** logical fermion per cell;
it is not evidence that c exceeds d-1. No thermodynamic CFT determination
or sufficient infrared smallness threshold for the chosen s is claimed.

The exported entries are floating-point evaluations of the displayed
mathematical formulas. The finite vacancy ground states determine the
compensators and transport coefficient. Their known reflection symmetry
is imposed after diagonalization to remove numerical mixing across a
small even/odd doublet; the eigenvector residual is checked independently.

## Verification

Run with Python and NumPy:

```
OPENBLAS_NUM_THREADS=1 python3 contrib/controller_instantiation/local/checks.py --export
```

The recorded [`results.json`](results.json) verifies:

- the complete three-site core spectrum and both encoded endpoint matrices;
- the wrapper code, spectrum and ports in the 256- and 640-dimensional
  selected occupation sectors, with code residuals below 4e-13;
- equalization at both cell lengths and the direct helix transport overlap;
- **all 57,600 columns** of the full microscopic bond, with Hermiticity
  error below 7e-15, exact U(1) and marker-number conservation, and agreement
  of its separately evaluated components on sampled columns (every 127th).

The finite wrapper spectrum check is restricted to the indicated occupation
sector; the earlier analytic selector/wrapper exclusion argument handles
other global sectors. The exhaustive two-site matrix check, in contrast,
does include all local illegal states. Neither is a finite-size extraction
of central charge.

The independent `port_projection_checks.py` also reads every actual hopping
column and contracts it through the endpoint matrix-unit compressions of
the full finite-cell codes. For all four selected length pairs, the result
matches the desired hopping to relative error below 5.8e-13. It checks X#
and R# contributions separately, and repeats the contraction with the
compensators removed to exhibit the unwanted length-dependent factors.
Run it with the same interpreter after `checks.py`; its output is
[`port_projection_results.json`](port_projection_results.json).

Provenance: Codex, GPT-6 family, 2026-10-28 UTC. The assembler was written
by the coordinating agent and receives an independent direct port-projection
check from a separate circuit agent. No human review is claimed.
