# Additional material from "Prove c d conjecture"

Imported 2026-10-27 after comparing the available chat with repository commit
`58ef8a9`. This folder records scientific details that were absent or only
briefly mentioned. It does not announce a proof or counterexample to c-d.

## Claim

The setting remains a Hermitian, one-site translation-invariant Hamiltonian
H_L = sum_j h_(j,j+1) on (C^d)^{tensor L}, with one fixed two-site bond and
couplings independent of L. The central charge is the total Virasoro c of the
unitary conformal sectors, not an entropy slope assigned to a nonconformal
state.

The chat contains an alternative first-order encoded-cell proposal with
r^2 logical complex-fermion channels and proposed physical dimension
O(r log r). Its additional symmetry and auxiliary-state constructions are
recorded below. The required microscopic compiler and infrared theorem are
not supplied, so the proposal remains conditional.

## Status

**Heuristic as a counterexample construction**, with elementary analytic
lemmas identified separately and proved in the accompanying notes. The
catalog is a collection of established examples with source and method
labels, not an exhaustive classification. No new lattice numerics or
independent verification of the original project's numerical runs is claimed.

| File | What is added | Relation to existing work |
|---|---|---|
| [torus_encoding.md](torus_encoding.md) | A species-torus charge assignment, its row/column-sum selection rule, exclusion of marginal Cooper transfer, and the conditional O(r log r) interpreter proposal | Refines the joint boundary compression already in C8; provides a symmetry-based alternative to the parity-only, shifted-momentum target |
| [massive_auxiliary.md](massive_auxiliary.md) | A direct-sum embedding of arbitrary finite excited-cell spectra into free massive auxiliary fermions, with exact algebra separated from the needed integration/RG estimates | Addresses the general-cell obstruction identified in `contrib/ir_control`, without declaring it solved |
| [qubit_catalog.md](qubit_catalog.md) | Explicit Hamiltonians, exact or numerical critical regimes, and range/Hilbert-space qualifications | Expands the brief family list in `specs/_d2_background.md`; supports the C15-C17 program |
| [proof_extensions.md](proof_extensions.md) | Bond averaging, the two-dimensional commutant criterion, the time-reversal-plus-reflection XYZ reduction, and closure under exact finite-depth circuits | Organizes consequences and extensions of the existing C1/C6 and permutation-chain results |

## What this does and does not add

The broad quadratic bound, the U(1)-qubit reduction, the SU(d) fundamental
case, abstract joint completely positive boundary compression, long-cell
selection, and the need for local control of leakage were already recorded.
They are referenced where needed rather than presented as discoveries here.

The species-torus proposal would make c/(d-1) unbounded **as d increases** if
its missing steps were completed. It does not establish unbounded c at one
fixed d. The latter is the stronger objective of the repository's second-order
Route A. A finite-block isometry or a fixed-alphabet clock alone proves neither
claim.

The chat's later statement that only a constructive-RG extension remains is
not adopted: a complete nearest-neighbour interpreter, its alphabet count,
initialization and illegal-sector analysis, and its boundary observables
remain obligations. The frozen-cell construction also retains the compatible-
ring/translated-vacuum qualification. The repository's newer mobile-cell
analysis in C13 is a separate conditional repair, not an automatic consequence
of the chat construction.

## How it was established and checked

- Source: the user's October 2026 chat titled "Prove c d conjecture", imported
  at the user's request. The original assistant's model identity was not
  available in the retrieved record.
- Import and review: Codex (GPT-6), with parallel reviews of the torus,
  auxiliary-state and catalog material on 2026-10-27. Each note identifies its
  analytic checks or primary-source verification. No human referee is claimed.
- Source limitation: the available retrieval truncated two long construction
  replies and the end of the catalog. Only visible material is used; unseen
  conclusions and unavailable citation targets have not been reconstructed.
  The chat's PDF attachment is not imported.
- Private chat links, account identifiers and attachment URLs are omitted;
  these notes contain the scientific content needed by collaborators.

## Reproducing

The new algebraic arguments are written in the notes and require no numerical
run. The catalog links primary sources for its literature claims. Existing
C1, C6, C8 and infrared-control scripts remain in their original locations.
This import adds no numerical evidence for a c-d violation.

## Relation to the rest of the repository

The import supplements [STATUS](../../STATUS.md),
[REPORT](../../REPORT.md), [C8](../../codex/work_C8/RESULT.md),
[infrared control](../ir_control/README.md), and the
[d=2 background](../../specs/_d2_background.md). The separate C14 claim against
the velocity-weighted Casimir proposal was already recorded in
[C14's progress log](../../codex/work_C14/PROGRESS.md); it is not a new result
of this import, and it is not a counterexample to c <= d-1.
