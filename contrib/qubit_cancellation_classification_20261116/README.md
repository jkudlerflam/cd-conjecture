# Cancellation near the qubit blockade: two classical branches

Date: 2026-11-16. AI assistance: OpenAI Codex, GPT-6 family. No independent
human or agent review. This contribution keeps every particle-number sector.

## Claim and scope

The microscopic model is a fixed finite Hermitian nearest-neighbour bond on
qubits, repeated with one-site translation symmetry, independent of system
size. We study analytic perturbations of `H0 = sum n_j n_(j+1)`, on periodic
rings and in the bulk. Any reference to c means the sum of the usual
nonchiral Virasoro charges of the **complete** unitary conformal bulk theory.
No central charge is measured or bounded here.

**[Proved algebra]** The complete first-order-invisible space has seven real
dimensions before removing scale and common onsite rotations. The complete
second-order operator, including complex amplitudes and number-changing
terms, is given in [CLASSIFICATION.md](CLASSIFICATION.md). A nonzero leading
adjacent-pair amplitude prevents a classical cancellation. With it set to
zero, all second-order diagonal cancellations reduce to a chemical potential
plus a nonnegative distance-two repulsion.

There are **two** exponentially degenerate classical coexistence loci:

* The dilute locus forbids particle separations one and two.
* The dense locus forbids both `11` and `000`: consecutive particles have
  spacings two or three. This locus is missed by restricting the chemical
  counterterm to the vacuum-side value used in the earlier nested calculation.

**[Derived local algebra; full-matrix checks]** On the dilute locus, the
fourth-order operator has only a nearest-neighbour hop and diagonal terms.
After cancelling the hop and chemical term, its two penalties are
`(2/3)|q-3|^2 D3 + |q-1|^2 D4`. They cannot both vanish. All those fourth-order
counterterms are tangent to a three-parameter rank-one-parent family, not
just the two previously exhibited curves. The exact parent identity is
proved; it does not constrain arbitrary detunings or supply a c bound.

**[Derived local algebra; full-matrix checks]** The dense locus has an
additional invisible third-order source `f`. Its fourth-order scalar locus
is `|f|^2 = 3|q-1|^2`. It is not tangent to a positive rank-one parent near H0.
Fifth-order conversion leaves either `(q,f)=(1,0)` or a circle with a fixed
phase relation. A sixth-order correlated-move witness is nonzero on the
remaining locus. The precise scope and verification status of this
finite-order obstruction, including higher microscopic jets, are in
[HIGHER_ORDER.md](HIGHER_ORDER.md).

**[Derived sixth-order slices; exact and numerical full-matrix checks]**
The dilute exceptional values q=1 and q=3 generate distance-four and
distance-three hopping after the correct further projection. At q=0,-1,
the next tuned constraint excludes separations through eight. The complete
sixth operators on these four slices, their relation to transfer-matrix
resonances, and the limits of a valley argument are in
[DILUTE_SIXTH.md](DILUTE_SIXTH.md). In particular, degree two is not a
permanent ceiling on the effective one-particle band after further tuning.

## Status and limitations

This is a cancellation classification at stated perturbative orders, not
an exhaustion theorem for qubit critical points. In particular:

* Classical fourth-order detunings can select other constrained subspaces.
* The dilute branch, including its exceptional `q=1` and `q=3` kernels,
  still needs an all-orders dynamical analysis.
* The dense sixth-order witness rules out making the whole stated space
  scalar at that order; it does not bound the c of its surviving dynamics.
* All finite-ring perturbation series are taken at fixed L. A
  thermodynamically uniform elimination and control of neutral/composite
  sectors are additional, unproved requirements.

The mathematical local-boundedness lemma still needed by compactness is
stated in [HIGHER_ORDER.md](HIGHER_ORDER.md). Neither finite parameter
dimension, finite-size continuity, nor a one-particle hopping range replaces
that lemma.

## Establishment and reproduction

The derivation uses local transition paths and nested Schur complements.
The numerical implementation independently embeds the literal 4x4 bond;
it does not import the C185/C186 or earlier nested-cancellation code.
The first two orders reproduce and extend the interpretation of C185/C186;
their general second-order formula is prior work, not a new discovery here.

From the repository root, with the dependencies in `requirements.txt`:

```bash
python3 contrib/qubit_cancellation_classification_20261116/check.py
python3 contrib/qubit_cancellation_classification_20261116/symbolic_check.py --sixth --generic-sixth --dilute-sixth
python3 contrib/qubit_cancellation_classification_20261116/check_dilute_sixth.py
python3 contrib/qubit_cancellation_classification_20261116/check_jets.py
```

The recorded runtimes are approximately 33, 74, 72, and 2 seconds,
respectively, with one BLAS thread on this machine. They write
`checks.json`, `symbolic_results.json`, `dilute_sixth_results.json`, and
`jet_results.json`. The largest ring has L=16. No raw eigenvector arrays
are committed.

There are 18 general complex second-order comparisons, 30 fourth-order
comparisons, 16 dense fifth-order comparisons, 40 dense sixth-order
witness comparisons, 24 complete dilute sixth-order comparisons, 20
higher-jet checks, and 80 exact-parent bond checks. Symbolic tests verify
the full fourth-order matrices with free complex amplitudes, the full
dense fifth operator on its cone, the generic sixth witness after exact
polynomial reduction, and all four dilute sixth matrices with rational
arithmetic. The largest numerical identity error is below 1e-10.
Direct spin-chain diagonalization separately checks the leading fifth-
and sixth-order scaling on finite rings. These are distinct implementations
within this contribution, not independent scientific review.

## Relation to the project

Inputs: [compactness note](../smaller_d_unbounded_20261116/COMPACTNESS.md),
[revised qubit conjecture](../d2_updated_conjecture_20261116/README.md),
[C16](../../codex/work_C16/RESULT.md), [C185](../../codex/work_C185/PROOF.md),
[C186](../../codex/work_C186/PROOF.md), and the
[facilitated-term audit](../cusp_validation_20261109/FACILITATED_REPAIR.md).

The older `STATUS.md` headlines do not turn conditional infrared arguments
into lattice-to-CFT theorems. The later qubit note correctly leaves even a
finite universal B(2) open. Our calculations neither confirm microscopic
c=2 nor justify c<=2. They address one singular neighborhood needed for a
possible finite bound, and expose a further branch that must be controlled.

The rank-one ground-space facts were checked against Bravyi and Gosset,
[Gapped and gapless phases of frustration-free spin-1/2 chains](https://arxiv.org/html/1503.04035),
Theorem 1 and Propositions 1--3. Their open-chain gap theorem is not a
central-charge theorem or a statement about every nearby Hamiltonian.
