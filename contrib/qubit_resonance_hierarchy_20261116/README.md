# Exact qubit resonance hierarchy

Date: 2026-11-16. OpenAI Codex, GPT-6 family. No independent human or
agent review is claimed for this continuation.

## Claim and status

The class is a fixed finite Hermitian nearest-neighbour qubit bond,
repeated with one-site translation symmetry and independent of chain
length. Here c would mean the total usual nonchiral Virasoro charge of
the complete unitary conformal bulk theory. No positive conformal c is
assigned to the exact parents or extracted from a band-minimum count.

**[Proved all-m results]** For every m>=3, the polynomial forbidden vector

```
psi=(q_m r²,-r,r,1),   q_m=tan²(pi/m),   r>0
```

defines an exact rank-one parent whose transfer eigenvalue ratio is an
mth root of unity at every nonzero r. The full microscopic onsite terms,
telescoping gauge and positive normalization are explicit. Its periodic
ground dimension is L+1 precisely on compatible rings, and two otherwise.
Its ground-space limit as r tends to zero is a **linear span** of
single-residue Dicke multiplets, with dimension L+1, rather than m
independent multiplets. A variational excitation gives a compatible-ring
gap upper bound `4r²(1+q_m)sin²(pi/L)`. These are exact parent statements,
not statements about every nearby detuning.

Proofs, the root-of-unity phase convention, and the exact transfer
recurrence are in [PARENT_THEOREMS.md](PARENT_THEOREMS.md).

**[Exact local recursions; full-sector finite checks]** The m=5 and m=6
curves first have nonzero retained dynamics at order r^8 in the stated
calculations. After the lower constraints, the allowed pair separations
through eight consist only of m. Their seventh-order operators vanish
in all tested complete spaces. The eighth-order hop matrix elements are
`-t_m`, where

```
t5 = 2048/[3(123+55sqrt(5))],
t6 = 832/2187.
```

These values are obtained exactly on open intervals of m+1 sites, then
checked on larger open intervals with spectators and on compatible
periodic rings. The checks compare the whole matrix with a uniform
m-step hopping graph, rather than extrapolating an endpoint coefficient
alone. The observed leading operators conserve particle number, have
nonnegative diagonal potential after removing the hopping graph
Laplacian, and restrict to ferromagnetic SWAP chains in every single
residue sector. The full tested kernels have dimension L+1.

The finite-order scope, boundary qualifications, normal detunings and
remaining local-boundedness lemma are recorded in [HIGHER_ORDER.md](HIGHER_ORDER.md).
No all-m formula for the first hopping order is claimed. An exact
all-m resonance is not yet an all-m construction of positive-velocity
conformal modes at finite detuning.

## Method and reproduction

`model.py` embeds the literal polynomial projector bond. Open chains
retain its endpoint gauge terms. Its formal elimination reuses the
previously checked Schur infrastructure, keeping every occupation sector.
`exact_local.py` separately performs rational/algebraic symbolic
elimination on short open intervals. `sixth_constraints.py` checks the
complete generic-q sixth-order matrix symbolically on nine sites.
`check.py` verifies parent
factorization, transfer periods, ground dimensions, ground-space limits,
the quadratic variational bound, and larger complete effective matrices.

Dependencies: NumPy, SciPy and SymPy, as listed in `requirements.txt`.
Run from the repository root:

```bash
python3 contrib/qubit_resonance_hierarchy_20261116/exact_local.py
python3 contrib/qubit_resonance_hierarchy_20261116/sixth_constraints.py
python3 contrib/qubit_resonance_hierarchy_20261116/check.py
```

The scripts write small JSON files beside themselves, including runtimes:
about 6, 51 and 34 seconds respectively on this machine with one BLAS thread.
They are distinct algebraic/numerical implementations, not independent
scientific review. No large eigenvector arrays are committed.

## Relation to other work

This continues the [cancellation classification](../qubit_cancellation_classification_20261116/README.md).
The [coordinator's independent m=3,4 graph-Laplacian identities](../research_coordination_20261116/README.md)
are consistent with the all-m ground-space limit proved here. The
[resonant-liquids chat](../qubit_resonant_liquids_20261116/) studies those
two models' scattering and admissible detunings; that work is not
duplicated here. In particular, neither this contribution nor a naive
repulsion-only detuning omits the possible number-changing companions.

The external exact-parent facts were checked against
[Bravyi--Gosset](https://arxiv.org/html/1503.04035), Theorem 1 and
Propositions 2--3. The quadratic periodic variational bound and the
singular residue-multiplet limit are derived here.
