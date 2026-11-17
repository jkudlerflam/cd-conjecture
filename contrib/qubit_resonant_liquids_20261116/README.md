# Resonant qubit branches: threshold obstruction and detunings

Date: 2026-11-16.

## Claim

For the q=1 and q=3 exceptional dilute paths of the
[cancellation classification](../qubit_cancellation_classification_20261116/README.md),
the displayed sixth-order effective operators at their tuned chemical
threshold are positive graph Laplacians plus positive potentials. On
compatible periodic rings L=mM, m=4 or 3, their complete ground spaces
have dimension L+1 and consist of vacuum and single-residue ferromagnetic
states. Their low-momentum dispersion is quadratic. Three/four band
minima therefore do not establish three/four positive-velocity conformal
modes.

The microscopic models are Hermitian, one-site translation-invariant
nearest-neighbour qubit chains. Their perturbative coefficients and
epsilon are fixed independently of L. The exact statements here concern
the specified effective operators; finite-order fixed-L elimination does
not establish a uniform microscopic thermodynamic limit. No Virasoro c,
or c_eff, is assigned to a new microscopic phase.

## Status

- **[Proved for the displayed effective models]** Full cyclic residue-word
  conservation, graph-Laplacian identities, complete compatible-ring
  kernels, threshold density jump, and the undeformed two-body solution.
  Same-residue finite-momentum bound strings exist, but do not lie below
  the global two-free-particle threshold. Mixed-residue channels have
  repulsive or absent potentials and the scattering lengths listed in
  [the derivation](DERIVATION.md).
- **[Derived coefficients; finite-matrix checked]** The pair-amplitude
  jet alpha contributes both repulsion and facilitated particle creation
  at order six. It is not a number-conserving XXZ detuning. Conditional on
  this local formula, its positive-parent representation and exact
  alpha!=0 threshold kernel, vacuum plus m one-particle W states, are
  proved. This rules out the proposed finite-density cure at that
  threshold, not at every chemical potential.
- **[Numerical, fixed effective-basis convention]** A seventh-order
  residue-changing creation witness is nonzero. Its projection on the
  entire undeformed zero manifold vanishes on the tested compatible
  rings. A higher-order energy-dependent calculation is needed to infer
  actual splitting or dressed conservation laws.
- **[Numerical/exploratory]** A freely adjustable order-six nearest-neighbour
  hop mixes residue labels while preserving the earlier blockades. Two
  selected two-valley slices show no binding in the reported few-particle
  tests. [Detuning analysis](DETUNINGS.md) gives a fixed microscopic q=1
  test point, its coefficients, and the unresolved stiffness, pairing,
  momentum-mismatch and infrared tests. It is not a counterexample.

## How it was established and checked

- Author: OpenAI Codex, GPT-6 family, acting as an AI research assistant;
  no human scientific verification is claimed.
- `leading.py` constructs complete constrained occupation-space matrices
  independently of the microscopic elimination. The coordinating chat
  independently verified the undeformed graph identities and kernel
  counts with integer arithmetic in its
  [residue-sector check](../research_coordination_20261116/README.md).
- `higher_jets.py` uses the sibling's generic Schur recursion but constructs
  the microscopic Taylor coefficients here. It checks the entire sixth
  matrix, including complex alpha and beta, rather than only fixed-N
  blocks. This is not an independent implementation of that recursion.
- The alpha kernel proof applies the periodic rank-one-parent result in
  Bravyi--Gosset Proposition 3, opened and read, with the transfer matrix
  in its nontrivial Jordan case. The formula and source are in
  [DERIVATION.md](DERIVATION.md).
- The two-body checks use both a relative-coordinate reduction and a
  separate finite-ring occupation implementation. Numerical bounds are
  explicitly distinguished from analytic no-binding statements.

## Reproducing

From the repository root, with Python, numpy, scipy and threadpoolctl:

```bash
python3 -m pip install -r contrib/qubit_resonant_liquids_20261116/requirements.txt
python3 -B contrib/qubit_resonant_liquids_20261116/check.py
python3 -B contrib/qubit_resonant_liquids_20261116/scattering.py
python3 -B contrib/qubit_resonant_liquids_20261116/relative_detuning.py
python3 -B contrib/qubit_resonant_liquids_20261116/candidate.py
```

Each script writes its matching JSON file. The matching `.txt` files are
captured stdout. BLAS threads are limited to one inside the scripts.
Observed runtimes on the available local machine were about 42, 8, 15,
and 4 seconds respectively; they are recorded in the JSON outputs.

`check.py` covers graph identities through L=24 for q=1 and L=21 for
q=3, complete alpha spectra through the listed L<=16, 24 full sixth-order
jet comparisons on L=10,12,13, and seventh witnesses on L=12,...,15.
The largest sixth jet discrepancy was 6.40e-14. The graph zero counts are
exact integers and agree with L+1; matrix comparisons use floating-point
tolerances. Sibling source hashes and the baseline commit are recorded in
`checks.json`.

`scattering.py` uses 2000-site relative half-lines, N=2 rings L=48,72,96,
and N=3 rings L=24,36. `relative_detuning.py` uses relative cutoffs
512,1024,2048 at same/opposite-valley total momenta and gamma=0.02,0.1,0.3.
`candidate.py` checks six full detuned Schur matrices and exports the
fixed bond and low-degree bare momentum mismatches. No production
finite-entanglement calculation is claimed.

## Relation to the repository

This follows the new exceptional dilute branches without modifying the
earlier [large-chi registration](../qubit_convergence_20261116/README.md)
or claiming that its convergence blocker is resolved. It supplies a
many-body obstruction to naive valley counting and a concrete next
detuning to investigate. The general resonance hierarchy and possible
finite upper bound on interacting central charge remain separate open
questions. The result is summarized in the repository README and REPORT.
