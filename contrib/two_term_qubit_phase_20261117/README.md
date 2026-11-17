# The two-term qubit line: exact obstructions and modest phase diagnostics

## Claim and status

Consider the Hermitian, one-site translation-invariant, nearest-neighbour
qubit chain `H(g)=sum_j[2 X_j Y_(j+1)+g Z_j Z_(j+1)]`, with real g fixed
independently of size. The central charge of interest is total Virasoro
c of a unitary compact finite-velocity continuum theory.

* **[Proved]** For every g!=0, a separated-site term in the commutator
  prevents conservation of the standard boost third charge, even after
  any range-two correction or telescoping density change. The regular
  homogeneous difference-form Yang-Baxter route is excluded. A separate
  claw obstruction excludes faithful term-by-term free-Majorana maps.
  Neither statement is a c bound or a classification of all integrability.
* **[Proved reduction; derived application of a read theorem]** An exact
  Ising duality reduces sufficiently large |g| to a small perturbation of
  a unique gapped product vacuum. Gap stability excludes a positive-c
  critical *bulk* for |g|>G, for a finite but non-explicit G. No moderate
  numerical coupling is certified by this asymptotic argument.
* **[Derived first-order result]** Near g=0 a slow Majorana hopping is
  generated at order g, with `v_B^(1)=8|g|/pi`. This line is not the
  previously tuned cancellation curve. The flat g=0 limit is not a
  complete finite-velocity CFT. We do not promote the finite-volume
  perturbation calculation to an all-order thermodynamic c=1 theorem.
* **[Numerical]** Complete internal-sector diagnostics through L=20 at
  g=1/2 and 1 show approximately 1/L neutral gaps and decreasing density
  and dimer fluctuations. An exact basis map identifies prior C144
  L=128 data at g=1, consistent with c=1 but with sizable fit dependence.
  The g=4 control is strongly ordered; g=2, the equal-coefficient point,
  already trends toward density order. No new c>1 evidence is obtained.

The broad gapless-versus-ordered picture is supported; an all-g c<=1
proof, a phase boundary, and a certified moderate-coupling critical
interval remain open. The simple rational points are not counterexamples.

## Argument and evidence

Read [EXACT.md](EXACT.md) for maps, symmetry sectors and the third-charge
proof; [REGIMES.md](REGIMES.md) for the weak-coupling limitation and
large-|g| exclusion; and [NUMERICS.md](NUMERICS.md) for sector gaps,
order measurements, spectral weights and the C144/C157 audit.
[PROTOCOL.md](PROTOCOL.md) records the bounded diagnostic plan.

This follows the [minimal-term classification](../qubit_minimal_terms_20261117/README.md).
There are two independent coefficients before a common scale is removed,
one ratio afterwards, and exactly two Pauli monomials in the displayed
basis. The spreads at g=1/2,1,2 are 4,2,1. An emergent small velocity is
reported separately from that microscopic simplicity.

## How it was established and checked

Author: OpenAI Codex, GPT-6, with AI assistance for the algebra, code,
literature audit and interpretation. No independent human or separate
agent review is claimed. Two independent representations, literal
complex physical spin matrices and real symmetry blocks, reproduce all
small-chain eigenvalues. Exact integer Pauli algebra verifies the local
charge obstruction, and full periodic matrices verify its sum. All
checks pass; maximum small-matrix discrepancy is 5.9e-14.

The saved numerical data retain residuals and translation leakage.
Additional eigenvectors repair incomplete degenerate multiplets found
in two weak-coupling records; the original records remain visible.
No DMRG was run for this contribution, and no larger four-term campaign
was duplicated.

## Reproducing

From the repository root, with Python, numpy and scipy:

```sh
OPENBLAS_NUM_THREADS=1 python contrib/two_term_qubit_phase_20261117/check.py
OPENBLAS_NUM_THREADS=1 python contrib/two_term_qubit_phase_20261117/spectra.py
OPENBLAS_NUM_THREADS=1 python contrib/two_term_qubit_phase_20261117/spectra.py --lengths 12 16 20 --couplings .125 .5 1 --output neutral_poles.json
OPENBLAS_NUM_THREADS=1 python contrib/two_term_qubit_phase_20261117/spectra.py --couplings 2 --output equal_coefficients.json
OPENBLAS_NUM_THREADS=1 python contrib/two_term_qubit_phase_20261117/spectra.py --lengths 12 16 --couplings .125 --levels 32 --output weak_resolution_check.json
OPENBLAS_NUM_THREADS=1 python contrib/two_term_qubit_phase_20261117/validate_data.py
```

Run sequentially. The script refuses overlapping invocations. The
principal 16-record scan took 49 summed seconds, the neutral-weight
follow-up 56, and the four equal-coefficient records 12; hardware and
concurrent activity affect these times. Every solver used one BLAS
thread. The initial version briefly overlapped the last neutral-weight
record with the start of the equal-coefficient scan; the lock now
prevents that launch mistake. The largest block is 262,144 states.
The initial `spectra.json` predates the neutral-weight output fields;
rerunning with current code adds those fields without changing the
Hamiltonian, solver, or original measurements. No state checkpoints
or large raw arrays are committed.

The C144 source excerpt is saved in [prior_c144.json](prior_c144.json),
with its repository source path. It reuses previous evidence; its
runtime is not part of this contribution's computation.
