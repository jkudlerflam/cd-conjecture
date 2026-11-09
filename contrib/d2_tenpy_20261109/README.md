# Independent full-qubit DMRG of d=2 candidates

OpenAI Codex, GPT-6 family, 2026-11-09. Numerical investigation, not a
counterexample claim. Uses TeNPy 1.1.1, NumPy and SciPy in an isolated environment.

## Protocol recorded before running

Run the literal nearest-neighbour d=2 Pauli Hamiltonian, with uniform fields
and open ends, without symmetry projection. Prioritize C183's full-chain
valley-retuned V=16 bond and a simple rational approximation supplied by the
independent simplification agent. First validate the model and solver against
an independent exact small-chain diagonalization and the existing Davidson
implementation. Scan the field to avoid confusing the dressed vacuum with a
finite-density gas. Then compare L=32,48,64 (and larger if practical) with
increasing bond dimension and multiple initial states at promising points.

Report the full entropy and density profiles, final sweep changes, discarded
weights, energy, X/Y expectations and one-body correlations. A fit c>1 alone
is not a pass: require convergence with bond dimension/initialization and
growth with size calibrated by c=1 and c=2 controls. The two-leg XX control
has range two and is only a calibration, never the proposed microscopic model.
Several bulk cut windows are fixed here (1/8, 1/5, 1/4, 1/3 edge exclusion).
The main two-size estimator uses the average of the three central cuts.

For open chains, a boundary particle sees a different coordination number;
the chosen convention is exactly the uniform Pauli field in the user's
Hamiltonian. A scalar energy shift sets the all-up product energy to zero,
and both shifted and literal energies are saved.

The C183 antiunitary symmetry forbids uniform X but does not by itself forbid
uniform Y; both are measured. No projection onto fixed particle number is made.

## Reproduction

Install `physics-tenpy==1.1.1 numpy scipy`. For example:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python run.py --name c183 --L 32 --hz -7.9463 --chi 96 --sweeps 16 --out results/c183_L32.json
```

The simpler full-spin parity-family test uses `--name c181153`.
Run **both** `--parity even` and `--parity odd` at identical couplings
before selecting the lower energy; parity projection is rejected for
the transverse-DM family. See `RESULT.md` and the independent `REVIEW.md`
for convergence limits and the exact finite-size control values.
Its c≈2 reading at L=32→64 is not sustained at L=128: all parity
pairings give a 64→128 slope below one. The larger states remain
incompletely optimized. Neither c>1 nor asymptotic c=1 is established.

For expensive future runs, use `--save-state /private/tmp/my_state.pkl`.
Continuation with `--resume-state /private/tmp/my_state.pkl` verifies
the length, Hamiltonian and parity, preserves the existing bond dimension,
and can be combined with `--strict-local` for a longer reorthogonalized
local eigensolve. Checkpoints are local working files, not Git artifacts.
The save/resume path was checked at L=8 to an energy difference below
4e-14. Only load trusted checkpoints produced by this code.

`analyse.py` saves the finite-chain summary and scientific plots. An
infinite-chain branch and its method checks are described separately in
`INFINITE_VALIDATION.md` and `INFINITE_LIQUID_BRANCH.md`; failed attempts
are explicitly tagged and excluded from physical estimates.

Sources: [TeNPy model documentation](https://tenpy.readthedocs.io/en/stable/intro/model_details.html),
[Calabrese and Cardy](https://arxiv.org/abs/hep-th/0405152),
and repository `codex/work_C183/RESULT.md`, `codex/work_C186/RESULT.md`.
