# Independent infinite-chain diagnostic, 2026-11-09

**No infinite-chain central-charge or ground-state correlation-length result is established here.**
The original `infinite_c183_cell4_chi32.json` is rejected. Its shifted energy density
`+1.3929312281` exceeds the all-up trial value zero; late sweep energies fluctuate
by order one and sweep norm errors are approximately 0.4--0.75. The tiny final
`norm_error` follows post-run canonicalization and does not validate optimization.
That JSON now explicitly carries the status `failed_unconverged_do_not_use_for_physics`.

## Reproducible method issue

The literal bonds in `run.py` are Hermitian. This audit did not identify a bond,
Pauli-normalization, scalar-offset, or infinite unit-cell convention error.

In the installed TeNPy 1.1.1, `MPOEnvironment.init_first_LP_last_RP` defaults to
`force_init_method='iter'`. For the random complex product initialization used here,
the iterative/GMRES initializer reports residuals of order one against a requested
`1e-11`. The resulting **actual effective Hamiltonian** can violate Hermiticity:
in the C183 initial 4-dimensional problem,

```
H_eff - H_eff^dagger = (0.5498674953929421 i) I + residual,
||residual||_F = 4.01e-16.
```

Thus the observed defect is an imaginary scalar environment contribution, rather
than a missing conjugate microscopic bond. Such a scalar leaves exact mathematical
eigenvectors unchanged, but is outside the Hermitian input contract of the local
solvers. This check alone does not prove it caused every observed failure.

The documented finite-range initialization bypass is

```python
resume_data = {'init_env_data': {'start_env_sites': cell}}
engine = dmrg.TwoSiteDMRGEngine(psi, model, options, resume_data=resume_data)
```

The independent diagnostic uses this route by default. Initial effective-matrix
Hermiticity defects then fall to approximately `5e-16`. It does **not** alter the
model. No installed TeNPy source or original root script was changed.

## Checks and limits

| Run, cell 4 | Direct local energy/site | Additional result | Interpretation |
|---|---:|---|---|
| XX, chi 16, ARPACK, iterative init | -0.636586792877 | Initial relative Hermiticity defect up to 0.00429 | Finite-chi control survives this initialization defect |
| XX, chi 16, ARPACK, contracted init | -0.636586792877 | Defect <4.9e-16; exact energy is -2/pi = -0.636619772368 | Correct control energy to 3.3e-5; not a precision finite-entanglement study |
| C183, chi 16, reorthogonalized Lanczos, contracted init | -0.060089875195 | Late sweep norm errors 0.2--0.6; returned/local energy mismatch 6.4e-4 | Failed optimization; discard phase diagnostics |
| C183, chi 16, contracted init, 20 environment sweeps per two optimization sweeps | -0.062528750972 | Density approximately 0.01104; returned/local difference 2.5e-14; final sweep norm 1.9e-11 | Plausible low-density competing variational state; energy not yet converged |

The last run's final several energy densities still differ by approximately
`6e-6`. Its entropy is approximately 0.064 and its reported correlation length
approximately 4.16, but these are **not** measurements of the global ground-state
phase. A dilute, weakly entangled branch can be metastable next to the denser
state found in finite-chain calculations. Convergence and initialization dependence
must be resolved before assigning physical significance.

The XX control took 2.9--3.1 seconds; the two completed C183 contracted runs took
22.1 and 24.4 seconds. C183 ARPACK trials were stopped after a bounded diagnostic
attempt and produced no accepted data. Their absence must not be treated as a
physical result.

## Reproduction

The script is `independent_infinite_check.py`; results are the
`results/independent_infinite_*.json` files. For example:

```sh
/private/tmp/cd-conjecture-venv/bin/python contrib/d2_tenpy_20261109/independent_infinite_check.py --name xx --chi 16 --cell 4 --sweeps 24 --out /private/tmp/xx_check.json
/private/tmp/cd-conjecture-venv/bin/python contrib/d2_tenpy_20261109/independent_infinite_check.py --name c183 --chi 16 --cell 4 --sweeps 40 --diag default --update-env 20 --out /private/tmp/c183_check.json
```

The script caps numerical threads at one and measures the energy by explicitly
contracting each literal local term in the final canonicalized MPS. The returned
iDMRG energy density is checked separately. `--env-init iter` reproduces the
initialization problem; the small first-matrix C183 audit is preserved as
`results/independent_infinite_initialization_audit.json`.

For future infinite calculations, require local/returned energy agreement, stable
energy and observables over further sweeps, independently controlled environments,
multiple initial densities/phases, and increasing bond dimensions. A final canonical
norm alone is insufficient. Until these checks pass, use the independently validated
finite-chain results for the scientific status report.
