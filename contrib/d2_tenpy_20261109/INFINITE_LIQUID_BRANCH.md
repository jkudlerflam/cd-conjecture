# C183 infinite-chain liquid-branch search, 2026-11-09

This extends `INFINITE_VALIDATION.md`. It finds a substantially better infinite
variational branch than the low-density branch in the first audit. It does not
establish a central charge or a counterexample.

## Initialization that preserves the dilute gas

An uncorrelated product state at density 0.2 pays an order-V energy for neighboring
particles. It repeatedly evolves toward the low-density branch. Simply closing a
central block of a finite-chain MPS also fails to preserve its physical state:
the two cut virtual bases need not match. The tested central-block initialization
had energy density +0.0384 after canonicalization.

Instead initialize the exactly excluded coherent gas with bond dimension two:

```
A[0] = [[1,0],[1,0]]
A[1] = [[0,z exp(i phi)],[0,0]]
z = sqrt(((1-2 rho)^(-2)-1)/4),  phi = -pi/2.
```

The matrix identity `A[1] A[1] = 0` forbids neighboring occupied sites. The transfer
matrix gives density `rho`; this was independently checked at rho=0.12, 0.2, 0.3,
with adjacent-occupation expectations below 2e-16 in magnitude. This initial state
is not an eigenstate of the full spin Hamiltonian. At rho=0.2 its literal C183
energy density is -0.02268102, so optimization can build the virtual dressing
while retaining a dense gas. The implementation is `hardcore_coherent_state` in
`independent_infinite_check.py`, also exposed by `infinite.py --hardcore`.

## Completed chi=32 experiment

The exact C183 bond parameters from `run.py` are used without modification:
Jx=0.4979983, Jy=-0.5020017, Jz=4, Dz=0.0152122, Dy=0.13693064,
hz=-7.9463. The unit cell contains eight sites. There is no imposed number,
parity, or spin conservation. This initial run uses no mixer, chi=32, 80 optimization
sweeps, 20 environment sweeps per two optimization sweeps, and Nmax=80 local
Lanczos vectors with full reorthogonalization.

| Quantity | Initial completed liquid run |
|---|---:|
| Direct shifted energy per site | -0.0627222840416 |
| Returned minus direct energy | -3.14e-13 |
| Mean density | 0.18877569 |
| Density range over unit cell | 0.18806051--0.18949466 |
| Mean bond entropy | 1.01963201 |
| Transfer correlation length, sites | 59.57854 |
| Last five sweep-energy span | 4.57e-9 |
| Final sweep norm error | approximately 2e-5 |
| Runtime, seconds | 154.7 |

The mean density and energy agree qualitatively with the large finite-chain
liquid branch. The root finite-chain energies at L=64 and 128 give a bulk energy
increment of about -0.06272676, several micro-units below this finite-chi
variational energy. The energy found here is also well below the earlier
low-density branch (-0.06252875).

These are properties of an approximate variational state, not yet of a validated
global ground state. Small energy drift does not by itself bound observable errors.
In particular, later local-solver checks expose residuals around 7e-5 at Nmax=80
and around 1e-5 at Nmax=160, with the local iteration cap reached. More accurate
solver and bond-dimension checks are required before using entropy/correlation
length for a central-charge estimate.

## Restarting the state

The first product-like state can use finite-range contracted initial environments.
For a saved long-correlated MPS, short contracted environments produce a large
transient and can destroy the branch. The transfer-matrix initializer preserves it:

```python
resume_data = {'init_env_data': {'force_init_method': 'TM'}}
```

This differs from the defective default iterative/GMRES initializer discussed in
the audit. The diagnostic exposes the choice as `--env-init TM`. The original
good state is preserved in the temporary file
`/private/tmp/cd_c183_i32_hardcore.pkl`; it is not tracked in Git. New runs checkpoint
their MPS every ten sweeps. Loaded checkpoints are canonicalized before initializing
environments.

## Reproduction

```sh
/private/tmp/cd-conjecture-venv/bin/python contrib/d2_tenpy_20261109/independent_infinite_check.py --name c183 --chi 32 --cell 8 --sweeps 80 --diag default --update-env 20 --density .2 --hardcore --save-state /private/tmp/cd_c183_i32_hardcore.pkl --out /private/tmp/liquid32.json
```

The corresponding saved repository result is
`results/infinite_c183_c8_chi32_hardcore.json`. Product-start failures, central-block
warm-start failures, and short continuation diagnostics are saved separately and
must not be pooled into a central-charge fit.

## Capped chi=64 continuation: provisional c approximately one

The preserved chi=32 state was continued to chi=64 using transfer-matrix
environments, Nmax=160, and the same eight-site cell. The run was intentionally
stopped after preserving the sweep-10 checkpoint, so computational effort could
move to the newly promising parity-preserving microscopic family. It was not
stopped by a numerical convergence test.

| Quantity | chi=32 initial liquid | chi=64, continuation sweep 10 |
|---|---:|---:|
| Direct shifted energy/site | -0.062722284042 | -0.062724543941 |
| Mean density | 0.18877569 | 0.18872430 |
| Mean bond entropy | 1.01963201 | 1.16265513 |
| Correlation length, sites | 59.57854 | 135.29003 |
| Density span over cell | 0.00143415 | 0.00037362 |

The two-point expression `6 (S64-S32) / log(xi64/xi32)` is **1.04635**. This is a
provisional consistency check for the finite-chain c approximately one
interpretation. It is not a converged finite-entanglement estimate, and its small
excess above one is not evidence for a counterexample.

The chi=64 final local-solver residual is 1.48e-5, the sweep canonical norm error
is 7.94e-5, and the mean energy still changes by 1.46e-9 between sweeps 8 and 10.
The entropy span over the cell is 0.00153075. There is only one higher-chi
checkpoint; entropy/correlation-length stability under further sweeps and chi
increases has not been established.

The state is preserved in `/private/tmp/cd_c183_i64_hardcore.pkl`; its independently
canonicalized measurements and progress/residual checks are saved in
`results/infinite_c183_c8_chi64_checkpoint10.json`. The formula and limitations
are also saved in `results/infinite_c183_provisional_fes.json`. The script
`inspect_infinite_checkpoint.py` reproduces the checkpoint measurements without
changing the stored state.

The original `infinite.py` now labels every run experimental, checks direct local
energy against the returned energy, records actual local-solver residuals, and
keeps a separate `ground_state_or_phase_validated=False` flag even if its numerical
stability thresholds happen to pass. It supports the exact-exclusion initial state
through `--hardcore`; short XX and C183 smoke runs confirmed the updated interface.

## Independent local-eigensolver audit

A six-sweep chi=32 continuation with ARPACK and transfer-matrix environments
finished in 373.8 seconds. Its maximum local eigenvector residual was 4.92e-10,
compared with order 1e-5--1e-4 for the capped Lanczos solves. After canonicalization:

```
direct energy/site = -0.062722248850
mean density      =  0.1887713077
mean entropy      =  1.0195131299
correlation length= 59.55275588 sites
```

Using this baseline instead changes the two-point chi=32-to-64 expression to
**1.04667**. Thus improving the local solver does not reveal a hidden c approximately
two signal in this diagnostic. However, its sweep norm drift grew to 0.0028 even
though the returned/direct energy difference is only 7.9e-10. The full infinite
optimization is not demonstrably converged. The result is saved as
`results/infinite_c183_c8_chi32_arpack_checked.json`; it must be read with this
limitation. More environment convergence and further bond dimensions would be
needed for an asymptotic statement.

No further transverse-DM infinite jobs remain running. This point now has useful
provisional evidence consistent with c approximately one; the newly promising
parity family is the higher-priority c>1 target.
