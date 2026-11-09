# Full microscopic transverse-DM chain: DMRG checkpoint

AI assistance: GPT-6 Codex, 2026-11-09. Standard: numerical physics. No counterexample is established by this contribution.

## Model and independent validation

We use C183's full-chain-retuned V/Delta=16 bond, with Jx=0.4979983, Jy=-0.5020017, Dz=0.0152122, Jz=4 and Gamma_xz=-Gamma_zx=sqrt(0.3)/4. The longitudinal field is scanned. The actual two-state Pauli chain is used with uniform field on every open-chain site; no effective-space projection or particle-number constraint is imposed.

The numerical solver is C183's two-site Davidson DMRG unchanged. `full_chain.py` adds physical X and Y profiles to the existing Z, entropy and energy outputs. `source_manifest.json` records the source files and hashes. All runs use two BLAS threads.

[numerical] At L=12, hz=-7.9463, independently assembled sparse Pauli ED and DMRG agree in energy to 8.3e-12, all cut entropies to 2.6e-10, and all Y expectations to 4.6e-10. The open-chain entropy is Sbar=0.105790 and central raw density is 0.062451. The central Y response is -0.323598, while X is zero within numerical error. This Y response is not merely an open-boundary effect: periodic ED at the same parameters has uniform Y=-0.400024, X approximately zero, and N=2.99039 on twelve sites. The periodic gap is 0.00143059. A second periodic ED at L=14 gives Y=-0.408604, N=2.99030 and gap 0.00191003, so substantial Y persists with a different allowed-momentum grid. These are finite-system statements and do not establish the relevance of any induced bulk source.

Files: `results/DM16_h9463_12.json`, `results/ED_DM16_h9463_12.json`, `results/ED_PBC_DM16_h9463_12.json`.

## Field ladder

[numerical] Complete L=24 sweeps, chi=16,32,48,48,64,64:

| hz | raw central density | Sbar, three central cuts | Y central | broad profile slope |
|---:|---:|---:|---:|---:|
| -7.9500 | 0.01248 | 0.04663 | -0.05995 | -0.199 |
| -7.9480 | 0.03156 | 0.06346 | -0.13393 | -0.378 |
| -7.9470 | 0.07734 | 0.14690 | -0.22240 | 0.029 |
| -7.9463 | 0.11717 | 0.23593 | -0.28308 | 0.519 |
| -7.9455 | 0.15286 | 0.29069 | -0.32129 | 0.819 |
| -7.9450 | 0.17276 | 0.32399 | -0.33383 | 1.016 |
| -7.9440 | 0.20579 | 0.39942 | -0.34472 | 1.422 |
| -7.9420 | 0.24574 | 0.44786 | -0.33507 | 1.427 |
| -7.9400 | 0.27193 | 0.47030 | -0.31271 | 1.261 |

All final local residuals are below 1e-9, last-sweep energy movements below 2.5e-12, and discarded weights below 1.3e-12. Density and entropy profiles, rather than just their central values, are saved in `results/lad_*.json`.

The single-profile estimates are sensitive to the fitted interval. For example, at hz=-7.9440 the fitted slopes are 1.42, 2.12, 1.54 and 2.00 when excluding respectively the outer eighth, quarter, third and three eighths of the chain. This is not a stable estimate of c. Nor does the small absolute entropy by itself exclude criticality; an additive nonuniversal constant matters.

The density in this table is the raw spin-down density. It includes virtual dressing and does not independently identify the density of a gas of dressed quasiparticles. The uniform Y response is permitted by the surviving antiunitary symmetry. It does not prove that the chain is gapped.

## Controls

[numerical, checked against exact free-fermion correlation matrices] The same Davidson solver and central-three-cut estimator distinguish XX (c=1) from two independent XX legs (c=2) on the chosen 24-to48 size pair:

| control | measured size-pair slope | exact finite-size slope | maximum entropy error |
|---|---:|---:|---:|
| XX | 1.0637380 | 1.0637380 | 4.9e-12 |
| two XX legs | 2.1199804 | 2.1203001 | 3.8e-5 |

The largest energy error is 1.2e-6 for the two-leg L=48 control, whose discarded weight is 2.8e-8. The exact low-density 48-to64 control slopes are instead 1.756 and 4.239 because of particle-number crossings, so those short adjacent sizes would be a poor uncalibrated test of c.

## Matched size tests

Selected before their output: hz=-7.9450, which has raw central density near 0.17 at L=24, and hz=-7.9440, the first ladder point with the broad profile estimate above 1.2. Both are being continued to L=48 with chi=16,32,48,64,64,96,96. At hz=-7.9450, L=48 has Sbar=0.4636847 and the matched 24-to48 slope is 1.20918. Its central density rises from 0.17276 to 0.23968 and Y=-0.42421. Fitted-window slopes are 1.353, 1.302, 0.257, -0.535, so this is not a stable c>1 estimate. The final local residual is 2.14e-9, discarded weight 8.94e-12 and last-sweep energy movement 1.1e-10; the run took 358 seconds. The other field remains pending. No positive conclusion is drawn from the L=24 profile readings alone.

## Reproduction

Use Python with numpy and scipy, from the repository root:

```bash
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 python3 contrib/d2_dmrg_priority_20261109/run_campaign.py --stage ladder
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 python3 contrib/d2_dmrg_priority_20261109/controls.py --lengths 24 48 --chis 16 32 64 96 96
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 python3 contrib/d2_dmrg_priority_20261109/run_campaign.py --stage sizes
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 python3 contrib/d2_dmrg_priority_20261109/analyse.py
```

`check_ed.py` validates a saved DMRG result using an independent sparse Pauli assembly; `--periodic` instead computes the same bond on a ring. `summary.json` collects the completed results and fitted-window sensitivity.
