# Continuing when suitable compute and checkpoints become available

No cluster job was launched by this contribution. The paths below come from
historical records and must be checked on the actual compute host. Do not
copy 30-100 GB checkpoints to the current 48 GiB Mac to try production.
Run from a clone of this repository on a suitable existing allocation.

## 1. Authenticate code, capacity and source files

Install `requirements.txt` into an isolated Python environment. CUDA mode
also needs a CuPy build matching the compute host's CUDA runtime. No specific
Slurm account, partition, host, or old allocation is assumed available.

Check every production/measurement source against `source_manifest.json`.
Reserve approximately 256 GB host memory and one 80 GB GPU for a single
sym40 run. Keep about 300 GB scratch free for source/output states and atomic
checkpoint replacement. Do not run two sectors concurrently within that
budget. The benchmark's memory flag is a **declared available budget**,
not a request for memory from a scheduler or a detector of competing jobs.

Known historical candidates for source checkpoints:

| Stage | Committed source run JSON | Historical checkpoint |
|---|---|---|
| L512 odd 2048 | `codex/work_C195/typhon_out_x2048/sym40_n316_L512_P1_rand1.json` | `/scratch/lustre/jkudlerflam/c19x/c195x2048/out2048/ckpt/sym40_n316_L512_P1_rand1.npz` |
| L512 even 1536 | `codex/work_C195/typhon_out_b/out512/sym40_n316_L512_P0_rand1.json` | `/scratch/lustre/jkudlerflam/c19x/c195b/out512/ckpt/sym40_n316_L512_P0_rand1.npz` |
| L256 even 1024 | `codex/work_C190/gpu/typhon_out/outX1024/sym40_n316_L256_P0_rand1.json` | `/scratch/lustre/jkudlerflam/c18x/c190g/outX1024/ckpt/sym40_n316_L256_P0_rand1.npz` |

`continuation.py` checks tensor count, parity blocks, sweep metadata, bond
dimensions, finiteness, energy offset and freshly recomputed energy/entropy/
density. It hashes the source/checkpoint and refuses mismatches or overwrite
of a configuration. The full-state CPU checks can take substantial time;
their cost is not the tiny local validation runtime. Missing files halt the
preparation. Source completion is not assumed to prove convergence.

## 2. Benchmark before committing to production

The local smoke command, already tested:

```sh
python3 contrib/qubit_convergence_20261116/benchmark.py \
  --numpy --L 32 --chi 64 --memory-budget-gb 2 \
  --out /tmp/qubit-benchmark-local.json
```

On an actually allocated GPU host, in an exclusive output directory:

```sh
python3 contrib/qubit_convergence_20261116/benchmark.py \
  --L 512 --chi 2560 --memory-budget-gb 256 \
  --out /tmp/qubit-benchmark-2560.json
```

This benchmarks the **real sym40** tensors and C195 update routines, not
the complex p153 bond. It deliberately uses a fixed 100 matvecs, so its
residuals are diagnostics, not convergence evidence. It allocates the full
random MPS and host environments and reports build, solve, split and entropy
costs. GPU mode is supplied but could not be exercised on this Mac.
If the real run fails memory or timing gates, record that failure; do not
relax convergence criteria or silently change chi.

## 3. Prepare and run the odd 2560 continuation

After verifying the checkpoint is present, substitute its actual path if it
has moved. The following creates only a new configuration:

```sh
python3 contrib/qubit_convergence_20261116/continuation.py \
  --source codex/work_C195/typhon_out_x2048/sym40_n316_L512_P1_rand1.json \
  --checkpoint /scratch/lustre/jkudlerflam/c19x/c195x2048/out2048/ckpt/sym40_n316_L512_P1_rand1.npz \
  --chi 2560 --memory-budget-gb 256 \
  --out contrib/qubit_convergence_20261116/runtime/campaign/cfg/odd2560.json

C190G_BACKEND=cupy C195_OFFLOAD=1 C190G_GPU_MIN_CHI=128 C190G_SVD_DEV=1 \
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 \
python3 codex/work_C195/run195.py \
  contrib/qubit_convergence_20261116/runtime/campaign/cfg/odd2560.json \
  contrib/qubit_convergence_20261116/runtime/campaign/odd2560
```

The new tag is `sym40_n316_L512_P1_rand1_nextchi2560`. Source tensors are
read only; output tensors are written into the new directory. Repeating the
driver with the **same configuration and output directory** resumes after
an interrupted sweep and skips a finished run. A `done` output can mean the
sweep cap was reached; inspect the actual stop and convergence checks.
Never reuse an output directory after modifying a configuration.

For even, first call `continuation.py` with its L512 even source above and
`--chi 2048`; run it in `runtime/campaign/even2048`. Then use that new run JSON
and its new `ckpt/*_nextchi2048.npz` as the source for `--chi 2560` in a separate
configuration/output directory. Preparation rejects a direct 1536->2560 jump.
Use L256 even's source with `--chi 1536` to match the small-size reference's
parity accuracy. These three branches were configuration/resume tested at
small L through the same preparation function, not run at production size.

## 4. Controls and measurements

`control_configs/xx{1,2}_L{256,512}.json` are explicit matched C195 configs.
Run each with the same driver into separate fresh output directories.
They preserve free-level midpoint filling, both entropy windows and the
registered chi/stop schedule. Their exact finite-size expectations are in
`controls.json:planned_matched_controls`; compare to those rather than forcing
the finite-size result to the limiting integer central charge.

For full correlations, C193's interface is:

```sh
python3 codex/work_C193/ckpt_meas/meas193.py run \
  RUN_JSON CHECKPOINT_NPZ NEW_MEASUREMENT_JSON --numpy
```

Those three capitalized arguments are placeholders for newly verified files.
The original measurement script does not reject every recorded discrepancy;
explicitly apply the registered energy, entropy, density and metadata checks
to its `check` output before fitting. At chi=2560 the unchanged GPU route
cannot fit the complex MPS on an 80 GB GPU. CPU measurement also needs a
separate stack-memory/time budget. Do not start it merely because DMRG fit.

Run `ana193.py states` from its own directory (so it finds `frozen.json`) on
the resulting full measurement files. Apply the registered 128/256/512
drift diagnostics with both K_s estimators and matched controls. Keep the
boundary-corrected fit secondary. The short `gwin=48` production measurements
do not replace these full correlators.

## 5. Evaluation and preservation

Retain all sweeps and unrounded values. Apply every numbered gate in
`REGISTRATION_20261116.md`, including both sectors, stationary previous
stages, C4, numerical energy-order margin, density and exact controls.
Use canonical final entropies for c_two and last-sweep entropies for C4;
publish both. The current `assessment.py` deliberately reports only the
**present** fallback assessment; it must not be used to label a future run a
pass without implementing/reviewing all the new-run gates.

Keep checkpoints out of Git. Commit the small raw outputs, source hashes,
resource measurements, pass/fail/undecided table and any dated amendment.
Do not revise the frozen registration after seeing new outputs. New evidence
may change the numerical verdict; a CFT interpretation still requires its
separate infrared assumptions.
