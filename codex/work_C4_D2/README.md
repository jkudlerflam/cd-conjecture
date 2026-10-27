# D2, nearest-neighbor qubit-chain search

Run `python3 d2.py --help`. This run used Python 3.12, NumPy 2.4.6,
SciPy 1.17.1, and Numba 0.65.1.
Matplotlib is optional for the summary figure. No external data or network
access is needed. Use the working directory as the output directory.

The scientific result is in `report.md`, with derivations and estimator
limitations in `methods.md`. `calibration.md` and `candidates.md` contain the
numerical tables. JSON files in `results/` retain parameters, spectra,
entropies, residuals, and individual Fourier-density overlap clusters.

## Reproduce

```sh
python3 d2.py selftest
python3 d2.py calibrate --sizes 8 10 12 14 16 18 20 22 --output results/calibration_fixed.json
python3 d2.py points --random 1200 --output results/global_points.json
python3 d2.py scan --points results/global_points.json --output results/global_scan.json
python3 make_boundary_points.py
python3 d2.py scan --points results/boundary_points.json --sizes 12 14 --output results/boundary_scan.json
python3 d2.py optimize --input results/scan.json --maxiter 40 --seeds 5
python3 d2.py refine --input results/global_scan.json --ids 1256 --sizes 10 12 14 16 18 20 22 --output results/example1256.json
python3 exact_controls.py
python3 free_fermion_control.py
python3 density_check.py
python3 summarize.py
```

The initial `results/scan.json` is the 238-point pilot, generated with
`python3 d2.py points --random 100`, then `python3 d2.py scan`.
It is a subset of the 1,338-point global scan and must not be counted twice.
The first optimizer uses five family-distinct seeds from that pilot.

Screening at L=14,16 uses `results/shortlist_points.json` (49 pilot points)
and `results/global_shortlist_points.json` (163 further points). Their
selection thresholds and the second optimizer are stated in the report.
`screen_velocity.py` adds the finite-momentum velocity-drift diagnostic.
These screening stages can be reproduced from the saved point files.

`refine` searches all momenta by default. An explicit `--momentum` instead
assumes that ground momentum. `--af-momentum` assumes m=0 for L divisible by
4 and m=L/2 otherwise, and is intended only for the even-length, zero-field
AF-XY benchmark convention. `--half-filling` restricts to N=L/2 and is checked
for U(1) conservation. It is not a generic-chain option. `--field-split 1`
uses the literal left-site field density for the Fourier estimator; the
default is the symmetric density, `--field-split .5`. Both produce the same H.

Scans and refinements resume from their existing JSON outputs. Calibration
also resumes. Optimization traces start a new run and overwrite their output,
so use a new filename to retain earlier trials. Output writes are atomic in
the final code. Each new CLI run writes its PID to `<output>.pid`.

## Cluster

Transfer this directory to Typhon, install the dependencies in a suitable
Python environment, and submit from the directory with, for example,

```sh
bash submit_cluster.sh --partition=YOUR_AUTHORIZED_PARTITION
```

The wrapper splits the 1,338 tasks into batches of at most 1,000 array indices.
There is one task per saved global point, requesting one CPU
thread and 24 GB per task. It runs coarse ED at L=14,16,18 and all three
diagnostics at L=14,16,18,20,22 for every point. This is CPU code; it does not
request or use a GPU. No partition, Python module, network access, or cluster
availability was assumed. The script was prepared but not submitted or tested
on Typhon. Adjust the array range when changing the point file; use
`D2_POINTS=...` to specify that file.

For the boundary file use `D2_POINTS=results/boundary_points.json` and override
the wrapper automatically adjusts the array range. In particular, point b43 has only local
results through L=18. Its slower, nearly degenerate L=20 run was interrupted
to reserve the local budget for stronger candidates; the array can finish it.

For more parameter samples, generate a larger file with
`python3 d2.py points --random 10000 --output results/cluster_points.json`,
then set `D2_POINTS=results/cluster_points.json` when running the wrapper.
Adjust `D2_BATCH_SIZE` and `D2_CONCURRENCY` to the site's array limits.

## Canonical and superseded outputs

Use `results/calibration_fixed.json`, not `results/calibration.json`.
The latter records the first arbitrary-start eigensolver attempt, which
missed part of an exactly degenerate XX overlap at L=8. The corrected solver
starts from h_n|0> and sums the weight within each degenerate eigenspace.

Use `results/ferromagnet_half.json` for the near-ferromagnetic XXZ controls.
`results/ferromagnet.json` was an interrupted trial with an incorrect fixed
ground momentum at L=2 mod 4 and is invalid at those lengths.
`results/ferromagnet_correct.json` is an interrupted, slower, unrestricted
trial through L=20 with the correct momentum convention. Its valid rows agree
with the final charge-resolved calculation.

The closed-book instruction was respected. No pre-existing project files were
read or listed, and no project memory outside this directory was accessed.
`STATE.md` is a local handoff, not a modification of the parent project's memory.
