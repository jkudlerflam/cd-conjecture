# D3 closed-book ED campaign

Read `REPORT.md` for the verdict, calibration numbers, partial phase catalogue, and unresolved candidates. `METHODS.md` gives the derivations and flags the continuum assumptions. `catalogue.csv` contains one row per retained finite-size calculation, including raw entropy fits, regression errors, cut dependence, symmetry sector, and velocity qualifications. `survey.png` and `survey.svg` summarize the survey.

The code requires Python 3, NumPy and SciPy. Matplotlib is only needed for the figure. The tested versions are NumPy 2.4.6, SciPy 1.17.1, and Matplotlib 3.10.9. This is CPU sparse exact diagonalization. No GPU library is used.

```sh
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
python3 d3.py --selftest --output selftest.json
python3 d3.py --family uls3 --sizes 6 9 12 --levels 8 --output new_uls3.jsonl
python3 d3.py --family sp4 --n 96 --sizes 6 8 10 --levels 4 --output new_sp4.jsonl
python3 d3.py --family su2 --n 192 --sizes 6 8 10 --levels 4 --output new_su2.jsonl
python3 d3.py --family orbital --sizes 6 8 10 --output new_orbital.jsonl
python3 d3.py --family spin1 --n 96 --sizes 6 9 12 --output new_spin1.jsonl
python3 d3.py --family sp4 --sizes 12 --sp4-cartan --no-cache-momenta --coeff 0.9978589232386035 0.06540312923014306 --label candidate_pi48 --output new_candidate.jsonl
```

The d3.py measurement outputs append. Most scan helpers do the same. The Sp(4) Casimir script regenerates its own output. Use a fresh output filename for a fresh d3.py run. `--coeff` specifies the operator coefficients described in `METHODS.md`. `--charge` overrides the additive charge sector and can exclude the global ground state. `--sp4-cartan` uses the two Sp(4) charges, whose sufficiency at larger sizes is qualified in the report. `--no-cache-momenta` saves memory by releasing each sparse momentum block after diagonalization.

For U(1)-only spin-1 points, `d3.py` by itself computes the requested magnetization sector, which need not contain the global ground state. `extra_scan.py u1`, `refinement.py u1`, and `near_saturation.py` audit all magnetizations. `global_velocity.py` additionally minimizes the adjacent-momentum excitation over all magnetizations. This distinction matters in the saved data.

Other scripts reproduce specific parts of the campaign.

* `extra_scan.py` runs the denser spin-orbital grid, U(1) grid, or initial-deformation charge audit.
* `refine_sp4.py` and `refinement.py` run the larger selected sizes and the SU(4) neighborhood.
* `sp4_followup.py fine` scans 768 angles at L=8. `sp4_followup.py crossing` locates the first finite-size 0/pi crossing at L=6,8,10,12.
* `sp4_towers.py` labels current-momentum eigenstates by the global Sp(4) Casimir. `sp4_tower_analysis.py` tests the proposed split-velocity description.
* `diagnostics.py` measures total SU(2) spin in selected states. `u1_occupations.py` measures the conserved nonzero-color occupation at beta=1.
* `fix_charge.py` creates the global-ground replacements for the initial U(1) deformations. `analyze.py` merges them. `invalid_spin1_L9.jsonl` is a quarantined, incorrect-sector run and is never loaded by analysis.
* `catalogue.py`, `plots.py`, and `crosschecks.py` regenerate the catalogue, figure, and checks from the saved outputs. `small_sizes.py` supplies L=4 points for commensurate energy fits.

The initial odd-L spin-1 sector bug and the empty-momentum-block bug were corrected. The corrected odd-L ULS calculation agrees in both energy and entropy with the independent color-sector implementation. These failures are recorded in the report rather than hidden.

## Typhon continuation

`cluster_manifest.json` and `typhon.sbatch` are prepared, not submitted. `cluster_job.py` processes one manifest entry per SLURM array task. It reads only the manifest, constructs the Hamiltonian, and writes `cluster_results/task_NNNN.jsonl`.

```sh
python3 make_cluster_manifest.py
python3 cluster_job.py --index 0 --dry-run
# After copying the scripts and manifest to a Typhon working directory:
sbatch --array=0-153%4 [allowed partition/account options] typhon.sbatch
```

Use the task count printed by `make_cluster_manifest.py` if the manifest has changed. The numerical payload is ready. Typhon's current partition/account policy was not inspected in this closed-book run, so those scheduler options remain environment-specific. The script requests one CPU and 128 GB per task, with a 24-hour limit, and never requests a GPU. Start with the L=14 Sp(4), L=12 spin-3/2, and L=16 spin-1 tasks before the most expensive entries. The supplied memory and time are requests, not measured requirements or guarantees for the largest sizes.

The manifest includes larger compatible calibration lengths, the small-angle Sp(4) candidates, high-entropy SU(2) points, the spin-orbital neighborhood, and full magnetization sweeps for selected U(1) points. For a U(1) sweep, compare all charge-task energies before choosing a ground-state entropy. Some large tasks omit entropy to keep the first spectral pass affordable. This local campaign did not submit any cluster job.
