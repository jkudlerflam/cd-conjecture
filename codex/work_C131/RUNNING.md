# Running the C131 numerical search

The code uses Python 3, NumPy, and SciPy. It writes only inside the working directory. No cluster connection or submission was made in this run.

Small local checks:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 checks.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 pilot.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 refine_pilot.py
```

The first pilot uses 12 starts and a cap of 1000 objective evaluations per start, separately at L=2,3,4. The refined pilot reruns the same starts with a cap of 4000. The full search is intended for Typhon, from a directory containing optimize_q.py and submit_c131.slurm and an environment providing NumPy/SciPy:

```sh
sbatch submit_c131.slurm
python3 summarize.py results
```

The SLURM array has 224 tasks, 32 starts for each L in 2,4,6,8,3,5,7, ordered in that sequence. It runs at most 8 tasks concurrently. Each task has 4000 objective evaluations available, 12 GB memory, one CPU, and a 24-hour limit. It targets the apollo partition named in the workspace instructions. It uses CPU sparse diagonalization, so it does not request a GPU. Cluster account/QoS requirements and Python environment availability cannot be tested from this sandbox.

One individual task can be reproduced with:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 optimize_q.py --L 8 --seed 0 --budget 4000 --out L8_seed0.json
```

Seeds 0 through 4 are Heisenberg, XX, an even-ring staggered-rotation image of Heisenberg, XXZ near Delta=-1, and a critical anisotropic XY point. Among the remaining seeds, multiples of three perturb Heisenberg and the others use generic nine-coefficient directions. A Householder rotation moves each starting direction away from the coordinate singularities of the eight-angle parametrization. The search covers scale directions through normalization; no sign quotient is taken.

Every final JSON includes the full nine coefficients normalized to the numerical M, both energies, both eigenvector residuals, a primal/dual numerical bracket for M, termination status, seed, and wall time. A budget-limited run is recorded as such. The summary never identifies a local-search value with a certified global supremum.

For larger-ring finalists, rerun from additional independent eigensolver starts and enlarge ncv before interpreting small differences as drift. A residual checks proximity to an eigenvalue; it is not a proof that an eigenpair is the lowest. This code uses a fixed generic starting vector for reproducibility and checks the full spin space, with no imposed symmetry sector.
