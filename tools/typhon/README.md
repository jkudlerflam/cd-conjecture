# IAS Typhon CPU jobs

Run commands from the repository root on the Mac. The driver uses Python 3.11 or newer, Git, OpenSSH and rsync. Numerical packages live on IAS.

## Access and Linux environment

Both SSH commands below authenticated as `pace` on 2026-11-24. The login alias uses `ias-gateway` as its SSH jump host. Both hosts accepted a key offered by the Mac's local agent. The workflow disables agent forwarding.

```bash
ssh -o BatchMode=yes -o ForwardAgent=no ias-gateway 'id -un && hostname -f'
ssh -o BatchMode=yes -o ForwardAgent=no ias-cluster 'id -un && hostname -f'
python3 tools/typhon/typhon.py setup
```

The gateway is `ssh1.sns.ias.edu` and the login node is `typhon-login1.sns.ias.edu`. Setup is repeatable. It loads the IAS Linux module `anaconda3/2025.06` and creates `/home/pace/cd-conjecture/envs/numerics-v1` using that module's site packages. It installs Linux wheels for TeNPy 1.1.1 and CVXPY 1.9.2 in the project environment and checks that CVXPY detects Clarabel and SCS. Every batch job sources the module initialization script, loads that module, activates the environment and sets numerical thread counts from the CPU allocation.

The validation used Python 3.13.5, NumPy 2.1.3 and SciPy 1.16.0. The shared home directory was writable from both the login and compute nodes. The login node also reported write access to `/mnt/lustre/flash`. This workflow stores its environment and runs in shared home. [Batch environment record](../../runs/typhon/heisenberg4-20261124/logs/runtime.json)

## Submit, inspect and retrieve

Choose a fresh run name for each submission. This example launches one exact diagonalization.

```bash
run="heisenberg4-$(date -u +%Y%m%dT%H%M%SZ)"
python3 tools/typhon/typhon.py submit "$run" tools/typhon/heisenberg4.job.json
python3 tools/typhon/typhon.py status "$run"
python3 tools/typhon/typhon.py logs "$run"
```

`submit` reports the Slurm ID and saves it as `job-id` on both hosts. `status` queries the live queue and the allocation's accounting record. `logs` shows the last 60 lines of standard output and standard error. If a connection drops during submission, inspect the same run name. Its exclusive remote directory prevents duplicate submission under that name. An active job can also be located with `ssh ias-cluster 'squeue --me'`.

Once the accounting row reads `COMPLETED` with `ExitCode=0:0`, fetch it and check the example. To retrieve the already verified job, set `run=heisenberg4-20261124` and start here.

```bash
python3 tools/typhon/typhon.py fetch "$run"
python3 tools/typhon/check_heisenberg4.py "runs/typhon/$run"
```

Fetching also works for a failed or active job. A live fetch records the accounting state sampled before the copy. Fetch again after completion for final outputs. A fresh checkout can recover the run by name.

| Location | Contents |
| --- | --- |
| `ias-cluster:/home/pace/cd-conjecture/runs/<run>/src/` | Working copy with repository-relative source and input paths |
| Same remote run directory | Original `snapshot.tar.gz` and SHA-256, manifest, job ID, `logs/` and `results/` |
| Local `runs/typhon/<run>/` | Snapshot, manifest, job ID, accounting, compute-node versions, stdout, stderr and results |

New downloads are ignored by Git. The small four-site validation directory is tracked in full.

## Run a different calculation

Copy the example job recipe and edit these fields. [Example recipe](heisenberg4.job.json)

| Field | What to supply |
| --- | --- |
| `files` | Enumerate the solver, its project-local imports and each input using repository-relative paths. The recipe and batch helpers are included automatically. |
| `command` | Give the executable and each argument as separate strings. The command runs in `src/`. Direct outputs to `../results/`. A staged shell wrapper can use the absolute `$TYPHON_RESULTS` when a solver has a fixed output path. |
| `seed` | Record the actual seed, or the seeds for each sector. Set the solver's seed through its own CLI, input or existing deterministic formula. Use JSON `null` for a calculation without randomness. |
| `resources` | Optionally override `cpus-per-task`, `mem`, `time` or `qos`. Defaults are one CPU, `2G`, `00:05:00` and `short`. Jobs use one node and one task on partition `typhon` with account `physics`. |

The upload reads the listed files from the working tree, including edits and new inputs. It accepts individual files and excludes environment directories, symlinks, compiled packages and prior runs. The original archive preserves the submitted bytes even if a solver writes beside its source. The manifest records the Git revision, status of selected paths, command, seed and Slurm resource request. The input files contain the model parameters. The driver checks the archive's SHA-256 on upload and fetch. Place large restart files outside `results/` to keep them on IAS.

## Verified four-site spectrum

On 2026-11-24 the project's sparse eigensolver diagonalized both spin-parity sectors of the periodic spin-½ Heisenberg ring at `L=4`, with `Jx=Jy=Jz=1/4` in the Pauli convention and zero field. Its starting-vector seeds were `20265098` and `20265099`, from `20261110 + 997*L + parity`. Slurm job **17755440** completed on **typhon-node62** in **17 seconds** with exit **0:0**, one CPU and 2 GB. [Solver](../../contrib/d2_order_one_20261110/ed_screen.py) · [Input](heisenberg4.json) · [Slurm accounting](../../runs/typhon/heisenberg4-20261124/accounting.psv)

For four periodic sites the Hamiltonian is `H = (S₀ + S₂) · (S₁ + S₃)`. Coupling the two sublattice triplets to total spin 0 or 1 gives energies −2 and −1. A sublattice singlet gives energy zero. The ground energy and first excitation gap therefore have independent exact references.

| Observable | Typhon | Exact | Absolute error |
| --- | ---: | ---: | ---: |
| Ground energy | −1.9999999999999998 | −2 | 2.22 × 10⁻¹⁶ |
| First gap | 0.9999999999999998 | 1 | 2.22 × 10⁻¹⁶ |

The checker uses absolute tolerance `1e-11`, ten times the requested solver tolerance for these order-one energies. The largest eigenpair residual was `5.74e-16`. The checker also matches the job and node to Slurm accounting and the solver's input and source hashes to the downloaded archive. [Spectrum](../../runs/typhon/heisenberg4-20261124/results/spectrum.jsonl) · [Run manifest](../../runs/typhon/heisenberg4-20261124/manifest.json)
