# Run ledger

All submissions use `python3 tools/typhon/typhon.py submit RUN RECIPE`.
Source bytes, SHA-256 hashes, seed, resources and Git revision are retained
per run. A successful allocation is not a convergence certificate.

| Run | Slurm ID | Recipe | Resources | State / numerical outcome |
|---|---|---|---|---|
| quartet-pilot-20261125a | 17764244 | pilot.job.json | 4 CPU, 16 GB, 1 h | Running; dense checks passed; L128 even sector, chi ladder through 512 requested |

Remote states: `/home/pace/cd-conjecture/runs/quartet-pilot-20261125a/states/`.
Fresh initial state, seed 1. Submission uses working-tree source at initial
upstream 389eadff; exact source hashes are in the uploaded snapshot and
`results/provenance.json`. The new files were uncommitted at submission.
No large MPS files are fetched or committed.
