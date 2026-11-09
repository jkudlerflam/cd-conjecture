This log was written by a Claude Opus 5.5 agent (C166, 2026-11-08). Astra was not used (out of credits).

# C166 progress log

- Copied unchanged from codex/work_C163: dmrg.py, dmrg_c.py, est.py, vumps.py, fesx.py, fesfit_c158.py, c163lib.py, tk.py. inputs/cand_{1,2,59,60,61,63}.json are copies of codex/work_C155/typhon_out/refine/. inputs/ref_open.json collects the earlier open-chain S_mid, c_prof, E (C159 L = 32 to 128, C163 L = 192 and 256, C155 L = 32 to 128), with the source file of each entry.
- New: c166lib.py (bonds, ED, DMRG with arbitrary edge fields, MPS observables), ed_check.py (task 1), vumps_cells.py (task 2), open_chain.py (tasks 3 and 4), summarize.py (collection), typhon/ (job inputs).
- Pass criterion written into RESULT.md before the first run.
- ed_check.py → ed_check.out, ed_check.json. All seven bonds: MPO and bond matrix agree to 1e-16, open DMRG = ED to 1e-12 at L = 8, 12, 16, VUMPS reproduces C163's χ = 16 energies to 1e-10.
- Local VUMPS scan (local/scan_tasks.txt, χ 8 to 32, seeds 0 and 1, n = 1 to 8, 10 to 14, 20, 21, both XZ bonds) → local/vc/, local/scan.log. Every even cell converges to the same Néel state as n = 2 (energy identical to 1e-12). Odd cells give non-injective cats at higher energy.
- Local open chains (local/oc/): hi_T_40 at L = 32, 64, 128 without field and with edge fields; dist_T_50 at L = 64.
- Prediction for the Typhon open chains registered in RESULT.md before the job started.
- Typhon job 17602270 submitted 2026-11-08 (PENDING, Priority). Directory /scratch/lustre/jkudlerflam/c166/, inputs typhon/tasks_heavy.txt (6 tasks, 6 threads each) and typhon/tasks_light.txt (46 tasks, 1 thread each), 40 h limit. Outputs oc/*.json, vc/*.json (+ .npz checkpoints), logs/heavy.log, logs/light.log, summary.txt (written by summarize.py at the end).
- To collect: `scp -r typhon:/scratch/lustre/jkudlerflam/c166/{oc,vc,summary.txt,logs} codex/work_C166/typhon_out/` (skip vc/*.npz if large), then `cd codex/work_C166 && python3 summarize.py typhon_out > typhon_out/summary_local.txt`. summarize.py works on partial output too (files are written per χ for VUMPS and per run for open chains).
