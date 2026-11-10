C193 log (Claude Opus 5.5 agent, 2026-11-10)

1. Read the spec, z2dmrg.measure, run190.py, models.py, C190 PROOF/RESULT, C192 CHECKPOINT, onebody.py, offline.py, pairconf.py, REPORT 9.160, 9.161, 9.178. Found that CP is the P_j = XX + YY bond correlator, not a pair function.
2. Wrote the registration and the addendum into RESULT.md before opening any candidate Gb, Gf or CP. File listing and energies only were read before that point.
3. Controls (c193.py controls): trial ranges R1 to R4, froze R3, R4 for E1 and ranges 2 to 4 for E2. Q0 holds. Window check on C176 runs (effcheck.py). Calibration section written.
4. Candidates (c193.py candidates), drift (drift.py). No range changed after step 3.
5. Wrote PROOF.md and the results. No Typhon job, no subagent, local only.
