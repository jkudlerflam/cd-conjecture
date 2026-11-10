C193b log (Claude Opus 5.5 agent, 2026-11-10). Typhon times from Typhon's date (EST).

1. Read the C193b brief, stage one RESULT/PROOF and scripts, C190 and C190g RESULT, z2dmrg.measure and the checkpoint code, run190(g).py, models.py, C176 fibdmrg.measure (Dpair, Cnn), C192 sq2.py, pairconf.py, CHECKPOINT.
2. Copied the registration B0 to B3 into RESULT.md before any output existed.
3. 07:25 Listed checkpoints (read only). Every GPU run of outG has its checkpoint (52 files, χ = 384, 512, 768 at L = 64, 128, 256), and xx1_L256 is in the CPU tree outB/ckpt. No final state is missing. lib/ holds unchanged copies of z2dmrg.py, z2dmrg_xp.py, models.py (md5 equal to the Typhon copies).
4. Wrote meas193.py, val193.py. Local numpy validation at L = 10 and 12 (PROOF §6), contraction and operators exact to 4e-15.
5. 07:32 Job 17611514 (c193, smoke, three controls). xx2 L = 256 in 30 s including the 17 s load. Checks to 1e-12.
6. 07:34 Job 17611534 (c193, production, 42 states: controls, then L = 256, 128, 64 candidates). Skips finished outputs.
7. Controls collected, trial ranges, frozen.json, calibration section of RESULT.md written before any candidate output was opened. B0 holds.
8. 07:46 Job 17611534 finished, all 42 outputs done, every state check passes (E to 2.7e-12, S̄₆ to 2e-14, stored Gb to 1.2e-14). Collected to typhon_out/out and typhon_out/logs.
9. ana193.py all (analysis.txt, analysis.json). Registered letter B3, carried by the structure-factor estimate only; with the pair estimate B2. B1 fails at sym40_n18 (η_b off by 0.054).
10. Unregistered diagnostics (diag193.py): the qubit chain's local η_b falls with r/L, the controls and the effective model do not. A fit with a pinned charge phase at the ends (θ_c Dirichlet) gives η = 0.325 to 0.330 (3/16) and 0.353 (1/8), equal to the effective model, and K_s ≈ 1.47 and 1.30 at L = 256. Written as R4.
11. Rerun commands for new states (χ = 1024) in RESULT.md. No job of anyone else touched, one c193 job at a time, nothing outside ckpt_meas edited.
