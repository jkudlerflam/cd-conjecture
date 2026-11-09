L = 192 stage of C181's partially polarized (FFLO) point (coordinator, Claude Opus 5.5, 2026-11-09).
Effective model of C173 (PROOF eq. 3.2) at bare couplings (the opposite-valley pair is bound), |w| = 1, n = 3/16,
phi = 1.53 and 1.555 (just off the symmetric line, valley imbalance 0.07 to 0.13), with N and N +- 1 for the parity gap,
plus c = 1 and c = 2 controls. C176's blockade-exact DMRG (fibdmrg.py, unmodified), chi to 512.
ref/ holds C181's converged L = 64 and 128 runs (Typhon backfill), so summarize.py gives 64 -> 128 and 128 -> 192.
At L = 128 these points read c_two = 1.66 to 2.18 against the c = 2 control's 1.86 to 2.15 (REPORT section 9.150).
Registered before output: consistent with c = 2 if c_two(128 -> 192) >= 1.6 in at least two windows at both phi with
L = 192 truncation below 1e-7 and controls within 0.15 of 1 and 0.25 of 2, and L*gap not rising by more than 50 percent
from 128; a crossover toward c = 1 if c_two(128 -> 192) <= 1.3 in two windows at either phi or L*gap doubles.
Needs python3 with numpy and scipy. Run: bash run_desktop.sh     then python3 summarize.py
This is the effective model, not the full qubit chain.

Amendment (2026-11-09, coordinator, registered before any ctl1fix output). cfg/ctl1_L192_N36.json was written with
"t2": 0.0, while C181's reference ctl1 runs in ref/ (L = 64, 128) keep the fibdmrg default t2 = -1. So ctl1_L192 is a
different Hamiltonian (energy per site -0.342 against -0.614) and its 128 -> 192 reading is meaningless. This was a
config error by the coordinator. The calibration uses ctl1fix_L192_N36 instead, identical to ref's ctl1 except for L and N.
ctl1_L192_N36 is kept and reported as a config error. ctl2 and the six candidate configs match ref/ exactly.
