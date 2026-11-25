# Validation ledger

| Check | Status | Evidence / limits |
|---|---|---|
| Six-file historical table reproduction | Pass | Byte-identical numeric table, analysis/historical_table_reproduced.txt |
| Filter widths, harmonics, references and fit windows | Completed; unstable | 12,672 deterministic variants and 480 additional envelope/momentum variants |
| Known-power synthetic extraction | Historical estimator fails | p=4 becomes 5.376; forward-filtered pure-power calibration passes; neutral identification still fails |
| Literal Pauli Hamiltonian, end fields and gauge | Pass | L6 error <2e-15, local and IAS |
| Quartet and enlarged operator matrix | Pass | L12 random complex state <3e-16; production offsets L24 parity cat <7e-16 |
| Checkpoint and restart | Pass for tested cases | Forced deadline after one L8 chi4 sweep; reload at chi16 improves E by .0177 and converges to ED within 4e-15; L128 E/Z/S roundtrip unchanged; all three pilot stages retained |
| L=128 reference | Control passes; target does not match | Control Smid agrees with rounded same-chi references through chi256; target chi128 differs by .16107 despite sweep convergence, and chi256 is still unconverged |
| CPU / memory pilot | Completed | Pilot 50:58, TotalCPU 3:21:35, Slurm MaxRSS 509872 KB; timing experiment's MKL count stayed one despite 1/4/8 requests, so no speedup inferred |
| Target chi256 sweep convergence | Fails | Shelved at 24 sweeps; abs(Delta E)=1.473e-6, abs(Delta S)=.006754, above 1e-9 and 1e-5 tolerances |
| Final saved-state measurement | Pass | Both chi256 states have parity +1 and norm errors <1e-9; repeated full control pair/quartet arrays agree exactly with the original measurement |
| Target quartet chi convergence | Fails in tested window | chi128 to partial chi256: 280% and 93% relative L2 changes in Q(2,5), Q(3,6), ref47, width32, r24–44; essentially unchanged at fixed momentum |
| Full L256/384 production gate | Not passed | Filter bias repaired only for calibrated controls; neutral identification and target convergence remain unresolved |
| Final publication integrity | Pass | All five allocations COMPLETED 0:0; archive and per-file hashes verified; final source/input provenance checked; Python compilation and git diff --check pass; largest published file 2090849 bytes; no HDF5 checkpoints included |

The smaller quadrature eigenvalue has no certified emergent-charge label.
Finite data cannot establish the asymptotic unitary compact Virasoro c=2
CFT, and a relevant matched quartet would still need a nonzero coupling.
