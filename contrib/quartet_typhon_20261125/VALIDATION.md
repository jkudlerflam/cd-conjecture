# Validation ledger

| Check | Status | Evidence / limits |
|---|---|---|
| Six-file historical table reproduction | Pass | Byte-identical numeric table, analysis/historical_table_reproduced.txt |
| Filter widths, harmonics, references and fit windows | Completed; unstable | 12,672 deterministic variants and 480 additional envelope/momentum variants |
| Known-power synthetic extraction | Historical estimator fails | p=4 becomes 5.376; forward-filtered pure-power calibration passes; neutral identification still fails |
| Literal Pauli Hamiltonian, end fields and gauge | Pass | L6 error <2e-15, local and IAS |
| Quartet and enlarged operator matrix | Pass | L12 random complex state <3e-16; production offsets L24 parity cat <7e-16 |
| Checkpoint and restart | Pass at small L; pilot warm restart active | Forced deadline after one L8 sweep; reload converges to ED within 2e-15; L128 E/Z/S roundtrip unchanged |
| L=128 reference | Pending | Match parity and chi; report residual convergence |
| CPU / memory pilot | Running; thread diagnostic complete | Timing test reports only one MKL thread despite 1/4/8 requests; no parallel speedup inferred |

The smaller quadrature eigenvalue has no certified emergent-charge label.
Finite data cannot establish the asymptotic unitary compact Virasoro c=2
CFT, and a relevant matched quartet would still need a nonzero coupling.
