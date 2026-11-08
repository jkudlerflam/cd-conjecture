# C155: An adversarial numerical search for c > 1 in B_2

Your folder is `../work_C155/`. Read `specs/_qubit_c_context.md` first, then `codex/work_C131/` (its adversarial optimizer and Typhon packing, `typhon_packed.slurm`) and the survey sections of REPORT.md.

## Goal

The earlier evidence for c ≤ 1 at d = 2 is random sampling. Replace it with a search that tries to make c large.

## Tasks

1. Build a fast, calibrated central-charge estimator for a qubit bond. Candidates: periodic ED at L = 10 to 20 using the Casimir coefficient together with a velocity from the lowest momentum-2π/L excitation; the half-ring entanglement entropy at two sizes; open-chain DMRG entropy profiles with Calabrese–Cardy and parity-oscillation corrections for refinement. Calibrate it on known points: transverse-field Ising (1/2), XX and Heisenberg (1), tricritical Ising if a qubit realization is available, and two decoupled Ising chains blocked into a qubit-like test where possible. Report its bias and spread, and how it behaves at gapped points and at z = 2 points (it must not report large c there).
2. Maximize the estimator over B_2 (8 parameters, use the SO(3)-reduced parametrization) with many restarts: Nelder–Mead or CMA-style search, seeded randomly and also at the known strata, multicritical intersections, and the XZ neighbourhood. Run it on Typhon as one packed 48-core node job, with checkpointing, and refine every candidate above 1.1 by DMRG at L up to 128 or more.
3. Report the largest c found, where, and whether it survives refinement. If any candidate survives, characterize it fully (spectrum, entanglement, operator content). If none does, state the search's coverage and the estimator's calibrated error.
