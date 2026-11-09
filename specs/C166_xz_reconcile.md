# C166: Reconcile open-chain DMRG and infinite-chain VUMPS in the XZ neighbourhood

Your folder is `../work_C166/`. Read `specs/_qubit_c_context.md` first, then REPORT.md §9.125 (the conflict), `codex/work_C163/RESULT.md` and PROOF.md §4.2 (the pre-fixed criterion), `codex/work_C163/typhon_out/` (summary.txt, `big/*.json`, `fes/T40.json`, `fes/T50.json`), and `codex/work_C155/typhon_out/` (summary.txt and `refine/cand_*.json`, which hold the bonds of the rising candidates). The validated tools are in `codex/work_C163/` (`dmrg.py`, `vumps.py`, `fesx.py`, `fesfit_c158.py`, `est.py`, `c163lib.py`). Read them before writing new code, and reuse them.

## Why

Every d = 2 nearest-neighbour bond that still reads c > 1 sits in one region, the XZ neighbourhood: J diagonal about (−0.58, 0, 0.55), DM along y about −0.5 to −0.6. The two methods disagree there.
- Open-chain DMRG (truncation below 1e-9) on hi_T_40 and dist_T_50 reads c_prof 1.27 and 1.31 at L = 256, with two-size slopes 6 dS/d ln L rising to 1.44 and 1.41. S(x) at fixed x still grows with L, and the profile is CFT-shaped with c ≈ 1.25. C155 candidates 1, 2, 59, 60, 61, 63 (same neighbourhood) rise to between 1.37 and 1.91 at L = 128.
- Two-site VUMPS on the same bonds saturates S at 0.70 and 0.68 from χ ≈ 96, ξ stays near 46 and 35, the FES slope is about 0.01, and the variational energy per site lies below every open-chain estimate, which approach it from above.
- The transfer-matrix spectrum has complex leading pairs (arg ±0.38π, ±0.60π), so correlations are incommensurate. The open chains have S(1) = ln 2 exactly, so the DMRG state is time-reversal symmetric.

A gapped bulk with ξ ≈ 46 should saturate a 256-site open chain, and it does not. The job is to find out which method is wrong, and so whether this region contains a c > 1 bond (a counterexample to c ≤ 1 at d = 2) or not.

## Hypotheses to separate

(i) The bulk breaks time reversal (or another discrete symmetry) weakly. The open chain is a symmetric cat that is still crossing over, so its S picks up ln 2 and slow edge effects, while VUMPS picks one symmetry-broken sector.
(ii) VUMPS is trapped. A two-site cell incommensurate with the wavevector, or χ ≤ 192, gives a gapped variational state that is not the ground state, and the true state is critical.
(iii) The bulk is gapped but with a crossover length well above 256 sites (for example a nearby critical point with a small gap), so open chains have not yet seen saturation.
(iv) Something else, for example a bug or convention mismatch between `dmrg.py` and `vumps.py` (same bond, same orientation of DM, same sign of B). Check this first, cheaply, on a bond with a known answer and on the XZ bonds at small sizes (compare open-chain DMRG energies with exact diagonalization at L ≤ 20, and the VUMPS e₀ with the L → ∞ extrapolation of periodic ED or DMRG).

## Tasks

1. Consistency of the two codes. Same Hamiltonian in both, verified by ED at small L. State what you checked.
2. VUMPS beyond the two-site cell. Unit cells n = 1 to 8 (and n commensurate with the observed wavevectors, e.g. from arg 0.38π and 0.60π), χ up to 512, several random initializations per cell. Report e₀(n, χ), S, ξ, the FES slope, and whether the lowest-energy state for each n breaks time reversal (measure ⟨σ⟩ per site, the spin current or chirality, and S on each bond of the cell).
3. Open chains to L = 512 (χ up to 384 or what truncation below 1e-8 requires) on hi_T_40 and dist_T_50 and on two of the C155 candidates (60 and 61). Record c_prof, the two-size slopes, S(x) at fixed x versus L, and the mid-chain bond energy. Compare that bond energy (averaged over the incommensurate period) with the VUMPS e₀. If the open-chain bulk energy is above the VUMPS e₀ by more than the open-chain finite-size correction allows, VUMPS has found a lower state the open chain is missing, and vice versa.
4. Pinning test. Add a weak field at one or both edges along the direction of any symmetry-broken order parameter that VUMPS finds (and, separately, along the generators of time reversal-odd operators) and redo the open chains at L = 128 and 256. A cat collapses and S drops by about ln 2 in the bulk, while a critical state's profile is insensitive at leading order.
5. Also run the three C163 path first-stage passes (KI22, KI26, N34, bonds in `codex/work_C163/` path output) through FES at χ ≥ 256, the second stage of the pre-fixed criterion that was never run.
6. Verdict. Fix the criterion before looking (state it in RESULT.md before running): a c > 1 candidate needs open-chain fits above 1.2 that do not fall from L = 256 to 512 with truncation below 1e-8, a VUMPS FES slope above 1.2 at χ ≥ 256 for the lowest-energy cell, and agreement within error between the open-chain bulk bond energy and the VUMPS e₀. Then either an explicit c > 1 bond with its evidence, or an explanation of the conflict (which hypothesis, with the measurement that decides it) and the statement that the XZ neighbourhood has no c > 1 bond within the region tested. Tag every claim proved, derived, numerical, heuristic or assumed.

## Rules

Large runs go to Typhon as one packed 48-core node job (`ssh -o LogLevel=ERROR typhon`, SLURM only, never on the login node, scratch /scratch/lustre/jkudlerflam/c166/). Submitting to Typhon is always fine. The bastion rate-limits, so keep ssh calls few, retry with waits, and copy outputs with scp -r (not tar over ssh). Small checks (ED, χ ≤ 64) can run locally. Record job ids, expected outputs and the pass criterion in RESULT.md. Begin RESULT.md with a line saying which model did the work. Astra is out of credits, so do the work yourself and never call codex directly. No numerical claim without the file it came from. Do not edit STATUS.md, REPORT.md or files outside `codex/work_C166/`, do not commit or push (the coordinator does), never stash, reset or switch branches, and do not pip install. Rules otherwise as in the shared context.
