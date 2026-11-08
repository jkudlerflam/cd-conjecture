# Shared context for C152 to C156: is c ≤ 1 for every nearest-neighbour qubit chain?

Paths are repository-relative. Read `STATUS.md` (items 1, 4 and the sub-blocks of item 5 dated 2026-11-07), `specs/_Cd_target.md`, `contrib/d2_gap_reduction/` (README.md and its four .md files), `codex/work_C141/` (RESULT.md, PROOF.md sections 5 and 6), `codex/work_C143/` and `codex/work_C144/` (RESULT.md, PROOF.md), and the REPORT.md sections on the d = 2 survey (search REPORT.md for "C6_d2", "103,000" and "9.37").

## The question

Every translation-invariant nearest-neighbour qubit chain can be written, up to a telescoping term, as

    H_L = Σ_j [ σ_j^T J σ_{j+1} + B · σ_j ],

with J a real 3 × 3 matrix and B ∈ ℝ³. Modulo a common SO(3) rotation and positive scale this is the 8-dimensional compact space B_2. The claim to prove or refute: at every point of B_2 where the chain is gapless with a conformal infrared, the total central charge satisfies c ≤ 1. This is the d = 2 case of the c ≤ d − 1 conjecture, which is false at large fixed d (STATUS item 1). We want a derivation at the theoretical-physics standard (real derivations, not just evidence: RG, bosonization, exact solutions, controlled perturbation theory, calibrated numerics), and rigorous proofs wherever reachable.

## What is established

- Uniform U(1) class (XXZ + z-field + z-DM): c = 1 at every gapless point (classical Bethe ansatz and bosonization; rigorous finite-size asymptotics at fixed parameters by Kozlowski below half filling). The campaign proved uniform finite-size bounds there (C129 to C149).
- Free-fermion locus (C143 §2.4): one complex fermion per site with a single-harmonic band, so c ≤ 1.
- Zero-field XYZ (Baxter): critical only where two couplings coincide, c = 1.
- Symmetric ferromagnetic rank-one exchange with any uniform field: c ≤ 1/2, via the Bravyi–Gosset–Liu–Wong gap theorem (arXiv:2607.14401), and a larger dominant-axis criterion (contrib/d2_gap_reduction results 1 and 2).
- Every bond with an exact additive charge Σ n_j·σ_j reduces to XXZ, commuting Ising, or the staggered-charge class a(XX − YY) + uZZ + bZ (equivalently XXZ in a staggered field), whose open region is −1 < u < √(1 + b²/4), b > 0, u ≠ 0 (d2_gap_reduction results 3 and 4).
- Singular strata (C141, C143, C144): Fibonacci point, frustration-free family, band-onset surfaces, and the XZ point, each with emergent c ≤ 1 in every direction checked. Near XZ the tuned Ashkin–Teller cone has a gapped wedge with a c = 1/2 line. The YX and ZY directions at XZ show a slow sector with c ≈ 1 at L = 64, not converged (Typhon jobs 17580192 and 17580432 pending, outputs in /scratch/lustre/jkudlerflam/c144/).
- Numerical: about 103,000 sampled bonds, the staggered branch to L = 256 and the floating phases found no c > 1. Known qubit models with c = 3/2 (O'Brien–Fendley, Rahmani–Zhu–Franz–Affleck) need three-site terms, which is d = 4 after blocking.

## What is missing

(i) The parity-symmetric class (C152), (ii) the staggered-charge box (C153), (iii) bonds without any symmetry (C154), (iv) an adversarial numerical search aimed at c > 1 rather than random sampling (C155), and (v) a literature check for any known nearest-neighbour qubit chain with c > 1 (C156). The structural heuristic so far is "one species per site": a qubit carries one fermion after Jordan–Wigner, and a second gapless species would need bound states that form a separate fluid, which C141 argues against (cluster energies) only heuristically.

## Rules

Rules of `specs/_bound_background.md` apply (tag claims proved, derived, numerical, heuristic or assumed; no probabilities in files; Python only in your folder; plain prose without em-dashes; PROGRESS.md, RESULT.md, PROOF.md). No package installs (NumPy, SciPy, numba and mpmath are available; C144's NumPy DMRG in `codex/work_C144/dmrg.py` is validated and may be copied). Typhon is available via `ssh -o LogLevel=ERROR typhon`: SLURM only, never compute on the login node, pack tasks into one node job (partition typhon, --cpus-per-task=48) because the job limit is often full, scratch under /scratch/lustre/jkudlerflam/<your run>, and never block on a queued job (record the job id and output path in PROGRESS.md). Literature lookups by web search and fetch are fine.
