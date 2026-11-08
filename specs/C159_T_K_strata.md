# C159: c ≤ 1 on the time-reversal strata T and K

Your folder is `../work_C159/`. Read `specs/_qubit_c_context.md` first, then `codex/work_C154/` (RESULT.md and PROOF.md, especially sections 1, 4 and 7, where the complete list of symmetric strata is proved and the conditional theorem is stated), `codex/work_C68/RESULT.md` (where strata T, K and N are defined), `codex/work_C156/RESULT.md` (literature: Kitaev–Γ, Kitaev–Heisenberg–Γ, compass and related chains, which lie in these strata after site-dependent rotations), and `contrib/d2_classification/`, `contrib/d2_floating_ising/`, `contrib/d2_two_fluid/` READMEs.

## Setting

Stratum T: B = 0 with an arbitrary real J (symmetric exchange including Γ-type off-diagonal terms, plus DM), whose only internal symmetry is Kramers time reversal; dimension 5. Ogata–Tasaki (the LSM theorem for time reversal with half-integer spin per unit cell) forbids a unique gapped ground state, so every point of T is gapless or has degenerate ground states. Stratum K: bonds whose only internal symmetry is one square-positive antiunitary (time reversal composed with a π rotation about a principal axis), with B perpendicular to that axis; dimension 5. Neither is a solved class. The literature models in them (Kitaev–Γ chains with emergent SU(2)₁, the Kitaev point with c = 1/2, compass chains) all report c ≤ 1.

## Tasks

1. T: use the Ogata–Tasaki constraint and the LSM anomaly to fix which infrared theories are allowed (gapless with the anomaly matched, or symmetry-broken degenerate, including dimerized phases), and derive that the gapless ones have c ≤ 1 at the theoretical-physics standard. Natural tools: the anomaly-matching statement that a minimal theory carrying the time-reversal LSM anomaly is the c = 1 SU(2)₁-type boson or a c = 1/2-type structure at special points; the site-dependent rotations that map families in T to known chains (Kitaev–Γ, XYZ, Heisenberg, compass); bosonization around the U(1) and XYZ points it contains; numerics to locate any region of T with an unexplained critical phase. Decide whether an anomaly-matching theory with c > 1 (for example SU(2)₂, c = 3/2, or two bosons) can be reached in T, and if it can, where.
2. K: same programme with one antiunitary. Map its gapless loci, connect them to the U(1), parity and XYZ classes where possible, and derive c ≤ 1 or find an obstruction.
3. Calibrate with DMRG and finite-entanglement scaling (C144's `dmrg.py`, C153's `vumps.py` and `fes.py`, C158's `fesfit.py`, all validated), on Typhon. Use entanglement, not Casimir energies, for c.
4. Theorem-style summary, every claim tagged.

The local machine is heavily loaded by other agents; run anything beyond a few seconds on Typhon (one packed 48-core node job, `ssh -o LogLevel=ERROR typhon`, scratch /scratch/lustre/jkudlerflam/c159/, never block on a queued job).
