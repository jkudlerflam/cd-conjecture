# C158: Infinite-size verification near the XZ point

Your folder is `../work_C158/`. Read `specs/_qubit_c_context.md` and `specs/C157_xz_counterexample.md` first, then `codex/work_C144/` (RESULT.md, PROOF.md, dmrg.py), `codex/work_C156/RESULT.md` (the Kitaev chain in a [110] field is the XZ point with a field along the bisector of its two axes, reported critical for 0 < h < 0.511K with no c measured), and the VUMPS and finite-entanglement-scaling code left by C153 (`codex/work_C153/vumps.py`, `fes.py`, `udmrg.py` and PROGRESS.md). You may read `codex/work_C157/` and `codex/work_C155/` as they appear, but never write there.

## Goal

An independent, effective-model-free determination of the total central charge near XZ, by methods suited to the long crossover there (emergent velocities of order δ or δ², energy scales near 1e-3 of the bond). The test that matters: is there a nearest-neighbour qubit bond with c > 1?

## Tasks

1. Validate the infinite-size machinery (VUMPS or iDMRG with finite-entanglement scaling, S = (c/6) log ξ, and the Pollmann–Moore scaling of ξ with bond dimension) on known points: critical Ising (1/2), XX and Heisenberg (1), the XZ + ZX point (two decoupled Majoranas, 1), and XZ itself (where the exact degeneracy must be handled, for example by a small symmetry-breaking regulator or by working in a fixed sector of the B algebra).
2. Measure c along (a) the bisector-field line h = XZ-type bond + B (x̂ + ẑ)/√2 for several B in (0, 0.5) in the units of arXiv:2209.06221, (b) the YX and ZY directions at XZ with moderate δ, and (c) points on the 7-dimensional tuned hyperplane dJ_zx + ⟨a⟩ dJ_yy = 0 proposed by C157 when available. Use bond dimensions large enough that ξ exceeds the crossover length, and report c with error bars from fits over several bond dimensions, cross-checked by finite open chains (L = 128 to 512) where feasible.
3. For any point with c clearly above 1, characterise it: the entanglement spectrum, the transfer-matrix spectrum (scaling dimensions from ratios of correlation lengths), and the momentum structure, enough to identify the CFT.

Rules as in the shared context. Run the heavy work on Typhon (`ssh -o LogLevel=ERROR typhon`, scratch /scratch/lustre/jkudlerflam/c158/, SLURM only, one packed 48-core node job, never block on a queued job; record job ids in PROGRESS.md).
