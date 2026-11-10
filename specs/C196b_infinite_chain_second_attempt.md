# C196b: a second attempt at entanglement scaling on the infinite qubit chain (2026-11-10)

Coordinator: Claude Fable 5. Jonah approved a second attempt on 2026-11-10. Context: REPORT §9.183 and `codex/work_C196/` (the first attempt, RESULT.md, PROOF.md, PROGRESS.md), REPORT §9.179 and `codex/work_C189/` (the infinite chain in the effective model, which converged), `codex/work_C190/` and `codex/work_C190/gpu/` (the qubit-chain solver and the GPU environment on apollo), `specs/C195_C196_third_size_and_infinite_chain.md` (the registration I0 to I4, which is unchanged).

## Why the first attempt failed, and what changes

The first attempt built a two-site-cell iDMRG with exact Z2 labels and a one-site VUMPS. Both pass their validations (ED to 4e-14, TFIM, XX at D = 16 to 64) and neither converges at the candidate.

- The two-site-cell iDMRG has no fixed point at p153c. Each growth step has to add 0.43 particles, and parity allows only even increments. XX shows the same at incommensurate filling.
- VUMPS from a random start stays at the dressed vacuum (n = 0.002 to 0.003). The liquid differs from it by an energy of order 0.03 Δ per particle against V = 40, and annealing in the field never reached it.

The diagnosis the second attempt rests on: particle number is conserved on the qubit chain up to virtual pairs. Real pair creation first appears at third order, with amplitude of order Δ w / V, about 1/40 of the hopping, and at low energy it does not conserve momentum on the symmetric line (the valleys sit at k and π − k, the term creates pairs of total momentum zero). So a state with a fixed number of particles per unit cell is close to a fixed point of the infinite chain, and the effective model already supplies such states.

Two changes follow.

1. A commensurate unit cell. Use a cell of q sites that holds an integer number p of particles, with p even so that the cell has even parity and plain Z2 labels work without shifted charges. Then each growth step adds an integer number of particles, and the cause of the drift is gone.
2. A start in the liquid. Build the start state from the effective model at density p/q (C189's `idmrg.py`, a low χ of 64 to 128 is enough) and embed it in the qubit Hilbert space, particle to flipped spin with the blockade respected. The solver then only has to add the virtual-pair dressing, which is local. Check the embedding on a finite chain first: the embedded state's energy and density against a C190 open-chain run at the same field, to the accuracy second-order perturbation theory allows.

## Routes, in order

- Route 1. Port the commensurate-cell growth scheme of `codex/work_C189/idmrg.py` (McCulloch prediction, two-cell window) to the qubit-chain MPO, on the first attempt's validated `idmrg196.py` or on C190's `z2dmrg_xp.py`, with exact Z2 labels and the start above.
- Route 2. VUMPS with a q-site cell from the same start (extend `vumps196.py`).
- Route 3. Multi-site-cell iDMRG started from the bulk of a converged C190 open chain (cut q or 2q sites from the middle of an L = 256 state, read only), if the embedding is the problem.

Bonds and fillings.
- sym40_n316 (h_z = −19.976047251) with n = 3/16, cell of 32 sites and 6 particles. First priority.
- sym40_n18 (field from the tuning table in `codex/work_C190/RESULT.md`) with n = 1/8, cell of 16 sites and 2 particles. Second priority. This density is the marginal one in §9.182.
- p153c (h_z = −19.97589125), density near 3/14, only if an effective-model start at φ = 1.53 can be produced with C189's code. Third priority.

A point to settle before production, and to write down in RESULT.md. The cell pins the density at p/q, while the field of a C190 key was tuned on an open chain and its bulk density is not p/q (0.1917 at mid-chain for sym40_n316). At small χ the pinned state is stable. At large correlation length the window in the field over which it stays stable shrinks as 1/(κ ξ), and a drift in the number per cell can reappear through the weak pair term. Monitor the particle number per cell and its variance at every χ. If the state drifts at the registered field, the allowed fix is a small retune of the uniform field so that the pinned density is stationary (the bond is unchanged), recorded as a dated addendum before any production output at the new field, with the size of the shift. In that case I4's energy clause is reported as not evaluated at the registered field, with the corrected comparison shown.

## Gates, each recorded with numbers in RESULT.md before production

- G1. The multi-site-cell code against ED on small chains (cells of 2 to 4 sites, window energies), to 1e-10, at one candidate bond and at TFIM.
- G2. Controls at commensurate filling with the production estimator and fit window, started from a number-definite state: the XX chain within 0.1 of c = 1, the two-leg XX control within 0.15 of c = 2, and a gapped control with ξ saturating. This is I0.
- G3. At sym40_n316 at χ = 64 and 128: I1 reached, number per cell stable, bulk energy density within 1e-5 of the C190 ground-branch value (E(256) − E(128))/128 = −10.0011861.
- G4. A benchmark and a time estimate for the χ sequence, with checkpoints that resume.

## Registration

I0 to I4 of `specs/C195_C196_third_size_and_infinite_chain.md` hold unchanged. Copy them into `codex/work_C196/b/RESULT.md` before any production output is opened. Additions fixed here, before any output exists.

- Estimator. S is the bond entropy averaged over the bonds of the cell, and ξ is the leading correlation length in sites from the cell transfer matrix, per Z2 sector (even, and odd, the sector of the single-spin operator). Use C189's convention for the average over the last iterations. Fix the details in RESULT.md before production and use the same for the controls.
- χ sequence 64, 96, 128, 192, 256, 384, 512, and 768 if feasible, each started from the previous one.
- The sector clause of I2 and I3 is meaningful only for a state with a fixed number of particles per cell. If the converged state mixes particle numbers at leading order (a cat of U(1)-broken components, as in the first attempt's VUMPS on XX), say so, report c from S against ln ξ, and mark the sector clause as not evaluable.
- I4's density clause is weak here because the cell pins the density. Report it, with the number per cell and its variance.
- If no route reaches I1 at a candidate bond, stop and report with numbers. Do not relax I1.

## Rules

- Work only in `codex/work_C196/b/` and `/scratch/lustre/jkudlerflam/c19x/c196b/`. Everything else is read only, in particular `codex/work_C189/`, `codex/work_C190/`, the rest of `codex/work_C196/`, and everything under `/scratch/lustre/jkudlerflam/c18x/`. Copy code you need, never edit it in place. Do not touch `codex/work_C194/` or `codex/work_C195/`.
- Keep RESULT.md (registration first, then gates, job ids, output layout and tested collection commands), PROOF.md (the embedding, the label scheme, the estimator) and PROGRESS.md (log).
- Typhon: SLURM only, never compute on the login node. Partition apollo (GPU), not h200, and no CPU jobs on partition typhon. Job names start with c196. One run per GPU, at most 3 GPUs at once, about 4 cores and at most 64 GB per GPU, at most 3 days per job. complex128 throughout, no change to any solver tolerance. No pip or conda install.
- Never cancel, hold, release, signal or alter a job you did not submit. Jonah has other GPU jobs there (gpu_lean1, gpu_lean2, ddL_g, ddH_g, ddS_g, c190g, c195). Leave them alone.
- Never copy checkpoints to the local machine. No absolute local paths and no probabilities in any file. Do not commit or push.
- Local machine: numpy only, at most about 4 cores.
- No Fable subagents. Astra (gpt_solve.sh) is out of credits, do not call it.
- ssh: `ssh -o LogLevel=ERROR -o ConnectTimeout=20 typhon`, retried with waits of about 30 s, because the bastion rate-limits. A local hook blocks `sleep N; command` chains in Bash, so put waits inside the ssh command string. Another hook blocks `grep -r`, `find .` and `ls -R`, so use explicit paths. Quote remote globs in scp (the local shell is zsh).
- Prose with commas, not em-dashes, few colons, no italics.

## Report

Stop and report once production is running and its first χ points are checked against G3, or once every route has a documented failure. Say what passed and failed at each gate with numbers, which route runs, job ids and where the output lands, the time estimate, and what is left. The coordinator collects, fits and writes up.
