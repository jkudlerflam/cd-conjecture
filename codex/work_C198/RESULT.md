# C198: the open-chain size sequence of the effective model, L = 128 to 1024 (registration, 2026-11-10 21:10 EST, before any output)

Coordinator, Claude Fable 5. Jonah asked on 2026-11-10 for everything that could strengthen the d = 2 candidate while L = 512 runs.

## Question

On the qubit chain the two-size readings fall with size, c_two = 2.17 to 2.36 (64 → 128), 2.00 (128 → 256, converged, REPORT §9.192), and an interim 1.58 at matched χ = 1024 for 256 → 512 (C195, unconverged). The effective model of the blockaded particles (C173 eq. 3.2) reads c = 2.0 on the infinite chain out to ξ ≈ 600 (§9.179). Its open chains at n = 3/16 already read c_two(W = 6) = 1.82 (64 → 128), 1.87 (128 → 256) and 1.89 (192 → 256) at χ ≤ 1024 (C176, `codex/work_C176/typhon_out/summary.txt`), below 2 and rising. C198 extends that sequence to L = 512 and 1024 with larger χ, so that the shape of the open-chain sequence of a model whose infinite chain is c = 2 can be set beside the qubit chain's.

## Runs (C176's run.py, fibdmrg.py and analyze.py unchanged; Typhon folder /scratch/lustre/jkudlerflam/c19x/c198/)

- sym316 (symmetric line, |w| = 1, φ = π/2, bare couplings, N = 3L/16) at L = 128, 192, 256, continued from C176's checkpoints (copied, the source untouched) with two more sweeps at χ = 1536.
- sym316 at L = 384 (χ to 1536), L = 512 and L = 1024 (χ to 2048), with C176's schedule below that.
- Controls with the same estimator: hard rods c1 (c = 1) at L = 512, two XX legs c2 (c = 2) at L = 512 and 1024, χ to 1536.

## What is reported (no pass or fail is registered; this is a comparison)

- c_two(L → 2L) for 64 → 128, 128 → 256, 256 → 512 and 512 → 1024 at W = 1, 6, 12, each size at its last χ, with the entropy step of its last χ doubling and the final truncation, and the same at matched χ.
- The controls' readings at 256 → 512 and 512 → 1024, as the calibration of the estimator at these sizes.
- Read as follows, stated before the output: if the effective-model sequence approaches 2 (non-decreasing, or decreasing toward 2 from above) at the sizes where the qubit chain is measured, the qubit chain's sizes can be compared with it directly. If the effective model also falls clearly below 2 at 512 → 1024 while its infinite chain reads 2.0, the open-chain two-size estimator is not reliable at these sizes for this model and the qubit-chain decline cannot be read as a flow by itself.

Job: see below once submitted.

Job 17615895 (c198, partition typhon, one 48-core node, 3 days), submitted 2026-11-10 21:12 EST. Collect with `scp -r typhon:/scratch/lustre/jkudlerflam/c19x/c198/summary.txt typhon:/scratch/lustre/jkudlerflam/c19x/c198/logs typhon_out/` and `scp 'typhon:/scratch/lustre/jkudlerflam/c19x/c198/out/*.jsonl' typhon_out/out/` (never out/ckpt), then `python3 -E -s analyze.py typhon_out/out` with C176's analyze.py.

The single-node job 17615895 was cancelled before it started and replaced by nine one-run jobs (typhon/c198_one.slurm, typhon/submit_split.sh), 17615896 to 17615904, named c198_<tag>, so that each can start in a gap of the queue. Start estimate at submission 2026-11-12 09:52 EST.
