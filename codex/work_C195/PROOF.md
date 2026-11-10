This work was done by a Claude Opus 5.5 agent (C195, 2026-11-10), not by Codex.

# C195 PROOF: storage of the offloaded sweep and memory estimates

## 1. Why the offloaded sweep computes the same numbers

A two-site sweep of `z2dmrg_xp.dmrg` (C190g) at bond j reads four arrays, the left environment L_j, the right
environment R_{j+2} and the tensors M_j, M_{j+1}, plus the MPO tensors W_j, W_{j+1}. It writes M_j, M_{j+1}, the
bond label n0_{j+1}, and one new environment (L_{j+1} when moving right, R_{j+1} when moving left). It sets the
environment that can never be read again in this sweep to None (R_{j+2} moving right, L_{j+1} moving left).
Nothing else is touched. So a sweep only needs, on the device, the four arrays of the current bond, and of
these two are carried from the previous bond (moving right, L_j is the environment just built and M_j is the B
tensor of the previous split; moving left, R_{j+2} and M_{j+1} likewise). Per bond the offloaded sweep therefore
uploads R_{j+2} and M_{j+1} (right) or L_j and M_j (left), and downloads the new environment and the finished
tensor. Every arithmetic call (TwoSite, kron4, davidson, split, env_left, env_right, schmidt_probs_center, the
<Z> profile) receives the same arrays in the same order as in the C190g code, so the results are the same up to
the memory layout of the uploaded copies, which can change the summation order inside BLAS or cuBLAS at the
level of rounding. The local numpy check (validate_off.py) shows bit-identical energies, entropies and <Z> in
the mixer sweeps and differences of 1e-14 in E and 2e-12 in S in the mixer-off sweeps, the latter also with one
BLAS thread, which is consistent with layout effects only. The resume path (a segment that starts on the GPU)
does the right canonicalization and the right environments site by site, again with the same calls in the same
order.

## 2. Memory

Complex128, MPO bond dimension w = 5 (channels start, Z, end, S+, S-; models.candidate). An environment is a dense
D × w × D array, 80 D² bytes. An MPS tensor is D × 2 × D, 32 D² bytes. At L = 512 the bond dimension is
D_j = min(χ, 2^min(j, L−j)), so it is at χ for 493 bonds at χ = 768 and 1024 and for 491 bonds at χ = 1536.
One environment per bond is stored at any time (left of the moving bond the left ones, right of it the right
ones).

| χ | one environment | one tensor | all environments | all tensors | host total |
|---|---|---|---|---|---|
| 768 | 47 MB | 19 MB | 23 GB | 9.3 GB | 33 GB |
| 1024 | 84 MB | 34 MB | 41 GB | 16.5 GB | 58 GB |
| 1536 | 189 MB | 75 MB | 93 GB | 37 GB | 130 GB |

Without offloading all of this sits on the 40 GB card, which is why L = 512 at χ ≥ 1024 does not fit (the L = 256,
χ = 768 runs of C190g held 27 GB). With offloading the device holds the arrays of one bond, the TwoSite blocks
(about one more environment), the Davidson basis (2 × 32 vectors of 2D² complex numbers, 2.1 GB at χ = 1024 and
4.8 GB at χ = 1536) and the step-1 intermediate (20 D² complex numbers, 0.75 GB at χ = 1536). The measured values
are in RESULT.md (benchmark).

The final measurements (entropies, one-body matrices, C_P) run on the host after dmrg() has returned, so the
environments are freed by then. They copy the MPS once (z2dmrg.measure and z2dmrg.entropies each work on a copy,
one after the other), so the peak is about two MPS, 75 GB at χ = 1536.
