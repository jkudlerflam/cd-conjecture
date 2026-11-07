# C136: The attractive U(1) spin-1/2 class at every filling, and through the ferromagnetic crossover

Your folder is `../work_C136/`. Paths below are repository-relative, so prefix them with `../../` from your folder. Read `specs/_Cd_target.md`, `specs/C129_d2_u1_bethe.md`, all of `codex/work_C129/` (RESULT.md and PROOF.md, in particular Steps 1 to 3, the dilute bound and the ferromagnetic band), all of `codex/work_C134/` (RESULT.md, PROOF.md, and REFEREE.md, an independent referee report that found the proof correct and lists presentation repairs), and `codex/work_C134/LITERATURE.md`. Rules of `_bound_background.md` apply (tag claims proved, numerical, heuristic or assumed; no probabilities in files; Python only in your folder; browsing tool only for literature; plain prose without em-dashes; PROGRESS.md, RESULT.md, PROOF.md). Your sandbox cannot reach Typhon; write any large computation as a SLURM-ready script and say so.

## Setting

C134 proved, for the half-filled XXZ ring H = J sum (Sx Sx + Sy Sy + Delta Sz Sz) at the minimizing flux, and every Delta in [-1 + delta, 0] and L >= 2,

    F_L >= L e_0 - 1200 J/(delta^3 L),   L (F_{2L} - 2 F_L) <= 3000 J/delta^3,

by a real-axis argument that uses the negativity of the scattering kernel (L1 norm 2 gamma/pi - 1 < 1 with Delta = cos gamma). The referee notes the dominant small-delta term is of order delta^(-5/2). C129 reduced the whole U(1) class (XXZ plus DM_z plus field, all lengths) to twist-optimized sectorwise defects Lambda_L^(n) = F_{2L}(2n) - 2 F_L(n) plus pi^2 M/(4L), and proved the dilute bound L Lambda/M_0 <= 2 pi^2 n^3/L, the ferromagnetic band L Lambda/M_0 <= (Delta + 1) L^2, and zero defect for Delta <= -1.

## Goal

Prove, with explicit constants,

    L (E_{2L} - 2 E_L) <= C M

for every bond of the U(1) class with J_z/J in [-1, 0] (that is, Delta in [-1, 0] after the C129 reductions), every field and DM angle, and every L >= 2. This would be the first complete rigorous theorem on an open set of interacting spin-1/2 chains, and a full stratum of the d = 2 cover (C131).

## Tasks

1. Every filling. Extend the C134 argument from n = L/2 to every sector 0 < n < L/2, with constants uniform in the density rho = n/L on [rho_0, 1/2] for fixed rho_0 > 0, and on all of (0, 1/2] after combining with the C129 dilute bound. With a Fermi boundary the root density is supported on an interval; decide whether negativity of the kernel still gives a monotone counting function and a Neumann inverse with norm uniform in the interval, and how the midpoint quadrature and the clipping argument change at the Fermi points. Identify every place where the half-filled structure (roots on the whole line, no boundary) was used in C134. Kozlowski arXiv:1508.05741 Corollary 4.3 gives the pointwise 1/L term below half filling and may supply the boundary analysis; check what is uniform.
2. Through the crossover. Replace the delta-dependent constant by a bound uniform for Delta in (-1, 0]. The C129 band covers (Delta + 1) L^2 <= A, so it suffices to prove a constant uniform in delta for L >= sqrt(A/delta) (lengths beyond the crossover length of order delta^(-1/2)). Decide whether the C134 constants, tracked with their dependence on both delta and L, give this, for example because the dangerous terms are small once L exceeds the crossover length. If not, state the exact estimate that is missing.
3. Assemble. Combine Tasks 1 and 2 with C129 Steps 1 to 3 (field, DM, twist) into a single theorem on the stated class with one explicit C, and say precisely where any convention (odd L, the two-bond L = 2 ring, M versus M_0 with a field) enters.
4. Numerics. Check every intermediate inequality and the final constant with the Bethe solvers of C129 and C134 on a grid Delta in {-0.9999, -0.999, -0.99, -0.9, -0.5, 0}, densities rho in {1/L, ..., 1/2}, and L up to a few hundred, with ED cross-checks at L <= 12.
5. Optional, only after 1 to 3: say what replaces the kernel-sign argument for 0 < Delta < 1 (the kernel is positive there, still of L1 norm below one), and whether any part of the argument survives.
