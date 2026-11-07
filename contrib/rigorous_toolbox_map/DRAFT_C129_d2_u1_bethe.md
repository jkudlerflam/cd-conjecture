# C129 (draft): C(2) on the U(1) spin-1/2 class by rigorous Bethe ansatz

Draft spec from `contrib/rigorous_toolbox_map/`, for coordinator review. Paths are repository-relative. Your folder would be `codex/work_C129/`. Read `specs/_Cd_target.md`, `specs/_bound_background.md` with all updates, REPORT.md sections 9.75 and 9.85 to 9.88, `contrib/rigorous_toolbox_map/TOOLBOX.md` (Sections 4, 7 and 8), and `codex/work_C125/RESULT.md`. Rules of `_bound_background.md` apply (tags proved, numerical, heuristic or assumed; no probabilities in files; Python only in your folder; browsing tool only for literature; plain prose without em-dashes; PROGRESS.md, RESULT.md, PROOF.md). Large computations go in a SLURM-ready script, not submitted from the sandbox.

## Setting

Jonah has said that C(2) < infinity alone is a result he wants. The general spin-1/2 bond is, up to on-site rotation and gauge, XYZ plus a Dzyaloshinskii–Moriya vector plus a uniform field (TOOLBOX Section 8.1). Its U(1)-symmetric part is

    h = J_perp (Sx Sx + Sy Sy) + D (Sx Sy - Sy Sx) + J_z Sz Sz - (B/2)(Sz (x) 1 + 1 (x) Sz),

with S = sigma/2. Every point of this class is integrable. The target of this run is

    C_U1 = sup over (J_perp, D, J_z, B) with M > 0, and over all L >= 2, of L (E_{2L} - 2 E_L)/M < infinity,

with an explicit constant. This is the first stratum of the d = 2 cover in `DRAFT_C131_d2_cover.md`. No vacuum-compatibility restriction should be needed in this class (see step 2), which is itself worth checking.

## Literature to use (statements as read in TOOLBOX)

- Yang–Yang, Phys. Rev. 150 (1966) 321: the ground state of the periodic XXZ ring in each sector of fixed down-spin number is unique and given by Bethe's ansatz (Perron–Frobenius plus continuity in Delta). Check the exact Delta range covered and whether twisted boundary conditions are included.
- Kozlowski, arXiv:1508.05741: for H = J sum(sigma^x sigma^x + sigma^y sigma^y + Delta sigma^z sigma^z), every Delta > -1 and density D in ]0, 1/2], existence and condensation of ground-state and particle-hole Bethe roots, with an all-order asymptotic expansion of the counting function (Theorems 3.4, 3.5, 4.1). The run must decide whether the remainders are locally uniform in (Delta, D) and whether the 1/L term of the sector energy follows rigorously from the counting-function expansion.
- Duminil-Copin–Kozlowski–Krachun–Manolescu–Tikhonovskaia, arXiv:2012.11675: six-vertex free energy, Bethe-root condensation for Delta < 1, sector free energies to O(1/N) per site (Theorem 2) and O(1/(zeta (N - 2n) N)) (Theorem 22). Not sharp enough for the Casimir term alone, but its condensation proof may give uniformity in parameters.
- Tasaki arXiv:1708.05186, Corollary 1a, and the Kohn budget of REPORT §9.75: twist energies cost O(M/L).
- Aizenman–Duminil-Copin–Warzel arXiv:2002.02543 for the Néel regime Delta > 1.

## Elementary reductions (the run should write these out as proofs first)

Step 1, the field drops out. [Assessment, elementary] S^z_tot is conserved, so E_L = min_n [E_L^(n) - B(L/2 - n)] with E_L^(n) the lowest energy at B = 0 with n down spins. If n_* attains the minimum on the L ring, the 2L ring sector 2 n_* has field energy exactly twice that of the L ring, so

    E_{2L} - 2 E_L <= Lambda_L^(n_*),   Lambda_L^(n) := E_{2L}^(2n) - 2 E_L^(n).

The bound is therefore uniform in B once the sectorwise defect Lambda is controlled for all n <= L/2 (and n >= L/2 by spin flip).

Step 2, the DM term becomes a twist. [Assessment, elementary] With J = sqrt(J_perp^2 + D^2) and theta = atan2(D, J_perp), the rotation exp(-i theta sum_x x S^z_x) maps every bulk bond to J (Sx Sx + Sy Sy) + J_z Sz Sz and puts a flux Phi = L theta on the closing bond. So E_L^(n) is the sector ground energy of XXZ with Delta = J_z/J and twist Phi, and the 2L ring carries twist 2 Phi. The sign of J_perp is absorbed into theta, so odd L is allowed and no sublattice rotation is used. Check that M = J/2 for |Delta| <= 1 and record M(Delta) in general (the sanity script in `contrib/rigorous_toolbox_map/d2_sanity.py` gives M = 0.5, 0.522, 0.625 for (Delta, D) = (0, 0), (0, 0.3), (1.5, 0) at J_perp = 1).

Step 3, twists cost O(M/L). [Assessment, sketch] Let Phi_m minimize E_L^(n)(Phi) and psi its ground state. The trial exp(i (Phi' - Phi_m) sum_x x S^z_x / L) psi at twist Phi' changes each bond by a rotation of angle alpha = (Phi' - Phi_m)/L. The first-order term is alpha times the total current in psi, which vanishes by minimality, and the second-order term is at most (alpha^2/2) J ||Sx Sx + Sy Sy|| = J alpha^2/4 per bond. Hence 0 <= E_L^(n)(Phi) - min E_L^(n) <= J pi^2/(4L) up to O(L^{-3}), uniformly in Delta and n. The run must make the higher-order remainder explicit. Consequently

    E_{2L} - 2 E_L <= [min_Phi E_{2L}^(2n_*) - 2 min_Phi E_L^(n_*)] + pi^2 M/(4L),

and the problem is the twist-optimized, field-free sectorwise doubling defect of XXZ.

Free-fermion check (Delta = 0). [Numerical check done in this survey] With closed shells, E_L^(n) = -J sin(pi n/L)/sin(pi/L), so L Lambda_L^(n) = -J L sin(pi rho)[1/sin(pi/(2L)) - 2/sin(pi/L)] = (pi J/4) sin(pi rho)(1 + O(L^{-2})), with rho = n/L. It is uniform in rho, at most 0.862 J at L = 4 and tending to pi J/4. Reproduce this in PROOF.md, including the open-shell parities, which the twist optimization handles.

## The core lemma

Prove, with explicit constants,

    min_Phi E_{2L}^(2n) - 2 min_Phi E_L^(n) <= C_1 M / L

for all L >= 2, 0 <= n <= L/2 and all real Delta. Suggested decomposition, each piece being a separate statement in PROOF.md:

(a) Compact interior. For Delta in [-1 + delta, Delta_max] and rho in [delta, 1/2], use Yang–Yang identification and Kozlowski's expansion to get E_L^(n)(Phi) = L e(rho) + (2 pi v/L)[-1/12 + quadratic charge and twist terms] + r_L with |r_L| <= epsilon(L)/L, epsilon -> 0 uniformly on the compact set. Then L Lambda is bounded by a constant times v plus a finite-L_0 term, and v <= 4eM (C109). The decisive question is local uniformity of Kozlowski's remainders. If the paper only gives pointwise expansions, identify the step where uniformity fails and whether compactness in (Delta, rho) repairs it.

(b) Low density rho -> 0 at fixed Delta. Dilute magnons on a ferromagnetic background. Here a direct Bethe analysis (n fixed, L -> infinity, then n/L small) should give Lambda = O(J sin(pi rho)/L) as at Delta = 0, with the two-body phase shift as the only interaction input.

(c) The ferromagnetic end Delta -> -1+. Here v -> 0 and the crossover length diverges, the d = 2 face of lesson 3. At Delta = -1 every sector ground energy equals -J L/4 and Lambda = 0 exactly. The required statement is that Lambda stays O(M/L) through the crossover. Kozlowski's theorem excludes Delta = -1, so this piece needs its own argument (for example a variational trial built from the isotropic multiplet plus the Delta + 1 perturbation, or a Perron–Frobenius monotonicity in Delta within each sector).

(d) Delta >= 1. Heisenberg has logarithmic corrections, but only an upper bound on L Lambda is needed, and the data (Q_L = 4, 2.79, 2.62, 2.56 at L = 2, 4, 6, 8) decrease. For Delta > 1 at rho = 1/2 the chain is gapped with Néel order (ADW), and Lambda should be exponentially small. Away from rho = 1/2 the sector is a Luttinger liquid of domain walls and falls under (a).

(e) Delta < -1. Sector ground states phase separate, and Lambda should be negative of order the domain-wall energy. A rigorous upper bound by a trial state (one domain on the 2L ring) suffices.

## First concrete step

Write Steps 1 to 3 as complete proofs, with the explicit remainder in Step 3. Then read Kozlowski arXiv:1508.05741, Sections 3 and 4, and state in PROOF.md exactly which constants in Theorems 3.4 and 3.5 depend on (Delta, D), whether they are continuous, and whether the 1/L energy term is derived there or must be derived from the counting function. That reading decides whether piece (a) is a corollary or a new paper.

## Falsification tests

1. Solve the logarithmic Bethe equations with twist for the sector ground states (Yang–Yang quantum numbers), for L up to 4096 and 2L up to 8192, on a grid Delta in {-0.999, -0.99, -0.9, -0.5, 0, 0.5, 0.9, 1, 1.5, 3}, rho in {1/L, 2/L, ..., 1/2} (subsampled), Phi optimized. Record L Lambda_L^(n)/M. The lemma is falsified if any fixed (Delta, rho) family, or the families Delta -> -1 or rho -> 0 at fixed L/n relations, shows growth without bound. Cross-check against exact diagonalization at L <= 10 sector by sector.
2. Check Step 3 numerically: E_L^(n)(Phi) - min E_L^(n) <= J pi^2/(4L) on the same grid at small L by exact diagonalization. [Numerical, done in the survey, `contrib/rigorous_toolbox_map/twist_check.py`] For L = 4 to 8, all sectors n <= L/2, Delta in {-0.9, 0, 0.5, 1, 2} and 121 twists, the largest ratio of the twist spread to J pi^2/(4L) is 0.671. Extend to larger L and to Delta near -1.
3. Check Step 1 end to end: for the 36 points of `contrib/rigorous_toolbox_map/d2_sanity_output.txt`, the bound E_{2L} - 2E_L <= Lambda^(n_*) + pi^2 M/(4L) must hold at every listed L.

## What success means

A theorem that L(E_{2L} - 2E_L) <= C_U1 M on the whole U(1) spin-1/2 class, for every L, with explicit C_U1, and a list of which inputs are literature theorems and which are new. Partial success, (a) plus (b) plus (e), would already give a rigorous bound on the class minus neighbourhoods of Delta = -1 and of Delta = 1 at half filling, which is the shape of result the cover in C131 needs.
