# C155 derivations behind the estimators and labels

This run was done by a Claude Opus agent (2026-11-07), not by Codex. C155 is a numerical search, so this file only records the derivations its code relies on, with the checks that confirm them.

## 1. Momentum-sector matrix elements (est.py)

[derived, checked] Let T be the one-site translation and, for a representative r with orbit period p_r and a momentum k with k p_r ∈ 2πℤ, set |r, k⟩ = p_r^{-1/2} Σ_{l=0}^{p_r - 1} e^{-ikl} T^l |r⟩. Write every basis state as y = T^{s_y} r_y with r_y its representative. For a translation-invariant H,

    ⟨r_y, k| H |r, k⟩ = Σ_y H_{y r} e^{i k s_y} (p_r / p_{r_y})^{1/2},

where the sum runs over the computational states y reached from r by one bond term, and terms with r_y incompatible with k are dropped. Proof: insert the definitions. The l-th term of the ket sum, after translating the bra-side orbit sum by T^{-l}, equals the l = 0 term, because H commutes with T and the phase e^{-ikl} is absorbed into the bra phase. That gives p_r copies of the l = 0 term, and dividing by the normalizations gives the formula. Incompatible r_y have orbit sums that vanish identically in sector k. Since H_{yr} = h[out, in] for the bond that produced y, every matrix element is a fixed coefficient times one entry of h, which is why est.py stores (row, column, entry index of h, coefficient) once per ring and assembles a new bond by a gather and a bincount. Check: the union of all sector spectra agrees with dense diagonalization to 1e-14 at L = 6, 8 for random B_2 bonds and to 2e-13 for a random range-3 bond at L = 8.

## 2. Screening estimators

[derived from standard CFT finite-size formulas]
- Periodic ring: S_L(l) = (c/3) log[(L/π) sin(πl/L)] + s_1 + (non-universal oscillations). At l = L/2 the difference between two sizes is c_S = 3 [S_{L2}(L2/2) - S_{L1}(L1/2)] / log(L2/L1), independent of the velocity and of s_1.
- Casimir: E_0(L)/L = e_0 - π v c / (6 L^2) + O(L^{-2-x}), so πvc/6 = (E_0(L2)/L2 - E_0(L1)/L1) / (L1^{-2} - L2^{-2}).
- Velocity: a level at total momentum shifted by 2π/L from the ground state, with scaling dimension x, sits at 2πv x/L. The estimator assumes x = 1 for the lowest such level. If the lowest level at that momentum has x < 1 (a low primary that happens to carry momentum 2π/L on a given ring, the sliding-mode situation), v is underestimated and c_cas overestimated; if it has x > 1 (Ising, L_{-1}σ with x = 9/8) c_cas is underestimated, which is the 0.448 at the Ising point. This bias is why the search uses min(c_S, c_cas) and why the refinement uses only velocity-free entanglement estimators.
- Gap exponent: z_eff = -log[Δ(L2)/Δ(L1)] / log(L2/L1). A CFT gives z_eff → 1, a gapped phase 0, a z = 2 point 2, and an exponentially split cat or dimer doublet a large z_eff for the first gap, in which case the second gap carries the physical scaling.
- Open chain: S(l) = (c/6) log[(2L/π) sin(πl/L)] + s_0 + b (-1)^l f(l), fitted on L/8 ≤ l ≤ 7L/8 either with an explicit (-1)^l term or after three-point smoothing; the two-size version uses smoothed midpoints, 6 ΔS / log(L2/L1).

## 3. Additive-charge test (crit.py, labels)

[derived] A uniform charge Q = Σ_j n·σ_j commutes with H_L = Σ_j h_{j,j+1} on every ring iff the bond commutator is telescoping, [h, n·σ ⊗ 1 + 1 ⊗ n·σ] = A ⊗ 1 - 1 ⊗ A for some single-site A. Sufficiency: the sum over bonds telescopes on a ring. Necessity: [H_L, Q] is a sum of two-site terms C_{j,j+1}, and on a ring of length at least 3 a sum of translates of one two-site operator vanishes only if the operator is of telescoping form (the standard argument comparing the parts of C that act nontrivially on both sites, which must cancel because each such part appears on one bond only). The staggered charge Σ_j (-1)^j n·σ_j gives the condition with alternating sign, [h, n·σ ⊗ 1 - 1 ⊗ n·σ] = A ⊗ 1 + 1 ⊗ A, on even rings. Both are linear in (n, A), six real unknowns, and the label is set when the null space is nonzero. Single-site A can be taken traceless since identity parts cancel. Check: every XXZ, XXZ + DM_z, field-along-axis, and staggered-box seed gets the expected label, and XZ + YX gets none (labels_known.out).

[derived] Free-fermion locus: after a rotation taking n to ẑ, the bond is quadratic in Jordan–Wigner fermions iff it contains no ZZ, XZ, ZX, YZ, ZY terms and no transverse field, that is J n = 0, J^T n = 0 and f ∥ n. The label tests this with the common null vector of J.

## 4. What these estimators cannot see

[heuristic] Small rings do not resolve sectors with velocity below about 2π/(L2 Δ_UV); near the XZ point the emergent velocities are O(δ^2) (C144), so the ring estimators there are a lower bound on any slow sector's contribution only if the sector is already in its scaling regime. A gapless sector with crossover length beyond the largest DMRG chain (256 sites) is invisible to every estimator used here. In particular the Kitaev chain in a [110] field has a reported low-energy scale near 1e-3 K (C156), well beyond these sizes, so its central charge is not decided by this run.
