This red-team review was done by Claude Opus 5.5 (agent C174, 2026-11-08), not by Codex. Astra was not used (out of credits).

# C174 result: red team of the conditional theorem c ≤ 1 on B_2

Full report in REDTEAM.md, short derivations in PROOF.md, log in PROGRESS.md. One local script, rank_check.py (output rank_check.out, seconds). No Typhon job.

## Verdict

[derived] No critical hole. C154 PROOF 7 is a valid implication. It is not sound modulo only (A_W), (H)/(A6) and A_S. And (H) is not a mild topological assumption. Given (A2), the closure of a c > 1 component cannot contain any controlled conformal or gapped point (Lemma R1, proved). So (H), restricted to the components that matter, says that every c > 1 component accumulates on a controlled singular point. With (A5) this is equivalent to "no c > 1 in N". The theorem reduces B_2 to N and proves nothing inside N.

## Holes, with severity and narrowest fix

- R1 [high] (H) is equivalent to the conclusion on N. Fix: restate PROOF 7 as a reduction, with (H*) or (A6) as the stated input. Connectivity arguments (C168) cannot prove it.
- R2 [high] Read literally, (A5) covers every singular point. C143's list includes band-onset hypersurfaces in N, so (A5) would contain the P3 exclusion (second onset inside a liquid) and the Lifshitz junctions of C154 PROOF 5.4. Fix: an explicit controlled singular list S_ctrl.
- R3 [high] Direction-wise control is not neighbourhood control (a quantifier interchange). This bites at XZ (curved approaches, and the unresolved Néel-cat region), at h_F (C161: ε → 0 and L → ∞ do not commute at soliton onsets) and along γ (soliton-onset layers). Fix: weighted blow-up charts covering each punctured neighbourhood, or remove the layers from S_ctrl.
- R4 [medium-high] Stratum results hold at generic points. A c > 1 component can attach only at a c > 1 point of a stratum, an isolated multicritical point that numerics miss. Uncovered: C152 tricritical Ising points and the Peotta region, C153 large b, the cyclic SU(2)_1/ordered boundary and its coarse grid. Fix: a junction inventory per stratum.
- R5 [medium-high] Irreducible c > 1 theories have only c ≤ 1 or gapped neighbours and need no second species. Examples are Z_N parafermions (N ≥ 5, c = 8/7 at N = 5, k = 2 in N), N = 1 superconformal minimal models with 1 < c < 3/2, W_3 minimal models and SU(2)_3. They are natural at the melting of Z_m-locked tongues with m ≥ 5. Fix: add a pattern P4 with its codimension count.
- R6 [medium] The K_II sheet is under no assumption. It has a repulsive c > 1 window (K < 2/15, or 1/(4K_1) + 1/(4K_2) > 2). Fix: assign it to (A6) or state an A_W-type assumption.
- R7 [medium] Unlisted inputs: (A1), (A2), (SC), (A_sol′), the asymptotic Schrieffer-Wolff step of C163 Proposition G, Lemma 3 at spin 1/2, and pending or uncollected jobs (17597010, 17596442, 17594765, 17602270, the C163 FES stage). Fix: list them in STATUS.
- R8 [medium] Translation broken with a larger cell. Partial fix supplied [proved]: only one bond crosses a cell boundary, so det h(K) is a degree-2 trigonometric polynomial for every period. That gives c_MF ≤ 1 at any period in the parity class (extending C152 P2 from p ≤ 2) and at most two soft momenta at onset (checked to 2e-14, p ≤ 7). Beyond mean field it is open. c must also be defined per superselection sector, or cats and first-order surfaces read c = 2.
- R9 [low-medium] Lemma 3 is linear spin waves about a product state, uncontrolled at onset for dressed polarized phases at spin 1/2. In N its conclusion survives by codimension. Fix: re-derive it that way.
- R10 [low] The strata list is complete for on-site product symmetries (the degenerate-singular-value case closes, PROOF 3). Spatial-reflection strata (dimension ≤ 5) are unanalysed.
- R11 [low, checked sound] Gap stability, "c constant on Σ", redundancy of the λ-neighbourhood, the incommensurate count, and (A2) in the form actually used.
- R12 [low] The only symmetry-free connectivity test (the 1282 path) was never completed.

## Left open

A dynamical exclusion of c > 1 fixed points in N away from S_ctrl, including P4. The junction inventories of R4. Weighted-chart coverage at XZ, h_F and γ.
