This work was done by a Claude Opus 5.5 agent (C170, 2026-11-08), not by Codex. Astra was not used (out of credits).

# C170 result: referee of C165 §1.2, §2 and §4 (the SU(2)_1 identification)

Report in REFEREE.md, derivations in PROOF.md, log in PROGRESS.md. Scripts written here, all independent of C165's code: lr2.py (exact second-order classical relaxation), check_c1.py (corrected quartic on random harmonics), clmin.py (nonlinear classical minimization, outputs cl_g02*.out), qed.py (sparse ED), frame.py (exact frame and translation identities), heff.py (generated uniform field), plateau.py and plateau2.py (Néel plateau edge), pfull.py (ED of C165's plateau-edge bond). No Typhon use.

## Verdict

C165's conclusion survives: for weak harmonics there is no staggered Ising transition inside the liquid at spin 1/2, and c ≤ 1 near the switch. At the hidden Heisenberg point the argument is stronger than C165 states, because it rests on an exact symmetry. Three statements are wrong and one scope gap needs a different argument, which I supply.

## Claims

- [proved, correct] §1.2 exact map of the parent to the spin-1/2 AF XXZ chain in a longitudinal field (spectra equal to 5.5e-14). The full bond commutes exactly with T_μ exp(−i(a − π)S^z) in the quantum frame (4e-15).
- [proved, correct] §2.1 rigid-cone linearity, for every q ∉ {0, π}.
- [derived, wrong as stated] §2.2 "DM gives no quartic term". DM gives a concave u² coefficient −D² sin²(q/2), as large as or larger than Γ, and Γ_xz interferes with D_x. Corrected: N⁴ coefficient = c_2 − c_1 with c_2 = (g² + g_xy²) cot²q/4 and c_1 = sin²(q/2)|tan(q/2)Γ_c − iD_c|², where Γ_c = Γ_xz − iΓ_yz and D_c = D_y + iD_x. Checked to all digits on random harmonics, and by nonlinear minimization (D_x = 0.25 removes the g = 0.2 window, as predicted). Qualitative rule kept: charge 2 is the only convex source.
- [derived, correct] §4.1 at h_eff = 0, Δ_eff = 1: one SU(2)_1 primary contains the staggered and spiral orders, and the coexistence (spin-4) operator has dimension ≥ 8.
- [derived, correct and strengthened] §4.2 at the KT point. The selection rule P − Qα ∈ 2πZ is exact. Momentum-π neutral operators (N^z, dimerization) are forbidden at any order and for any α. Charged relevant or marginal operators need q ∈ {0, π}. Answer to the spec's question: no relevant operator at momentum π is ever generated. At momentum 0 the uniform field is generated at second order (hx·Γ, hx·D; checked by ED, quadratic scaling). It is the parent's tuning parameter, so the SU(2)_1 point sits at h_eff = 0, not at hz = 0.
- [derived, gap closed] §4.2 on the spin-flop line. The quantum switch at hz ≠ 0 is the Néel plateau edge hz_c1(Δ), a PT transition, not SU(2)_1, so the dimension-8 argument does not apply there. Replacement argument: an off-resonance U(1)-symmetric effective chain, one boson, and coexistence needing a bound magnon, which weak harmonics cannot create because the parent's spinons have S(0) = −1 and no zero-energy resonance.
- [derived, wrong] "Quantum shift of the Néel onset is O(1)". At hz = 0 it is O(harmonics²). At hz ≠ 0 the O(1) separation comes from the field (hz = 0.05: classical Δ = 1.0003, quantum 1.321).
- [numerical, wrong] "Gapped Néel up to hz = 2.2 at Dz = 2". The parent edge is hz_c1 = 0.78. The VUMPS energy at hz = 2.2 is 0.1 above two independent ED estimates of the ground state, so the VUMPS state is metastable.
- [numerical, not decisive] Line Q2 at χ ≤ 40 cannot place a KT onset, and its Dz = 1.25 and 1.3 points lie on the parent's liquid side of the PT edge (Dz = 1.32 at hz = 0.05).

## Consequences

- REPORT §9.124 and C165 RESULT and PROOF §6 should drop "DM and fields do not enter" and use the corrected c_1.
- A_S should be stated along the quantum switch set (the KT point and the plateau edge h_eff = hz_c1(Δ_eff) for Δ_eff > 1) as the absence of a second low-energy sector (bound magnon) at order-one amplitude. "Parent Δ near 1, fields up to the plateau edge" describes an exponentially small field range.
- Job 17601089, family A2: its fields are mostly above the parent plateau edges (0.27 at Dz = 1.6, 0.78 at Dz = 2.0). Any A2 "gapped Néel" reading above the edge must be compared in energy with a magnetized seed before it is used. The edges themselves are sampled only by the lowest one or two fields.

## Open

- The bound-magnon exclusion at the plateau edge is derived for weak harmonics only. At order-one amplitude it is exactly what A_S assumes, and C172's continuity route would need it.
- §4.4 and §4.5 were checked only for consistency, not refereed.
