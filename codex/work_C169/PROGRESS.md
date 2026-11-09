This run was done by Claude Opus 5.5 (claude-opus-5-5), not by Codex. Astra was not used (out of credits).

# C169 log

2026-11-08
- Read specs/C169_C174_swarm.md, specs/_qubit_c_context.md, STATUS.md, REPORT.md §9.110 to §9.126, work_C154/RESULT.md, work_C165/RESULT.md, work_C157/RESULT.md, work_C156/RESULT.md, work_C163/PROOF.md §4 and §5.
- Framing: c is invariant under every exact map considered, so the hunt reduces to recognizing a known c > 1 model inside α^{-1}(B_2), or finding a separating invariant per class.
- frustration.py: offset-diameter and component invariants for B_2 and the targets. The first component count used rings, where wrap-around merged cluster-Ising's three chains. Switched to open chains: B_2 maximum 2 over 4095 subsets, targets 3 to 4.
- freecut.py: free chains of any period. The first version sampled on radius 1.7 (noisy negative-power coefficients) and had the XX sign wrong. Fixed both: span ≤ 4, at most 4 unit roots in 3800 chains, rank-1 bond gives span 2.
- Coordinator note (C174 R5): Z_5 parafermion at a melting tongue. Literature search (Rydberg period-p melting). Derived: chirality is never symmetry-forbidden for translation-Z_m, the Lagrange bound, and the spin-1/2 specifics. An early draft wrongly had reflection forbidding chirality. The kink argument shows reflection maps kinks to kinks (it acts as CP), so this was corrected before write-up.

2026-11-09
- Coordinator note (C171 R1, R4): Potts × half-floating boson and the repulsive-host c = 2. Added PROOF §9.2 and §9.3.
- claw.py: claw-free Pauli-aligned B_2 subsets, 20 classes, all ordinary JW-free.
- interface.py: qubit-interface test at b = 2, 3. Every target fails, and the controls pass. The "Ising × Ising ladder" line is not a control for passing (it is not blocked B_2 either).
- Literature: Yamaguchi–Chiba–Shiraishi classification (symmetric nearest-neighbour spin-1/2, with fields) used in PROOF §4.
- No Typhon job, no DMRG. Wrote PROOF.md and RESULT.md.
