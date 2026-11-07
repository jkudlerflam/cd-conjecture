This run was done by a Claude Opus agent (2026-11-07), not by Codex.

# C144 result: the XZ point and its tuned Ashkin-Teller cone

Derivations are in PROOF.md, the log is in PROGRESS.md. Scripts (NumPy/SciPy/numba, all in this folder): xzlib.py, kernel.py, bmodel.py, validate_bmodel.py, bscan.py, bdiag.py, first_order.py, ff_open.py, ff_mf.py, qfree.py, qscan.py, dmrg.py, test_dmrg.py, dmrg_run.py, cfit.py, corr_run.py, yx_ed.py, line_ed.py. Typhon inputs are in typhon/, copied outputs in typhon_out/.

## Answer to (Q1): is Q_L bounded near XZ?

Yes, at the physics level, uniformly in the direction of approach, and quantitatively close to the XZ value.

- [proved] Closed form at XZ on even rings: Q_L(XZ) = 2L tan(pi/(4L))/cos(pi/(2L)), decreasing from 2.3431 (L = 2) to pi/2. Odd rings give negative values.
- [proved, C131] Q_L(h) <= Q_L(XZ) + 4L^2 delta for all L, with delta = M[h - XZ]. This covers L <= (kappa/4 delta)^{1/2}.
- [derived, conditional] Beyond that, Q_L = L(D^A_L + D^B_L)/M with an A sector of velocity 4(1 + O(delta)) (or gapped) and a B sector of energy scale O(delta). The bound is sup_L Q_L <= max{2.343 + kappa, (pi/2)(1 + C delta) + (pi/8) v_B^max/M + O(kappa^{-1})}, with v_B^max = O(delta). Assumptions: (S1) the A/B separation of scales (derived in section 1 and 3.1); (S2) no Casimir overshoot inside a sector, assumed, and checked for massive free Majoranas, whose Q_L only falls below the critical value; (S3) c_B bounded, with the B energy-density drift lowering Q, not raising it (derived; the effective B chain at the tuned point has L D^B_12 = -21 beta^2). Q1 does not need c_B small, only bounded.
- [numerical] ED of Q_L for L = 2 to 10 (rings up to 20 sites; Typhon job 17580166, 86 runs over 22 direction families: 8 normal directions, the 2 rotations, 6 random directions, and the cone with and without mixing, B_x and B_y added, at amplitudes 0.05 to 0.6). The largest value anywhere is 2.3476 at L = 2 (random direction r4, s = 0.05), against Q_2(XZ) = 2.3431. Over L >= 4 the maximum is 1.847 (YY, s = 0.4). Gapped directions decay to 0 (B_x, s = 0.1: 2.24, 1.15, 0.61, 0.33, 0.18 at L = 2 to 10). The cone gives 2.17, 1.71, 1.59, 1.58, 1.53 at s = 0.05, and smaller values at larger s.

## Answer to (Q2): c_B on the tuned cone

c_B <= 1/2 on the whole Ashkin-Teller tuned cone. The first-order cone sits inside a gapped wedge (c_tot = 0), and a single line in the wedge carries c_tot = c_B = 1/2 with A gapped. No c_tot > 1.

- [proved] On the cone the second-order B Hamiltonian has no hopping at all. The kernel K_jl = <da_j (H_A - E_0)^{-1} da_l> is real symmetric, so the b_j b_{j+1} terms cancel. H_B = t sum b_j - 2 beta^2 sum_{r >= 2} K(r) b_j b_{j+r}, a purely quartic Majorana chain at t = 0.
- [derived and checked] The tail is K(r) = (-1)^r/(2 pi r) + O(r^{-3}). Its staggered sum is (1/pi) log L - 0.0577, so the kernel is marginal (an eps_B eps_B/|x - y| interaction).
- [numerical] This effective chain reproduces the full qubit ring at beta = 0.004. Its own ED to L = 40 at t = 0 gives a twofold ground state with <b_j> = +-(2/pi)(-1)^j. Mean field explains this: B self-generates a staggered-sign hopping and becomes a critical Majorana chain in a period-4 gauge. Its velocity grows like log L. Any t != 0 gaps it.
- [derived] Feedback on A: the staggered <b_j> dimerizes A, gaps it, and cuts the logarithm at log(1/beta).
- [numerical] DMRG of the full qubit chain (open chains, own NumPy DMRG, validated against ED and exact free fermions to 1e-8 or better). At L = 128 and beta = 0.15, 0.3, 0.5 the averaged c fits are about 1 for alpha/beta <= 0.3 and >= 1.0. They fall below 0.1 in a wedge around alpha/beta = 0.55 to 0.8. Its width in alpha is 1.17, 1.02, 0.98 times beta^2 log(1/beta). Inside the wedge A is strongly dimerized (a_stag up to -0.35) and b alternates with amplitude near 2/pi. One line inside it (b_mid = 0, alpha/beta = 0.64 to 0.68) is critical. Two observations support c = 1/2 there. At L = 64 the a-a correlator decays exponentially (xi about 1.7) and b-b decays as r^{-1.7 to -1.9}, the Ising energy correlator, so only one Majorana is critical. At beta = 0.5, L = 128 the parity fits give 0.46 / 0.44. Single-chain profile fits on the beta = 0.3 line range from 0.3 to 1.0, because the B mass is pinned differently near the open ends. On the ordered-B side, entropies near the line also carry a spurious log 2 from a cat of the B edge zero modes: S(L/2) sits at 0.87 for alpha = 0.199 to 0.203 and drops to 0.17 at 0.205 (L = 64). So entropy-based estimates there are unreliable, and the correlation structure is the cleaner evidence. The L = 128 to 256 finite-size scaling of this line is still running (see below). No fit anywhere exceeds 1.09.
- Convergence: truncation <= 1e-9 and energy changes <= 1e-9 in the last sweep for the L = 128 grid. Two points have dE about 1e-5 to 1e-6 and are not used for conclusions: beta = 0.15 at alpha/beta = 1.0, and the L = 64 point alpha = 0.203. The beta^2 log(1/beta) width of the wedge and the c = 1/2 line come from three values of beta only. The line's c needs the pending L = 128 to 256 runs to be called converged.

## The full tuned locus in B_2

- [proved] At XZ the first-order effective Hamiltonian on the degenerate manifold, for a general tangent vector, is H_eff^(1) = dJ_xz E_0^A + (dJ_zx + <a> dJ_yy) sum_j b_j, with <a> = E_0^A/L -> -2/pi. All nine other directions project to zero; this was checked to 1e-13 at L = 10.
- So the tuned locus is the hyperplane T = {dJ_zx + <a> dJ_yy = 0}, dJ_zx = (2/pi) dJ_yy as L -> infinity. It has codimension one (7-dimensional) in the 8-dimensional normal space. The rest of the normal space is (XX + ZZ, B_y, YX, ZY, B_x, B_z).
- The exact AT tuned curve (t_eff = 0) is shifted to alpha/beta = 0.64 to 0.68 at beta = 0.15 to 0.5.
- Directions that break the A/B factorization:
  - XX + ZZ and B_y hybridize the Majoranas. They stay free and were solved exactly: they gap A at first order and leave a z = 2 B touching on T, of Pokrovsky-Talapov type, single species, with c in {0, 1} nearby.
  - B_x and B_z are Jordan-Wigner strings of the form sigma_A times a classical B string. They are strongly relevant and gap everything. Numerically S_mid = 0.26 and Q_L -> 0.
  - YX and ZY are also P-odd, but their B parts do not commute. [numerical, not converged] They gap A and leave a slow sector that looks like c of about 1 at L = 64, with a period-3 commensuration on rings. L = 128 and 256 runs are pending.
  - XY and YZ are rotations.

## Task 4: c_tot > 1 at a qubit bond? (EM-XZ) and (NC)

- [numerical, with derivation] The tuned cone of XZ does not host c_tot > 1. On the AT plane c_tot is 1 outside the wedge, 0 inside it, and 1/2 on the line. Elsewhere on T the A sector is gapped by the factorization-breaking component.
- (EM-XZ) holds, with three amendments. First, v_B = O(delta^2 log(1/delta)) on T. Second, A can be gapped, by hybridization or by spontaneous dimerization. Third, c_B is bounded in every direction: <= 1/2 on the AT plane, and about 1 at most on YX/ZY, with A gapped there.
- (NC) as stated (c_B <= 1/2 on the tuned cone) holds on the AT plane. It should be restated as c_tot <= 1, because on YX/ZY the emergent sector may carry c of about 1 while A is gapped. With that restatement it holds everywhere tested.

## Pending Typhon jobs

- 17580192 (c144_d2). Tuned-line scaling at beta = 0.3: alpha = 0.194 to 0.1985 at L = 128, 192, 256 (the L = 64 part is done). Also YX (L = 128, 256), ZY (L = 128), and cone + YX, cone + (XX + ZZ), cone + B_x (L = 128). Outputs: /scratch/lustre/jkudlerflam/c144/t_L*_a*.jsonl, yx_L*_e*.jsonl, zy_L128_e0.1.jsonl, cone{yx,mix,Bx}_L128.jsonl, slurm_17580192.out.
- 17580432 (c144_d3). alpha = 0.199 to 0.205 at L = 128, 192, 256, where the L = 64 to 96 entropy maximum lies. Outputs: /scratch/lustre/jkudlerflam/c144/u_L*_a*.jsonl, slurm_17580432.out.
- Analysis for both: python3 cfit.py <files>. The decisive quantity is max over alpha of S(L/2) at each L, whose slope in log L is c/6.

## References

Kohmoto, den Nijs and Kadanoff, PRB 24, 5229 (1981) (Ashkin-Teller chain). Rahmani, Zhu, Franz and Affleck, PRL 115, 166401 (2015) and PRB 92, 235123 (2015) (interacting Majorana chains). Brzezicki, Dziarmaga and Oles, PRB 75, 134415 (2007), and arXiv:0903.1682 (exact solution of the one-dimensional compass model, macroscopic degeneracy at the isotropic point; found by search, not re-read in full). Calabrese and Cardy, J. Stat. Mech. P06002 (2004) (open-chain entanglement). C131 cover lemma, C141 sections 5 and 6, C143 sections 2.2, 4.2 and 5.
