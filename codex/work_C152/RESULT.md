This run was done by Claude Opus agents (2026-11-07), not by Codex. It was interrupted and then resumed by a second agent, which continued from PROGRESS.md and the kept Typhon job.

# C152 result: c ≤ 1 on the parity-symmetric qubit class

Derivations are in PROOF.md (sections 1 to 7 from the first agent, 8 to 11 from the resumed run), the log is in PROGRESS.md. Scripts in this folder: model.py, check_free.py, ring_ceff.py, period2_det.py, period2_mf.py, period2_search.py, dmrg.py (C144, with optional on-site terms), dmrg_pt.py, dual_check.py, free_open.py, summ_dmrg.py, fss.py, xzj_table.py, make_tasks2.py, typhon_scan.py. Typhon inputs in typhon/, copied outputs in typhon_out/. Tags: [proved], [derived] (physics standard), [literature], [numerical], [heuristic].

## Verdict

No point of the class with c > 1 was found or derived. c ≤ 1 is proved on the free locus and for every period ≤ 2 mean-field state, exact on the U(1) and zero-field XYZ loci, and derived in neighbourhoods of all of these, on and near the Ashkin-Teller plane, at the Ising-chiral junctions away from XZ, and at strong repulsion. Strong attraction rests on a cluster-energy heuristic. The only places where the open-chain numerics read above 1 are the XZ window at h_z ≈ 0 (φ = 1.45, L = 64 fits 1.13 to 1.55), where the free chain reads 1.0, the derivation gives one boson, and C158 has shown that the same region produces crossover readings above 1 that converge to 1. The finite-size scaling there (L = 32, 64, 128, 192) is queued in Typhon job 17597010. As input A4 of the C154 conditional theorem, c ≤ 1 on the parity class holds at the derived standard away from the XZ window and at intermediate couplings far from all controlled regions. Those intermediate couplings rest on numerics.

## The class

[proved] Bonds commuting with Π = ∏σ^z are, after a z rotation, J_x XX + J_y YY + D(XY − YX) + J_z ZZ + h_z Z (4 parameters modulo scale). Jordan-Wigner gives the interacting Kitaev chain w c†_{j+1}c_j + h.c. − Δ(c_j c_{j+1} + h.c.) + V n_j n_{j+1} − μ n_j with w = J_x + J_y + 2iD = t e^{iφ}, Δ = J_x − J_y, V = 4J_z, μ = 4J_z + 2h_z (checked against BdG to 1e-5). The XZ (compass) point of C143/C144 lies in the class (parity axis along y) at J_x = −J_y = D, J_z = h_z = 0, i.e. cos φ = 0, Δ = t, μ = V = 0.

## Gapless loci and their central charges (PROOF section 10)

| locus | dimension in the class | c | status |
|---|---|---|---|
| G1: chiral Luttinger liquid (v_L ≠ v_R), including the U(1) liquid and, at D = 0, the floating phase (1/4 < K < 1/2) | open (4) | 1 | proved at V = 0, exact at Δ = 0, derived elsewhere, numerical at 9.2, 9.3 points |
| Ising surfaces (topological/trivial; CDW/disordered; z-FM/disordered/x-FM pair at attraction) | 3 | 1/2 | proved at V = 0, derived (mass argument), literature (Chepiga-Mila, Peotta et al.), numerical |
| Pokrovsky-Talapov and band-edge surfaces bounding G1 (filling, CDW edge, pairing threshold Δ_c(φ, K), Δ_c = t sin φ at V = 0) | 3 | z = 2 boundary of a c = 1 phase | proved at V = 0, derived |
| KT boundaries of G1 inside the symmetric hyperplanes D = 0, h_z = 0 and the Ashkin-Teller plane | 2 | 1 | literature (D = 0), derived |
| First-order surfaces (Luttinger liquid to ferromagnet at attraction; phase separation) | 3 | none | literature (Peotta et al.), heuristic |
| Ashkin-Teller locus: crossing of two Ising surfaces on cos φ = 0 = h_z, s > 1 | 2 | 1 | derived (PROOF 8.3, 8.4) |
| Zero-field XYZ critical lines, incl. the 8-vertex point where KT and PT boundaries meet | 1 | 1 | exact (Baxter), literature (Chepiga-Mila) |
| Junctions J: Ising surface meets the boundary of G1 (one-sided cubic onset) | 2 | ≤ 1 | proved at V = 0 and at mean-field level, derived beyond away from XZ |
| Tricritical Ising points where Ising surfaces turn first order (attraction) | ≤ 2 | 7/10 | possible by analogy with interacting Majorana chains; not located here (Set 3 of job 17597010) |
| XZ point and its window | 0, window O(V) | ≤ 1 emergent | C143, C144, C154, C158; h_z ≈ 0 window pending here |
| Ferromagnetic Heisenberg point (endpoint of G1 at attraction) | 0 | z = 2 | exact |

The frustration-free (Peschel-Emery, Katsura-Schuricht-Takahashi) surface is a disorder surface inside the gapped Z2 phase, not critical.

## Main results

1. [proved] Free locus (Prop. F) and Lemma P2: for any nearest-neighbour BdG chain with one orbital per site and period ≤ 2, det H_BdG(K) is a quadratic polynomial in cos K, so at most four Majorana-Weyl zeros and c ≤ 1 (500-model numerical check of the Fourier content to 4e-16; 20,000 random and 60 adversarial models give max c_free = 1). Hence every Hartree-Fock-Bogoliubov state of the class, uniform or period 2, has c_MF ≤ 1. This is where the class differs from the interacting Majorana chains of Rahmani et al., whose c = 3/2 phase already appears at mean-field level through a three-site term.

2. [proved] Ashkin-Teller plane and duality (PROOF 8.1, 8.2). At cos φ = 0, h_z = 0 the Majoranas A = (a + b)/√2, B = (a − b)/√2 decouple at the quadratic level: H = Σ[−t_A ε^A_j − t_B ε^B_j − (V/4) ε^A_j ε^B_j] with t_{A,B} = (t ± Δ)/2, the self-dual Ashkin-Teller chain with unequal couplings. B_j → (−1)^j B_j gives the exact duality (t, Δ, V) → (Δ, t, −V) (ring spectra agree to 2e-14 at L = 8, dual_check.py). This is C144's plane XZ + α ZX + β YY, and the duality makes C144's half-plane cover all of it.

3. [derived] Two-Majorana picture (PROOF 8.3). For s = Δ/(t sin φ) < 1 the two Majoranas form one Dirac fermion with unequal velocities. The velocity anisotropy has dimension K + 1/K (irrelevant for K ≠ 1, exactly marginal at the free point), the single-species masses are forbidden by translation (Kramers-Wannier on each species), and the hopping phase and field enter as the current and density of the Dirac fermion, so G1 is one boson with c = 1. For s > 1 the same bilinears are the masses μ ± 2t cos φ of the two Majorana combinations, which gives the Ising surfaces (c = 1/2) and their codimension-2 crossing, the Ashkin-Teller line (c = 1), not Ising × boson. This continues Prop. F to V ≠ 0 and explains the DM-induced splitting of a c = 1 line into two Ising lines seen by Peotta et al.

4. [derived] The spec's multicritical worry (Ising line meeting a Luttinger region, potential c = 3/2) resolves in this class into (i) Ising-Ising crossings, which are Ashkin-Teller c = 1, and (ii) Ising-chiral junctions, which are one-sided cubic Lifshitz points with c ≤ 1, protected beyond mean field by the margin c_3 = t² cos² φ/(4Δ) except within O(|V|) of XZ. A c = 3/2 point needs six low-energy Majorana-Weyl components from two Majorana species per site. The quadratic theory allows at most four, so only a strongly coupled slow sector, which exists only in the XZ window, could supply the rest.

5. [derived] Around the U(1) plane (bosonization, PROOF 4), at strong repulsion (hard-core rods, CDW domain walls; one species) and at strong attraction (cluster energies (n − 1)V favour phase separation over a pair-plus-particle two-fluid phase; [heuristic] for Δ ≠ 0, Bethe ansatz for Δ = 0): one species, c ≤ 1.

## Numerics [numerical]

- Calibration at L = 64 (open chains, central-window fits): 8-vertex point 0.98, free chiral 0.97 (exact free profile 0.96), free Ising 0.56 (open-boundary excess), free near-XZ chains 0.99 to 1.01 at L = 64 and 0.99 at L = 256.
- Ring-ED candidates at φ = 0.9 (c_eff 2 to 6 at L = 12 → 16) are c = 1 by DMRG at L = 64, 96 (fits 0.86 to 1.13, wide window 0.96 to 1.00 at L = 96): level-crossing artifacts.
- Near-XZ family (264 runs at L = 64): repulsive lines c = 1, attractive lines gapped with one transition each (Ising fits 0.2 to 0.53 at V = −0.8; 0.81 to 0.90 at V = −0.4, a narrow c = 1 sliver or crossover, under refinement). XZ window at small h_z and φ = 1.45: fits 1.13 to 1.55 (largest at s = 1.1, V = 0.4, h_z = 0), unresolved at L = 64 (PROOF 8.5, 9.3).
- Small-ring slices (V = 6): Luttinger region 0.7 to 1.4 (shell oscillations), CDW 0, isolated spikes at single level crossings, no extended region above 1.4.

## Candidates with c > 1

None established. The only open candidates are the XZ-window points above (φ = 1.45, h_z ≤ 0.05, s = 0.9 to 1.1, |V| = 0.4 to 0.8), with L = 64 fits up to 1.55. They match the crossover pattern already diagnosed by C158 (a slow strongly coupled sector whose entanglement has not saturated), and the free chain at the same points reads 1.0. They are counted as unresolved, not as c > 1, until job 17597010 returns.

## Pending Typhon job

- 17597010 (c152_j2), one 48-core node, 14 h limit, 267 DMRG runs; PENDING (Priority) at write-up; the scheduler estimate moved from 2026-11-08 09:06 at submission to 2026-11-12 13:59. Inputs typhon/c152_j2.slurm, typhon/run_d2.sh, typhon/dtasks2.txt (make_tasks2.py). Outputs /scratch/lustre/jkudlerflam/c152/j2/dmrg_{32,48,64,96,128,192}.jsonl, dmrg_err.log, slurm_17597010.out. Contents: Set 1, h_z = 0 near XZ (φ = 1.45 and π/2, s = 0.9, 1.0, and 1.1 at 1.45, V = ±0.4, ±0.8) at L = 32, 64, 128; Set 2, fine h_z scans at φ = 1.3, V = −0.4, s = 1.0, 1.1; Set 3, strong attraction V = −1.5, Δ = 0.5, φ = 0, 0.8 (tricritical Ising and pair-particle coexistence); Set 4, the unresolved meeting region of Peotta et al. (φ = 1, h_z = 0, V = −2.3 to −1.6, Δ = 0.25 to 0.4, L = 48, 96); Set 5, XZ-window points at h_z = 0.05, 0.1, the φ = 1.45 attractive transition, and L = 192 at the two largest h_z = 0 fits. Analysis: mkdir typhon_out/j2 && scp 'typhon:/scratch/lustre/jkudlerflam/c152/j2/dmrg_*.jsonl' typhon_out/j2/, then python3 fss.py typhon_out/j2/dmrg_*.jsonl. The decisive numbers are c(L, 2L) = 6[S̄(2L) − S̄(L)]/ln 2 at the XZ-window points.
- 17594701 (c152_dmrg) COMPLETED (33 min, all 283 runs, empty error log). Outputs /scratch/lustre/jkudlerflam/c152/dmrg_{64,96}.jsonl, copied to typhon_out/.
- 17594696 was cancelled; only its extra.jsonl slices exist (typhon_out/extra.jsonl).

## References used

Chepiga and Mila, PRB 107, L081106 (2023), arXiv:2206.11754. Peotta, Mazza, Vicari, Polini, Fazio, Rossini, J. Stat. Mech. (2014) P09005, arXiv:1403.4568. Mahyaeh and Ardonne, PRB 101, 085125 (2020), arXiv:1911.03156. Rahmani, Zhu, Franz, Affleck, PRL 115, 166401 and PRB 92, 235123 (2015). Kohmoto, den Nijs and Kadanoff, PRB 24, 5229 (1981). Miao et al., PRL 118, 267701 (2017). Katsura, Schuricht and Takahashi, PRB 92, 115137 (2015). C141, C143, C144, C154, C158 of this campaign.
