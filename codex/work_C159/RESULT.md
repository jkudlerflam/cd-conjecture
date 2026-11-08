This run was done by a Claude Opus agent (2026-11-07), not by Codex.

# C159 result: c <= 1 on the time-reversal strata T and K

Derivations are in PROOF.md, the log in PROGRESS.md. Scripts in this folder: tk.py (normal forms, spiral parent, symmetry test), scan_tk.py, refine_tk.py, summarize_tk.py, halffloat_test.py, khg_check.py, vumps_tk.sh, typhon/c159.slurm, with C155's est.py, dmrg.py, dmrg_c.py, C153's vumps.py and C158's fesx.py, fesfit_c158.py, all unmodified. Local pilot data in scan_local/.

## Verdict

No c > 1 was found or derived on T or K, and the two strata turn out to be more constrained than the symmetry-free stratum N. Two new structural facts do the work. (1) An exact antiunitary symmetry commuting with translation makes one of the two fields of any floating boson commensurate ("half-floating", Lemma A). (2) On T, the Kramers LSM anomaly rules out every c < 1 critical point that preserves the symmetry (Lemma B), so symmetric critical points of T have c >= 1, and Ising points occur only inside a spontaneously broken background. With these, every c > 1 candidate that could fill a hypersurface or an open set of T or K must contain an attractive Luttinger sector (Proposition C). SU(2)_2 and SU(2)_3 cannot be reached in T. The upper bound c <= 1 is derived on controlled neighbourhoods (spiral frame at zero field, AF Heisenberg, cyclic SU(2)_1, ferromagnetic onsets, free locus, literature families), and on all of T and K it is conditional on one assumption narrower than C154's (A6). The numerical tests are queued on Typhon and have not run.

## Task 1: stratum T (B = 0, Kramers Theta only, dimension 5)

Allowed infrared theories [derived, PROOF 2]:
- Symmetric gapless phases carrying the anomaly. A floating liquid in T is half-floating. The spin density along the spiral axis has a commensurate k_F = pi/2 (Theta forces zero magnetization along every axis), and only the transverse angle is incommensurate. Its only allowed relevant operator is the umklapp cos 4phi, so it is stable if and only if K > 1/2. It ends by a KT transition into Neel order along the spiral axis (Theta broken) or into a dimerized phase (translation broken), with a Gaussian c = 1 line between the two. The SU(2)_1 regions (AF Heisenberg, cyclic stratum) belong here.
- Spontaneously broken phases. These are Theta-broken magnetic phases (ferromagnet, Neel, canted, commensurately locked spirals) and translation-broken ones (dimerized, period p).
- Transitions. At a symmetric critical point, c >= 1 (Lemma B). Ising (c = 1/2) points occur only between two phases that both break Theta or both break translation. There are also commensurate-incommensurate (z = 2) and first-order transitions. Singular points (the Kitaev point XZ, the ferromagnetic Heisenberg point) are outside the CFT class.
- [proved] Lemma B on the free locus cap T: every gapless point has c = 1 and never 1/2, from the band.
- [derived] Lemma B is consistent with every known point. C144's c = 1/2 line near XZ sits inside a dimerized A background. In Peotta et al. at h_z = 0, the Ising lines separate phases that break Theta. The c = 1/2 phase of Luo et al. is in a bond-alternating chain outside B_2, which has no Kramers LSM.

c > 1 in T [derived, PROOF 2.4 and 4]:
- SU(2)_2 does not carry the anomaly (k even), and its nine spin-1 operators would need 9 > 5 tunings. SU(2)_3 carries the anomaly but also needs 9 tunings, except on a curve of the cyclic stratum, where C154 found nothing. SU(2)_1 x Ising needs at least 2 tunings.
- Ising x boson with period-2 order. The commensurate vertex e^{2 i phi} at momentum pi neutralizes sigma, so sigma e^{2 i phi} (dimension 1/8 + K) is allowed. The candidate needs 2 tunings unless K > 15/8, and its ordered neighbour exists only if K > 2. In N the same candidate needed only one tuning at every K.
- Ising x boson with uniform Theta-odd order needs 2 tunings (epsilon and sigma d_x phi).
- Two half-floating bosons. Since the anomaly is Z_2-valued, exactly one boson carries it. This candidate needs no tuning only if K1 + K2 > 2 and both K > 1/2.
- So an anomaly-matching c > 1 theory is not excluded by symmetry. On loci of dimension >= 4 it can only sit next to attractive liquids: K > 1, which for a zero-field parent means Delta < 0, and for Ising x boson K > 15/8, that is Delta < -0.67. Such liquids lie on the ferromagnetic side of the XY regime, including the easy-plane neighbourhood of the ferromagnetic Heisenberg point. That is where the Typhon scan TFM looks.

c <= 1 in T. [derived] It holds near the U(1) and XYZ strata, because the spiral parent of every T bond is at zero field and the harmonics resonate only at spiral angles 0, pi/2 and pi. It holds near the AF Heisenberg point, where SU(2)_1 plus DM leaves c in {0, 1} and reproduces Baxter at D = 0. It holds near the cyclic SU(2)_1 region, where the translation action -P forbids every relevant operator. It holds at ferromagnetic onsets (one soft species, C154 Lemma 3) and on the free locus cap T. [literature] It holds on the Kitaev-Gamma (with and without DM) and Kitaev-Heisenberg-Gamma families. khg_check.py shows that KHG is a 2-dimensional family in the interior of T (spectra agree to 2e-14, only Theta survives), with four c = 1 liquids and Theta-broken ordered phases (Yang, Nocera and Affleck 2020).

## Task 2: stratum K (one antiunitary A = Theta R_pi(z), B in the plane, dimension 5)

- [proved] Normal form (Lemma 1): K is XYZ + D_z(XY - YX) + a field (B_x, B_y, 0), which is C155's RA coordinates. It has two sheets. On K_I (J_x J_y + D_z^2 > 0) the spiral axis is z and the parent is a zero-field XXZ chain. On K_II (J_x J_y + D_z^2 < 0) the spiral is a staggered pi rotation about an in-plane axis u, and the parent is XXZ in the field u.B.
- [derived] Allowed infrared theories. There is no anomaly, so unique gapped phases exist (field-polarized). Ising transitions are codimension 1 when the order parameter is odd under A (moment along z, uniform or staggered) or translation-odd. There are half-floating c = 1 liquids. On K_I a liquid connected to the parent has a commensurate density and is stable for K > 1/2, as in T. On K_II the transverse correlations are commensurate and the magnetization along u is incommensurate. Transverse anisotropy then gaps the liquid unless K < 1/2, which is the XYZ-in-a-field picture of Chepiga and Mila. KT, commensurate-incommensurate and Gaussian transitions also occur.
- [derived] c > 1 in K. The counting is the same as in T without the anomaly: two half-floating bosons with x1 + x2 > 2, or Ising (sigma at pi) x boson with x > 15/8, in terms of the dimension x of the commensurate vertex. Everything else needs at least 2 tunings.
- [derived and literature] c <= 1 holds near U(1) cap K and in the spiral regime of K_I. In K_II it holds in the field-induced-liquid regime. It also holds at ferromagnetic and field-polarized onsets (Lemma 3), on the Garate-Affleck and Chan et al. family (XXZ + D_z + transverse field: one c = 1 liquid, two Ising-ordered phases), on tilted-field XXZ, and on the Kitaev chain in a [110] field (rank-one boundary of K, c = 1.01 +- 0.03 from C158).

## Task 3: numerics

- [numerical, local, seconds per run] Ring structure factor at L = 12: the spin correlations along the spiral axis peak at k = pi at the C68 T and K anchors, and away from pi at the N anchor. Every random T or K pilot point has exactly one nontrivial uniform symmetry. A pilot near the ferromagnet in T found an attractive point (parent K = 2.35) that reads c_S = 1.48 at rings 12/16 and 0.71 at rings 14/18, an incommensurate shell effect. With a one- or three-site cell, VUMPS cannot converge to a Theta-symmetric state in T (the MPS form of LSM), so two-site cells are used.
- [pending] Typhon job 17597033 runs the scans of T, K and the ferromagnetic and Heisenberg neighbourhoods (about 64,000 points at rings 12/16), FES at the C68 anchors, and refinement by rings 16/20 and DMRG at L = 32, 64, 128. Section "Pending" below lists its outputs.

## Task 4: theorem-style summary

- Lemma 1 [proved]: normal form and sheets of K.
- Lemma A [derived]: an exact antiunitary commuting with translation makes the field it reverses commensurate, with a shift of 0 or half a period.
- Lemma B [derived for Ising, tricritical Ising and Potts; assumed for the other minimal models via Ruelle-Verhoeven and the anomaly-freeness of their symmetries]: no symmetric c < 1 critical point in T.
- Proposition C [derived, kinematic]: in T and K, a c > 1 candidate of codimension <= 1 needs an attractive half-floating sector, x > 1 (two bosons) or x > 15/8 (Ising x boson). SU(2)_2 and SU(2)_3 are excluded in the interior of T.
- Theorem [derived, conditional on A1 to A5' of C154, the pending parity class (C152) and the XZ hyperplane (C157), and H_T, H_K]: every conformal point of T and K has c <= 1. In T every conformal point that breaks no symmetry has c = 1 exactly, and c = 1/2 occurs only inside a broken background. H_T and H_K can be replaced by (A6_TK) [assumed]: no point of T or K hosts two independent gapless sectors one of which is an attractive half-floating boson.

## Open regions

- No unexplained critical phase was found. The region still uncontrolled is the attractive (K > 1) half-floating liquids of T and K_I far from the U(1) stratum, in particular the easy-plane, DM-twisted neighbourhood of the ferromagnetic Heisenberg point. That is the only place where Proposition C leaves a codimension-0 or codimension-1 window for c > 1. No candidate there survives the pilot.
- K_II liquids (K < 1/2) and the K_I branch with zero commensurate shift (stable only for K > 2) have not been mapped.

## Pending Typhon job

- 17597033 (c159, partition typhon, one node, 48 cores, 8 h limit). Status PENDING (Priority), scheduler estimate 2026-11-12 13:59. Directory /scratch/lustre/jkudlerflam/c159/.
- Outputs:
  - scan/{T,K,TFM,KFM,TAF}_*.jsonl (stage A)
  - fes/{T,K,N}_C68.json (stage V)
  - refine/*.json (stage B)
  - summary.txt, written by python3 summarize_tk.py at the end of the job
  - logs/c159_17597033.out and logs/*.log
- What to read:
  - In summary.txt, the count of "symmetric Ising-like" points in T against K, which tests Lemma B.
  - c_obj by parent Delta, which tests Proposition C.
  - The DMRG c_two and profile fits at L = 128 for every candidate.
- Criterion fixed in advance: a c > 1 candidate needs fits above 1.2 at L = 128 with truncation below 1e-7, and then FES at chi >= 256.
