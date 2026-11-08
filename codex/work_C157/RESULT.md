This run was done by a Claude Opus agent (2026-11-07), not by Codex.

# C157 result: no c > 1 from the emergent Majorana chain at XZ

Derivations are in PROOF.md and the log is in PROGRESS.md. All scripts are in this folder: kernel_q.py, hfk.py, hf_bmodel.py, bsmall.py, block_check.py, bis_ed.py, corr_bis.py, cscale.py, casimir.py, rzfa_dmrg.py, dmrg_d.py. C144's xzlib.py, kernel.py, bmodel.py, dmrg.py, dmrg_run.py and cfit.py are copied unmodified. The Typhon inputs are in typhon/. C156's Typhon outputs (read only) are copied to c156_copy/.

## Verdict

No counterexample. Within second-order perturbation theory around XZ, and in every direction tested numerically, no nearest-neighbour qubit bond near XZ reaches c_tot > 1. The XZ loophole, a second Majorana species B, is real, but it closes for a structural reason. The only self-interaction of B is induced by virtual A excitations, and it has the sign that dimerizes B rather than the sign that would split its Fermi points.

## Effective B model and its reachable region (Tasks 1 to 3)

- [proved] The whole P-even normal space at XZ is an exact two-leg Majorana ladder (PROOF 1.1). Its only quartic term is the YY plaquette a_j b_j. The A/B-mixing directions (XX + ZZ, B_y) are pure bilinears.
- [derived; C144 validation at beta = 0.004] On the A/B-preserving part of the tuned hyperplane the B chain is h_B = lam sum b_j - 2 sum_{r >= 2} K(r) b_j b_{j+r}, with lam = t/beta^2. Explicit couplings: 2K(2) = +0.1496 (an RZFA term), -2K(3) = +0.1019 (an O'Brien-Fendley term), and K(r) = (-1)^r/(2 pi r) + O(r^-3). The tail becomes exponential when A is gapped (K(2) = 0.056, 0.025, 0.005 for A dimerizations 0.05, 0.2, 0.5).
- [proved] The kernel K is positive semidefinite (Lehmann), with K~(0) = 0 exactly at XZ. [numerical] K~(q) rises monotonically to its maximum at q = pi. So the induced interaction is -beta^2 sum_q K~(q) b_q^dag b_q, attractive in every channel and most attractive in the staggered (Ising energy) channel.
- [derived] Matched to Rahmani-Zhu-Franz-Affleck (PRB 92, 235123), g_RZFA = 2 beta^2 K(2) > 0, with |t_RZFA| = |t|. The reachable region is the attractive half-line g/t in (0, infinity), covered as lam runs from infinity to 0, plus a positive-semidefinite long-range tail. RZFA's c = 3/2 Ising + Luttinger phase sits at -2.86 < g/t < -0.285, on the repulsive side, and is not reachable. Their tricritical point (g/t of about 250) corresponds to lam of about 6e-4, but the log-divergent staggered tail dimerizes B already at lam of order log(1/beta) (C144's wedge, g/t of about 0.1).
- [derived plus numerical] Hartree-Fock (hfk.py, calibrated on RZFA's mean-field Lifshitz point -0.29, which it reproduces) gives h3/h1 = (2/pi)|U_2|/|h1| < 3/7 for the attractive RZFA term (analytic; numerics 0.4267). The full kernel gives at most 0.94 (critical A), 0.85 and 0.64 (gapped A), at every lam. The Lifshitz route to c = 3/2 needs h3/h1 < -1/3 or > 1. Zero count 2 and c_MF = 1/2 everywhere.
- [numerical] ED of the B chain (L = 24 to 40 Majoranas): lam = 0.1, 0.25 (and C144's 0.5 to 1.5) are twofold, dimerized and gapped. lam >= 2 is Ising. lam = 0 is a staggered critical chain with c_B = 1/2. No phase has c_B > 1/2.
- [derived] With mixing, the hybridized band's zero equation (PROOF 2.8) has at most one root in sin^2 k unless h3/h1 > 1 or < -1/3, so c_tot <= 1. P-odd directions are relevant couplings of the A order and disorder fields (dimension 1/8), so there is no second-order B Hamiltonian for them. YX/ZY are a period-3 twisted XYZ chain (C154's finding, left to C154).

## Does any qubit direction reach a c > 1/2 phase of B?

No, at second order. The two candidate routes are c_B = 3/2 (Lifshitz, repulsive side, excluded by sign) and c_B = 7/10 with A critical, i.e. c_tot = 6/5 on the Ising-to-dimer boundary. The second is pre-empted. The A-mediated eps_B eps_B/|x - y| interaction is marginal at Ising and strongly relevant (RG eigenvalue 8/5) at a putative B tricritical point, with the dimerizing sign. The entry into C144's wedge is the Ashkin-Teller/XXZ Kosterlitz-Thouless route, c_tot = 1 to 0 [heuristic, consistent with C144's L = 128 fits, which fall smoothly from 1 to 0 and never exceed 1.09].

## Best explicit candidate: the bisector field (Kitaev chain in a [110] field, from C156)

- Bond: h = X (x) Z + B (s^x + s^z)/sqrt2 per site, B = 2h/K. This is the unique field direction invariant under (pi rotation about x + z) x (site reflection).
- [proved, block_check.out, 1e-14] Exact two-qubit-cell form: H = sum_i [ sigma^z_i + g^z_i sigma^x_i sigma^x_{i+1} + (B/sqrt2)( sigma^x_i (1 + g^z_i) + g^x_i (1 + sigma^z_i) ) ]. Here sigma is the critical A Ising chain, g are the flat B qubits (the Gauss operators of XZ) acting as its bond signs, and the field couples the A order field to (1 + g^z) and the A disorder field to g^x. [heuristic] This frustration of two strongly relevant couplings explains the tiny emergent scale.
- [numerical] c estimates. Central-entropy pair estimates c = 6 dS-bar/d log L (C156 job 17594698): at h = 0.4 they are 0.44, 0.94, 0.98 for L = 48 to 64, 64 to 96 and 96 to 128, with smoothed profile fit 0.97 at L = 128. At h = 0.45 they are 1.26, 0.89, 0.97, with smoothed fit 1.07. At h = 0.3 they are 1.05, 0.84 (L to 96). In the soliton phase (h = 0.6) they are 0.43, 0.30, 0.21 (gapped). Convergence: truncation 1e-9 or below, last-chi energy change 1e-7 or below. The pair estimate settles near 1 from L = 64 on.
- [numerical] Calibration of the estimators on a known c = 3/2 phase (RZFA, g = -1, rzfa_dmrg.py with a validated d = 4 DMRG): smoothed fit 0.91 at L = 64; 1.50 / 1.58 profile and 1.69 smoothed at L = 128; c_pair(64, 128) = 1.24. The bisector line is clearly below that at the same L. So c = 1 on the bisector line, not above, with the total central charge reading 1 whether or not A is separately gapped. That point was not settled: Casimir fits are dominated by incommensurate oscillations.
- Other candidates (Typhon job 17595396, pending): the Ashkin-Teller wedge boundary at beta = 0.5 (alpha = 0.21 to 0.27, L = 64 to 192), and mixing plus plaquette, XZ + alpha ZX + 0.3 YY + 0.1 (XX + ZZ) with alpha from -0.15 to 0.19 (L = 64, 128).

## Why the XZ loophole closes

B's dynamics on the tuned hyperplane come only from virtual A excitations. Second order makes them a positive-semidefinite (attractive) density-density interaction peaked at the staggered channel q = pi. The c = 3/2 phase of interacting Majorana chains requires repulsion there (g_RZFA < 0, Fock h3/h1 < -1/3). Attraction instead orders the staggered bond density, which is C144's dimerized wedge with c_tot = 0 and a c_tot = 1/2 line. The tricritical route is cut by the same log tail. Mixing directions gap A and add only bilinears. Caveats: (a) higher orders at finite beta could in principle change the sign of the induced r = 2 coupling, but C144's full-chain DMRG to beta = 0.5 found no c > 1.09; (b) the P-odd directions are non-perturbative. Of those, the bisector reads c = 1 at L = 128, and YX/ZY are with C154.

## Typhon

- C157 job 17595396 (packed 48-core node, 27 tasks), PENDING (Priority) at the end of this session. Outputs: /scratch/lustre/jkudlerflam/c157/{bis_h0.3,bis_h0.45}_L{192,256}.jsonl, at_L{64,128,192}_a*.jsonl, mix_L{64,128}_a*.jsonl, with .out files and slurm_17595396.out. Analysis: python3 cscale.py label <files> (pair estimates), python3 cfit.py <files> (window fits).
- C144 jobs 17580192 and 17580432 were still RUNNING (7.3 h and 6.5 h) at the end of this session. None of the requested L = 128 to 256 YX/ZY or tuned-line outputs existed yet in /scratch/lustre/jkudlerflam/c144/. Nothing was collected and nothing was deleted.
- C156 job 17594698 was RUNNING. Its h = 0.2 and 0.3 L = 128 points were not yet written when last copied.

## References

Rahmani, Zhu, Franz and Affleck, PRL 115, 166401 and PRB 92, 235123 (2015), arXiv:1505.03966 (Hamiltonian convention Eq. 4, spin form Eq. 7, phase boundaries -2.86, -0.285, about 250). O'Brien and Fendley, PRL 120, 206403 (2018). Sorensen, Gordon, Riddell, Wang and Kee, PRR 5, L012027 (2023), arXiv:2209.06221. arXiv:2211.15598 (revisits the RZFA C-IC endpoint). C143 section 4.2, C144, C154 (twisted XYZ form of YX/ZY), C156 (bisector identification and job 17594698).
