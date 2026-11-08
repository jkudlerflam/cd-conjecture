This run was done by a Claude Opus agent (2026-11-08), not by Codex. Astra was not used (out of credits).

# C160 result: c <= 1 near the ferromagnetic Heisenberg point

Derivations are in PROOF.md, the log in PROGRESS.md. Scripts in this folder: lsw_minima.py, twomag.py, csl_tilt.py, csl_fine.py, texture_periodic.py and texture_refine.py (pilot), typhon/ (dmrg_task.py, texture_task.py, summarize.py, c160.slurm), with C159's est.py, dmrg.py and dmrg_c.py copied unmodified.

## Verdict

Outcome (a), at the physics standard. On a full punctured neighbourhood of h_F every conformal bond has c in {1/2, 1}, and every attractive (large-K) liquid there is a single boson. This removes the window that C159 Proposition C left around h_F, in T, K and N alike, because the argument uses no symmetry. The one non-standard input is the semiclassical correspondence (SC) at large spin per coherence length, the same input that gives K ~ 1/sqrt(1 + Delta) for XXZ near Delta = -1. No c > 1 bond was found.

## Mechanism

Near h_F the chain is a quantum Landau-Lifshitz ferromagnet with spin s = 1/2 per site and a single unit-vector field. Its effective Planck constant is hbar_eff = 1/(s l) = 2 eps^{1/2}, or 2 eps in the wedge around the pure-DM directions. Here l is the coherence length (PROOF 1.3). So the infrared is fixed by the classical ground state and its soft modes. Three facts about the classical problem, a Tonelli Lagrangian on S^2 with the chain coordinate as time, then limit the soft sector to one canonical pair.
- (M1) [proved, given Mather's graph theorem] Classical ground states are uniform or singly periodic. Minimizing measures sit on a Lipschitz graph over S^2, and by Poincare-Bendixson the only recurrent orbits of a flow on S^2 are fixed points and closed orbits. So there is never a second phason.
- (M2) [derived; non-conjugacy of minimizers, Poincare's instability of minimizing closed orbits] A periodic ground texture has real positive Floquet multipliers, so its Hessian has no zero mode at nonzero Bloch momentum. This excludes period-doubling or incommensurate Ising order inside a soliton lattice.
- (M3) [derived] An extra zero mode at Bloch momentum 0 has a nodeless normal component. Its spin-symplectic pairing with the phason is the change of the solid angle enclosed by the orbit, which is nonzero, so the two form one canonical pair with omega ~ k^2 (z = 2). The other possibility (a mode along the family of orbits) is not a minimizer, because the compressibility diverges first.
Uniform states obey the continuum form of C154 Lemma 3 (one soft point). So the possible infrareds are c = 0, one compact boson (c = 1), one Ising mode (c = 1/2), or z = 2 points. The patterns P1, P2 and P3 of C154 and C159 each need a soft mode independent of the phason, which M1 to M3 forbid.

## Task by task

- Task 1, blow-up [derived]. J = -1 + eps (S0 + A(D)), B = eps b gives E[n] = int [ |n'|^2/2 + eps D.(n x n') + eps (n^T S0 n + b.n) ]. Leading phases come from the trust-region minimization of V = n^T S0 n + b.n on S^2. A unique minimizer is a gapped ferromagnet (open set). The hard case is a first-order surface between two canted states (codim 1). Its edge is an Ising pitchfork, c = 1/2 (codim 2, codim 1 inside K). A circle of minima (codim 4) and the pure-DM sphere carry wedges of relative width eps with c = 1 floating liquids (chiral soliton lattices, K ~ 1/hbar_eff), commensurate pinned phases, Pokrovsky-Talapov transitions and saturation onsets. In T near h_F there is no Ising point (b = 0 forces the Theta-broken pair), consistent with C159 Lemma B.
- Task 2, magnon kinematics [proved]. Lemma 3': with no DM along the moment, a stable polarized state has its magnon minima only at k = 0 or pi. So two degenerate minima at +-k_0 never occur for a one-site nearest-neighbour qubit bond. C154 Lemma 3 covers the case, and its "+-k_0 pair" item is empty for stable references. The J1-J2 mechanism needs a second harmonic (cos 2k) or a vanishing stiffness, and near h_F the stiffness is 1 - O(eps). Numerical check: 5916 stable random references, no interior minimum (lsw_minima.out, BdG check 7e-14).
- Task 3, bound states [proved for the number-conserving part, derived for the full bond]. The exact two-magnon bound state is E_b(K) = 2c0 + V + tau(K)^2/V, with V = 4J_z and tau = 2|a| cos(K/2 - phi). This is the Bethe 2-string, shifted by the DM twist (ED agreement 1e-4 at L = 32, twomag.out). The effective coupling g_eff = 8|J'_perp|(|J'_perp| - |J_z|)/|J_z| is repulsive on the easy-plane side, which gives one Luttinger liquid. On the easy-axis side it is attractive, and the exact string energies per magnon decrease monotonically with N (McGuire collapse), so magnetization reverses first order. Pair creation sources the pair channel directly. No regime with pairs and single magnons both soft at finite density. This reproduces field-induced XXZ near saturation (one species).
- Task 4, numerics. [numerical] On the K-stratum chiral soliton lattice, the most dangerous P1 candidate (a planar CSL that might tilt into "floating + SSB"), the tilt Hessian has its minimum at Bloch k = 0 with a nodeless mode. Its overlap with the phason is 0.995 to 1.000. A first-order transition to the canted ferromagnet always comes before the tilt instability (kappa_1 - kappa_c = 0.0002 to 0.057 for h = 0.15 to 0.55, csl_tilt.out, csl_fine.out, csl_firstorder.out). So the floating + SSB phase does not exist in that family. [pilot, inconclusive] The lattice texture periodicity test at N = 400 to 600 is not used as evidence. [pending] The DMRG and texture tests are queued on Typhon (below).
- Task 5, theorem [derived, conditional on (SC)]. There is eps0 > 0 such that every conformal bond within eps0 of h_F (h_F excluded) has c in {1/2, 1}. This also gives C159's (A6_TK) on that neighbourhood. PROOF section 8 lists the lemma tags.

## Narrowest remaining assumption

(SC): as hbar_eff -> 0, each linear classical canonical pair becomes one compact boson, a tuned real linear mode becomes Ising, quadratic pairs stay non-conformal, and gapped modes stay gapped. Secondary inputs: Mather's graph theorem and non-conjugacy of minimizers, applied to the continuum functional. Lattice corrections enter by persistence of hyperbolic invariant circles, a heuristic step. eps0 is not quantified.

## Pending Typhon job

- 17599583 (c160, partition typhon, one node, 48 cores, 8 h limit). Submitted 2026-11-08, PENDING (Priority). Directory /scratch/lustre/jkudlerflam/c160/.
- Outputs:
  - dmrg/{label}_L{64,128,192,256}.json. Open-chain DMRG of the CSL family v = (-1, -1, -1 + kappa eps^2/2, 0, 0, -eps, -h eps^2, 0, 0) at eps = 0.5, chi up to 160. Labels h30_k040, h30_k070, h30_k085, h30_k095, h30_k105, h30_k130 (h = 0.3), h00_k060 (U(1) calibration, XXZ with Delta_AF = -0.96, K of about 5.6), h45_k060 and h45_k085.
  - texture/tex_{0..11}.jsonl: 48 random pure-DM-wedge bonds, N = 2000, 16 starts each, and the spectral test of single periodicity.
  - summary.txt, written by summarize.py at the end of the job.
  - logs/c160_17599583.out and logs/*.log.
- Collect with `scp -r typhon:/scratch/lustre/jkudlerflam/c160/{summary.txt,dmrg,texture} .`
- Criteria, fixed in advance:
  - On the CSL side, c_prof, c_sm and the two-size c_two should approach 1 from below, with no window above 1.2 at L = 256 and truncation below 1e-7. The slow approach is expected because K is of order 1/hbar_eff.
  - The canted side should read near 0.
  - In the texture stage, any bond with unexplained (non-harmonic) peak weight above 1e-2 should be re-minimized before it is read as a counterexample to M1 on the lattice.
  - These sizes probe 15 to 20 soliton periods (period about 13 to 17 sites at eps = 0.5).
