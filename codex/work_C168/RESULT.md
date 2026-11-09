This work was done by Claude Opus 5.5 (agent C168, 2026-11-08 to 2026-11-09), not by Codex. Astra was not used (out of credits).

# C168 result: the symmetry-free stratum N

Derivations are in PROOF.md and the log is in PROGRESS.md.

Scripts in this folder:
- stage-A analysis and refinement: stageA_analyze.py, refine_stageA.py;
- screen calibration: calib_c2.py;
- Lemma D checks: imp_mix.py, imp_swap.py, imp_swap2.py.

Typhon packages are typhon/ (job 17602943) and typhon2/ (job 17603367). Copied outputs are in typhon_out/.

C155's est.py, dmrg.py and dmrg_c.py are copied unmodified as est_c155.py, dmrg_c155.py and dmrg_c_c155.py. Only the import lines of dmrg_c_c155.py are renamed. C158's fesx.py and fesfit_c158.py and C153's vumps.py (via C165) are copied unmodified in typhon2/.

## Verdict

Outcome: a strictly narrower residual. Neither (H) nor (A6) is derived on N, and no c > 1 bond was found. Two stage-A points await the second stage of the pre-fixed rule (section 4).

The dynamical content of (A6) on N now reduces to five named softening events (R_N), none of which has a small parameter. What is newly derived:
- At a generic double onset out of any gapped phase of N, periodic and entangled references included, the two dilute species are immiscible or alternate, so c ≤ 1 on a full neighbourhood. This closes the entry of P3 from gapped phases at codimension 2 (Lemma D).
- Onsets at codimension 1 carry one species, for any reference state. This replaces C154 Lemma 3, as C174 R9 asks (Lemma O).
- Clock tongues of every order melt through c ≤ 1 at weak locking. This is C174's P4 in that regime (Lemma P).
- Three-state Potts × floating boson (C171 R1) has codimension 2, not 1, in N.
- c ≤ 1 holds on a full neighbourhood of the AF Heisenberg point in N (Lemma H).

The sampling evidence was much weaker than stated before. C155's screen detects only about 1 percent of c = 2 chains, and its 28 stage-A flags had never been refined. With a calibrated two-stage screen and DMRG to L = 128, the uniform sample bounds the volume fraction of any open c > 1 region by p < 7e-4 to 1.4e-3 at 95 percent confidence, provided the two pending points fail.

## Theorem-style summary

1. [proved, standard] Lemma 1(a). Translation alone has no anomaly (H^3(Z, U(1)) = 0), and a filling constraint needs a microscopic U(1). So no LSM- or filling-type constraint acts on the emergent charges of a floating liquid in N.
2. [proved for H0, derived with W; numerical support] Lemma 1(b, c). A symmetry-free qubit chain with range 3 (Jordan-Wigner fermions with complex NN and NNN hopping, plus a generic NN term) has two floating bosons (c = 2) on an open set. On proxies at rings 12/16 the median c_S is 1.99 for two Fermi seas and 1.00 for one. So both routes suggested in lead 1 are impossible in principle, since they use only one qubit per site and translation. Any exclusion must use the nearest-neighbour range and dynamics.
3. [derived] Lemma O. In N, at a generic point of the boundary of any gapped phase, one branch softens at one momentum: one floating boson (c = 1) or Ising (c = 1/2). The codimension-2 points are:
   - (O1) two incommensurate minima;
   - (O2) a real mode at 0 or pi together with an incommensurate one;
   - (O3) real modes at 0 and pi (two Majoranas, c ≤ 1);
   - (O4) Lifshitz points.
   Charge-conjugation actions of translation, the "C row" of C171 R3, and swapped boson pairs arise only from (O3) pairs, which give one boson.
4. [derived, numerical checks] Lemma D. At an (O1) point with unequal effective masses (relative imbalance eta ≠ 0), no bound states and no zero-energy transmission resonance, the dilute species are impenetrable. Their order costs 0.31 eta^2 hbar^2 n^2/M per unlike interface (imp_mix.out: E(ABAB) - E(AABB) = 19.6 eta^2/L^2, L = 16 to 32; C164 had 0.3). Exchange is O(hbar^2 n^3 a_t/M) (imp_swap2.out: dE_x L^3 → -120). So the pseudospin anisotropy |Delta_ps| is of order eta^2/(n a_t), which diverges at the onset: segregation (first order) or alternation (pseudospin gap). Hence c ≤ 1 on a full neighbourhood. In the swap model at eta = 0.5 the alternating weight falls from 0.33 to 0.15 (L = 12 to 32), and at eta = 0 it stays at 0.52. Equal masses, a resonance, or bound-state degeneracy costs codimension 3, and then only thin wedges |eta| ≲ (n a_t)^{1/2} with an easy-plane pseudospin could carry c = 2. Lemma D fails in K, where the antiunitary equates masses, and that is C162's domain.
5. [derived, literature inputs] Lemma P. While a tongue's neighbourhood is one compact boson with near-Gaussian vertex dimensions, the locking vertex V_{m,0} and the dual vertex V_{0,1} have dimension product m^2/4. They can be relevant together only for m ≤ 3, where the self-dual points are Ising and Potts, with c ≤ 1. For m ≥ 5 a c = 1 window intervenes, and chiral drift gives Pokrovsky-Talapov transitions. So the Z_m parafermion with m ≥ 5 (c = 2(m-1)/(m+2) > 1), and P4 in general, needs clock amplitudes of order the exchange, and it has codimension at least 2 in N.
6. [derived] Potts × floating boson (C171 R1, c = 9/5). The chiral Huse-Fisher operator (spin 1, dimension 9/5) is translation invariant and allowed in N, so the conformal point has codimension at least 2. The generic transition at the would-be one-tuning surface is chiral, non-conformal or split by a floating phase.
7. [derived] Lemma H. Near the AF Heisenberg point in N, translation forbids tr g and N^a to all orders. The chiral fields h_{L,R} = c1 B ± c2 D shift the two chiral bosons by different incommensurate wavevectors when |h_L| ≠ |h_R|. All marginal terms except J^z_L J^z_R then oscillate out, giving c = 1, and the special loci give XXZ-type c in {0, 1}. This is uniform in direction, which is the C174 R3 concern for this point, and it extends C159 3.2 and C165 4 there.
8. [derived] Proposition F, frontier classification. This uses C174 Lemma R1: lead 3 is dropped in its connectivity form, because bounds only travel from higher c to lower c. Given (A1) with c per sector, (A2), S_ctrl (C174 R2) and Lemmas O, D, P and H, the outer frontier of any cluster of c > 1 components of N consists of the events in the residual below, or points of S_ctrl. So c ≤ 1 on N follows from (R_N).
9. [assumed] Residual (R_N), on N outside S_ctrl and outside the neighbourhoods of Lemmas D, H and P:
   - (i) In a c = 1 liquid of N (any translation action), no sharp excitation independent of the liquid's boson reaches zero energy while the liquid survives. This applies at momentum 0, pi, 2 pi/3 or incommensurate, and includes composites of host quanta. It is A_S and A_W restated for N, and it covers P1 from the floating side, P2, P3 by guest onset, Potts × boson, the C row and C163's Lemma-L gaps.
   - (ii) No codimension-2 point where a c ≤ 1 transition of a gapped phase (Ising, two-Majorana, Potts or chiral clock) coincides with the onset of an independent branch. This is the P1 corner (O2).
   - (iii) No codimension-3 equal-mass, resonant or bound-state-degenerate double onset with an easy-plane pseudospin.
   - (iv) No irreducible c > 1 fixed point (P4: Z_m parafermion with m ≥ 5, N = 1 superconformal or W_3 minimal model, SU(2)_3) at clock amplitudes of order the exchange.
   - (v) No c > 1 component whose outer frontier is entirely first order.
   (R_N) is implied by (A6) and is strictly narrower, because items 3 to 7 are now derived.

   Where it bites:
   - (i) and (ii) need a lattice-scale branch to come down at harmonics of order the exchange. No small parameter exists there, because in the U(1) class, near h_F and near the AF and hidden Heisenberg points the infrared has one sector. The dilute method of Lemma D does not reach a guest inside a dense liquid. The guest is transparent to the host, because host backscattering needs a momentum 2 k_F that the guest cannot absorb on shell. So a two-component Luttinger liquid is self-consistent there. A Wentzel-Bardeen instability argument was tried and rejected, since the exact compressibility matrix stays positive semidefinite, as the Bethe-ansatz Hubbard chain near saturation shows. At the corner (ii) the dilute particle couples to the critical Ising sector only through an exactly marginal energy defect, so a P1 surface t + g n = 0 emanates without obstruction.
   - (iii) and (iv) need extra tunings and coupling-sized amplitudes.
   - (v) is invisible to local softening analysis, and only sampling addresses it.
10. Scope. K_II is outside this run and stays under (A6)/(H*), as C174 R6 asks. Lemma D does not apply there.

## Numerics

### Stage-A gap (lead 1)

[numerical] C155 stage A (12,000 uniform points, rings 12/16) has 636 points with crit < 0.15, and 28 of them have c_obj > 1.1. These 28 had never been refined (C155 refine.py reads only search and gap-closing outputs).

### Screen calibration

[numerical, Typhon job 17602943, typhon_out/summary.txt] The proxies are 1600 range-3 qubit chains with known c. Fraction flagged:

| family, c | F1: c_S(12,16) > 1.5 | F1 and F2: also c_S(16,20) > 1.5 | C155 rule |
| :-- | --: | --: | --: |
| near-free, c = 1 | 0.170 | 0.016 | 0.008 |
| near-free, c = 2 | 0.641 | 0.341 | 0.018 |
| interacting, c = 1 | 0.107 | 0.012 | 0.006 |
| interacting, c = 2 | 0.331 | 0.178 | 0.006 |

### Refinement

[numerical, same job] Stage S1 ran all 416 stage-A points with c_S(12,16) > 1.5 or a C155 flag at rings 16/20. Nine of them have c_S(16,20) > 1.5. Stage S2 ran open-chain DMRG at L = 32, 64, 128 (chi up to 128, truncation ≤ 3e-10) on 45 points (all of those, plus c_obj(16,20) > 1.1 and all 28 flags). 43 of them read c_prof(128) ≤ 1.16.

Two meet the first stage of the pre-fixed rule (c_prof(128) and c_two(64 -> 128) above 1.2, truncation below 1e-7):
- idx 10193 (v = (-0.307, -0.004, -0.114, -0.712, 0.246, 0.176, -0.291, 0.222, 0.400)). c_prof is 0.51, 0.61, 1.94. Its L = 128 energy lies 0.06 above the continuation of L = 32 and 64, and its 64 -> 128 energy slope is 9e-4 above the 32 -> 64 slope, the largest anomaly in the set. My reading is an unconverged L = 128 state [heuristic until the reseeds report].
- idx 450 (v = (-0.083, -0.366, -0.186, -0.045, 0.640, 0.249, 0.183, 0.337, -0.451)). c_prof is 0.09, 1.62, 1.21, falling from L = 64, with energy slopes consistent to 3e-5. This looks like an incommensurate crossover [heuristic].

### Bound

[numerical, statistics] If both points fail the second stage, zero detections give:
- p < 7.3e-4 (95%) and 1.1e-3 (99%), using the near-free efficiency eta = 0.341;
- p < 1.4e-3 (95%) and 2.2e-3 (99%), using the interacting efficiency eta = 0.178.

Here p is the normalized S^8 volume of any c > 1 set, and the DMRG stage is assumed to pass a true c = 2 point (C164's two-XX control reads 1.73 at L = 64). With the C155 rule alone the 95 percent bound would have been about 4e-2. Caveats: the proxies have range 3 and are close to free, the bound uses C155's coordinate measure, and it says nothing about small-volume sets.

### Side remark

[numerical] crit < 0.15 passes only 6 to 17 percent of c = 1 proxies, and the median stage-A c_S is 0.93. So the stage-A "gapless fraction" of 5.3 percent likely undercounts the floating volume.

### Local runs

[numerical] refineA_table.out covers 18 of the 28 flags at L ≤ 64 locally, and agrees with Typhon to 1e-6 in energy. calib_c2_local.out is a 324-sample pilot of the proxies. imp_mix.out, imp_swap.out and imp_swap2.out are the Lemma D checks. typhon2/loc/ holds small-chi VUMPS pilots on 10193 and 450 (two-site cell, chi 8 to 24, damped, err 0.08 to 0.2, so not converged).
- On 10193 the pilot sits at e = -0.96827, well above the DMRG bulk slope of about -0.977. Its transfer spectrum has an eigenvalue at -1 (xi of order 3e4). That suggests the ground state needs a cell longer than two sites, and the four-site FES in job 17603367 tests this.
- On 450 the pilot sits at e = -0.87335 against -0.87869, and the pair slope over chi 8 to 24 is about 0.8.
These pilots carry no weight either way.

## Typhon jobs

- 17602943 (c168): COMPLETED 2026-11-09 01:17, 1 h 06 min, 48 cores. Outputs are in /scratch/lustre/jkudlerflam/c168/{out,logs,summary.txt}, copied to typhon_out/.
- 17603367 (c168b): PENDING (Priority), scheduler estimate 2026-11-12 14:26, one 48-core node, 10 h limit. It runs FES for 10193, 450 and the control 10540 (two-site and four-site cells, chi 16 to 256) and DMRG reseeds (L = 128 with seeds 1 and 2, L = 192). Outputs are in /scratch/lustre/jkudlerflam/c168/b/{fes,dm,logs,fesfit.out}.
  - Collect with `scp -r 'typhon:/scratch/lustre/jkudlerflam/c168/b/{fes,dm,logs,fesfit.out}' typhon_out/b/`.
  - Pre-fixed reading: a point is a c > 1 candidate only if the FES slope is above 1.2 at chi ≥ 256 and the reseeded DMRG at L = 128 and 192 keeps c_prof and the two-size slope above 1.2 at energies no higher than the L = 128 run of 17602943. Otherwise the bound above stands.
