This run was done by Claude Opus 5.5 (claude-opus-5-5), agent C184, not by Codex. Astra was not used (out of credits).

# C184 result: multicritical points with c > 1

(Registration first. The verdict and results are appended below the registration once the data are in.)

## Registration (written 2026-11-09 before any new numerical output was looked at)

The C152 job 17597010 outputs were copied and summarized before this registration (c152j2_summary.out). That was existing data from another agent, and it is reported separately, not counted as a test of anything registered here.

Target crossing. The region of the parity class at φ = 1, h_z = 0, V ≈ −2t, Δ ≈ 0.3t to 0.5t (fermion units of C152, t = 1), where the two Ising surfaces into which DM splits the attractive x-FM/z-FM transition meet the chiral Luttinger region and its first-order boundary (Peotta et al., arXiv:1403.4568, left unresolved there and by C152). C174 R4 names it among the uninventoried junctions.

Predictions (derived in PROOF.md §3 from Majorana-species counting and the AT/Lifshitz structure of the parity class):
- P1. Every junction found in this slice has c ≤ 1. At the best-tuned junction point, the open-chain pair estimate c_pair(96, 192), with S̄ the mean of the four central cuts, is at most 1.1.
- P2. No grid point has c_pair above 1.2 at both pairs (48, 96) and (96, 192).
- P3. The two Ising lines meet the Luttinger region either at a Lifshitz junction (c read between 0.5 and 1, falling or flat with L) or along the first-order boundary (entanglement jumps by about ln 2 at one parameter value, cat states).
- P4. Controls at the same sizes read close to their values: XX chain 1.0, critical transverse-field Ising 0.5, and the decoupled Ising × XX zigzag (outside B_2, c = 3/2) above 1.3.

Pass criterion for a c > 1 candidate: the adapted C163 §4.2 rule of the spec (open-chain fits above 1.2 not falling with L at truncation below 1e-7, ring-spectroscopy confirmation, controls at the same sizes).

## Verdict

No candidate. No multicritical point with c > 1 was found or constructed. The inventory (PROOF §1 to §4) shows that a crossing changes only the codimension count. Every c > 1 crossing still needs two independent soft sectors at one bond, or a single-species irreducible theory, and both are already the content of (A6)/(R_N). Lemma M (PROOF §3) makes this concrete. Near any quadratic reference a nearest-neighbour qubit chain has at most four soft Majorana–Weyl components at any crossing, and every c > 1 candidate needs at least six, or a strongly interacting species. So a c > 1 crossing can only sit in the four regions where a non-quadratic species exists: the XZ flat band, the Fibonacci cusp, strong attraction, and near the ferromagnet.

Best candidate bond: none. The most promising crossing outside the cusp was chosen and tested. That is the Peotta et al. meeting region of the parity class, which sits in the kinematic window where Ising × Gaussian has its lowest commensurate codimension (2). On rings up to L = 20 it reads c ≈ 1 on the band and a level crossing at the meeting point. Example bond on the band: h = 0.520 XX + 0.020 YY + 0.421(XY − YX) − 0.55 ZZ (h_z = 0), with ring c_S(16, 20) = 0.98. The registered open-chain test is in Typhon job 17604462, pending.

## Key claims

- [derived] Crossing framework (PROOF §1). At a crossing the theory is the same field (c ≤ 1: AT, Heisenberg point, eight-vertex, Lifshitz), a decoupled product with codimension (1.1), or a coupled fixed point with smaller c.
- [derived; agrees with C174 R1] The closure lemma does not exclude multicritical c > 1 (PROOF §2). A c > 1 multicritical locus is consistent with it and is excluded only by (H*)/(A6).
- [proved for quadratic chains, derived for weak interactions] Lemma M (PROOF §3). At most four soft Majorana–Weyl components at any crossing near the free locus, the AT plane, the U(1) class and the slow-B XZ directions. So c ≤ 1 at every perturbatively reachable crossing.
- [derived] Ising × Gaussian (PROOF 4.1).
  - At a commensurate Gaussian surface σ·V couplings are relevant (k ≥ 3) and the crossing locks, except in the window K > 15/8 with σ even under Π. There k = 2, and this coincides with C167's (A_W1)(a).
  - With a floating boson it is P1 (k = 1), open only at order-one harmonics.
  - The pair-unbinding mechanism (Ruhman–Altman, arXiv:1705.00635) needs a pair liquid, which nearest-neighbour attraction preempts by cluster binding [heuristic].
- [derived] TCI × Ising (PROOF 4.3). ε_Tε_I (dimension 6/5) is always allowed, so k ≥ 4. In the parity class, the only interspecies term Vε^Aε^B is also the only source of the TCI-driving self-interaction, which needs the other species gapped. So the crossing flows away [perturbative in V/t; heuristic at V ~ t].
- [derived] Ising² × boson needs eight Majorana–Weyl components against Lemma M's four. Potts × boson and Z_5 reduce to C168, C169 and C171. SU(2)_2 reduces to P1. Two Gaussian surfaces of one boson give c = 1, and of two bosons need the cusp. Lifshitz points are not conformal. (PROOF 4.2 and 4.4 to 4.10, with the summary table in 4.10.)
- [numerical, new collection of C152's job 17597010, never collected before; c152_j2_copy/, c152j2_summary.out]
  - The worst XZ-window point (φ = 1.45, s = 1.1, V = 0.4, h_z = 0) reads c_pair = 1.46, 1.24, 1.13 for (32, 64), (64, 128), (128, 192), falling toward 1.
  - The Peotta region (C152 Set 4) reads at most 1.15 at (48, 96).
  - Strong attraction (V = −1.5t) is gapped, consistent with first order. No TCI point was located.
- [numerical, new collection of C154's job 17594765, never collected before; c154_cyc_copy/] The cyclic stratum's SU(2)_1/ordered boundary looks like a level crossing. The triplet collapses (L·gap 2.03 → 0.107 along g/b = −0.6) while c_cas falls (0.93 → 0.69). No SU(2)_3 tower at any grid point.
- [numerical, ring ED, ring_phi1_table.out, ring_cS.out] Peotta slice at φ = 1, h_z = 0. The band reads c_S(16, 20) = 0.96 to 0.99, and the LL control reads 0.97. The meeting point at Δ = 0.3, V = −1.6 is a level crossing (entropy jump, then flat).
- [numerical, exact] Estimator calibration (ff_controls.out). With S̄ the mean of the four central cuts, c_pair(48, 96) and (96, 192) read: XX 0.989, 0.993; TFI 0.515, 0.507; XX×XX (c = 2) 1.977, 1.987; TFI×XX (c = 3/2) 1.504, 1.501. The estimator sees 3/2 and 2 at these sizes.

## Against the registration

P1 to P3 are not yet testable as registered, because the open-chain pairs at (48, 96) and (96, 192) are in the pending job. Ring level is consistent with P1 and P3. P4 holds for the exact free-fermion controls, and the DMRG qubit controls are in the job.

## What is left open

PROOF §8: the job's open-chain test; the K > 15/8 window of 4.1(a), which is (A_W1)(a); a fine scan of the cyclic SU(2)_1/ordered boundary for SU(2)_3; the TCI × Ising obstruction at V ~ t; and order-one harmonics under (R_N).

## Typhon job

Job 17604462 (c184), one 48-core node, 14 h limit, submitted 2026-11-09, PENDING (Priority) at write-up with estimated start 2026-11-09 15:29.
- Directory: /scratch/lustre/jkudlerflam/c18x/c184/. Inputs: c184.slurm, run_task.sh, tasks.txt (418 tasks from make_tasks.py), and the copied dmrg_pt.py, model.py, dmrg.py, dmrg_d.py, controls.py, summ.py.
- Expected outputs:
  - out/dmrg_48.jsonl, out/dmrg_96.jsonl, out/dmrg_192.jsonl (labels ph{φ}_D{Δ}_V{V});
  - out/ctrl_48.jsonl, out/ctrl_96.jsonl, out/ctrl_192.jsonl;
  - out/summary.txt (python3 summ.py out), out/err.log, slurm_17604462.out.
- Collection, from the cd_conjecture directory:

      mkdir -p "codex/work_C184/typhon_out"
      scp -r -o LogLevel=ERROR typhon:/scratch/lustre/jkudlerflam/c18x/c184/out "codex/work_C184/typhon_out/"
      scp -o LogLevel=ERROR typhon:/scratch/lustre/jkudlerflam/c18x/c184/slurm_17604462.out "codex/work_C184/typhon_out/"
      cd codex/work_C184 && python3 summ.py typhon_out/out

  Then read the result against P1 to P4 above. P2's threshold is c_pair > 1.2 at both (48, 96) and (96, 192), applicable only on the L = 192 grid (Δ = 0.4 to 0.7, V = −2.6 to −1.2). Compare with the exact controls in ff_controls.out.

Literature opened in this run: arXiv:1403.4568 (abstract) and arXiv:1705.00635 (abstract). Search snippets only, not opened: arXiv:2211.15598, arXiv:1012.4518, arXiv:2510.05988.
