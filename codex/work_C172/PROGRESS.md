# C172 progress log

Model: Claude Opus 5.5 agent (2026-11-08). Astra not used (out of credits). No Typhon job.

- Read specs/C169_C174_swarm.md (shared part and section C172), specs/_qubit_c_context.md, STATUS.md, REPORT.md 9.110 to 9.126, codex/work_C154/RESULT.md, codex/work_C165/ (RESULT.md, PROOF.md, PROGRESS.md, classical.py conventions, typhon/make_tasks.py families).
- Plan: (1) exact hidden-symmetry submanifolds of the A_S region in the polar frame; (2) a locking lemma for an independent momentum-pi Ising sector against the liquid's own density vertex; (3) the kink-parity dichotomy at the edge of a period-2 phase, reducing the edge part of A_S to two-kink binding; (4) small ring ED (L <= 20) testing only the derived two-kink statement at Neel points of C165's families.
- Wrote ed_kink.py (matrix-free ring ED, momenta from <psi|T|psi>); calibrated on Heisenberg L = 10 (E0 = -18.06 in Pauli units, matches the known -4.5154 x 4).
- Batch 1 (run_batch.sh): parent P0 and P7 (Dz = 2, hz = 0, 0.6) and P1 to P8 with C165-type harmonics, L = 10 to 18. Parents give identical b(L), as field-independence of binding in U(1) requires. Harmonics shift b by at most 0.05. Past the parent edge (hz = 1.0, 1.4) the kink gap closes first.
- Coordinator note received: C170 (codex/work_C170) corrected A_S to the quantum switch set and named the bound-magnon exclusion at order-one amplitude as the target, the plateau edge at Dz = 2 being 0.78 and C165's "Neel to hz = 2.2" metastable. Re-aimed the work at that statement. Read C170 RESULT.md and PROOF R4.2.
- wall_pt.py: exact first-order projection onto one- and two-wall sectors. Found the charge-2 contact hop of a reversed spin and the threshold |P| = |A| (sharp at a = pi, E = -(w + t^2/w)). Proved Lemma W (half-line Dirichlet inequality) for charge-0/2 plus longitudinal field.
- wall_scan.py (120 random bonds) and wall_target.py (72 targeted bonds at the sharp pitch): charge-1 terms and fields never bind for |P| < |A|; they only raise the two-wall energy.
- Batch 2 (run_batch2.sh): Q1, Q2 (amplitudes 1.4, 1.28, ratio 0.9), R0 to R2 (Dz = 1.6), S0, S1 (pitch a = 2.2). Shifts of b at most 0.105 against 2D_K = 0.9 to 1.9. Summary in ed_summary.txt.
- Wrote PROOF.md and RESULT.md. No Typhon job, no Astra, no commits.
