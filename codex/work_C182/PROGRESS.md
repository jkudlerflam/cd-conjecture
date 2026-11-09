# C182 progress log

This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

- Read specs/C181_C185_counterexample_push.md, specs/_qubit_c_context.md, STATUS.md, specs/C175_C180_fibonacci_swarm.md, REPORT §9.137 to §9.140, RESULT.md of C173, C175, C179, C180, C176 (registration only, job pending), C181 (registration only), contrib/cusp_scattering_20261109 (README, DERIVATION). Local load average 94 to 250 on 10 cores at start, so heavy work goes to Typhon.
- Registered Q1 to Q5 in RESULT.md before writing any code.
- pairband.py: exact two-body bands on the half-line (banded Dirichlet restriction, R up to 4000). Self-check: E(K = π) at |w| = 1, 2 equals the collaborator's closed form to all printed digits (−4.519526534666399, −6.122547914236252).
- sym_pairs.py → sym_pairs.out: φ = π/2 pair bands for |w| = 0.3 to 3.9 (run in background, slow under load).
- nbody.py (C175's nbody_box.py generalized to any φ with pairband's hop table); check at |w| = 2: trimer 0.0183 below dimer + atom, matching C175's 0.018.
- clusters.py (cluster ladder N ≤ 5 at valley-sum momenta with golden refinement), ladder6.py (balanced hexamer at K = π), make_tasks.py, select_stage2.py, typhon/c182.slurm. C176's fibdmrg.py, run.py, analyze.py copied unmodified. Smoke test at L = 24 passed (local/smoke).
- Registered Q6 and Q7 (stage-2 rule) before submission.
- Typhon job 17604443 submitted (partition typhon, one node, 48 cores, 320 GB, 3 days; 904 stage-1 tasks then stage 2), scratch /scratch/lustre/jkudlerflam/c182/. Queue was long at submission (C176's 17604009 still pending).
- sym_pairs.out done (|w| = 0.3 to 3.9). pairmass_exact.py: exact 1/m_pair positive and monotone on (0, 4), 4 → 0.643.
- Local small-box ladders at |w| = 3 (local/cl_w3_d0_small.json, local/cl_w3_d0.75_small.json, local/l6_small.jsonl, summary local/ladder_w3_summary.out): per-particle energies at φ = π/2 e2..e6 = −0.0975, −0.0877, −0.1107 (R 26), −0.1087 (R 18), −0.1194 (R 14, hexamer at K = π); at δ = 0.75 E_b: −0.0247, −0.0324, −0.0388, −0.0450. Hexamer at |w| = 2, R = 14: −0.0600 (box too small, Typhon has R = 16, 18).
- local/mini (L = 32, 64, χ ≤ 96; 2 points plus controls, Q8 not formally covering it): controls c1 1.05, c2 2.12 (W = 1); |w| = 3, π/2, N/L = 1/4 forms a droplet of density 0.29; |w| = 2, δ = E_b, N/L = 1/6 forms a droplet of density 0.23 with valley ratio 0.77 and c(32 → 64) = 2.5 to 2.6. mini_an.py summarizes.
- Registered Q9 and launched local/b2 (liquid side, L = 40, 80). Registered Q10 and launched local/b3 (follow-up of the |w| = 2, δ = E_b lead). A first b3 launch at n = 0.2 was killed within a minute, before any output was read, because 0.2 is below the droplet density; relaunched at n = 1/4 and Q10 corrected accordingly.
- local/b2 done (summary.out): |w| = 3 liquids at n = 0.3, 0.4 read c ≈ 0.87, homogeneous; |w| = 2 at n = 0.3 is a droplet (n_mid 0.35); |w| ≤ 1 at n = 0.2 inhomogeneous with high readings. Controls 1.10 and 2.17 (W = 1).
- local/b3 done (summary.out): the |w| = 2, δ = E_b lead is a droplet (n_mid 0.328, 0.313, 0.297 at L = 40, 80, 120), parity gap −0.010 (L = 80) and +0.014 (L = 120), so not a lead by Q10.
- Registered Q11, ran local/b4: |w| = 2 at n = 0.4 homogeneous, c ≈ 0.85 to 0.95; |w| = 1 at n = 0.35 edge-depleted (n_mid 0.414, 0.388), c_prof(80) = 2.0; |w| = 0.5 at n = 0.35 is a droplet of density 0.47.
- Registered Q12, ran local/eos (L = 40, N = 1 to 20, φ = π/2): |w| = 1 monotonic μ, |w| = 3 self-bound, |w| = 2 clusters of 4 to 6, |w| = 0.5 flat μ for n 0.2 to 0.3 (coexistence).
- Typhon 17604443 still pending at 12:30 (estimated start 20:08). RESULT.md written with the local results and the collection command. Not waiting on the job.
- local runs finished; eos.out and b3/summary.out refreshed. select_stage2.py amended (homogeneity filter) and re-uploaded while job 17604443 was still pending; amendment recorded in RESULT.md Q7.
