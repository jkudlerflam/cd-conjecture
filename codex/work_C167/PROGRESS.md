# C167 progress log (Claude Opus 5.5, 2026-11-08)

- Read spec C167, shared context, C159/C163/C165 RESULT and PROOF, C161 REFEREE, C160 PROOF 1 to 3, C162 and C164
  RESULT, C164 typhon summary, REPORT 9.120 to 9.126, C166 RESULT (for the KI22/KI26 status).
- Typhon queue checked: job 17602270 (C166) and about 80 other jobs PENDING (Priority). Decided to keep numerics local.
- charges.py: charge-1 / charge-2 decomposition of the harmonics. Lemma R0 checked (removing the charge-1 part of
  200 random bonds leaves a pi_m-symmetric bond to 1.2e-13; XYZ + D_z bonds have no charge-1 part). The 16
  attractive C159 points have charge-1 fractions 0.09 to 0.48 (all exchange, B = 0) and charge-2 0.03 to 0.54.
  hi_K_1 (K_I) has field fraction 1.02 and charge-2 0.77 at parent Delta = -0.15; it is the one window point whose
  open-chain c rose above 1 at L = 256 (1.45).
- hi_K_1, two-site and four-site VUMPS to chi 96 / 64 (loc/hiK1_n2.json, loc/hiK1_n4.json): period-2 in-plane
  order, staggered moment 0.183 converged, S saturates (c_pair 0.29, 0.17, 0.08 over chi 32 to 96), xi 22 still
  creeping. No period-4 component. Open-chain two-size energy slopes approach the VUMPS e0 = -0.6562052.
- path_hiK1.py: harmonics scaled by lam. lam <= 0.98: c_pair = 1.0 liquid, staggered moment decays as xi^(-x_pi)
  with x_pi = 1.03 to 1.1 for lam <= 0.85 (parent K = 1.06), then 0.98, 0.81, 0.72, 0.62, 0.56, 0.47 at
  lam = 0.9, 0.95, 0.96, 0.97, 0.975, 0.98. lam = 1: ordered. lam = 1.05, 1.1: polarized (xi = 2). So the
  period-2 order enters where the boson's own V_c reaches dimension 1/2 (KT), and hi_K_1 sits just inside it.
- A hysteresis restart from the polarized state (loc/hyst) failed to converge (non-injective), discarded.
- Running: fine scan lam = 0.985 to 1.0 (chi to 48), open chains of hi_K_1 at L = 64 unpinned and with two edge
  pinnings (oc/), and the charge-separated runs q1only / q2only at lam = 1 (Proposition R test).
- Fine scan (chi to 48): x_pi = 0.36, 0.23, 0.12 at lam = 0.985, 0.99, 0.995 with c_pair 1.00, 0.85, 0.58
  (ordered side, finite-chi crossover). Last-pair c_pair never above 1.06 on the path.
- Charge-separated bonds at lam = 1: charge-2 only is a c = 1 liquid (Theta, Pi exact, no moment); field only is a
  c = 1 liquid with x_pi = 0.75 to 0.81. Only the combination reaches x_pi = 1/2.
- Pinned open chain L = 64 (oc/hi_K_1_LRst_64.json): flat profile at the VUMPS entropy, c_prof 0.62 (unpinned C159
  value 1.37). Unpinned L = 64 rerun and the swapped-registry pinning were stopped for time (C159 has the unpinned).
- T points (loc/tpts/): ising_T_43 c_pair 1.00 (Theta-symmetric), hi_T_34 c_pair 0.91 to 0.98 with an in-plane
  finite-chi moment decaying as xi^-0.32 (unidentified). hi_T_39 stopped for time.
- Typhon job 17602938 submitted (c167, /scratch/lustre/jkudlerflam/c167/), PENDING. Contents in RESULT.md.
- Wrote PROOF.md (Lemma R0, Proposition R, Lemma C, Lemma Theta2, Lemma O, negative results, numerics, residual A_W1)
  and RESULT.md.
