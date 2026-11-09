This work was done by a Claude Opus 5.5 agent (C172, 2026-11-08/09), not by Codex. Astra was not used (out of credits). No Typhon job.

# C172 result: closing A_S analytically

Derivations in PROOF.md, log in PROGRESS.md. Scripts: ed_kink.py (matrix-free ring ED with momenta, `fam` copied from C165), run_batch.sh and run_batch2.sh (drivers), wall_pt.py (exact first-order projection onto one- and two-wall configurations), wall_scan.py and wall_target.py (scans). Outputs: ed_*.out, ed_*.json, ed_summary.txt, wall_scan_1.out, wall_scan_2.out, wall_target.out. The target is A_S in C170's corrected form: along the quantum switch set (KT point, Néel plateau edges), no second low-energy sector, in particular no bound magnon, at order-one harmonic amplitude.

## Verdict

A_S is not derived in full. It is reduced to a strictly narrower residual, and it is derived at order-one amplitude at both ends of the switch set. The key step is a frame-free form of the bound-magnon question. At the edge of any period-2 phase in N, a route to c = 3/2 exists only if a kink-even excitation softens before the kink, which needs a kink-even state below the two-kink continuum. In the Ising regime this is excluded for every harmonic amplitude below the point where the transverse exchange block becomes singular, with a finite margin. Near the KT end it is excluded by the universal spinon S-matrix. In between it holds numerically at amplitudes up to 1.4 with a resolution of about 0.1. No violating bond was found.

## Key claims

- [proved] Exact on-site SU(2) symmetry exists only on the twisted Heisenberg curve J ∈ ℝ·SO(3), B = 0 (PROOF 1.2). So no exact identity keeps a hidden SU(2) at order-one amplitude, and since ring momentum is the only exact quantum number in N, the SU(2)_1 towers have no non-crossing protection. Leads 1 and 2 of the spec close negatively.
- [derived] Locking lemma (PROOF 2). An independent momentum-π Ising sector couples to every liquid vertex V with σ-frame momentum π + δ_V. If x_V < 15/8 and |δ_V| ℓ_g ≲ 1, the two lock into a gapped period-2 state or a single-boson soliton liquid. On T (δ = 0) this recovers C159's second tuning. Near T, any c = 3/2 hypersurface of N is pinched onto C159's two-tuning locus. Deep in the liquid the lemma only constrains, since the bare couplings are not computable.
- [derived] Kink-parity dichotomy (PROOF 3). Kink parity is exact on rings (odd versus even length). At an edge of a period-2 phase either kinks soften first, and the period-2 order becomes a vertex of the kink liquid (c = 1), or a kink-even state softens first, which opens floating + SSB, the ordered side of C154's P1. The second case needs a kink-even bound state below 2D_K. Any binding is fatal, so this has to be excluded exactly.
- [proved within first-order perturbation theory, derived for the chain] Lemma W (PROOF 4.3). With transverse block A σ⁺σ⁻ + P σ⁺σ⁺ + h.c. and det M = (|A|² − |P|²)/4, charge-0 and charge-2 harmonics of any size and any longitudinal field, in the Ising regime Dz ≫ |A|, |P|: no kink-even state below the two-kink continuum when |P| ≤ |A|, with margin 2(|A| − |P|) against extra contact attraction. The proof is a half-line Dirichlet inequality. The bound is sharp: at a = π a bound magnon appears for every |P| > |A| with energy −(w + t²/w), reproduced to four digits. For C165's family |P| ≤ |A| means g² + g_xy² ≤ 1, the whole range where its parent is the polar parent, and at |P| = |A| the transverse block is singular. Charge-2 contact hopping is the quantum counterpart of C165's convex quartic, and the spin-1/2 hard core beats it until the block becomes singular.
- [numerical] Charge-1 terms in the Ising regime: 120 random and 72 targeted bonds (wall_scan_*.out, wall_target.out). Charge-1 terms and fields only raise the two-wall energy, and the threshold never moves below |P| = |A|.
- [derived, conditional on the KT point lying in the SU(2)_1 basin at that amplitude] Near the KT end the kink-even channel at energy D_K is governed by the massive SU(2)-Thirring spinons, which have no bound states. Irrelevant corrections vanish as D_K → 0 (PROOF 4.4). The width of this neighbourhood is not computed.
- [numerical] Full ring ED, L = 10 to 18, at Dz = 1.6 and 2 (PROOF 4.5, ed_summary.txt). At harmonic amplitudes up to ‖·‖_F = 1.4 (charge-2 up to 0.9 of the threshold, all nine harmonic types), the shift of D_2 − 2D_K relative to the parent is at most 0.105, against 2D_K = 0.9 to 1.9. Along field paths through the edge the kink gap closes first while D_2 stays near 0.6. A weak bound state below about 0.1 is not resolved.

## Residual A_S′ (strictly narrower)

On the plateau-edge part of the switch set at intermediate anisotropy, between the KT crossover and the Ising regime (roughly 1.25 ≲ Δ_eff ≲ a few), at order-one harmonic amplitude with non-singular transverse block (|P| < |A|), the two-kink channel has no bound state below the two-kink continuum. It is a two-body statement on a gapped background, so it can be tested directly by odd/even-sector DMRG or a two-kink ansatz on VUMPS, with resolution well below 0.1. The original A_S was a many-body statement about the liquid. Away from the switch set the liquid interior stays under C154 (A6), now with Lemma L's detuning condition.

## Open

- A proof at intermediate anisotropy. The natural route is to extend Lemma W's Dirichlet margin beyond first order, for example by a variational bound on the two-kink sector.
- Lemma W with charge-1 terms is numerical only.
- The KT-end statement assumes the SU(2)_1 basin at that amplitude (C170 derived it for weak harmonics).
- C165's Typhon screen (job 17601089) was not used or redone. C170 notes that its A2 fields mostly lie above the parent plateau edges.
