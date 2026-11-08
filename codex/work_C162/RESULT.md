This work was done by a Claude Opus agent (2026-11-08), not by Codex. Astra was not used (out of credits).

# C162 result: two-valley soliton onsets

Derivations are in PROOF.md, the log in PROGRESS.md. Scripts in this folder: kfun.py (general K functional), search_caseB.py, probe_caseB.py, opt_caseB.py, bound_scan.py (bound_scan.out), lemma_check.py (lemma_check.out), scatlen.py (scatlen.out). All local runs took seconds to about a minute. No Typhon job was submitted.

## Verdict

At a soliton onset over a symmetry-protected two-valley band (the K stratum, and the U = I R_pi(n_bg) class, which are the only two strata with such protection) the dilute two-valley soliton gas is immiscible. One valley wins, the Z2 breaks spontaneously, and c = 1. In K this rests on a proved classical rearrangement inequality and is derived for the quantum gas. In U it is derived at the tail level. No explicit c = 2 bond was found. The narrowest remaining assumption is about accidental (unprotected) two-valley soliton onsets in N and T far from K and U, at eps of order one.

## Mechanism

- The dilute gas of hard-core solitons maps to a pseudospin XXZ chain on the squeezed lattice, J_z = c (a_{++} + a_{--} - a_{+-} - a_{-+}) with c = pi^2 hbar^2 n^3/(3M), and J_perp from valley exchange (pass-through or flip-flop), exponentially small in 1/hbar_eff [derived]. c = 2 needs |J_z| < |J_perp|. S^z is conserved at every order except on measure-zero density lines, because a valley flip carries the incommensurate momentum 2 delta [derived].
- Rearrangement lemma [proved, continuum, and lattice for J_z <= 0]: in K the DM density (n x n')_z does not involve n_z and W is even in n_z, so replacing n_z by |n_z| preserves the energy. Every minimizer has single-signed n_z, and V_{++}(R) < V_{+-}(R) at every separation. This is the classical shadow of Eisenberg-Lieb, strengthened from the SU(2) point Delta = -1 to Delta < -1.
- Scattering-length calculus [proved: monotonicity identity da/dlambda = k2 int psi^2 dV/dlambda, and additivity for a single exponential channel. derived: sign rule J_z/c = -c_2 I (X^+ - X^-)(Y^+ - Y^-) from the subleading channel]. The slowest channel gives additive hard rods and no pseudospin coupling. So in K, J_z < 0 by a power of eps, while |J_perp| is exponentially small, giving Delta -> -infinity [derived].
- Tail ordering [proved operator inequality lam0 >= m_w^2 - m_u^2 + min(2 phi' + 2(L - L_bg)), numerical for the bracket]: in 4719 random K shapes every tilt-unstable onset had m_u^2/m_w^2 >= 4.6. So the tilt tail is slow, like valleys attract and bind, and the onset is first order into the same-tilt CSL, as in C161 [derived].
- Beyond the dilute limit [proved classically, derived quantum]: no alternating or mixed soliton lattice is a classical minimizer at any density, and a valley flip is gapped at every spacing.

## Numerical claims

- C161's tilt eigenvalue reproduced (-0.0670).
- Off-axis K bond (psi = 0.3), lemma_check.out: n_z single-signed, tails B_R = 0.36 and B_L = 1.03, both positive, U_alt - U_same > 0 at T = 10 to 22, like binding.
- scatlen.out: the K_0 formula is reproduced to 1e-6, and the sign rule holds in the K-like, U-like and no-symmetry cases.

## Construction side (coordinator's request)

PROOF section 10 covers this. Accidental degeneracy in N is allowed but must lie a finite distance from K and U. Oscillatory tails do flip the sign of the valley-odd part, but semiclassically they bind every channel. Easy-plane pseudospin requires J_z to change sign inside a non-binding two-valley region, and the window is exponentially thin. J1-J2 has Delta of order one only because magnons are point-like. The pseudospin scale at eps = 0.1 and L <= 256 is about 1e-7, too small to test. The remaining target is a direct lattice search at eps of order one.

## Narrowest remaining assumption (A_sol) [assumed]

No accidental two-valley soliton onset in N or T far from K and U, in particular at eps of order one, combines non-binding interactions with J_z > -|J_perp| (including the undetermined mass-imbalance correction at order n^2).
