# Referee report on C75, the 15-letter route-A bond

## Claim

Astra C75 (`codex/work_C75/`) claims that the 23-letter route-A bond can be compressed to a 15-letter one-site translation-invariant nearest-neighbour bond on C^15 (alphabet L_a, R_a, D_(A,a), D_(B,a), X_(A,a), X_(B,a), M_0, M_A, M_B) with exact U(1). Its n = 15 member (cell length 16, D = 1000, t_d = beta 2^-84, g = 2^-256, K = 2^1152, lambda = 2^-640) has a quadratic second-order reference with 30 simple Fermi points, distinct sea centres and positive pair-sum margins, certified including finite-g and finite-K corrections. The conditional infrared prediction is c = 15 > d - 1 = 14 for the frozen chain, or c = 16 > 15 at d = 16 after a one-vacancy lift.

This report checks the finite-dimensional claims and asks what the compression costs.

## Status

The finite claims hold up. The bond, legal-cell structure, code, vanishing first order, the factor one half relative to C33, the Gaussian structure, the uniform zero-port selector at finite K, and the n = 15 certificate were all reproduced or independently checked. I found no fatal gap. The frozen d = 15 prediction stands on the same conditional footing as the 23-letter parent, with two fixable gaps (a nonzero-lambda sector exclusion and a many-body elimination of illegal words) that the 23-letter version did not need. The d = 16 vacancy version has one new serious gap. The exact grammar invariance that the 23-letter vacancy analysis leaned on is lost, and the domain-wall moves become unpenalised across a vacancy. The infrared conclusion remains conditional (heuristic) in both cases, as C75 itself says.

What would upgrade it is a written proof of the nonzero-lambda sector exclusion along the lines sketched in F3, a local Schrieffer-Wolff step for the K layer (F4), and a perturbative bound in J/t and 1/K for the vacancy-lift admixture of illegal material words (F5).

## How it was established and checked

- Author: Claude (Opus 5.5) referee subagent for Jonah Kudler-Flam, 2026-10-28. No Codex call was made. All computations are local, NumPy/SciPy only.
- The bond was rebuilt from the text of PROOF.md Section 1 alone (`my15.py`). `bond15.py` was imported only as the object under test, from a scratch copy.
- C75's own scripts were rerun in a scratch copy of `codex/work_C75` and `codex/work_C33` outside the repository, with the repository layout kept.
- Tags: [proved] means a finite argument I checked line by line or an exact identity; [numerical] means floating-point evidence at resolved parameters; [heuristic] means a physics argument not made rigorous; [assumed] means an input taken from elsewhere.

## Findings by task

### 1. The bond (`check_bond.py`, `check_bond.json`)

[numerical, exact in floating point] The rebuild matches `bond15.bond` entry by entry, maximum difference 0.0, at a generic n = 2 point (all phases, t_d, g, field and mu_sel nonzero, K = 100) and at an n = 15 point, for lambda = 1 and 0.37. Both have 374 nonzero entries. The matrix is exactly Hermitian and commutes exactly with the on-site charge sum. The comparison ties psi_A = n phi_A + pi/2 and psi_B = pi/2, as `bond18.py` does. PROOF Section 1 lists psi_A, psi_B as free, but the certified point satisfies this relation, so nothing changes.

[numerical] On a three-site open chain the total bad-pair number commutes with the 18-letter parent Hamiltonian (commutator norm 0) and fails to commute with the 15-letter one (norm sqrt 2). This is the structural price of the compression. C75 states it (PROOF Section 2, "proved limitation") and treats it correctly at fixed cell. Its consequences beyond a fixed cell are F3 to F5 below.

### 2. Legal cells, code, first order, second-order matrices (`check_cells.py`, `check_cells.json`)

[numerical, exhaustive] Enumerating every word of n data letters plus a marker, after a frozen marker, gives exactly (5n+3)2^n legal words for n = 1 to 5. Each has at most one excited letter, matching the list in PROOF Section 2.

[numerical] The projected zero-excitation cell Hamiltonian has ground energy -2cos(pi/(n+2)) and gap 2[cos(pi/(n+2)) - cos(2pi/(n+2))] to machine precision for n = 1 to 6. At finite K = 37 the full orientation space (all 2^n L/R strings, illegal ones included) has a spectrum independent of the data bits to 0.0 for every bit pattern. The code is therefore exactly 2^n-fold degenerate at every K, as claimed.

[proved, and numerical] First order vanishes at all K. Every term of H_0 preserves the number of excited letters in each cell separately (also off-grammar), and the port raises it by one in both incident cells. Numerically G^dagger V G = 0 exactly at K = infinity (n = 1, 2, 3) and at K = 100, 1000, 10^4 (n = 1, 2), and V G lies entirely in the (1,1) excitation sector.

[numerical] The g^2 coefficient of the two-cell second-order Hamiltonian, computed from the literal rebuilt bond on the legal two-cell space, was compared with PROOF eq. (3), with the folded graph transcribed independently. At the C75 phases (psi_A = pi, psi_B = pi/2, resolved D = 6, t_d = 0.07) and n = 1, 2, 3, the hopping T agrees with -2 alpha^4 Re C to relative error below 2e-15. The onsite matrix agrees with -2i alpha^4 Im A to 6e-11 (n = 2) and 1e-8 (n = 3) relative. That onsite matrix is only 6e-3 and 6e-5 of T there, and at n = 1 it vanishes, so the residual reflects cancellation in the diagonal subtraction rather than a discrepancy. The order-g^0 term is a scalar on the whole code. The Gaussian residual of the g^2 coefficient on the complete 4^n pair code (global Jordan-Wigner order) is below 7e-16.

[proved, and numerical] The factor one half follows algebraically. With a common carrier (W = U, Z = V in C33's `fast_loop.fast`), C33's hopping is -2 alpha^4 Re[(V^dagger R U)^T + U^dagger R V]. Since R is Hermitian, Re(V_j^dagger R U_i) = Re(U_i^dagger R V_j), so the two terms are equal and T_C33 = 2 T_C75. The onsite term doubles in the same way. C33's code can only realise a common loop with a real loop entry (loopflux = pi/2), so I checked that case numerically. With psi_A = pi/2, psi_B = 0 the literal T divided by C33's T is 0.5 entry by entry to 2e-15 at n = 1, 2, 3, with the same sign. The complex-loop case rests on the algebra plus the eq. (3) agreement above.

[numerical] At finite K, ||H_2(K) - H_2(infinity)|| times K is 0.446, 0.494, 0.499 (n = 1) and 0.346, 0.376, 0.379 (n = 2) for K = 10^2, 10^3, 10^4. This is the predicted O(1/K) with an O(1) constant, far inside PROOF's 2^64/K. Exact ED of H_0 + lambda V on the invariant closure of the code (dimensions 178 and 5338, illegal states included, K = 100) agrees with E_K + lambda^2 spec H_2(K) with error over lambda^4 constant at 0.001986 (n = 1) and 0.000726 (n = 2) for lambda = 0.2, 0.1, 0.05. My closures are smaller than C75's 436 and 15404 because I seed only the code, not every legal word. Both are invariant spaces.

### 3. Domain wall replacing the vacuum token: selection and code (`check_selector.py`, `check_selector.json`)

[proved, checked line by line] PROOF Section 8 is correct as far as I can tell. Token hops, conversions and the switch preserve every bad-pair projector. A vacuum move changes only the pair to its left, and only when the predecessor is R, D or X. These split off as partial isometries a_j with a_j a_j^dagger <= Pi_bad,(j-1) and at most k - 1 of them. A legal fragment has at most one head, with kinetic row sum at most 2. The Young inequality with eta = K - 4 then gives E_v(k) - (k-1)_+/(K-4) <= E_K(k) - mu_sel <= E_v(k). The cells decouple at lambda = 0 because the conversion's control letter is "any marker" and marker positions are exactly conserved.

[numerical, full Hilbert spaces] Single cells with k = 0 to 4 data sites (every word, all illegal and multi-token sectors, dimension up to 62208) at K = 6, 10, 100 all satisfy the two-sided bound. Example, k = 4, K = 10: E = -1.8443 within [-2.2321, -1.7321]. Marker-free rings of 2, 3 and 4 sites satisfy E >= -N/(K-4). At K = 100, N = 4, E = -0.039992 against the bound -0.041667, so the bound is nearly tight there. A run with Delta = 1.5, which violates the stated hypothesis Delta > 2 + g, still satisfied the bound at these sizes. The hypothesis is sufficient, not shown to be needed.

[numerical] The inherited midpoint selector gives mu_sel = 1.90158663762416, argmin of e_b at b = 16 over b <= 2000, and a minimum density gap of 2.077e-5 = 1.36 x 2^-16. This agrees with C75's rational certificate. Since 1/(K-4) < 2^-17, length 16 stays selected at lambda = 0 for all lengths and sectors. Code uniqueness survives as claimed.

### 4. Rerun of the certificates (`rerun_compare.py`, `rerun_compare.json`, `spot_n15.py`, `spot_n15.json`)

[numerical, exact reproduction] `run_certificate.py` (about 17 s), `checks.py` and `selector_checks.py` ran without error in the scratch copy. All n = 15 certificate files are bit-identical to the committed ones: certificate, root certificates, full-loop bands, rational bands, finite budget, and the five .npz integer files. `selector_certificate.json` and `small_checks.json` agree to floating-point noise.

[numerical, independent spot check] My floating-point evaluation of eq. (3) at n = 15, D = 1000, t_d = beta 2^-84 reproduces C75's 512-bit integer matrices F and M entry by entry up to a single common factor, beta^2/4 for both. The relative spread is 3e-15 for F and 6e-15 for M, over entries with a dynamic range of 7e32. So the hopping, the loop-induced onsite term and their relative normalisation are what eq. (3) says. A Bloch scan of C75's whitened matrices on a 4 x 10^5-point grid finds exactly 30 inertia changes, ascending 0 to 15 and back to 0. They sit within one grid step of the reported roots. Recomputing the margins from the reported roots gives delta_2 = 0.014554 and delta_4 = 1.85e-7, consistent with the certified lower bounds 0.0145540 and 1.1444e-7.

[proved, arithmetic checked] The finite-budget logs match the formulas. The penalty term is log10(16n) + 41.349 + 64 log10 2 - 1152 log10 2 + 255.902 = -27.89. The conversion term is -14.94. Both are far below 2^-35, and alpha^4 >= 64/17^6 > 2^-19 holds.

### 5. Does the huge K create problems the 23-letter version lacked?

[proved at fixed cell, numerical at small n] Not at the level of the finite certificate. K is the largest scale (K >> B_0 = 2^16 >> 1 >> g >> lambda), so every word with a bad pair sits at energy at least K - 3B_0 and no illegal state is near-resonant with anything retained. At lambda^0 the finite-K corrections are exactly bit-independent (spread 0.0 above), so they cannot compete with the logical field lambda^2 E_unit 2^-140. At lambda^2 they are O(lambda^2/K) (constant K times difference above). PROOF bounds them by 2^64/K at n = 15, which is 10^-27.9 in whitened certificate units. That covers the marker-port process through the illegal R M_A pair, which is kept rather than discarded. Higher orders in lambda through illegal words are smaller still. Retained scales are not undercut.

The real cost of large finite K is structural. See F3 and F4.

### 6. Comparison with the 23-letter route

- Frustration-freeness. [proved] Neither frozen route is frustration-free in the walker sector, since both have a sine-profile vacuum. In the 23-letter bond, however, the grammar penalty annihilates the legal ground state and the legal sector is exactly invariant. In the 15-letter bond the ground state carries an O(B_0^2/K^2) admixture of illegal words, and the vacuum energy shifts by O(1/K) (for example -0.0066 at n = 2, K = 37). The superselection by grammar is lost. Nothing in the frozen certificate needed it, because C75 replaced it by Sections 6 and 8.
- Mobility. [proved] Marker positions remain exactly conserved, so the chain is still frozen, with the same broken-translation reading as the parent. The mobile variant C32 was built on exact grammar invariance and would have to be redone. That is irrelevant at d = 15, since C32's alphabet is in the thousands.
- Vacancy lift. See F5. This is the one place where the compression loses something that matters.

## Gaps, classified

| # | Gap | Class |
| :-- | :-- | :-- |
| F1 | `root_certificate_*_n15.json` still records the inherited finite_g_power = -247, while the schedule uses g = 2^-256. C75's STATE flags this, but RESULT and the JSON can mislead a reader. PROOF lists psi_A, psi_B as free parameters, but the code ties them. | cosmetic |
| F2 | Section 6 (finite-K comparison) is proved only for an n = 15 two-cell problem with frozen markers, as C75 says. | cosmetic (correctly scoped) |
| F3 | Nonzero-lambda global sector exclusion. The parent has it (counterexample_campaign/route_a/README.md Sections 3 to 5), but its Section 5 needs every bad-pair projector to commute with H, which fails here (task 1). C75 proves selection only at lambda = 0 and calls the rest open "as in the parent", which understates the parent. A repair appears available [heuristic sketch, not written out in full]. Section 8 is an operator lower bound valid in every sector of a cell, so wrong cells obey E_l(K) >= l u/16 + gamma with gamma >= l(2^-16 - 1/(K-4)) >= 2^-17, using u <= E_16(K) <= E_16(infinity). The parent's boundary-energy lemma (2)-(4) then gives E_wrong - S_m >= k(gamma - 5B) > 0 with B = sqrt 2 lambda, illegal words included and with no splitting argument needed. Marker-free rings are covered by the -1/(K-4) density bound against e_16 < -2^-8. | fixable |
| F4 | Many-body elimination of illegal words. In the 23-letter chain the legal sector is exactly invariant. Here an extensive chain needs a local (Bravyi-DiVincenzo-Loss type) Schrieffer-Wolff layer at gap K with per-site perturbation B_0/K = 2^-1136 before the existing infrared argument starts. This is standard but not supplied, and it is one more layer in an already deep hierarchy. | fixable |
| F5 | Vacancy lift (d = 16). In the lift of COUNTEREXAMPLE_D24.md the bond acts only on adjacent occupied pairs. A move R R -> L R just to the right of a vacancy, or R M_0 -> L M_0, creates a material R-L pair with no penalty while the vacancy separates it, and the pair costs J K only when that gap closes. The 18-letter parent is immune, because its arrow-grammar moves preserve outer arrows. [numerical, `check_vacancy.py`] Ground states of small rings (2 cells, n = 3 with 1 vacancy, and n = 2 with 6 vacancies at half filling), started from legal states, show 0 illegal weight for the 18-letter lift at every K. For the 15-letter lift the illegal-material weight is 2.6e-2 (J = t) and 4e-6 (J = 0.01t), and these values are K-independent once JK >> t. For JK << t the weight scales like 1/K^2. The admixture is therefore controlled by J/t and 1/K, not by K alone. Under the D24-type hierarchy (J/t below 10^-500) it is quantitatively negligible [heuristic]. But the exact legal-material invariance used by the parent's sector proof and the vacancy matching premise does not transfer, and C75's "inherits the parent's vacancy matching problem" misses this new mechanism. | serious (d = 16 only) |
| F6 | Inherited inputs. These are bosonization and multichannel Luttinger stability of the 30 Fermi points (P), volume-uniform matching (M), the frozen-cell reading of "CFT", and an extreme hierarchy with no possible numerical test at the certified point. They are unchanged from the 23-letter candidate. | serious, inherited (not new) |

## Verdict

C75's finite-dimensional claims are correct as far as I can check. I rebuilt the bond independently and it matches entry by entry. The legal structure, the exact 2^n code at every K, the vanishing first order, the Gaussian g^2 structure and the factor one half relative to C33's common-loop matrices all reproduce. The uniform finite-K selector holds in full Hilbert spaces, and the n = 15 certificate reproduces bit for bit, with independent spot checks. The K = 2^1152 hierarchy creates no near-resonance with any retained scale. Nothing is fatal.

The compression is not free. It gives up exact grammar invariance, which costs the parent's nonzero-lambda sector proof (F3, fixable by the sketch above) and requires an extra Schrieffer-Wolff layer (F4, fixable). It also introduces an unpenalised wall-creation mechanism in the vacancy lift (F5, serious for the d = 16 claim only).

The frozen d = 15 chain is therefore a credible conditional candidate on the same footing as the 23-letter chain. It is not an established counterexample, because the infrared inputs remain heuristic. The d = 16 vacancy version needs F5 addressed before it can be quoted alongside the d = 24 lift.

## Reproducing

All scripts run from this folder with Python 3, NumPy and SciPy. `SCR` is a scratch directory holding copies `SCR/codex/work_C75` and `SCR/codex/work_C33` of the repository folders.

- `python3 check_bond.py SCR/codex/work_C75` takes about 3 minutes, mostly the dense three-site commutator.
- `python3 check_cells.py SCR/codex/work_C33` takes a few minutes.
- `python3 check_selector.py` takes a few minutes.
- `python3 spot_n15.py SCR/codex/work_C75` takes a few minutes. Run it after `run_certificate.py` there.
- `python3 check_vacancy.py` takes a few minutes.
- `python3 rerun_compare.py SCR/codex/work_C75 ../../codex/work_C75` must run after `run_certificate.py`, `checks.py` and `selector_checks.py` in the scratch copy.

`my15.py` is the referee rebuild of the bond and a generic closure builder. Each script writes the JSON file of the same name.

## Relation to the rest of the repository

This addresses the open item "The reduced alphabet has not undergone an independent referee check" in `codex/work_C75/RESULT.md`. It does not change the conditional status of the route-A candidates. It adds F5 as a new obstruction specific to lifting the 15-letter bond with a vacancy.
