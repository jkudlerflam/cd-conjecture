# Route A compile: a charge-conserving carrier interface and the same-schedule condition

## Claim

Chains are one-site translation-invariant, nearest-neighbour and Hermitian, built from spinless fermions on a fixed on-site alphabet (a spin chain after Jordan-Wigner), with exact fermion number and parity. The target is route A's second-order logical chain with c = n.

1. [proved] A frustration-free, charge-conserving interface (fermionic Bell-dimer ports whose excitations are mobile carriers) gives an exactly degenerate 2^n code, zero first order, exact fermion number and parity at all orders, and logical-operator-valued exchange, escaping C18's scalar-resolvent obstruction.
2. [proved structure, numerical scaling] The order-λ² bond is Gaussian at leading order in the conversion amplitude τ, with relative corrections of order τ².
3. [proved] Same-schedule lemma. A Hermitian definite junction hopping matrix T with on-site matrix M small in the T metric gives n partially filled bands, so diagonal (Jacobi) filters are not needed.
4. [numerical, as in C21 and C22] The plain cell P D^n P fails, with at most two gapless bands for every n tested.
5. [numerical] A central-switch cell gives a symmetric definite T and the Gaussian count c = n at n = 1 to 4, and a small arm asymmetry plus a phase displaces the Fermi seas (n = 2).
6. The frozen-pattern alphabet is d = 48. Equal arms and the cell length still need a selector.

## Status

Tags as stated. The infrared conclusion c = n is the Gaussian band count plus the weak-coupling input of `contrib/ir_control/` [assumed]. Upgrading it needs an equal-arm selector, control of the small scales of Section 6, and an IR check with displaced centres.

## How it was established and checked

- Author: Claude (Opus 5.5) subagent, 2026-10-28, for Jonah Kudler-Flam. No Codex.
- Checked by: internal cross-checks only (two extractions of the bond, truncation tests, τ and λ scaling).
- Inputs reused: C19's grading and Pauli-port bond, C18's obstruction, the target and pairing warning of `contrib/ir_control/` and `contrib/routeA_numerics/`, and the C7 selector idea. Astra C21 (`codex/work_C21/`, a Clifford cursor with J = −2RG_f) and Astra C22 (`codex/work_C22/`, a 17-state bus with T = [[a,b],[c,0]]) found the reversal failure of Section 2 independently. A folded-bus cell with an express track is being built in parallel in `contrib/routeA_spatial/`. Its switch idea is close to Section 4, and I have not compared numbers.

## Argument

### 1. The interface

Every site carries a dimer of two orbitals (β, α). Its four states are ∅ = (1,0), the code state, and three carriers T = (0,1) (neutral), H = (0,0) (charge −1) and D = (1,1) (charge +1). Data sites add one orbital c, the logical mode. At a port, β and α are the bonding and antibonding orbitals of a physical pair (p, s), so the code port state is a U(1) Bell pair. The junction bond is

    V = (p_R† p_L + h.c.) + κ (n_pR − ½)(n_pL − ½),     p = (β + α)/√2.

On the code p_L gives H at the port, p_R† gives D, and n_p − ½ gives T. All three have zero code expectation, so G†VG = 0 [proved].

Cell terms are on-site energies e_T, e_H, e_D, nearest-neighbour hops of each carrier type (T: β_x†α_x α_{x+1}†β_{x+1}, H: β_x†β_{x+1} on empty α's, D: α_{x+1}†α_x on full β's, amplitudes −w_X), and two conversions on data sites:

    deposit  −τ_A (c† α (1 − n_β) + h.c.)    T + c0 ↔ H + c1,
    pickup   −τ_B (β† c n_α + h.c.)          T + c1 ↔ D + c0.

Every term annihilates the code, so the cell is frustration free with an exactly degenerate 2^n code [proved]. C18's obstruction assumes that the only positive terms touching a port are static checks, so that a source F G evolves like the code. Here the port dimer trades carriers with its neighbour, so F G propagates into the cell and meets the data. The carrier gap is positive for e_X > 2|w_X| on a path, slightly less on the switch graph below (numerically 0.30 to 0.40).

In a second-order process V creates a carrier in each cell, the carriers move a fermion between a data site and the port, and V closes. For example n − ½ creates T in both cells, T deposits at a in cell j and picks up at b in cell j + 1, then p_R† p_L closes. Carrier hops are bilinears commuting with data operators, so at leading order each cell factor is a linear combination of c_a† or c_b [proved]. Corrections need three conversions, hence the τ² suppression. With w_D = −w_H the H and D branches cancel in the order-τ² one-cell terms, so M = O(τ⁴) [numerical: M/T falls by 4 when τ halves].

The poles are single-carrier eigenvalues on the cell graph, about three per site, so the pole bound 2n ≤ 4(k − 1) holds automatically, and the coefficients come from carrier propagators rather than Jacobi polynomials.

### 2. Why the plain cell fails

Write the bond as Σ_j c_j† T c_{j+1} + h.c. + c_j† M c_j. The right port of cell j sees data site x at depth n + 1 − x, the left port of cell j + 1 at depth x. Amplitudes decay with total depth, so T is a reversal times a decaying kernel and the logical network is one zigzag through the junctions. Random scans over the ten cell and bond couplings give c_max = 1 in 124 of 125 valid samples at n = 2 and 66 of 66 at n = 3, and the weak-decay family gives 1 for n = 2 to 5 except one value 2 [numerical]. C21 proves the same for its cursor (Fermi points only at 0 or π) and C22 for its bus.

### 3. The same-schedule lemma

[proved] Let T = T† be definite and M = M†. With h(k) = M + T e^{ik} + T† e^{−ik} = T^{1/2}(N + 2cos k)T^{1/2}, N = T^{−1/2} M T^{−1/2} (use −T if T < 0). If N has distinct eigenvalues ν_a with |ν_a| < 2, det h(k) has exactly 2n simple zeros at cos k = −ν_a/2, every band crosses zero twice, and the Gaussian count is c = n. The zeros of det h as a function of z = e^{ik} are symmetric under z → 1/z̄, so simple zeros on |z| = 1 survive any small perturbation of T and M. Grading of T (deep flavours exponentially weaker) does not matter, only the ratio M/T does.

[proved, and confirmed numerically to 1e-17] The hop-first half of the second-order sum is the adjoint of the density-first half with the ports exchanged. If the right-port and left-port source factors coincide, T = Z + Z†. A port automorphism therefore forces T Hermitian, h(k) = h(−k), and all Fermi seas are centred.

### 4. The central-switch cell and its rule table

The cell is P_L L^m M R^m P_R with n = m. Lane E (express) runs over every site. Lane U (data) exists only on R sites, each carrying one data orbital. The only E-U link is the bond (M, R_1). Both ports reach U only through M along arms of equal length, and swapping the arms is a graph automorphism fixing U, so both ports see the data with the same schedule.

Letters and local dimensions: P_L, P_R, L and M carry one dimer (4 states each), R carries an E dimer, a U dimer and c (32 states). The frozen-pattern alphabet is d = 4 + 4 + 4 + 4 + 32 = 48. Letter types are conserved.

Bond rules, by ordered letter pair:

- (P_L, L), (L, L), (L, M), (M, R), (R, R), (R, P_R): carrier hops on lane E with −w_T, −w_H, −w_D (T hops may carry a phase, and the two arms may differ slightly, see below).
- (R, R): the same three hops on lane U.
- (M, R): the same three hops from E(M) to U(R).
- (P_R, P_L): λV.
- On-site: e_T [T] + e_H [H] + e_D [D] on every dimer, and the two conversions on R sites.
- Optional non-FF term: μ_c n_c on R sites, which splits the code at order J and sets the filling.
- All other ordered pairs get a classical penalty U_pat > 0. This makes (P_L L^a M R^c P_R)* the only zero-penalty words. It does not enforce a = c, and a selector must do that (Section 6).

[numerical] T is symmetric to 1e-16 and negative definite. Its eigenvalues relative to the largest are 7.7e-4 at n = 2; 1.2e-3 and 1.0e-6 at n = 3; and 1.2e-3, 2.0e-6 and 1.5e-9 at n = 4. Since M/T ∝ τ², the lemma's condition |ν_a| < 2 holds at μ = 0 for τ ≤ 0.01 (n = 2), 0.001 (n = 3) and 3e-5 (n = 4), and c_max = n there. For larger τ only a fine-tuned window remains.

[numerical] Centres at n = 2, τ = 0.003. With symmetric arms every phase tried leaves T Hermitian, and the seas sit at 0 and π, as the lemma predicts. With the left arm's hopping scaled by 1.05 and a phase 1e-3 on all T hops, T is non-Hermitian (defect 7.6e-4), there are still four Fermi points, and the centres move to −π + 0.0138 and −0.0137. The inter-flavour pair-tunnelling momentum is then 0.055 mod 2π, not zero.

### 5. Checks by exact diagonalization

Rows diagonalize H0 + λV on open chains of 2 or 3 cells, keeping at most K carriers (H0 conserves carrier number and V changes it by 0 or ±2, so K = 2 is exact through order λ³). The effective Hamiltonian is the des Cloizeaux polar form on the code, H2 = −PVRVP comes from a separate resolvent solve, and E is the operator-norm error against λ²H2. Parameters: e = 1, w_T = w_H = 0.3, w_D = −0.3, τ = 0.3, κ = 1, with a phase 0.3 on T hops in the uniform rows and none in the switch rows.

| cell | n | cells | K | dim | gap | E/λ³ (λ = .08, .04, .02) | traceless error at λ = .02 |
| :-- | --: | --: | --: | --: | --: | :-- | --: |
| uniform | 1 | 2 | 2 | 616 | 0.400 | 0.1042, 0.1041, 0.1040 | 6.2e-10 |
| uniform | 1 | 2 | 4 | 7636 | 0.400 | 0.1043, 0.1042, 0.1040 | 6.3e-10 |
| uniform | 1 | 3 | 2 | 2816 | 0.400 | 0.1999, 0.2041, 0.2060 | 8.6e-10 |
| uniform | 2 | 2 | 2 | 4432 | 0.301 | 0.1042, 0.1041, 0.1040 | 6.6e-10 |
| uniform | 2 | 3 | 2 | 40384 | 0.301 | 0.2000, 0.2042, 0.2061 | 1.3e-9 |
| switch | 1 | 2 | 2 | 2524 | 0.366 | 0.1039, 0.1038, 0.1037 | 6.6e-14 |
| switch | 1 | 3 | 2 | 11456 | 0.366 | 0.1995, 0.2036, 0.2054 | 8.9e-14 |
| switch | 2 | 2 | 2 | 22912 | 0.296 | 0.1039, 0.1038, 0.1037 | 1.4e-15 |

In every row the H0 code residual and the first-order block vanish exactly, as do all number-changing matrix elements (checked on two cells). The traceless error falls by 6.2 to 7.9 per halving of λ (the smallest switch values approach the double-precision floor, since the switch cell's logical scale is only 1e-10 to 1e-8), the remainder after λ³H3 scales as λ⁴, and the next level stays 0.29 or more above the code band. The non-Gaussian part of H2 relative to its traceless part is 0.215, 0.055 and 0.014 at τ = 0.3, 0.15 and 0.075 (uniform, n = 1). The one-body bonds of Sections 2 to 4 use exact resolvents in the sector {code} ∪ {one carrier per cell} with at most three data fermions, which is exact for the one-body part. A second extraction that removes the conversion-free scalar analytically agrees to 1e-6 relative.

### 6. What is still missing

1. Equal arms. The automorphism needs a = c. Unequal lengths were not tested, while a 5 percent difference in arm hopping already gives a non-Hermitian part of 1e-3. Two identical C7-type walkers, one per arm on a track blind to the lanes, should select equal lengths [heuristic]. A three-state C7 track on L and R letters gives d ≈ 4 + 4 + 4 + 12 + 96 = 120 [heuristic accounting]. Cell-length selection and the mobile wrapper are not addressed.
2. Scales. The ratio of smallest to largest eigenvalue of T falls by about 10^{-3} per flavour, and the thresholds above follow τ ≲ (ratio)^{1/2}. Third-order terms presumably need λ/Δ below the ratio itself [heuristic]. This is existence at each fixed n, as for C19's geometric design, with exponentially small velocities.
3. Infrared. Exact U(1) removes the routeA_numerics pairing problem. Displaced centres are shown only at n = 2, and the third-order and τ² non-Gaussian terms have not been checked against the margins of ir_control.
4. The three-cell ED of the n = 2 switch cell (207k states) was stopped because of machine load. The uniform cell has all four (n, cells) cases.

## Reproducing

numpy and scipy only, single core, minutes per case on an unloaded machine. `model.py` (uniform cell) and `model2.py` (lane cells) build the Hamiltonians. `run_checks.py` and `run_switch.py ed|ed3` give the ED tables (`checks_results.json`, `switch_ed*_results.json`). `onebody.py` and `onebody2.py` give the exact one-body bond and the Fermi-point count. `scan_c.py`, `uniform_weak.py`, `switch_scale.py`, `switch_phase.py` and `switch_break.py` give the scans of Sections 2 and 4 (`scan_n*.json`, `uniform_weak_results.json`, `switch_scale_m*.json`, `switch_phase_*.json`, `switch_break_m2.json`).

## Relation to the rest of the repository

This addresses route-A bottlenecks 1 and 2 and relaxes C19's target. A spatial cell does not need orthogonal filters, only the same schedule at both ports plus a small breaking for displaced centres. It explains the C21 and C22 failures and gives a candidate cell meeting the condition at n ≤ 4. It does not change "Where things stand", since the selector, the scales and the IR remain open.
