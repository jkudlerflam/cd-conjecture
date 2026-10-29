# Referee report on C80 (12-letter bond), C82 (coloured vacancies, F3, F4, n = 16) and C83 (12-letter coloured lift, n = 13)

## Claims under review

- Astra C80 (`codex/work_C80/`) compresses C75's 15-letter route-A bond to 12 letters, L_0, L_1, R_0, R_1, D_A, D_B, X_A, X_B, E_0, E_1, M_A, M_B. Tokens carry no bit, the displaced bit sits on the next site and the last bit on the empty marker, the vacuum marker is identified with E_0 and the vacuum token with R_0. A certified n = 12 member (24 simple Fermi points, positive margins on all 300 pair sums) gives the conditional frozen prediction c = 12 > 11. C80 also states that one vacancy gives a 13-letter candidate with c = 13 > 12.
- Astra C82 (`codex/work_C82/`) repairs F5 of `contrib/referee_C75` with two vacancy colours V_0, V_1 carrying the grammar state (d = 17 on the 15-letter bond), proves F3 (nonzero-port sector selection) and F4 (a local Schrieffer-Wolff layer for the integer penalty), and certifies an n = 16 member, so c = 17 > 16.
- Astra C83 (`codex/work_C83/`, added to this review on request) applies the two-colour repair to the 12-letter bond (d = 14), shows that a single uncoloured vacancy fails because R_0 V L_0 hides an illegal junction, transfers F3 and F4, and certifies an n = 13 member, so c = 14 > 13.

## Status

All three finite constructions hold up. I found no fatal gap. The 12-letter bond rebuilt from C80's text matches `bond_shift.py` entry by entry, its legal cells, code, first and second order and small ED all reproduce, and C80's admitted payload-dependent splitting of the vacuum code turns out to be smaller and simpler than C80's bound allows. It is exactly a one-body field on the bits at order 1/K, of size about 0.5/K per bit, so it is not even an interaction at leading order. The colouring of C82 is correct, unique on the legal space, and gives K^-2 suppression in an independent test, including a test of the 12-letter version with nonzero bits that C83 did not run. C82's F3 proof checks line by line and also covers C80 and C83. Its F4 normal form is the standard local Schrieffer-Wolff construction and is adequate at the physics standard. The n = 12, 13 and 16 certificates reproduce bit for bit and pass independent floating-point spot checks.

One claim does not survive. C80's statement that one uncoloured vacancy gives a 13-letter candidate with c = 13 > 12 inherits F5 unchanged (illegal material weight 0.023 at J = t, independent of K). C83 itself proves the single-vacancy route fails for this grammar, and its d = 14 construction replaces it. The answer to the question whether C82's colouring applies to the 12-letter bond is yes, and C83 has done it with a certified n = 13 member, giving the conditional d = 14, c = 14 candidate. That is the strongest conditional statement in the campaign so far, and it rests on the same infrared assumptions as every earlier route-A candidate.

## How it was established and checked

- Author: Claude (Opus 5.5) referee subagent for Jonah Kudler-Flam, 2026-10-29. No Codex call was made. All computations are local, NumPy/SciPy only, on resolved diagnostic parameters except where integer certificates are reproduced.
- `my12.py` rebuilds the 12-letter bond from C80 PROOF Section 1 alone (letter names as strings, my own word decoder and closure builder). `check_bond14.py` builds the 14-letter lift from C83 PROOF eq. (3) on top of it. The Astra bond files were imported only as objects under test.
- All supplied scripts of C80, C82 and C83 were rerun in a scratch copy of the repository layout outside the repository, and the outputs compared with the committed files.
- Tags: [proved] means a finite argument I checked line by line or an exact identity, [numerical] means floating-point evidence at resolved parameters, [heuristic] means a physics argument not made rigorous, [assumed] means an input taken from elsewhere.

## Part I. C80, the 12-letter bond

### 1. The bond (`check_bond12.py`, `check_bond12.json`)

[numerical, exact in floating point] My rebuild matches `bond_shift.bond(p, 12)` entry by entry, maximum difference 0.0, at four generic points (n = 2 and 12, lambda = 1, 0.37, 0.21, K = 37 to 2^20, nonzero field, mu_sel, t_d, g and generic flux). Both have 205 nonzero entries, as C80 states. The 144 x 144 matrix is exactly Hermitian, commutes exactly with the on-site charge sum, and commutes exactly with the marker projector on each of the two sites. psi_A is tied to n phi_A + pi/2 in the code, as in C75. At the certified point this equals pi, as PROOF says.

### 2. Legal cells and the shifted bit (`check_cells12.py`, `check_cells12.json`)

[numerical, exhaustive] Brute-force enumeration over all 8^n data words and 4 terminal markers, after a frozen marker, gives exactly (5n + 3) 2^n zero-penalty words for n = 1 to 7. Each has at most one excited letter. My decoder (first non-L data letter is the head, R_0 there is the vacuum, later R letters carry bits, the terminal E carries the last bit) is a bijection onto (position, species, bits) and inverts PROOF eqs. (2)-(3) on every word. So the identifications E_0 = M_v and R_0 = D_v are injective on legal words, as claimed.

[proved] I checked rule by rule that every off-diagonal term changes only the bad-pair projector immediately to its left, and only for the vacuum move R_0 R_a -> L_a R_0 or R_0 E_a -> L_a E_0 with an R, D or X predecessor. All markers have the same set of allowed right neighbours, so the marker-to-first-letter penalty depends only on the right cell. This is what C80 Section 6 and C83 Section 6 use.

[numerical] Marker-free zero-penalty rings of 2 to 6 sites are only the all-L and all-R rings (`marker_free` in the JSON), as expected. No new zero-penalty configuration arises from R_0 doubling as the vacuum token or from the shared E_0 inside a frozen cell. The cost of the identifications is entirely off-grammar, as described next.

### 3. The payload-dependent vacuum splitting (`check_cells12.json`, `payload_splitting_*`)

C80 admits that at finite K the 2^n vacuum code is not degenerate, because the vacuum hop rule also acts on bit-carrying R_0 letters in the tail and creates a bad pair. I computed the order-1/K effective operator G^dagger P H Q H P G exactly for n = 1 to 11 and compared with ED of the complete zero-excitation cell space at n = 2, 3, 4.

[numerical, exact to rounding] Every intermediate state has exactly one bad pair, and the operator is diagonal in the bits. Its Walsh decomposition has no component of weight 2 or more (below 3e-16 for n up to 10). So at order 1/K the splitting is exactly a one-body field on the bits, sum_i w_i (1 - a_i)/K, which is a quadratic (chemical-potential) term in the fermion language, not an interaction. [proved, short argument] This follows because each bad intermediate is produced by one hop of one bit-0 letter, with amplitude set by the vacuum profile to its left, so the energy shift is linear in the bit indicators.

[numerical] The single-cell spread is K^-1 times 0, 0.25, 0.638, 1.083, ..., 4.509 for n = 1, ..., 11, growing roughly like (n - 1)/2. ED at K = 10^3 to 10^5 agrees with these coefficients to O(1/K^2). At n = 2 the two-cell spread is 0.00495 at K = 100, which is C80's number. At n = 12 the two-cell spread is about 10/K, far inside C80's bound 2^36/K.

[proved, arithmetic] At the certified K = 2^2048 the splitting is about 2^-2045. The physical band scale is lambda^2 E_unit, about 2^-1668, and the physical field is about 2^-1780. The splitting is therefore below the band scale by a factor of order 2^-377 and below the logical field by about 2^-265. C80's budget, which divides 2^80/K by lambda^2 E_unit, is correct in structure and very conservative. So the splitting is really below the lambda^2 band scale, as asked, and at leading order it does not even generate an interaction.

### 4. First order, second-order structure, ED (`check_cells12.json`)

[numerical] On the projected two-cell space at K = infinity, G^dagger V G = 0 exactly and V G lies entirely in the (1,1) excitation sector, for n = 1, 2, 3. At finite K = 100 the first-order block is also exactly zero. [proved] This holds at every K, because H_0 preserves each cell's excited-token number even off grammar.

[numerical] My H_2 at K = infinity agrees with the C75 18-letter parent `bond18.pair` H_2 in the bit basis to 8e-18, 4e-18 and 2e-18 absolute at n = 1, 2, 3 (entries of order 1e-2). The g^0 term is an exact scalar. The g^2 coefficient, extracted by Richardson extrapolation, is Hermitian and Gaussian on the full 4^n pair code to relative residual 1.6e-10, 1.1e-9, 3.4e-9. These residuals are the O(g^4) extrapolation error, not a non-Gaussian term. The exact parent agreement plus the C75 referee's 7e-16 Gaussian check of the parent are the stronger statement. This confirms C80's claim that the 12-letter g^2 operator is C75's, with no change of scale.

[numerical] ED of H_0 + lambda V on the invariant closure of the two-cell code at K = 100 (dimensions 160, 3302 and 54118 for n = 1, 2, 3, illegal words included) agrees with the quasidegenerate formula (7) with error over lambda^4 equal to 0.002075, 0.002073, 0.002073 (n = 1), 0.000731, 0.000730, 0.000729 (n = 2) and 0.000257, 0.000257 (n = 3) for lambda = 0.2, 0.1, 0.05. Residuals are below 4e-13. The n = 3 run is new. My closures are seeded with the code only, so they are smaller than C80's.

### 5. Selector at zero port (`check_selector12.py`, `check_selector12.json`)

[proved, checked line by line] C80 Section 6 carries over from C75. The vacuum moves that create a bad pair split off as partial isometries a_j with a_j a_j^dagger <= Pi_bad,(j-1) and at most k - 1 of them. Each legal fragment has at most one mobile head, because a later R_0 has an R predecessor and its move belongs to some a_j. The Young inequality with eta = K - 4 gives eq. (16).

[numerical, complete Hilbert spaces] Single cells with k = 0 to 5 data sites (every data letter at every site and every terminal marker, so all illegal and multi-token sectors, dimension up to 131072) at K = 6, 10, 100, with Delta = 2.2 > 2 + g, all satisfy E_v(k) - (k-1)_+/(K-4) <= E_K(k) - mu_sel <= E_v(k). Example, k = 5, K = 10: E = -1.9715 within [-2.4686, -1.8019]. Marker-free rings of 2 to 5 sites satisfy E >= -N/(K-4), nearly tight at K = 100 (N = 5: -0.049990 against -0.052083).

[numerical, reproduced] `selector_checks.py` reproduces mu_sel = 1.8561806413628503, density gap 5.47433e-5 and e_13 < -0.0072058 from rationals.

### 6. Reruns and certificate spot check (`spot_C80_n12.json`)

[numerical, exact reproduction] `run_certificate.py`, `checks.py`, `selector_checks.py` and `grammar_audit.py` ran without error in the scratch copy (about 9, 1, 5 and 0.3 s). All certificate files are bit-identical to the committed ones (certificate, exact root certificate, full-loop bands, finite budget, all five .npz integer files), except for timing fields. `small_checks.json` and `selector_certificate.json` differ only in floating-point noise at the 1e-14 level.

[numerical, independent] My floating-point evaluation of eq. (5) at n = 12, D = 1000, t_d = beta 2^-66 reproduces the 512-bit integer F and M matrices entry by entry up to one common factor beta^2/4 for both, with relative spread 2e-15 over a dynamic range of 9e25. A Bloch scan of the whitened matrices on a 2 x 10^5-point grid finds exactly 24 inertia changes, 0 up to 12 and back to 0, each within one grid step of a reported root. Recomputed margins from the reported roots are delta_2 = 0.015097 and delta_4 = 5.67e-7, consistent with the certified lower bounds 0.01509618 and 5.5894e-7.

[proved, arithmetic checked] The finite-budget logs match eq. (13). For example the penalty term is log10(192) + 33.011 + 80 log10 2 - 2048 log10 2 + 960 log10 2 + 213.208 = -54.94, as printed.

### 7. Do F3, F4, F5 recur for the 12-letter bond?

- F3 [proved]. It recurs as an open item in C80 (Section 9 calls nonzero-port selection open) and is closed by C82's argument, which I checked (Part II). The only bookkeeping change is that a marker followed by R_1 or E_1 is now a bad pair, but all markers share one outgoing state, so the term belongs to the right cell. C83 Section 6 states this correctly. With delta > 5.47e-5 at b = 13 and lambda = 2^-480, gamma > 5 sqrt(2) lambda by a huge margin.
- F4 [proved at the physics standard]. It recurs (the 12-letter H_0 does not preserve the legal space) and is covered by C82/C83's local Schrieffer-Wolff layer. The new feature is that the order nu^2/K term is not a scalar on the code at lambda^0. Its leading part is the one-body bit field of item 3 and is included in the budget.
- F5 [numerical, `check_vacancy12.py`]. It recurs unchanged in C80's uncoloured 13-letter lift. On two n = 3 cells with bits (1,0,0) and one vacancy, the illegal material weight at J = t is 0.0393, 0.0239, 0.0231, 0.0230 for K = 10 to 10^4, so K-independent, as for C75. At J = 0.01 t it is 3e-6 at K = 10^4. A typical hidden wall is R_0 V L_a, produced by a vacuum hop on a bit-carrying R_0 next to a vacancy, which is C83's example.

## Part II. C82, coloured vacancies, F3, F4, n = 16

### 8. The colouring (`check_vacancy12.py`, `check_vacancy12.json`)

[proved, numerical check] Both the 15-letter and 12-letter grammars factor exactly through a two-state automaton, allowed(a, b) iff o(a) is in I(b), for all 225 and 144 ordered pairs. For the 15-letter grammar my automaton equals C82's table up to relabelling the two states. For the 12-letter grammar it equals C83's table (R_1 and E_1 accept only state 1). C82's own table would not be correct for the 12-letter bond, which C83 notes.

[numerical, exhaustive] On all rings of 3 to 5 sites over the 14 letters, and 3 to 4 sites over the 17 letters, every zero-penalty coloured word has a legal squeezed material word, every uncoloured occupancy pattern of a legal material word has exactly one zero-penalty colouring, and N_material_bad <= N_bad holds without exception. [proved] The induction along the vacancy string in C82 Section 2 is correct and gives these statements for every length.

[numerical, independent of C82's code] My own coloured C75 lift reproduces C82's numbers to all printed digits: illegal weight 0.0100052, 1.07551e-4, 1.08213e-6, 1.08279e-8 at J = t, K = 10 to 10^4, and the uncoloured plateau 0.0256. At J = 0.01 t the coloured weight is 2.7e-6 at K = 10^3 and 1.4e-8 at K = 10^4, also falling like K^-2 once JK >> t. The 12-letter coloured lift with nonzero bits (1,0,0) gives 7.5e-3, 8.31e-5, 8.39e-7, 8.40e-9 at J = t, again K^-2. [proved] C82 eq. (4) gives the fixed-volume K^-2 bound and states correctly that it is not volume-uniform.

[proved, projected hopping] C82's argument that each projected swap has amplitude -t with no multiplicity factor is correct. The off-grammar reverse hops are kept in the bond, as they should be.

### 9. F3 (C82 Section 5)

[proved, checked line by line] The almost-additivity estimate |A_m - m u| <= B_port follows from cutting q - 1 ports in A_(qm). S_m <= m u + 2B_port uses one closing port and a product trial state. On a compatible ring with k wrong cells, removing at most 2k ports and bounding at most k runs gives -3k B_port, and E_l(K) >= l e_b + l gamma with u/b <= e_b gives (17). Since every cell length is at least one, gamma sum l >= k gamma. Marker-free rings are excluded by the -1/(K-4) density bound. The port norm sqrt(2)|lambda| is right because the stars for different input words are disjoint. The argument is complete for the frozen chain and makes no use of conserved bad-pair labels. It does not address the mobile chain, as C82 says.

### 10. F4 (C82 Section 4)

[proved at the physics standard] The integer Fourier decomposition A_m under exp(i theta N), the support enlargement by one site at each end, |m| <= |X| + 1, the generator S_1 = sum B_m/(kappa m), and the compression (10) with its negative sign are all correct. The volume-independent remainder argument is the standard connected-cluster bound of local Schrieffer-Wolff theory (Bravyi, DiVincenzo and Loss, and later work on local perturbation theory for commuting projector penalties). It is sketched rather than written out, which is acceptable here because the setting (commuting integer penalty, bounded-range perturbation, fixed order) is the textbook case. My order-1/K payload computation in item 3 is an instance of eq. (10) with m = 1 and matches ED to O(1/K^2). C82 is careful to claim neither an all-order block diagonalization nor a uniform gap of the N = 0 band. `check_sw.py` reproduces K^2 ||R_2|| = 7.34 to 7.42.

### 11. The n = 16 certificate (`spot_C82_n16.json`)

[numerical, exact reproduction] `certify_n16.py`, `budget_n16.py`, `check_vacancy.py`, `check_full.py`, `check_parent.py` and `check_sw.py` all ran. Every output is bit-identical to the committed one except timing fields.

[numerical, independent] At n = 16, t_d = beta 2^-92 my floating-point formula reproduces the integer F and M with common factor beta^2/4 and relative spread 3e-15 and 8e-15. The Bloch scan finds 32 inertia changes, 0 to 16 and back, within one grid step of the reported roots. Recomputed delta_2 = 0.014530 and delta_4 = 6.01e-7, consistent with the certified 0.0145298937 and 5.92937e-7.

[proved, arithmetic checked] budget_n16: log10(256) + 44.441 + 64 log10 2 - 1280 log10 2 + 281.473 = -37.73, as printed. C75's 2^64/K argument has no dependence on n except through B_0 = 2^16 and the vacuum gap, so it transfers to n = 16 as claimed.

### 12. Does C82's colouring apply to the 12-letter bond?

Yes. The 12-letter grammar is a two-state automaton (item 8), so two colours suffice and the same induction proves uniqueness and N_material_bad <= N_bad. This gives d = 14 with c = n + 1, and a violating count needs n >= 13. C83 has carried this out. Two differences from C82 matter.

- [proved] The incoming sets now depend on the bit (R_1 and E_1 accept only state 1). A reverse colour hop from R_0 or E_0 has a second, illegal channel, while one from R_1 or E_1 has none. This produces a bit-dependent virtual energy of order t^2/(JK) per vacancy contact, which couples the payload to the vacancy fluid. [numerical] In the coloured 12-letter ring, K times the energy difference between all-0 and all-1 tails is -0.997, -1.000, -1.000 with frozen vacancies (t = 0, C80's own splitting) and -1.290, -1.301, -1.302 with t = J (K = 10^2 to 10^4). The extra part scales as 1/K as predicted. It is covered by C83's eqs. (22)-(23), which choose K large compared with nu^2/(J^2 tau), so this is not a gap, but it is a new term that C83 does not mention explicitly.
- [proved] The finite-K budget must divide by lambda^2 E_unit (C80's structure), which C83 does.

## Part III. C83, the d = 14 lift

### 13. Bond and colouring (`check_bond14.py`, `check_bond14.json`)

[numerical, exact in floating point] My 14-letter bond built from C83 eq. (3) on my 12-letter rebuild matches `bond14.bond` entry by entry, maximum difference 0.0, at three generic points (n = 2, 3, 13, varied J, t, mu, lambda, K), for both the coloured 196 x 196 bond and the uncoloured 169 x 169 variant. C83's automaton table reproduces C80's complete allowed-pair table. The matrix is exactly Hermitian and commutes exactly with internal charge, occupied number and marker number.

[proved] The restricted no-go for one uncoloured vacancy is correct. R_0 V R_0 and L_0 V L_0 occur in legal embeddings, so R_0 V and V L_0 must both cost zero, and R_0 V L_0 then hides the forbidden R_0 L_0 pair. It excludes this repair, not every d = 13 design, as C83 says. The F5 numbers in item 7 show the uncoloured failure dynamically.

[numerical] K^-2 suppression: C83's zero-bit tests coincide with C82's (and with mine). My nonzero-bit 12-letter test in item 8 adds a case where R_1 and E_1 restrictions matter, with the same K^-2 law.

### 14. Reruns and the n = 13 certificate (`spot_C83_n13.json`)

[numerical, exact reproduction] `certify_n13.py`, `budget_n13.py`, `check_lift.py`, `check_frozen.py` and `check_sw.py` ran (about 11, 4, 1, 0.5 and 0.3 s). All outputs are bit-identical except timing fields.

[numerical, independent] At n = 13, t_d = beta 2^-74 my floating-point formula reproduces F and M with common factor beta^2/4, relative spread 2e-15 and 7e-15. The Bloch scan finds 26 inertia changes, 0 to 13 and back, within one grid step of the reported roots. Recomputed delta_2 = 0.014909 and delta_4 = 9.55e-8, consistent with the certified 0.0149087301 and 8.70283e-8. The pair-sum margin is the smallest of the four certified members, but it is a certified positive lower bound and the budget is applied to it.

[proved, arithmetic checked] budget_n13: penalty log10(208) + 35.635 + 84 log10 2 - 2304 log10 2 + 1152 log10 2 + 243.594 = -39.95, conversion -32.73 and local lambda term -49.59, as printed. The constant 2^84/K follows C80's pattern with 4^13 columns. Given item 3, the true constant is many orders smaller. The mobile requirement (22) is symbolic and correctly stated as such.

[not checked] C83 Section 9's harmonic bookkeeping (k_phys = Q_h/28, delta_phys > 3.10815e-9 = delta_4/28) depends on C81's Gaussian dictionary. I only checked the arithmetic 8.70283e-8/28 = 3.108e-9. It is an [assumed] input for this report.

## Gaps, classified

| # | Gap | Class |
| :-- | :-- | :-- |
| G1 | C80 RESULT and STATE claim a 13-letter vacancy candidate with c = 13 > 12. One uncoloured vacancy inherits F5 for the 12-letter bond (illegal weight 0.023 at J = t, K-independent), and C83 proves the single-vacancy repair impossible within this design. The claim should be withdrawn and replaced by the d = 14 lift of C83. | serious (C80's d = 13 claim only) |
| G2 | The C80 finite-K constants (2^36/K, 2^44/K, 2^80/K) are derived loosely in PROOF Section 5, by allowance of factors rather than line-by-line estimates. Direct computation shows the actual splitting is about 0.5(n - 1)/K per cell and purely one-body, so the bounds hold with enormous slack. A cleaner proof would use the exact order-1/K operator of item 3. | cosmetic |
| G3 | The root routines record inherited intermediate finite_g_power thresholds (-205, -264, -219) different from the final schedules (-216, -288, -256). All three folders say so, but the JSON fields can mislead a reader, as F1 did for C75. | cosmetic |
| G4 | C82/C83's local Schrieffer-Wolff remainder bound is a sketch of the standard connected-cluster estimate, not a full proof. Adequate at the physics standard; a citation would help. | cosmetic |
| G5 | In the 12-letter coloured lift the colour hopping produces a bit-dependent virtual energy of order t^2/(JK) that couples payload and vacancy fluid (item 12). It is controlled by C83's choice of K but not named there. | fixable (a sentence) |
| G6 | Mobile finite-J matching, the dressed-defect comparison with a gapless vacancy fluid, contact-dependent material onsite and field terms, and multichannel Luttinger stability. These are stated as assumptions by C82 and C83 and are inherited from the d = 24 lift. | serious, inherited (not new) |
| G7 | Frozen-chain infrared inputs (bosonization of the 2n Fermi points, volume-uniform matching, the frozen-cell reading of "CFT", no possible numerical test at the certified hierarchy). Unchanged from C75. | serious, inherited (not new) |

## Verdicts

- C80 (12-letter frozen, n = 12, c = 12 > 11): the finite construction is correct. The bond matches its text exactly, the legal-word map is a bijection, the second-order operator equals C75's, the certificate reproduces and spot-checks, and the admitted vacuum splitting is real but harmless (one-body, and below the lambda^2 band scale by about 2^-377). F3 and F4 recur and are closed by C82's arguments. Its d = 13 vacancy claim is not supported (G1, serious). The frozen d = 12 candidate stands on the same conditional footing as C75.
- C82 (d = 17, n = 16, c = 17 > 16): the colouring repair of F5 is correct, with unique legal colouring and K^-2 suppression reproduced independently. F3 is closed for the frozen chain. F4 is supplied at the physics standard. The n = 16 certificate reproduces bit for bit and spot-checks. No fatal or new serious gap.
- C83 (d = 14, n = 13, c = 14 > 13): the bond matches its text exactly, the two-colour automaton is the right one for the 12-letter grammar, the single-vacancy no-go is correct, leakage falls as K^-2 including with nonzero bits, and the n = 13 certificate reproduces and spot-checks. One new term (G5) needs a sentence. This is now the smallest-d conditional candidate, at d = 14 with predicted c = 14.

None of these is an established counterexample. All of them still need the inherited infrared inputs (G6 for the mobile lifts, G7 for everything).

## Reproducing

All scripts run from this folder with Python 3, NumPy and SciPy. `SCR` is a scratch directory with copies `SCR/codex/work_C80`, `SCR/codex/work_C82` and `SCR/codex/work_C83` of the repository folders, after their own certificate scripts have been run there.

- `python3 check_bond12.py SCR/codex/work_C80` takes a few seconds.
- `python3 check_cells12.py SCR/codex/work_C80` takes about a minute.
- `python3 check_selector12.py` takes about 25 minutes, mostly the complete five-site cell spaces.
- `python3 check_vacancy12.py` takes a few minutes.
- `python3 check_bond14.py SCR/codex/work_C83` takes a few seconds.
- `python3 spot_cert.py SCR/codex/work_C80 12 66 C80_n12`, `... work_C82 16 92 C82_n16` and `... work_C83 13 74 C83_n13` take a few minutes each.

`my12.py` is the referee rebuild of the 12-letter bond, the word decoder and the closure builder. Each script writes the JSON file of the same name (`spot_cert.py` writes `spot_<tag>.json`).

## Relation to the rest of the repository

This addresses the "no independent referee" notes in `codex/work_C80/RESULT.md`, `codex/work_C82/PROOF.md` and `codex/work_C83/RESULT.md`. It resolves F3, F4 and F5 of `contrib/referee_C75` in the sense stated above, and it adds G1 (C80's d = 13 claim) and G5 (a new bit-dependent colour term). It does not change the conditional status of any route-A candidate.
