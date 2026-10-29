# Referee report on C88 (two-colour mobile lift of the pinless 11-letter grammar, d = 13, n = 13)

## Claim under review

Astra C88 (`codex/work_C88/`) drops C84's boundary pin, keeps C84's eleven material letters L_0, L_1, R_0, R_1, D_A, D_B, X_A, X_B, E, M_A, M_B, and adds two vacancy colours V_0, V_1 through C83's automaton with E_1 deleted, so d = 13. It claims that the squeezed pin count F (L_1 followed by a marker in the material word) is conserved by every literal move, that on F = 0 the full finite-K material Hamiltonian equals C85's, so that C85's n = 13 material certificate applies unchanged, and that the conditional count is c = 14 > 12 = d - 1. Colouring is unique, the complete leading grammar correction is -(t^2 C + J^2 A + 2 J^2 lambda^2 Z)/(JK) with a halved colour-ambiguity counter C, illegal weight falls as K^-2, and frozen selection of F = 0 is proved by an extended F3 bound. Mobile selection of F = 0 is stated as an assumption S_pin (PROOF eq. (5)).

This is the construction suggested as finding H1 in `../referee_C84_C85/README.md`.

## Status

The construction is correct as stated. My rebuild of the bond from PROOF Section 1 matches `bond13.py` entry by entry with difference 0.0. F is conserved by every off-diagonal entry of the literal bond on every ring of 2 to 6 sites, legal or not, coloured or uncoloured, with port and vacancy hopping (12.7 million directed transitions at L = 6). On F = 0 the C88 and C85 ring Hamiltonians agree entry by entry at finite K and nonzero port on all 11^L words for L up to 5. Colouring is unique and the legal count matches the squeezed transfer-matrix count through L = 6. The grammar correction (10) is exact on every legal ring of 4, 5 and 6 sites. Leakage is K^-2 with nonzero bits, two vacancies and a nonzero port. The frozen-selection proof is correct. All supplied scripts rerun bit for bit, and the regenerated certificate is identical to C85's. I found no fatal or serious new gap.

The main new finding concerns S_pin. The pin part of S_pin does not need to be assumed in the form C88 gives it. A replacement argument (Part IV) shows that on any ring, at K = infinity and zero port, the lowest energy of a sector with pinned cells exceeds that of the pin-free sector with the same cell lengths (up to a port allowance of a few lambda per pin) by at least J times the top eigenvalue of a path matrix built from contact probabilities on each pinned cell. This holds for every vacancy number, every t/J and every cell length. Small-ring ED confirms the bound in every case and shows pin-pin repulsion, so pins do not proliferate. What remains of S_pin is a lower bound on the contact probability (a local, order-one quantity at fixed positive density) and the inherited wrong-length and marker-free comparison that every route-A mobile candidate already needs.

## How it was established and checked

- Author: Claude (Opus 5.5) referee subagent for Jonah Kudler-Flam, 2026-10-29. No Codex call was made. All computations are local, NumPy/SciPy only, at resolved diagnostic parameters.
- `my13.py` rebuilds the bond from C88 PROOF Section 1 alone (string letters, the o/I table, the move table, the diagonal u and eq. (1)), with my own pin counter, grammar counters and closure builder. It imports no Astra file. The previous referee's independent C84 rebuild (`../referee_C84_C85/my11.py`) serves as a second parent. Astra bond files (`bond13.py`, C84's `bond11.py`, C85's `bond14.py`) are imported only as objects under test.
- Supplied scripts were rerun in a scratch copy of `codex/work_C88` outside the repository.
- Tags: [proved] means a finite argument I checked line by line or an exact identity or exhaustive enumeration, [numerical] means floating-point evidence at resolved parameters, [heuristic] means a physics argument not made rigorous, [assumed] means an input taken from elsewhere.

## Part I. Bond, conservation, zero sector

### 1. Bond rebuild (`check_bond13.py`, `check_bond13.json`)

[numerical, exact in floating point] At four generic points (p = 2, 3, 5, 14, K from 11 to 2^20, generic lambda, port phase, flux, field, selector, J, t, mu), my coloured 169 x 169 and uncoloured 144 x 144 bonds equal `bond13.bond` with maximum entry difference 0.0. The coloured bond is exactly Hermitian, has 74 nonzero off-diagonal entries, and commutes exactly with N, M and Q. The allowed-pair tables agree on all 169 pairs.

[numerical, exact] At the bond level, h_11(C88) + K Pi_(L_1, marker) equals the Astra C84 bond (difference 0.0), the referee C84 rebuild (0.0), and the occupied block of C85's `bond14` divided by J (below 6e-14, rounding of J K/J).

### 2. Pin conservation (`check_pin13.py`, `check_pin13.json`)

[proved, exhaustive] I applied every off-diagonal entry of the coloured bond (material moves at nonzero g, t_d, port, vacancy hops with every colour, at a generic point) on every bond of every ring word over the 13 letters, legal or illegal, for L = 2, 3, 4, 5, 6. That is 148, 2886, 50024, 812890 and 12681084 directed transitions. F never changes. The same holds for the uncoloured control on 12 letters (up to 8709120 transitions at L = 6). C88's own `check_sector.py` covers only L = 4 (cosmetic). Since a move changes two adjacent sites and F depends only on squeezed neighbours, L = 4 already contains every local context, and L = 5, 6 add vacancy intervals between them.

[proved] The rule-table facts behind C88's general argument hold. No rule has (L_1, marker) as input or output. An L_1 output always has right partner in {R_0, D_A, D_B, X_A, X_B}, as does every L_1 input. No rule turns a right letter into L_1. A right letter changes marker type only with left input L_0, D_s or X_s. No rule turns a marker into a nonmarker. So [H, F] = 0 at every K, J, t, lambda, as PROOF Section 2 says.

### 3. Zero-sector equality with C85 (`check_bond13.json`)

[proved, exhaustive] On all 11^L full-occupancy ring words for L = 2 to 5 (up to 161051 words) at K = 37 and lambda = 0.17, both the C88 ring Hamiltonian and the C85 material ring Hamiltonian conserve F, their restrictions to F = 0 agree with difference 0.0, and H_C85 = H_C88 + K F holds identically (below 3e-14). This is PROOF eq. (3). The equality is of material operators at finite K and nonzero port, including all virtual paths, which is what the C85 material certificate needs. C88 correctly says that the two mobile lifts differ at finite K because their wrong-colour channels differ.

### 4. Colouring and N_matbad (`check_pin13.json`)

[proved, exhaustive] On every ring of 2 to 6 sites over the 13 letters, N_matbad <= N_bad holds for every word with N > 0, each occupancy pattern carries exactly one legal colouring, and the number of zero-penalty words with N > 0 (53, 376, 2685, 19267, 138488) equals sum_N binom(L, N) tr(A^N) with A the squeezed allowed matrix. So the zero-penalty words are in bijection with legal squeezed words, pinned ones included. The all-vacancy sector has two constant colourings, as stated.

[proved, by inspection] The two-colour minimality argument of PROOF Section 3 is correct within its class (one nonnegative pair penalty, every legal placement retained). A single vacancy letter accepting L_0 V, V L_0, R_0 V and V R_0 admits R_0 V L_0, whose squeezed junction is illegal.

## Part II. Grammar correction and suppression

### 5. The leading correction (`check_grammar13.py`, `check_grammar13.json`)

[numerical, exact to rounding] On every legal coloured ring of 4, 5 and 6 sites (2685, 19267, 138488 words, with 840, 7370, 63216 distinct illegal intermediates), the operator P B Q N_bad^-1 Q B P built from my literal bond at J = 0.7, t = 1.3, lambda = 0.37 and nonzero g, t_d, field and selector equals diag(t^2 C + J^2 A + 2 J^2 lambda^2 Z) to 4.4e-16. Every intermediate has one mismatch and the mixed hopping-material block vanishes. L = 6 goes beyond C88's own check (up to L = 5). The C coefficient is |I(a)| - 1, one on R_0 and E and zero elsewhere, as claimed.

[numerical] Energy test with nonzero port (ring of two p = 2 cells and one vacancy, J = 0.7, t = 1.3, lambda = 0.3, projected dimension 3528, full closure 45836). K(E_K - E_inf) is -1.74131, -1.75656, -1.76096, -1.76250, -1.76279 at K = 300, 10^3, 3 x 10^3, 10^4, 3 x 10^4, approaching the predicted -1.763171 with a 1/K remainder. The same geometry under C85's three-colour lift gave -3.107 in the previous report. The difference is the smaller colour-ambiguity counter, which confirms PROOF's remark that the tail-bit coefficient halves and the L/D/X coefficients disappear.

### 6. K^-2 suppression (`check_grammar13.json`)

[numerical] Ground-state illegal-material weight, my own geometries:

| Geometry | Dimension | JK = 10^2 | JK = 10^3 | JK = 10^4 |
| :-- | --: | --: | --: | --: |
| two p = 3 cells, bits (1,1),(0,1), one vacancy, lambda = 0, J = 1 | 3222 | 1.19e-5 | 1.21e-7 | 1.21e-9 |
| same, J = 0.1 | 3222 | 1.22e-7 | 1.21e-9 | 1.20e-11 |
| same, uncoloured | 1728 | 1.00e-3 | 9.41e-4 | 9.36e-4 |
| two p = 2 cells, bits (1),(1), two vacancies, lambda = 0.1, J = 1 | 103716 | 3.63e-7 | 4.27e-9 | 4.34e-11 |
| same, J = 0.1 | 103716 | 5.03e-9 | 5.02e-11 | 5.02e-13 |
| same, uncoloured | 40740 | 2.74e-5 | 2.68e-5 | 2.67e-5 |
| two p = 2 cells, bits (1),(1), one vacancy, lambda = 0.1, J = 1 | 16114 | 4.86e-7 | 5.79e-9 | 5.90e-11 |
| same, J = 0.1 | 16114 | 6.15e-9 | 6.21e-11 | 6.22e-13 |
| same, uncoloured | 10185 | 1.49e-5 | 1.43e-5 | 1.42e-5 |

Coloured weights fall as (JK)^-2 and uncoloured weights plateau. Residuals are below 7e-9. The pin weight is exactly zero in every closure, as conservation requires. The closures are smaller than in the C85 lift (for example 103716 against 192472), because there is one fewer colour.

## Part III. Frozen selection, certificate, deletion

### 7. Frozen F3 with pins (`check_frozen13.py`, `check_frozen13.json`)

[proved, checked line by line] PROOF Section 6 is correct. The all-string bound (19) is a compression statement for the pinless bond and so covers pinned cells of every length. Pinned cells of length 15 need the separate bound (22). Its quadratic-form step is right. With mu_sel < B_*, <H> - mu_sel >= -2B_* xy + (K - B_* - mu_sel) y^2 >= -B_*^2/(K - 2B_*). I checked this against the exact 2 x 2 minimum at three parameter sets. With B_* = 2^16 and K = 2^2560, epsilon_pin = 2^-2528 and a rational Taylor bound gives 2cos(pi/16) > 1.96153, so r_pin > 1.96 > 1. The cut-and-glue estimate (25) then gives a positive excess k[min(gamma, r_pin) - 5 B_port] for every defect configuration, with gamma > 2^-17 and B_port = sqrt(2) 2^-608.

[numerical] A crude norm estimate at schedule-like values (D = 1000) gives ||B_cell|| below 15 ||b_11|| of about 1.5 x 10^4, consistent with C85's B_* = 2^16.

[proved, numerical] The legal pinned words are frozen. For p = 2 to 6, every legal pinned cell L_a0 ... L_a(p-2) L_1 m inside a ring with a good neighbour cell (terminal L_0 E or vacuum R_0 ... R_0 E) keeps all of its data letters in every state of the invariant closure, at K = 50, with and without port, illegal states included. Only the marker changes, through the intercell port when the neighbour starts with R_0. So E_pin,15 = mu_sel exactly at zero field for the legal pinned word, and (22) is only needed for the illegal part of the F = 1 block. C88's caution about "every such word is frozen" concerns illegal words and is appropriate.

### 8. Reruns and certificate transfer (`rerun_C88.json`)

[numerical, exact reproduction] `run_all.py` (eight scripts) ran without error in about 17 s. Every output file is bit-identical to the committed one except timing fields in four JSON files. The ten certificate source files are byte-identical to C85's. The regenerated `root_certificate_exact_n13.json`, `finite_budget_exact.json`, `selector_exact.json`, `selector_certificate.json` and `finite_budget.json` are byte-identical to C85's files, and all five integer .npz files have identical arrays. `rational_bands_n13.json` differs from C85's only in its timing field. The material certificate is therefore C85's, which the previous referee spot-checked independently (26 roots, delta_2 = 0.0156745, delta_4 = 3.25949e-6). MEMBER.json carries the same values as C85's under renamed keys.

[proved, arithmetic] For the mobile K of PROOF (30), the third condition can be evaluated. log2(lambda^2 E_unit) is about -2074.6 and log10 tau_W is about 35.2, so tau is about 2^-2240 and K > 2 A_2/(J^2 tau) needs roughly K > 2^2244 at J = t = 1, below the material minimum 2^2560. Only the condition involving the unstated C_3 stays implicit, together with the infrared tolerance. C88 says this.

### 9. Letter deletion (`check_witness.py`, `check_witness.json`)

[proved] All thirteen creation witnesses in `sector_checks.json` check out with my bond. Each input is legal with F = 0 and lacks the letter, each output is legal with F = 0 and contains it, the move is local and the amplitude matches. C88 correctly limits the conclusion to this deletion mechanism and does not claim d >= 13.

## Part IV. S_pin: is the pin-free sector the mobile ground-state sector?

### 10. A replacement bound for pinned cells (referee result)

[proved at K = infinity, zero port, on finite rings; uses PROOF eq. (8)] On the legal subspace the unique-colour isometry gives H = H_hop + J sum_j 1(g_j = 0) b_j - mu N (PROOF eq. (8)). Here H_hop moves vacancies with amplitude -t independent of the letters, and the material letters are carried by hard-core particles that never pass one another. So in the particle-labelled picture the Hilbert space is (vacancy configuration) x (squeezed word), and the material term is a sum of bond operators each multiplied by a diagonal contact indicator c_j(x). Take the ground state psi of a sector in which some cell is pinned, L_a0 ... L_a(p-2) L_1 E. On that cell the material Hamiltonian is diagonal: field terms f a >= 0 and the marker contact mu_sel/2, with no move available (Section 7). Replace the cell's letters by the clock state phi = sum_j phi_j |L_0^j R_0^(p-j) E> on the all-zero vacuum code, keeping every vacancy amplitude of psi. The marker contacts are unchanged, the field energy does not increase, and the clock bonds add

    J <T_c> = -2J sum_(j=0)^(p-1) phi_j phi_(j+1) <c_j>_psi,

which is nonpositive for real nonnegative phi. Optimising phi gives

    E_(F=0, same lengths) <= E_pinned - J lambda_max(P(<c>)),               (R1)

where P(<c>) is the (p+1)-vertex path matrix with weights <c_j> on its edges. Each pinned cell can be replaced independently, so the bounds add, k pinned cells costing at least the sum of their terms. A pinned M_s cell maps to the legal F = 0 word L_0^p M_s with no greater diagonal energy (it is cheaper still by the Delta + D contacts). The pinned sector and the trial have the same N and the same cell lengths. Charge is not fixed, consistent with S_pin's "no artificial Q restriction". At nonzero port the trial loses the port matrix elements of the pinned marker in which its data letters would be R_0 (they are projected out), which changes the energy by at most a few B_port = sqrt(2) lambda per pinned cell. That is negligible at lambda = 2^-608. The ring-rotation (Z_N relabelling) bookkeeping of the particle picture is standard and does not affect the expectation values. I have not written it out line by line.

If the contact probabilities on the pinned cell are near the free-fermion half-density value alpha_contact = 1/2 - 2/pi^2, then lambda_max(P) is about 2 alpha_contact cos(pi/(p+2)). For b = 15 this is PROOF's reference cost (6), J alpha_contact 2cos(pi/16), about 0.58 J. So (6) is not only a reference number. Since the top eigenvalue of a nonnegative path matrix grows with its weights, (6) is a lower bound on the pin cost whenever every bond of the pinned cell has contact probability at least alpha_contact.

### 11. Small-ring ED (`check_spin.py`, `check_spin.json`, `check_spin_density.py`, `check_spin_density.json`)

[numerical] Complete legal sectors at fixed cell lengths for F = 0, 1, 2 (the union of closures of every legal squeezed word with that F, no charge or bit restriction) were diagonalised at J = 1, mu_sel = 0.5, D = 4, g = 0.07, t_d = 0.003.

| Ring | p | v | t/J | dim F=0 / F=1 | E_1 - E_0 | bound (R1) | mean contact | (E_1 - E_0)/(J 2cos(pi/(p+2))) | E_2 - E_0 - 2(E_1 - E_0) |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: |
| two cells | 2 | 0 | 1 | 676 / 312 | 1.4142 | 1.4142 | 1.000 | 1.000 | 0.0000 |
| two cells | 2 | 1 | 0.1 | 14196 / 6552 | 1.3244 | 1.0675 | 0.724 | 0.936 | 0.0642 |
| two cells | 2 | 1 | 1 | 14196 / 6552 | 1.1822 | 1.0836 | 0.766 | 0.836 | 0.0834 |
| two cells | 2 | 1 | 10 | 14196 / 6552 | 1.1773 | 1.1649 | 0.824 | 0.833 | 0.0139 |
| two cells | 2 | 2 | 0.1 | 56784 / 26208 | 1.2948 | 0.9529 | 0.647 | 0.916 | 0.0337 |
| two cells | 2 | 2 | 1 | 56784 / 26208 | 0.9984 | 0.8741 | 0.617 | 0.706 | 0.0887 |
| two cells | 2 | 2 | 10 | 56784 / 26208 | 0.9642 | 0.9494 | 0.671 | 0.682 | 0.0117 |
| two cells | 2 | 3 | 0.3 | 170352 / 78624 | 0.9444 | 0.6339 | 0.438 | 0.668 | 0.1603 |
| two cells | 2 | 3 | 1 | 170352 / 78624 | 0.8305 | 0.6973 | 0.492 | 0.587 | 0.0858 |
| two cells | 2 | 3 | 3 | 170352 / 78624 | 0.7950 | 0.7460 | 0.527 | 0.562 | 0.0321 |
| two cells | 3 | 0 | 1 | 5184 / 1728 | 1.6180 | 1.6180 | 1.000 | 1.000 | 0.0000 |
| two cells | 3 | 1 | 0.3 | 186624 / 62208 | 1.4901 | 1.3381 | 0.804 | 0.921 | 0.0774 |
| two cells | 3 | 1 | 1 | 186624 / 62208 | 1.4337 | 1.3294 | 0.820 | 0.886 | 0.0755 |
| two cells | 3 | 1 | 3 | 186624 / 62208 | 1.4183 | 1.3713 | 0.848 | 0.877 | 0.0428 |
| one cell | 2 | 3 | 0.3 | 1560 / 360 | 0.6886 | 0.4074 | 0.279 | 0.487 | n/a |
| one cell | 2 | 3 | 1 | 1560 / 360 | 0.4893 | 0.3942 | 0.278 | 0.346 | n/a |
| one cell | 2 | 3 | 3 | 1560 / 360 | 0.4244 | 0.3929 | 0.278 | 0.300 | n/a |
| one cell | 2 | 6 | 1 | 6552 / 1512 | 0.1897 | 0.1321 | 0.093 | 0.134 | n/a |
| one cell | 3 | 4 | 1 | 20160 / 3360 | 0.5921 | 0.4733 | 0.289 | 0.366 | n/a |
| one cell | 3 | 6 | 1 | 60480 / 10080 | 0.3636 | 0.2670 | 0.162 | 0.225 | n/a |
| one cell | 4 | 4 | 1 | 115920 / 15120 | 0.8001 | 0.6610 | 0.376 | 0.462 | n/a |

The bound column is J lambda_max(P(<c>)) from (R1), evaluated with the contacts measured in the F = 1 ground state. One-cell rings put a single cell of p data sites and its marker on a ring with v vacancies; p = 2, v = 3 and p = 3, v = 4 are half density. All rows are at lambda = 0. The lambda = 0.1 runs differ by less than 4e-4 J.

In every run (two-cell and one-cell rings, 0 to 6 vacancies, t/J from 0.1 to 10, lambda = 0 and 0.1) the pin cost is positive and at least the bound (R1). At lambda = 0.1 it needs the port allowance only at full occupancy, where the cost 1.41395 J sits 3e-4 J below the lambda = 0 bound 2cos(pi/4) J, well inside 3 B_port. At zero vacancies the cost equals the frozen clock energy 2cos(pi/(p+2)) J exactly. The cost falls with vacancy density in proportion to the contact probability, and in the one-cell rings the ratio of cost to J cbar 2cos(pi/(p+2)) stays between 1.0 and 3.3. It rises toward small t/J, because slow vacancies arrange themselves to favour the good cell's clock bonds. Two pins never cost less than twice one pin. E_2 - E_0 - 2(E_1 - E_0) vanishes at full occupancy and zero port (frozen pins are additive) and lies between 0.01 J and 0.16 J with vacancies, so pins repel and do not bind. At finite K (two p = 2 cells, one vacancy, t = J, lambda = 0.1, complete closures of dimension 295351 and 10766) the pin cost is 1.1969, 1.1867, 1.1825 J at K = 30, 100, 1000, approaching the K = infinity value 1.1820 J from above with K times the excess 0.45, 0.47, 0.48. So the order-1/K corrections favour the pin-free sector, as the counter argument in item 12 predicts.

### 12. Can pin defects proliferate?

[heuristic, physics standard] No, at least not in any regime that this construction uses. The pin cost is an order-J loss of clock kinetic energy, about J cbar 2cos(pi/16), where cbar is the contact probability. It does not compete with the vacancy kinetic energy t, because pinned and good cells have the same vacancy hopping and the bound (R1) already allows the vacancies to adjust to the pinned configuration. At half density with t of order J it is about 0.5 J to 0.6 J per pinned cell. That is many orders of magnitude above every small scale in the problem: the port band lambda^2 E_unit (about 2^-2075), the selector density margin gamma (about 2.8e-5 per site), the grammar corrections t^2/(JK) and J/K, and the 2^-35 normalized Bloch allowance. Entropy plays no role at T = 0, and the pin-pin interaction is repulsive or zero in every case computed, so there is no clustering route either. Finite-K corrections at order 1/K favour the good cell further, because the pinned cell contributes nothing to A or Z, and to C only through its marker, the same as a good cell. The good cell's R_0 head and R_0 tail bits add negative terms -t^2/(JK) and -J/K.

The one way the pin cost could become small is if the contact probability on pinned cells tends to zero, which would need vacancies to sit inside every material bond of that cell. At fixed positive material density this is a local expectation of order one in a gapless vacancy fluid, so it should stay bounded below uniformly in L. That uniform bound is what is still heuristic.

### 13. What remains of S_pin

[assessment] S_pin bundles three statements: (i) pinned length-15 cells lie above the pin-free sector, (ii) wrong-length cells and marker-free competitors lie above it, (iii) finite-defect positivity and the grand-canonical choice of N. Part (i) follows from (R1) at K = infinity with the single extra input of a uniform lower bound on contact probabilities. At finite K it follows at the physics standard. Parts (ii) and (iii) are the inherited dressed-defect comparison of C82/C83/C85 with no new pin content. So the pin adds no new kind of assumption to the campaign, and the d = 13 candidate stands on the same conditional footing as C85.

## Gaps, classified

| # | Gap | Class |
| :-- | :-- | :-- |
| R1 | S_pin's pin part is stated as an assumption. The replacement bound (R1) proves it at K = infinity on finite rings, up to a uniform lower bound on contact probabilities and a negligible port allowance. C88 should state (R1) and narrow S_pin to that contact bound plus the inherited wrong-length comparison. | fixable (C88 understates what it can prove) |
| R2 | PROOF eq. (6) presents J alpha_contact 2cos(pi/16) as a mean-field reference cost whose transfer to the mobile chain is "precisely the pin part of S_pin". By (R1) and monotonicity of the top eigenvalue of a nonnegative matrix, it is a lower bound on the K = infinity pin cost (up to the port allowance) once every bond of the pinned cell has contact probability at least alpha_contact in the pinned ground state. | cosmetic |
| R3 | `check_sector.py` checks conservation of F only on four-site rings. Exhaustion to L = 6 here, and the rule-table argument, settle it. | cosmetic |
| R4 | The ring-rotation (Z_N) bookkeeping in the particle-labelled form of (8), which (R1) uses, is not written out by C88 or here. Standard. | cosmetic |
| R5 | Mobile K in PROOF (30) depends on the unstated C_3 and the infrared tolerance. The tau condition is met by K = 2^2560 at J = t = 1 (needs about 2^2244). Inherited from C85 (previous report H4). | fixable, inherited |
| R6 | Section 9 rules out only this deletion mechanism, not a redesigned d = 12 lift. C88 says so. | cosmetic |
| R7 | Mobile infrared inputs: massive and finite-J vacancy matching, the wrong-length and marker-free dressed-defect comparison (S_pin parts (ii), (iii)), the t^2/(JK) payload-contact coupling in the stability problem, multichannel Luttinger stability and conformal-observable control. Stated as assumptions by C88. | serious, inherited (not new) |

## Verdict

C88 (d = 13, n = 13, conditional c = 14 > 12): the finite construction is correct. The bond matches its text exactly. The pin count is conserved by every literal move, including illegal configurations, all vacancy hops and the port. On F = 0 the material Hamiltonian is C85's at finite K, so C85's n = 13 certificate applies and in fact regenerates byte for byte. The two-colour lift gives a unique colouring, an exact named grammar correction with the halved colour counter, and K^-2 suppression with nonzero bits, two vacancies and a nonzero port. Frozen selection of F = 0 and length 15 is proved. The mobile assumption S_pin is weaker than C88 presents it. Pinned cells are excluded by a replacement bound up to a contact-probability lower bound, ED shows no sign of pin proliferation, and what remains is the inherited mobile defect comparison. This makes C88 the smallest-d mobile conditional candidate on the same footing as C85. It is not an established counterexample, because it still needs the inherited infrared inputs (R7).

## Reproducing

All scripts run from this folder with Python 3, NumPy and SciPy. `SCR` is a scratch directory with copies `SCR/work_C84`, `SCR/work_C85`, `SCR/work_C88` of the repository folders, and `SCR/orig_C88` an untouched copy, with `run_all.py` already run inside `SCR/work_C88`. The previous referee's `../referee_C84_C85/my11.py` must be present.

- `python3 compare_rerun.py SCR/orig_C88 SCR/work_C88 C88` writes `rerun_C88.json`.
- `python3 check_bond13.py SCR` takes about a second.
- `python3 check_pin13.py` takes about 6 s.
- `python3 check_frozen13.py` takes about 2 s.
- `python3 check_witness.py SCR/work_C88/sector_checks.json` is instant.
- `python3 check_grammar13.py` takes about 8 minutes, mostly the L = 6 enumeration and the two-vacancy closure.
- `python3 check_spin.py` and `python3 check_spin_density.py` take about 4 and 25 minutes. The second is dominated by the p = 3 and p = 4 one-cell rings at high vacancy number.

`my13.py` is the referee rebuild. Each check script writes the JSON file of the same name.

## Relation to the rest of the repository

This answers the "no independent referee" status of `codex/work_C88/`. It confirms that H1 of `../referee_C84_C85/README.md` is realised: the pinless 11-letter grammar lifts with two colours to d = 13 with the unchanged C85 material certificate. It adds the replacement bound (R1), which removes the pin from the list of new mobile assumptions. It does not change the conditional status of any route-A candidate.
