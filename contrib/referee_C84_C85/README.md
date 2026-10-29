# Referee report on C84 (11-letter frozen bond, n = 11) and C85 (three-colour mobile lift, d = 14, n = 13)

## Claims under review

- Astra C84 (`codex/work_C84/`) removes E_1 from C80's 12-letter bond, fixes the last displaced bit to zero, and adds a penalty K on L_1 followed by any marker (the "boundary pin"). The alphabet is L_0, L_1, R_0, R_1, D_A, D_B, X_A, X_B, E, M_A, M_B. A cell with p data sites carries n = p - 1 bits, the legal cell has (5p+3)2^(p-1) states, and the scalar carrier keeps 2p+1 vertices with an inactive last data endpoint. A certified frozen member at p = 12, n = 11, cell length 13 gives the conditional prediction c = 11 > 10 at d = 11.
- Astra C85 (`codex/work_C85/`) lifts the 11-letter grammar to a mobile chain. It claims that the grammar needs three vacancy colours (d = 14), that three are sufficient, minimal and give a unique legal colouring, that the complete leading grammar correction is -(t^2 C + J^2 A + 2 J^2 lambda^2 Z)/(JK), and it certifies an n = 13 member (p = 14, cell length 15), giving c = 14 > 13.

## Status

Both finite constructions are correct as stated. My bond rebuilds from the PROOF texts match `bond11.py` and `bond14.py` entry by entry, the legal-word bijection and the count (5p+3)2^(p-1) hold, the C84 cell graph is exactly C80's with its last conversion deleted and last bit fixed, the finite-K vacuum splitting is a pure one-body bit field and equals C85's closed formula (10) to rounding, first order vanishes, the g^2 port operator is Gaussian, and small ED agrees with the second-order formula at n = 1, 2, 3. C85's grammar correction (8) is exact on every legal ring of 4, 5 and 6 sites, its colouring is unique, and illegal weight falls as K^-2 with nonzero bits, two vacancies and a nonzero port. Both certificates rerun bit for bit and pass independent spot checks. I found no fatal gap.

The main new finding concerns the boundary pin. It is dynamically inert. No literal move or port transition of the C84 bond, legal or illegal, frozen or with vacancies, creates or destroys an (L_1, marker) adjacency, so the pin projector commutes with the pinless Hamiltonian and the pin only removes decoupled frozen sectors. C84 is unaffected (the pin is harmless there). For C85 it matters. The third colour exists only to encode the pin, and the pinless 11-letter grammar admits C83's two-colour automaton with E_1 deleted, which I checked passes the same uniqueness and K^-2 tests. That suggests a d = 13 lift with the same material certificate, so C85's conclusion that the 11-letter compression cannot go below d = 14 holds only within its restricted class.

## How it was established and checked

- Author: Claude (Opus 5.5) referee subagent for Jonah Kudler-Flam, 2026-10-29. No Codex call was made. All computations are local, NumPy/SciPy only, at resolved diagnostic parameters except where integer certificates are reproduced.
- `my11.py` rebuilds the C84 bond from C84 PROOF Section 1 alone (string letters, my own word enumerator, decoder and closure builder). `my14.py` builds the C85 lift from C85 PROOF Section 1 and eqs. (3)-(4) on top of it. The previous referee's C80 rebuild (`contrib/referee_C80_C82/my12.py`) serves as an independent parent. The Astra bond files were imported only as objects under test.
- All supplied scripts were rerun in a scratch copy of the repository layout outside the repository, and outputs compared with the committed files.
- Tags: [proved] means a finite argument I checked line by line or an exact identity, [numerical] means floating-point evidence at resolved parameters, [heuristic] means a physics argument not made rigorous, [assumed] means an input taken from elsewhere.

## Part I. C84, the 11-letter frozen bond

### 1. Bond, Hermiticity, charges (`check_bond11.py`, `check_bond11.json`)

[numerical, exact in floating point] My rebuild matches `bond11.bond` entry by entry, maximum difference 0.0, at four generic points (p = 2, 3, 5, 12, lambda from 0.1 to 1, K from 37 to 2^20, nonzero field, selector and generic flux, port phase -pi/4). It has 166 nonzero entries at a generic point, as C84 states. The 121 x 121 matrix is exactly Hermitian and commutes exactly with the on-site charge sum and with the marker projector on each site.

[numerical, exact] The compression identity C84 eq. (19), H_11 = J^dag H_12 J + K sum Pi_(L_1, marker), holds with error 0.0 against the independent C80 rebuild. This is the identity on which the C84 selector and F3 rest.

[numerical] On a bond whose left site is a marker, the zero-port operator is the identity on the marker type (deviation below 1.5e-14, rounding of the split onsite term). The pin sits on the bond (last data letter, own marker) and is intra-cell. So the cell factorisation used by F3 survives the pin.

[proved, `check_charges11.py`] The 24 additive-charge equations of the literal rules have rank 8, nullity 3, spanned by the identity, marker number and the charge (1), as C84 Section 8 says.

### 2. Legal words and the boundary pin (`check_cells11.py`, `check_cells11.json`)

[numerical, exhaustive] Brute force over all 8^p data words and 3 markers after a marker gives exactly (5p+3)2^(p-1) zero-penalty words for p = 1 to 7 (8, 26, 72, 184, 448, 1056, 2432). Each word has at most one excited letter. My decoder (first non-L letter is the head, R_0 there is the vacuum, later R letters carry bits, a terminal word has its last data letter L_0) is a bijection onto (head position, species, n bits) and inverts C84 eqs. (3)-(4). For p up to 6 the set equals C80's legal set with the last bit fixed to zero, relabelled E_0 -> E. Marker-free zero-penalty rings of 2 to 7 sites are only the all-L and all-R rings.

[numerical, exhaustive; proved] Applying every literal move to every legal word (p up to 6, every left marker) creates 18 types of bad pair, all of the predecessor type of the vacuum moves R_0 R_a -> L_a R_0 and R_0 E -> L_0 E, exactly as in C80. No move ever creates the pin pair (L_1, marker). The reason is structural. In every rule an output L letter is the left member of its pair, and its partner is R_0, D_s, X_s, E or M_s, with an L_1 output only in R_0 R_1 -> L_1 R_0. No rule has (L_1, marker) as input or an L letter as right input. So the pin creates no new zero-penalty illegal configuration and no new leakage channel. Its effect on the code is nil (Part III gives the stronger statement).

### 3. Finite-K vacuum-code splitting (`check_cells11.json`)

[numerical, exact to rounding] The order-1/K effective operator on the single-cell vacuum code, computed from the illegal intermediates, equals C85 eq. (10), -(1/K) sum_r (1 - a_r) sum_(j <= r) |v_j|^2, to 2e-15 for p = 1 to 9. Its Walsh components of weight 2 or more vanish (below 2e-16). So the splitting is a one-body bit field, as for C80, and the spread is K^-1 times 0, 0.25, 0.638, 1.083, 1.554, 2.037, 2.526, 3.019, 3.514 for p = 1 to 9, about (p - 1)/2. ED of the full zero-excitation cell space (dimension 4^p, p = 2 to 6, K = 10^3 to 10^5) agrees with these coefficients up to O(K^-2), K^2 times the error staying between 0.15 and 3.5 (the K = 10^5 entries for p >= 4 are at the rounding floor). At the certified K = 2^2048 the splitting is about 2^-2045, far below the physical band scale lambda^2 E_unit. C84's 2^80/K constant is correct in structure and very conservative.

### 4. Selector at zero port (`check_cells11.json`, `check_selector_numbers.py`)

[proved] C84 eq. (20) follows from the compression identity of item 1 and C80's all-string bound, with the C80 vacuum with last bit zero as upper trial state, as C84 says.

[numerical, complete Hilbert spaces] Single cells with k = 0 to 4 data sites, every data letter at every site and every marker (dimension up to 12288, all illegal and multi-token sectors), at K = 6, 10, 100, with Delta = 2.2, g = 0.1, t_d = 0.5, all satisfy E_v(k) - (k-1)_+/(K-4) <= E <= E_v(k). Example k = 4, K = 10: E = -1.8443 within [-2.2321, -1.7321].

[numerical] mu_sel = 1.85618064136285, density margin 5.474330e-5 (competitor b = 14), e_13 = -0.0072057833, vacuum gap 0.14791809, all as C84 states. log2 of gamma/(5 sqrt(2) lambda) is about 463.

### 5. First order, Gaussian structure, small ED (`check_second11.py`, `check_second11.json`)

[numerical, exact] For p = 2, 3, 4 (n = 1, 2, 3) the projected single-cell matrix of my 11-letter rebuild equals, entry for entry (difference 0.0), the C80 rebuild with its last conversion D_A E_1 -> D_B E_0 deleted, restricted to last bit zero. This is C84's claimed unitary equivalence of legal graphs (here literally an identity after relabelling). The constrained vacuum is an exact eigenvector (residual 2e-16).

[numerical, exact Taylor coefficients] G^dag V G = 0 exactly and V G lies entirely in the (1,1) excitation sector. Writing H_2(g) = -W^dag (H_0 - 2E_v)^-1 W with W = V G, I computed the Taylor coefficients in g exactly from the resolvent series. The g^0 term is a scalar (deviation below 4e-18), the g^1 term vanishes (below 4e-19), and the g^2 coefficient is Hermitian and lies in the span of 1 and c_i^dag c_j (Jordan-Wigner order, left-cell bits then right-cell bits) to relative residual 9e-16, 1.0e-15, 1.2e-15 of its norm at n = 1, 2, 3. Its traceless part is small compared with its scalar part at these couplings (ratio 7e-6, 8e-8, 2e-8), so the residual relative to the traceless part is 1e-10, 1e-8, 6e-8, which is the rounding floor of that subtraction, not a non-Gaussian term. The carrier argument in C84 Section 3 is the stronger statement.

[numerical] Projected two-cell ED of H_0 + lambda V (dimensions 436, 3392, 22336 for n = 1, 2, 3, split by charge sector) agrees with 2E_v + lambda^2 eig H_2 with error/lambda^4 equal to 6.77e-4, 6.76e-4, 6.75e-4 (n = 1), 2.08e-4, 2.07e-4, 2.07e-4 (n = 2) and 7.31e-5, 7.29e-5, 7.29e-5 (n = 3) at lambda = 0.2, 0.1, 0.05. The n = 3 run is beyond C84's checks.

[numerical] Finite-K ED of the invariant closure of the two-cell code (dimensions 803 and 14022 at n = 1, 2, illegal words included) agrees with the quasidegenerate formula C84 eq. (15) with error/lambda^4 of 7.3e-4 (n = 1) and 2.6e-4 (n = 2) at K = 100, constant in lambda. The zero-port code band width at n = 2 is 0.012741366 at K = 100 and 0.0012761712 at K = 1000, the same digits as C84 Section 7, from an independent build.

### 6. Reruns and the n = 11 certificate (`rerun_C84.json`, `spot_C84_n11.json`)

[numerical, exact reproduction] `run_certificate.py`, `checks.py`, `selector_checks.py`, `grammar_audit.py` and `alphabet_audit.py` ran without error (about 9, 7, 13, 1 and 0.3 s). All certificate files are bit-identical to the committed ones, including every .npz integer file, except for timing fields. `small_checks.json` and `selector_certificate.json` differ only by floating-point noise (below 1.5e-10 and 1.6e-14).

[numerical, independent] My floating-point evaluation of the carrier formula with 2p+1 = 25 vertices, sources only on D_0, ..., D_10, closing entry on the inactive D_11, D = 1000 and t_d = beta 2^-55 reproduces the 512-bit integer F and M matrices entry by entry up to the common factor beta^2/4, with relative spread 2.4e-15 and 2.0e-15. This confirms the inactive-endpoint structure. A Bloch scan of the whitened matrices on a 2 x 10^5-point grid finds 22 inertia changes, 0 up to 11 and back to 0, each within one grid step of a reported root. Recomputed margins are delta_2 = 0.0164205 and delta_4 = 9.2009e-7, consistent with the certified lower bounds 0.0164202883 and 9.11916e-7.

[proved, arithmetic checked] The budget logs match eq. (17). The penalty term is log10(176) + 29.885 + 80 log10 2 - 2048 log10 2 + 960 log10 2 + 213.208 = -58.10, as printed. The bound alpha^4 >= 64/14^6 > 2^-17 follows from sin x >= 2x/pi.

### 7. Alphabet accounting (ALPHABET.md, PROOF Section 8)

[not independently checked beyond item 1] The lower bounds a >= 7, 9, 11 are restricted to stated encoding classes, and C84 tags them that way. I checked only the charge rank. `alphabet_audit.py` reruns identically. Part III shows that one of C84's functions (the pin) is not needed dynamically, which does not change the letter count of C84 because the pin costs no letter.

## Part II. C85, the three-colour lift

### 8. Bond and charges (`check_lift14.py`, `check_lift14.json`)

[numerical, exact in floating point] My lift from C85 PROOF eq. (4) on my 11-letter rebuild matches `bond14.bond` entry by entry, maximum difference 0.0, at three generic points (p = 2, 3, 14, varied J, t, mu, lambda, K), for both the coloured 196 x 196 and the uncoloured 144 x 144 bonds. Both are exactly Hermitian and commute exactly with Q, occupied number N and marker number M. The (o, I) table reproduces all 121 ordered pairs of C84, and `bond14.allowed` agrees with mine on all 196.

### 9. Colouring, minimality, uniqueness (`check_lift14.json`)

[numerical, exhaustive] On all rings of 3, 4, 5 sites over the 14 letters, the zero-penalty words (313, 2127, 14497) are in bijection with the legal squeezed material words at every nonempty occupancy pattern, each with exactly one colouring, and N_material_bad <= N_bad holds without exception. [proved] C85's induction along the vacancy string gives this for every length.

[proved, by exhaustion over maximal rectangles] The fewest rectangles A x B inside the allowed relation that cover the required pairs is 3, whether one requires all 65 allowed pairs, the 56 pairs that occur in rings of cells with p >= 1, or those occurring in rings of cells of one fixed length. So three colours are necessary for this grammar with a pair penalty, also allowing nondeterministic colourings. C85's fooling set uses (M_A, M_A), which only occurs next to an empty (p = 0) cell. A fooling set inside fixed-length cells exists, for example (D_A, E), (L_0, M_A), (L_1, D_A), so the conclusion is unaffected (cosmetic). The allowed relation has exactly three distinct rows, {L_0, E, M_A, M_B}, {L_1}, {R, D, X}, so the deterministic three-colour table is unique up to relabelling.

### 10. The leading grammar correction (`check_grammar14.py`, `check_grammar14.json`)

[numerical, exact to rounding] On every legal coloured ring of 4, 5, 6 sites (2127, 14497, 98699 legal words, 1256, 10690, 87516 distinct illegal intermediates), the operator P B Q N_bad^-1 Q B P built from the literal bond at generic J = 0.7, t = 1.3, lambda = 0.37, nonzero g, t_d, field and selector, equals diag(t^2 C + J^2 A + 2 J^2 lambda^2 Z) to 8.9e-16. Every intermediate has exactly one mismatch, and the mixed hopping-material block is exactly zero. [proved, checked channel by channel] C85's derivation is correct. Only the vacuum move and the port narrow the incoming set of their first letter, the hopping intermediates have a vacancy as second site of the mismatch and the material ones a material site, and each illegal intermediate returns to P only by the inverse move.

[numerical] ED test on a ring of two p = 2 cells and one vacancy (projected dimension 3528, full closure 63189, J = 0.7, t = 1.3, lambda = 0.3): K (E_K - E_inf) is -3.0693, -3.0958, -3.1034, -3.1061 at K = 300, 10^3, 3 x 10^3, 10^4, approaching the predicted -3.10724 with a 1/K remainder. This includes the port term Z, which C85's own energy test (lambda = 0) does not probe.

### 11. K suppression with nonzero bits (`check_grammar14.json`)

[numerical] Ground-state illegal-material weight, my own geometries:

| Geometry | Dimension | JK = 10^2 | JK = 10^3 | JK = 10^4 |
| :-- | --: | --: | --: | --: |
| two p = 3 cells, bits (1,1),(0,1), one vacancy, lambda = 0, J = 1 | 4284 | 1.18e-5 | 1.21e-7 | 1.21e-9 |
| same, J = 0.1 | 4284 | 1.20e-7 | 1.20e-9 | 1.20e-11 |
| same, uncoloured | 1728 | 1.00e-3 | 9.41e-4 | 9.36e-4 |
| two p = 2 cells, bits (1),(1), two vacancies, lambda = 0.1, J = 1 | 192472 | 3.62e-7 | 4.27e-9 | 4.34e-11 |
| same, J = 0.1 | 192472 | 4.79e-9 | 4.99e-11 | 5.01e-13 |
| same, uncoloured | 40740 | 2.74e-5 | 2.68e-5 | 2.67e-5 |
| two p = 2 cells, bits (1),(1), one vacancy, lambda = 0.1, J = 1 | 21966 | 4.85e-7 | 5.79e-9 | 5.90e-11 |
| same, J = 0.1 | 21966 | 6.02e-9 | 6.20e-11 | 6.22e-13 |
| same, uncoloured | 10185 | 1.49e-5 | 1.43e-5 | 1.42e-5 |

Coloured weights fall as (JK)^-2, the uncoloured weights plateau. Residuals are below 1e-8. The weight on squeezed words that violate only the pin is exactly zero in every closure, which anticipates Part III.

### 12. F3 and F4 for C85

[proved, checked line by line] C85 Section 6 is C82's argument with the bookkeeping C80/C83 already needed, plus the pin. All three markers share an outgoing state, so the marker-to-first-letter penalty belongs to the right cell, and the pin is intra-cell (item 1). The lower bound (18) holds by compression. The almost-additivity and the k(gamma - 5 B_port) estimate are as in C82. With gamma > 2.8e-5 and lambda = 2^-608, the margin is about 2^590. `check_selector_numbers.py` reproduces mu_sel = 1.88918807655502, density margin 2.8116274e-5 (competitor b = 16), e_15 = -0.0048254990 and vacuum gap 0.11381150.

[proved at the physics standard] F4 (Section 7) is C82/C83's local Schrieffer-Wolff construction with kappa = JK, unchanged by a third colour. `check_sw.py` reruns identically.

### 13. Reruns and the n = 13 certificate (`rerun_C85.json`, `spot_C85_n13.json`)

[numerical, exact reproduction] `run_all.py` (all seven scripts) ran without error in about 19 s. Every output is bit-identical to the committed one except timing fields.

[numerical, independent] With 29 carrier vertices, sources on D_0, ..., D_12, inactive D_13 and t_d = beta 2^-66, my floating-point formula reproduces F and M up to beta^2/4 with relative spread 2.5e-15 and 3.1e-15. The Bloch scan finds 26 inertia changes, 0 to 13 and back, within one grid step of the reported roots. Recomputed delta_2 = 0.0156745 and delta_4 = 3.25949e-6, consistent with the certified 0.0156742281 and 3.2513114e-6.

[proved, arithmetic checked] Penalty term log10(208) + 35.229 + 84 log10 2 - 2560 log10 2 + 1216 log10 2 + 259.231 = -82.52, as printed, and alpha^4 >= 64/16^6 = 2^-18 is exact.

[proved, arithmetic] With tau_band of C85 eq. (28), log10 tau_band is about -674, so the third condition in (27) at J = t needs K above roughly 10^675, below 2^2560 (about 10^771). The fourth condition involves the unstated constant C_3. So the mobile member's K is defined only implicitly, which C85 says.

## Part III. The boundary pin is a superselection label (`check_pin.py`, `check_pin.json`)

[proved; numerical, exhaustive] For the C84 bond without the pin penalty, every move and port transition on every ring of 3, 4, 5 sites over the 11 letters (1584, 23232, 319440 transitions) conserves the number of adjacent (L_1, marker) pairs. With an uncoloured vacancy, every material move and every hop on every ring of 4 and 5 sites (40320 and 604800 transitions) conserves the number of squeezed (L_1, marker) pairs. The general argument is the one in item 2. An (L_1, marker) pair is never an input of any rule, no rule changes the right member of a pair into a marker next to an L_1 or turns an L into L_1 with a marker as right partner, and a marker changes only through rules whose left input is X_s, D_s or L_0, or through the port, which keeps it a marker. So Pi_pin commutes with the pinless Hamiltonian at every K and lambda. Consequently H_C84 = H_pinless + K Pi_pin has the same spectrum and eigenvectors on Pi_pin = 0, and the pin-violating sectors are frozen cells L...L_1 m with no clock kinetic energy, separated from the selected vacuum by about |E_v|, close to 2 at the certified length. The compression argument of item 4 gives the selector bound for the pinless bond directly, since H_pinless = J^dag H_12 J.

[proved, numerical] Without the pin, the allowed relation has two distinct rows, {L_0, L_1, E, M_A, M_B} and {R, D, X}, and is reproduced by C83's two-state automaton with E_1 deleted (o = 0 on L and markers, 1 on R, D, X; I = {0} on L, D, X, M_s, {0,1} on R_0 and E, {1} on R_1). On rings of 3, 4, 5 sites over these 13 letters the zero-penalty words (376, 2685, 19267) are in bijection with legal squeezed words, colourings are unique, and N_material_bad <= N_bad. In the same two-cell, nonzero-bit, one-vacancy geometry as item 11, the illegal weight is 1.19e-5, 1.21e-7, 1.21e-9 at JK = 10^2, 10^3, 10^4 (J = 1), and 1.22e-7, 1.21e-9, 1.20e-11 (J = 0.1). This is the same K^-2 law as the three-colour lift.

[heuristic, conditional] This suggests a d = 11 + 2 = 13 mobile lift. Its occupied material Hamiltonian on the physical sector Pi_pin = 0 is C85's, so C85's n = 13 certificate (or a new n = 12 one) would serve, and the conditional count would be c = 14 > 12 = d - 1. I have not rebuilt the whole C85 verification suite for it. The frozen pin-violating cells are gapped defects of O(J) energy and would need to be covered by the same (assumed) mobile dressed-defect comparison as every other wrong cell. The colour-dependent contact energy has coefficients c(a) = |I(a)| - 1 equal to 1 on R_0 and E and 0 elsewhere, still bit-dependent, so the t^2/(JK) field of C85 Section 4 recurs with smaller coefficients.

## Gaps, classified

| # | Gap | Class |
| :-- | :-- | :-- |
| H1 | C85 concludes that the 11-letter grammar needs three colours and that the compression does not reduce the coloured lift below d = 14. The minimality is correct for the pinned grammar, but the pin is an exact superselection label (Part III), and the pinless grammar admits a two-colour lift with unique colouring and K^-2 suppression. The headline conclusion holds only within C85's restricted class, and a d = 13 candidate appears available with the same material certificate. | serious (for C85's conclusion, not for its construction) |
| H2 | C84 states that the pin "enforces the fixed boundary bit". The fixed bit is enforced by deleting E_1 (the walker can only write 0). The pin only removes frozen, dynamically decoupled sectors. Harmless for every C84 result. | cosmetic |
| H3 | C85's fooling set uses (M_A, M_A), which needs an empty cell. A fooling set within fixed-length cells exists, so minimality stands. | cosmetic |
| H4 | The mobile K in C85 (27) involves an unstated constant C_3, so the mobile member is specified only implicitly. The "certified n = 13 member" is the material certificate. C85 says this, but MEMBER.json could state it more plainly. | fixable |
| H5 | C84's finite-K constants (2^36/K up to 2^80/K) are allowances, not line-by-line estimates. The splitting is exactly the one-body field of C85 eq. (10), about (p - 1)/(2K) per cell. Inherited from C80 (previous report G2). | cosmetic |
| H6 | Root routines record intermediate finite_g_power thresholds (-199, -228) different from the final g. Both folders say so. Inherited. | cosmetic |
| H7 | C84's alphabet lower bounds (a >= 7, 9, 11) are restricted to encoding classes. Correctly tagged; only the charge rank was checked here. | cosmetic |
| H8 | Mobile infrared inputs for C85 (finite-J matching, dressed-defect comparison with a gapless vacancy fluid, the t^2/(JK) payload-contact coupling inside the stability problem, multichannel Luttinger stability). Stated as assumptions by C85. | serious, inherited (not new) |
| H9 | Frozen-chain infrared inputs for C84 (bosonization of 22 Fermi points, volume-uniform matching, nonzero-port dressed-defect comparison). Unchanged from C75 and C80. | serious, inherited (not new) |

## Verdicts

- C84 (d = 11, n = 11, conditional c = 11 > 10): the finite construction is correct. The bond matches its text, the legal-word map is a bijection with the stated count, the cell graph equals C80's with the last conversion deleted, the vacuum splitting is the one-body field of C85 eq. (10), the g^2 port operator is Gaussian, ED agrees at n = 1, 2, 3, the selector transfers by compression, and the n = 11 certificate reproduces bit for bit and spot-checks. The pin creates no zero-penalty illegal configuration and does not interact with F3. It is in fact dynamically inert (H2, cosmetic). This is now the smallest-d frozen conditional candidate, on the same conditional footing as C75 and C80.
- C85 (d = 14, n = 13, conditional c = 14 > 13): the construction is correct. The bond matches its text, three colours are necessary and sufficient for the pinned grammar with a unique colouring, the named correction is exact including the port term, leakage falls as K^-2 with nonzero bits and nonzero port, F3 and F4 transfer, and the n = 13 certificate reproduces and spot-checks. Its interpretive conclusion that d = 14 is the floor for the 11-letter compression is not supported beyond its restricted class (H1). Dropping the pin appears to give a two-colour d = 13 lift, which the campaign should build and verify.

None of these is an established counterexample. Both still need the inherited infrared inputs (H8 for the mobile lift, H9 for both).

## Reproducing

All scripts run from this folder with Python 3, NumPy and SciPy. `SCR` is a scratch directory holding copies `SCR/codex/work_C84` and `SCR/codex/work_C85` of the repository folders after their own scripts have been run there, and untouched copies `SCR/codex/orig_C84`, `SCR/codex/orig_C85`. The previous referee's `../referee_C80_C82/my12.py` must be present.

- `python3 compare_rerun.py SCR/codex/orig_C84 SCR/codex/work_C84 C84` (and the same for C85) writes `rerun_<tag>.json`.
- `python3 check_bond11.py SCR/codex/work_C84` takes a few seconds.
- `python3 check_cells11.py` takes about 20 s.
- `python3 check_second11.py` takes about 15 minutes, mostly the n = 3 projected ED and the n = 2 finite-K closure.
- `python3 check_charges11.py` and `python3 check_selector_numbers.py` take a second.
- `python3 check_lift14.py SCR/codex/work_C85` takes a few seconds.
- `python3 check_grammar14.py` takes about 40 minutes, mostly the two-vacancy closure of dimension 192472.
- `python3 check_pin.py` takes about two minutes.
- `python3 spot_cert11.py SCR/codex/work_C84 11 55 C84_n11` and `... work_C85 13 66 C85_n13` take a few seconds each.

`my11.py` is the referee rebuild of the C84 bond, word enumerator, decoder and closure builder, and `my14.py` the rebuild of the C85 lift. Each check script writes the JSON file of the same name.

## Relation to the rest of the repository

This addresses the "no independent referee" notes in `codex/work_C84/RESULT.md` and `codex/work_C85/PROOF.md`. It confirms the C84 frozen d = 11 candidate and the C85 d = 14 lift as finite constructions, and adds H1, the observation that the pin is a superselection label, which reopens a d = 13 two-colour lift built on C84. It does not change the conditional status of any route-A candidate.
