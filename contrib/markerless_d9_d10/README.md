# Markerless cells: a nine-letter frozen chain and a ten-letter mobile chain

## Claim

The route-A mechanism of C84/C88 (`contrib/small_d_note`) does not need marker letters, and its mobile lift does not need coloured vacancies.

1. **Frozen, d = 9.** There is an explicit one-site translation-invariant nearest-neighbour bond `h9` on C^9 ⊗ C^9, with only the coefficients depending on an integer n, whose cells store n fermionic flavours. Its logical chain at second order in the port has 2n simple Fermi points with non-resonant pair sums at the members listed below, seven of them (n = 9, 10, 11, 12, 14, 16, 20) with an integer-interval certificate. With the infrared inputs (I1)–(I3) of the small-d note, the gapless sector above a fixed cell pattern has c = n. The n = 9 member gives c = 9 > 8 = d − 1. This is the same "frozen-cell reading" as the eleven-letter chain C84, so it replaces d = 11 by d = 9.
2. **Mobile, d = 10.** Adding one uncoloured, fully transparent vacancy letter V gives a bond `h10` on C^10 ⊗ C^10 of exactly C88's form, J embed(frozen bond) − t (vacancy hops) − μ (occupancy). With (I1)–(I4) and the vacancy dictionary of C88, which applies in the regime J ≪ t, c = n + 1. The n = 9 member gives c = 10 > 9, and the n = 10 member gives c = 11. This replaces d = 13 (C88) by d = 10. The first version of this bond used genuine on-site energies and was wrong: the referee showed that it does not select the cell pattern once vacancies are present. The definition is now C88's.
3. **Growing c at fixed d.** Members with 2n simple Fermi points were found at every n tried: 2, 5, 8, 9, 10, 11, 12, 14, 16 and 20, certified at 9, 10, 11, 12, 14, 16 and 20. No obstruction to larger n was seen, but no all-n statement is proved.
4. **Negative result.** The natural seven-letter version (no blind express letters) fails for a structural reason: its hopping is a reflection times a positive matrix, the flavour chains dimerize, and at zero field at most n mod 2 bands cross the Fermi level. A field does not rescue it at the points tested.

Nothing here measures a central charge. The conclusions are conditional in exactly the sense of the small-d note.

## Status

- [proved] The sector of words with no adjacent (L, R) pair is exactly invariant under every term of `h9` and `h10`, at every value of every coupling. The argument is three lines (below) and is confirmed exhaustively on all 9^L ring words for L ≤ 7 and all 10^L ring words for L ≤ 6 (`ring_checks.py`).
- [proved, machine-checked] Hermiticity, the charge Q, the head number and the occupied number.
- [numerical, exact linear algebra in double precision] On two cells with n = 2, 3: zero leakage, vanishing first order, a bit-independent λ² term, no g¹ term, and a λ²g² term that is a free-fermion bilinear (residual 5e-14 to 6e-12 of the hopping at well-conditioned points). `two_cell_checks.py`.
- [numerical] The arbitrary-precision scalar formula used for large n agrees with the brute-force letter model to 1e-15 (hopping) and 1e-12 (on-cell block) at n = 2, 3. `two_cell_checks.py`.
- [numerical] The two-cell calculation pads each end with one frozen letter so that the end-of-cell terms see a neighbour. Three cells on a genuine ring (n = 1, 2) reproduce the Bloch assembly μ^L + μ^R on the diagonal and T between neighbours to 4e-12 of the hopping. `ring_three_cells.py`.
- [numerical] Full diagonalization of two cells (3600 states, all head sectors) at finite couplings reproduces the free-fermion levels with an error that falls as g². `ed_two_cell.py`.
- [proved, computer-assisted] Integer-interval certificates of the logical Bloch data at n = 9, 10, 11, 12, 14, 16 and 20, with the C84 schedule (a = 4, D = 1000). The scalar graph of this bond is C84's with the marker vertex deleted, so C84's certificate pipeline applies after changing one line in each of three scripts (n + 2 express vertices become n + 1). `certificate/`.
- [numerical] The certified roots agree with the independent arbitrary-precision formula of this folder to 1.8e-15 (n = 9) and 8.9e-16 (n = 10), the precision of the certificate's own floating-point root list. `crosscheck_certificate.py`.
- [numerical, 100 to 160 significant digits] A second, better-separated family at a moderate express gap (a = 2, D = 100): 2n roots on the unit circle to 1e-67 or better, n right-movers and n left-movers, identical at two working precisions. These are not certified. `member.py`.
- [proved given C84's two-cell remainder] Exact-fraction budgets for g and λ at the certified members. The finite-g tail bound is re-derived for this bond. The two-cell port remainder 2^50 λ⁴ is C84's, with its two bond-dependent hypotheses checked here.
- [assumed] The infrared inputs (I1)–(I4) of `contrib/small_d_note`, and for the mobile chain the vacancy dictionary of C81/C83/C85.

What would upgrade it: (i) a line-by-line check of C84's two-cell contour argument for this bond; (ii) a vacancy computation for the mobile chain beyond one cell; (iii) an independent human read.

Referee. `REFEREE.md` is the report of an independent Claude subagent with an adversarial brief, with its own code in `referee/`. It rebuilt both bonds from this text, redid the invariance enumeration, the second-order calculation (to n = 4, and on a three-cell ring), the Fermi-point counts by a different method (n = 5, 9, 10 and the certified n = 9 member), and reran the n = 9 certificate. It confirmed the frozen d = 9 chain and the certificate, and it did not confirm the mobile chain as first written. The corrections it asked for are made below; where its report quotes the earlier text, the earlier text was wrong. Its addendum rebuilt the corrected mobile bond, reran the enumeration and the one-cell diagonalization, and judges the objection resolved at the level raised: the mobile chain now has C88's conditional status.

## How it was established and checked

- Author: Claude (Fable 5.2) for Jonah Kudler-Flam, 12 November 2026. No Codex call was made. `certificate/` reuses C84's scripts, which were written by Codex (Astra), with the change stated below.
- Checked by: the scripts in this folder, all run locally (Python 3.11, numpy, scipy, mpmath). No human has checked any of it.
- Read before starting: `STATUS.md`, `contrib/small_d_note/note.tex`, `contrib/mobile_below_13/README.md`, `codex/work_C84` (RESULT, ALPHABET, structured.py, full_loop_integer.py), `codex/work_C86/RESULT.md`, `codex/work_C87/RESULT.md`.

## The construction

### Letters and grammar

    L_0, L_1, R_0, R_1   data: bit a to the left (L) or right (R) of the head
    W                    vacuum head (walker)
    A, B                 active tokens, charges -1 and +1
    X_A, X_B             blind tokens on the express rail, charges -1 and +1
    V                    vacancy (mobile chain only)

Charges: q(L_a) = q(R_a) = 2a, q(W) = q(V) = 0, q(A) = q(X_A) = −1, q(B) = q(X_B) = +1.

A cell is a word L_{a_0} … L_{a_{j−1}} H R_{a_j} … R_{a_{n−1}} with one head H ∈ {W, A, B, X_A, X_B}, n data sites and physical length b = n + 1. A cell boundary is the junction between the R-run of one cell and the L-run or head of the next. There is no marker letter. The only forbidden ordered pairs are (L_a, R_b), which would be a cell with no head. A cell has 5(n + 1)2^n legal states.

Compared with C84: the three marker letters E, M_A, M_B are gone, and the walker is its own letter instead of "R_0 in context". The count is 4 + 1 + 2 + 2 = 9 against 4 + 4 + 3 = 11.

### The frozen bond

With a ∈ {0, 1}, s ∈ {A, B}, z_A(a) = (−1)^a, z_B(a) = 1, φ_B = 0, N_L = {R_0, R_1, W, A, B, X_A, X_B} and N_R = {L_0, L_1, W, A, B, X_A, X_B},

    h9 = w + (1/2)(o ⊗ 1 + 1 ⊗ o) + K Π_bad

    w = Σ_a [ −|L_a W><W R_a| − Σ_s z_s(a) ( |L_a s><s R_a| + e^{iφ_s} |L_a X_s><X_s R_a| ) ]      head motion
        − Σ_s Σ_{y ∈ N_L} |y X_s><y s|                                                             rail switch, far-left end
        − t_d Σ_s e^{iψ_s} Σ_{y ∈ N_R} |X_s y><s y|                                                closing link, far-right end
        + g |B R_0><A R_1|                                                                         conversion
        + λ ( e^{−iπ/4} |X_A B><W W| + e^{+iπ/4} |X_B A><W W| )  + h.c.                            port

    o = ν Σ_heads P_H + Δ Σ_s (P_s + P_{X_s}) + D Σ_s P_{X_s} + (f/2) Σ_ℓ q(ℓ) P_ℓ

Π_bad projects on the pairs (L_a, R_b), and "+ h.c." applies to all of w. At a generic point `w` has 82 nonzero off-diagonal entries. `bond9.py` is the definition: `bond(p)` returns (w, o, Π_bad) and `full_bond(p)` the assembled 81 × 81 matrix.

Every head moves by the one rule (L_a, H) ↔ (H, R_a). The port acts on two adjacent heads, which happens exactly when the left walker is at the right end of its cell and the right walker at the left end of its own. The left token is born blind on the express rail, crosses its cell without touching data, becomes active only at the far-left end, and returns the same way. The right token is born active at its left end. Both therefore see data site j at the same graph distance. This is C84's fold, with the entry vertex X_n in place of the marker.

### Exact invariance, and what it buys

[proved] Let F be the number of adjacent (L, R) pairs. Every term preserves F on every word.

- The head rule maps x L_a H y to x H R_a y. The new pairs are (x, H), (H, R_a) and (R_a, y), none of which has an L on the left and an R on the right. In the other direction the new pairs are (x, L_a), (L_a, H) and (H, y), again none of the form (L, R).
- The switch, the closing link and the port only replace heads by heads.
- The conversion changes bit values and a head, not orientations.

The cell pattern is conserved as well. Let B_i = [x_i is not an L]·[x_{i+1} is not an R], which is 1 exactly on a cell boundary. The head rule on sites (i, i + 1) keeps B_i = 0, leaves B_{i−1} = [x_{i−1} is not an L] unchanged because x_i goes from L_a to a head, and leaves B_{i+1} = [x_{i+2} is not an R] unchanged because x_{i+1} goes from a head to R_a. The other terms change no orientation and no head position. So every B_i commutes with `h9`, and for `h10` the cyclic boundary pattern of the squeezed word is conserved (the referee checked both exhaustively). In C84 the frozen pattern comes from on-site marker projectors. Here it comes from these two-site diagonal operators. In F = 0 every cell has exactly one head unless the whole ring is headless.

So the sector F = 0 is invariant at every K, including K = 0. Three consequences, none of which holds for C84:

1. There is no finite-K error and no need for K = 2^2048. K only raises the other sectors.
2. At zero port the 2^n code of each cell is exactly degenerate, and the ring spectrum is a sum of independent cell clocks.
3. The mobile lift needs no vacancy colours (next section).

### The mobile bond

    h10 = J embed(b9) − t Σ_{a material} ( |a V><V a| + h.c. ) − (μ/2)(n_occ ⊗ 1 + 1 ⊗ n_occ)
    b9  = w + (1/2)(o ⊗ 1 + 1 ⊗ o) + K Π_bad

The embedding vanishes outside the material-material block, so the material on-site energies are contact terms: a letter gets half of its on-site energy from each neighbouring bond whose other end is material. This is C88's bond h13 with b11 replaced by b9 and no colours. `bond10.py` is the definition. With vacancies removed, a physical word gives a material word, and F is defined on that squeezed cyclic word.

[proved] F is conserved by `h10` on every word. Vacancy hops do not change the squeezed word. Material terms act only on physically adjacent material pairs, and by the argument above they never create an (L, R) adjacency whatever the neighbours are. So the legal sector is (legal material word) × (any vacancy placement), which is the structure C88 obtains with two coloured vacancies and C82 showed fails for C75 with an uncoloured one. The difference is that C84's walker is identified by context, so a vacancy can hide a junction that its dynamics creates. Here no term creates a junction.

Why the on-site terms must be contact terms. For J ≪ t the vacancies form a fast fluid, and to first order in J every material bond operator is multiplied by the probability c that the two letters are in contact. With contact on-site terms the whole frozen Hamiltonian is multiplied by cJ, so the selector, the token energies a and D and the certified Bloch data are those of the frozen chain. With genuine on-site terms only the two-site part is multiplied by c, the cell energy becomes J(ν − 2c cos(π/(b + 1))), and this is positive once c < 0.92 at n = 9. Every head then costs energy and the headless ring lies below the cell pattern. This was the error in the first version, found by the referee.

[numerical, one cell] `mobile_selector.py` diagonalizes one cell plus N_v vacancies on a ring against the headless ring with the same N_v, at the frozen selector value of ν. The cell energy ε_cell must be negative. For n = 4 (frozen value −0.460):

| N_v | density | contact, t/J = 25 | contact, t/J = 1 | genuine, t/J = 25 | genuine, t/J = 1 |
| --: | --: | --: | --: | --: | --: |
| 0 | 0 | −0.460 | −0.460 | −0.460 | −0.460 |
| 1 | 0.17 | −0.369 | −0.396 | −0.119 | −0.204 |
| 2 | 0.29 | −0.288 | −0.324 | +0.189 | +0.072 |
| 3 | 0.38 | −0.222 | −0.260 | +0.436 | +0.305 |
| 5 | 0.50 | −0.135 | −0.166 | +0.764 | +0.636 |

With contact terms the ratio to the frozen value at t/J = 25 is 0.80, 0.63, 0.48 and 0.29, which is the contact probability of the vacancy fluid (0.291 at half occupancy on this ring, against 1/2 − 2/π² = 0.297 in the thermodynamic limit). At n = 9 (frozen value −0.1487) the contact convention gives −0.134, −0.120, −0.106 and −0.094 for N_v = 1 to 4 at t/J = 25, and the genuine convention gives +0.035, +0.221, +0.398 and +0.559 (`data/mobile_selector_n9.json`).

What this does and does not show. It shows the sign and the uniform rescaling for one cell with the port off. It does not test wrong-length competitors, several cells, tokens or the port. Away from J ≪ t the rescaling is no longer uniform (at t/J = 1 the ratio is 0.86 where the contact probability gives 0.80), so the selector midpoint shifts and the schedule needs J/t small, as C88's red team found for d = 13. How small is not quantified here: the bound has to be set against the selector window (ν within about ±0.020 at n = 9) and the certificate's 2^−35 allowance. At zero port and zero field the cell Hamiltonian is bit-independent for every J/t, so there is no analogue of C88's t²/(JK) bit field. With a field f the field is itself a contact term, and the referee's addendum measures a position-dependent bit energy of about 0.03 (J/t) f, which has to be included in that budget.

### Selector

At zero port a cell with n data sites has energy ν − 2cos(π/(n + 2)), so with b = n + 1 the density is e_b = (ν − 2 + a_b)/b with a_b = 2[1 − cos(π/(b + 1))]. This is the formula of the small-d note with μ_sel = ν, and the same midpoint ν = 2 − (δ_b + δ_{b−1})/2 selects b. For `h9`, sectors with headless stretches or no heads have zero or positive energy on the dead sites and lose to e_b < 0. For `h10` the same comparison is made with every material bond term multiplied by the contact factor, which is why the on-site terms have to be contact terms (previous section). `selector.py` tabulates ν, e_b, the density gap and the clock gap for every certified n. At cell lengths 13 and 15 it returns C84's and C88's selector values digit for digit (ν = 1.8561806413628503 with density gap 5.47433e-5, and ν = 1.889188076555019 with 2.8116274e-5), as it must, since the zero-port vacuum is the same clock.

| n | cell length | ν | e_b | density gap | clock gap |
| --: | --: | --: | --: | --: | --: |
| 9 | 10 | 1.7702928472898476 | −0.0148693 | 1.82e-4 | 0.2365 |
| 10 | 11 | 1.8059143706510356 | −0.0114488 | 1.18e-4 | 0.1998 |
| 11 | 12 | 1.8338586041390402 | −0.0090021 | 7.92e-5 | 0.1710 |
| 12 | 13 | 1.8561806413628503 | −0.0072058 | 5.47e-5 | 0.1479 |
| 14 | 15 | 1.8891880765550188 | −0.0048255 | 2.81e-5 | 0.1138 |
| 16 | 17 | 1.9120160428044660 | −0.0033882 | 1.56e-5 | 0.0902 |
| 20 | 21 | 1.9406853677266191 | −0.0018551 | 5.73e-6 | 0.0607 |

### Second order

After gauging the A carrier by the parity of crossed bits, each species sees a scalar graph with express vertices X_0 … X_n (diagonal a + D) and data vertices D_0 … D_n (diagonal a), where a = Δ + 2cos(π/(n + 2)) is the energy above the cell vacuum. The edges are H[X_{j+1}, X_j] = −e^{iφ_s} along the express rail (the amplitude for X_j → X_{j+1}), −1 along the data rail, −1 between X_0 and D_0, and H[X_n, D_n] = −t_d e^{iψ_s} (the amplitude for D_n → X_n). The left source is X_n and the right source is D_0. Conversions act at D_0 … D_{n−1}. This is C84's graph with its marker vertex M deleted.

With R_st = (H_s ⊗ 1 + 1 ⊗ H_t)^{−1} and Π_j = |D_j><D_j|, in units λ² α⁴ g² with α⁴ = (2/(n + 2))² sin⁴(π/(n + 2)),

    T_jj'   = i [ <X_n D_0| R_AB (Π_j ⊗ 1) R_BB (1 ⊗ Π_j') R_BA |X_n D_0>
                + <X_n D_0| R_AB (1 ⊗ Π_j') R_AA (Π_j ⊗ 1) R_BA |X_n D_0> ]
    μ^L_jj' = − [ <X_n D_0| R_AB (Π_j ⊗ 1) R_BB (Π_j' ⊗ 1) R_AB |X_n D_0>
                − <X_n D_0| R_BA (Π_j' ⊗ 1) R_AA (Π_j ⊗ 1) R_BA |X_n D_0> ]

    μ^R_jj' = − [ <X_n D_0| R_BA (1 ⊗ Π_j) R_BB (1 ⊗ Π_j') R_BA |X_n D_0>
                − <X_n D_0| R_AB (1 ⊗ Π_j') R_AA (1 ⊗ Π_j) R_AB |X_n D_0> ]

(μ^R is μ^L moved to the right factor with the species labels exchanged.) The logical Bloch matrix is h(k) = μ^L + μ^R + f + T e^{ik} + T† e^{−ik}. The signs were fixed by comparison with the brute-force letter model, not derived by hand. `scalar9.py`.

As in C84: at t_d = 0 the two terms of μ cancel (measured |μ|/|T| = 4e-11 in double precision, 3e-45 at 50 digits), and T is a phase times a matrix with positive Hermitian part. The closing link makes the bands non-reciprocal and separates the Fermi-sea centres.

The loop flux seen by species s is ψ_s − nφ_s. With φ_A = π/(2n), ψ_A = π and ψ_B = π/2 both species see flux π/2, the two scalar graphs are gauge-equivalent, the hopping is real up to an overall phase and μ is purely imaginary (measured: imaginary part of the gauged hopping below 1e-205, real part of μ below 1e-177). This is the matched choice that the certificate assumes. The uncertified family below was computed with φ_A = π/(2(n + 1)), where the fluxes differ slightly; at n = 9 the matched choice gives the same count (18 of 18, pair-sum separation 1.0e-4 at η = 1e-6).

## Certified members

[proved, computer-assisted] Schedule: a = 4 (so Δ = 4 − 2cos(π/(n + 2))), D = 1000, φ_A = π/(2n), φ_B = 0, ψ_A = π, ψ_B = π/2, t_d = β 2^(−e) with β = 1/d_{n+1}, d_0 = 1, d_1 = 1004, d_{k+1} = 1004 d_k − d_{k−1}, and a uniform field λ² α⁴ β² g² 2^(−F). Apart from the cell geometry this is C84's schedule.

| n | cell length | e | F | Fermi points | root separation > | pair-sum separation > | Hermitian hopping margin > | finite-g power | predicted c |
| --: | --: | --: | --: | --: | --: | --: | --: | --: | :-- |
| 9 | 10 | 35 | 83 | 18 | 0.0168079 | 8.8309e-5 (171 pairs) | 0.5769 | 2^−160 | 9 frozen (d = 9), 10 mobile (d = 10) |
| 10 | 11 | 42 | 92 | 20 | 0.0154175 | 8.7978e-5 (210 pairs) | 0.5733 | 2^−175 | 11 mobile |
| 11 | 12 | 48 | 101 | 22 | 0.0151683 | 4.8781e-5 (253 pairs) | 0.5695 | 2^−190 | 12 mobile |
| 12 | 13 | 56 | 111 | 24 | 0.0147956 | 3.0955e-6 (300 pairs) | 0.5653 | 2^−204 | 13 mobile |

For comparison C84's n = 11 member has pair-sum separation above 9.1e-7. Each certificate consists of a positivity certificate for the tree Gram matrix (integer Bareiss determinants), 512-bit integer enclosures of the loop matrices, rational endpoint inertia tests that bracket every root, and exact chord bounds on all pair sums. It allows a further normalized perturbation of 2^−35. The "finite-g power" column is the scalar threshold printed by the root routine. The couplings actually budgeted are in the next table.

[proved given C84's two-cell remainder] `certificate/finite_budget9.py` repeats C84's exact-fraction budget, PROOF.md (17), with the finite-K term deleted:

| n | g | λ | finite-g error | two-cell port remainder | allowance |
| --: | --: | --: | --: | --: | --: |
| 9 | 2^−163 | 2^−359 | 6.8e-13 | 1.4e-12 | 2^−35/8 = 3.6e-12 |
| 10 | 2^−177 | 2^−388 | 1.7e-12 | 1.7e-12 | 3.6e-12 |
| 11 | 2^−192 | 2^−419 | 9.6e-13 | 4.8e-13 | 3.6e-12 |
| 12 | 2^−207 | 2^−448 | 4.8e-13 | 9.7e-13 | 3.6e-12 |

For comparison C84 at n = 11 uses g = 2^−216, λ = 2^−480 and K = 2^2048. The two estimates are C84's, PROOF.md (10) and the paragraph before (17).

- Finite g. With γ the gap of the free excited pair, the Neumann series in the conversion gives ‖H_2(g) − H_2(0) − g² H_22‖ ≤ ‖VP‖² γ^(−1) (2g/γ)⁴ / [1 − (2g/γ)²], because odd orders change the charge of the data and vanish on the code. Here ‖VP‖ = √2 α² (two orthogonal port components), the pair conversion has norm at most 2g (at most one transition per basis state and cell), and γ ≥ 4 − 2t_d > 3 by Gershgorin on the path graph with diagonal 4. This gives C84's α⁴ g⁴/4 for g ≤ 1/2, so the bound is re-derived, not only inherited.
- Two-cell port remainder 2^50 λ⁴. C84's contour argument needs the vacuum band of two cells to be isolated by more than 1/32 at zero port and the port to have norm at most √2. Both hold: the nearest levels are clock excitations at 0.236, 0.200, 0.171, 0.148, 0.114, 0.090 and 0.061 for n = 9, 10, 11, 12, 14, 16 and 20, and the port maps |W W> to a vector of norm √2. I did not re-check C84's statement that the intertwiner has norm below two, so this constant is inherited.

Also used: α⁴ ≥ 64/(n + 2)^6.

Larger links fail, and the failures are kept in the record: e = 38 at n = 10, e = 42 at n = 11 and e = 50 at n = 12 each lose a root bracket. The same pipeline also certifies n = 14 (e = 68, pair-sum separation above 6.4e-7), n = 16 (e = 84, above 3.4e-6) and n = 20 (e = 112, 40 simple Fermi points, above 1.6e-6, so c = 21 for the mobile chain at d = 10). At n = 20, e = 122 fails because two pair-sum intervals overlap at the bracket width used. See `certificate/larger_n_attempts.log`.

What was changed in C84's code: `m=n+2;S=2*n+3` became `m=n+1;S=2*n+2` in `certify_tree.py`, `rational_bands.py` and `full_loop_integer.py` (plus a two-line header comment added afterwards), and one metadata line in `certify_roots.py`. `loop_scan.py` and `structured.py` are unchanged. `finite_budget.py` is replaced by `finite_budget9.py`, which drops the finite-K term (there is none) and uses this cell's α.

[numerical] The dictionary between the certificate's scalar graph and the nine-letter bond does not rest on C84's derivation. It was checked through the chain: brute-force letter model on two cells = scalar formula of `scalar9.py` (n = 2, 3, to 1e-15), and scalar formula = certificate roots (n = 9, 10, to 1.8e-15), at the certified parameters and 230 to 260 digits. In units λ² α⁴ g² the certificate's hopping is T/β² and its field is f/β².

## Uncertified members at a moderate express gap

This is a separate family and should not be read together with the certified one. Its lowest token level is only 0.07 to 0.08 above the cell vacuum, below the clock gap, so the finite-g bound above (which needs a pair gap of 3) does not cover it, and it uses a slightly unmatched express phase.

All at a = 2 (so Δ = 2 − 2cos(π/(n + 2))), D = 100, φ_A = π/(2(n + 1)), ψ_A = π, ψ_B = π/2, t_d = ηβ with β = 1/d_{n+1}, d_0 = 1, d_1 = a + D, d_{k+1} = (a + D)d_k − d_{k−1}. The uniform field is f = f_γ γ_min, with γ_min the smallest eigenvalue of the Hermitian part of the hopping.

| quantity | n = 9 (d = 9 frozen, d = 10 mobile) | n = 10 (d = 10 mobile) |
| :-- | :-- | :-- |
| cell length b | 10 | 11 |
| η, f_γ | 1e-5, 0.3 | 1e-6, 0.1 |
| Δ | 0.0810140527710052 | 0.0681483474218634 |
| log10 β, log10 t_d | −20.0856, −25.0856 | −22.0942, −28.0942 |
| simple Fermi points on the unit circle | 18 (9 right, 9 left) | 20 (10 right, 10 left) |
| largest distance of a root from the unit circle | 8e-68 | 8e-74 |
| root separation | 0.0180182283 | 0.0181039521 |
| pair-sum separation | 2.23995624e-4 over 171 pairs | 1.47281654e-4 over 210 pairs |
| smallest whitened slope at a root | 1.380 | 0.444 |
| first-order whitened error allowance | 7.7e-5 | 1.6e-5 |
| ν (selector midpoint) | 1.7702928472898476 | 1.8059143706510356 |
| selected density e_b | −0.0148693100 | −0.0114488438 |
| selector density gap | 1.82e-4 | 1.18e-4 |
| clock gap | 0.2365 | 0.1998 |
| predicted c | 9 frozen, 10 mobile | 11 mobile |

The roots agree to all 25 printed digits between 100 and 150 digits of working precision (n = 9) and between 110 and 160 (n = 10). The "error allowance" is d4 × (smallest slope)/4, the size of a perturbation in the whitened norm that could first make two pair sums meet. It is a first-order estimate, not a bound. In the whitened metric the closing-link block μ has norm 1.2 (n = 9) and 1.9 (n = 10), so these two members are not small deformations of the t_d = 0 point. A more conservative n = 9 member, η = 1e-6, has whitened ‖μ‖ = 0.12, the same 18 simple Fermi points, root separation 0.0183 and pair-sum separation 5.8e-5 (`data/member_n9_eta1e-6_f0.3_dps100.json`). The smallest eigenvalue of the Hermitian part of the hopping is 10^−51.6 (n = 9) and 10^−57.0 (n = 10) in units λ² α⁴ g², and the lowest token level lies 0.0809 and 0.0680 above the cell vacuum.

Other members with 2n simple Fermi points (`data/*.log`):

| n | D | η | f_γ | root separation | pair-sum separation |
| --: | --: | --: | --: | --: | --: |
| 5 (at a = 2.02) | 20 | 1e-2 | 0.3 | 0.151 | 8.1e-4 |
| 8 | 100 | 1e-4 | 0.3 | 0.0225 | 4.8e-5 |
| 11 | 100 | 1e-8 | 0 | 0.0156 | 6.0e-5 |
| 12 | 100 | 1e-9 | 0.1 | 0.0146 | 1.0e-4 |
| 14 | 100 | 1e-13 | 0.1 | 0.0121 | 2.9e-6 |

Two observations about the parameter choice.

- The link strength has to fall by a factor of about 10 to 30 per added flavour. At n = 9, η = 1e-4 loses two Fermi points, and at n = 10 so does η = 3e-6. The eigenvalues of the Hermitian part of the hopping span 11 to 18 orders of magnitude, so the velocities are extremely disparate, as in C84.
- A moderate express gap is better than C84's D = 1000. The root separation scales like 1/D (0.022 at D = 100 against 0.0022 at D = 1000 for n = 8), but D = 20 is too small from n = 8 on, where the Hermitian part stops being positive.

## Why seven letters fail, and why nine

[numerical, with a structural explanation] `seven_letter_negative/` is the bond with letters L_0, L_1, R_0, R_1, W, A, B and port (W, W) ↔ (A, B), (B, A). Its legal sector is invariant, the first order vanishes and the second order is exactly Gaussian. But the left token now enters from the right end of its cell and the right token from the left end of its own, so T = J S with J the flavour reflection and S a positive Gram matrix. Flavour j of one cell pairs with flavour n − 1 − j of the next, with hopping S_jj on one side and S_{n−1−j,n−1−j} on the other. Each pair is a dimerized chain, gapped at zero energy unless j is the middle flavour. Measured at zero field: 0 Fermi points of 4 at n = 2 and 2 of 6 at n = 3. A second conversion channel with a relative flux makes T complex but does not change the count. The referee scanned a uniform field with the same T and μ: n = 2 then has 2 Fermi points for fields outside the dimerization gap and n = 3 stays at 2 of 6, so the design still falls short of 2n, but "n mod 2" is a zero-field statement tested at one parameter point and n ≤ 3. It is not a no-go theorem.

[heuristic, restricted] Inside the pair-port mechanism I do not see how to go below nine, for these reasons.

1. Four data letters. A nearest-neighbour rule can force exactly one head per cell only if every data site carries the side of the head it is on. Energetics alone does not do it: free walkers between unoriented walls phase-separate.
2. One vacuum head, distinct from the data. Without a marker, a contextual walker (R_0 after an L) makes the boundary invisible.
3. Both tokens must reach flavour j at the same distance, which needs a carrier that crosses the cell without converting. Hermiticity means a letter produced by a conversion can convert back, so blindness needs its own letters. Encoding the rail in the orientation of the data instead needs a second head rule, and any second rule lets a token cross an R|L boundary into the neighbouring vacuum cell.
4. Blind letters of both charges. With one blind species the on-cell bounce of the directly created active token is not cancelled by its A ↔ B partner, and it is larger than the hopping by 1/β².

This is not a lower bound for other mechanisms.

## Checks

| script | output | what it shows | tag |
| :-- | :-- | :-- | :-- |
| `bond9.py`, `bond10.py` | stdout | Hermiticity, Q, head number, occupied number | proved |
| `ring_checks.py 7 6` | `data/ring_checks*.json`, `data/ring_checks_contact.log` (corrected mobile bond) | F conserved by every literal transition on all 9^L words (L ≤ 7, 33.9 million transitions at L = 7) and all 10^L words (L ≤ 6, 6 million at L = 6): 0 violations | proved |
| `two_cell_checks.py` | `data/two_cell_checks.json` | leakage, orders g⁰, g¹, g², Gaussianity, scalar formula against brute force | numerical |
| `ring_three_cells.py` | `data/ring_three_cells.json` | three cells on a genuine ring (n = 1, 2): one-particle block equals the Bloch assembly of the padded two-cell blocks to 4e-12 of the hopping | numerical |
| `ed_two_cell.py` | `data/ed_two_cell.json` | finite-coupling ED against the free-fermion levels | numerical |
| `selector.py` | `data/selector.json` | selector midpoint, density gap and clock gap | proved (closed form, evaluated at 40 digits) |
| `mobile_selector.py n` | `data/mobile_selector_n*.json` | cell energy with vacancies for both on-site conventions: negative with contact terms, positive with genuine ones | numerical, one cell |
| `certificate/run9.py n e` | `certificate/*_n{9,10,11,12}.json` | integer-interval certificate of 2n simple Fermi points and the separations | proved, computer-assisted |
| `certificate/finite_budget9.py n` | `certificate/finite_budget_n*.json` | exact-fraction budget for g and λ | proved given C84's two-cell remainder |
| `crosscheck_certificate.py` | `data/crosscheck_certificate_n*.json` | certificate roots against the independent scalar formula | numerical |
| `member.py n a D η f_γ digits` | `data/member_*.json` | uncertified members: roots, separations, slopes, selector | numerical |
| `scan_mp.py` | `data/scans_n5_n8.log` (also holds n = 2), `data/scans_n9_to_n14.log` | parameter scans and the other members, including the failing points | numerical |
| `seven_letter_negative/seven_check.py` | `data/seven_check.json` | the seven-letter failure | numerical |

## What is not done

1. The certificate covers the scalar Bloch data at second order in the port. The two-cell port budget reuses a constant from C84 whose hypotheses were checked but whose derivation was not repeated in full. The a = 2, D = 100 members are floating-point only.
2. No many-cell statement, no dressed-defect comparison at nonzero port, and no infrared control. These are (I1)–(I3), unchanged.
3. For the mobile chain, the vacancy dictionary and (I4) are taken over from C88, whose bond has the same form. The only vacancy numerics for this bond are the one-cell runs of `mobile_selector.py`. Wrong-length competitors at finite vacancy density, several cells, and the port with vacancies are not computed, and no bound on J/t is given.
4. Three-cell processes at order λ⁴ are not computed (a blind token that has crossed its cell can annihilate with the next cell's token at the other boundary, as in C84). The non-Gaussian terms at order λ² g⁴ are bounded by the finite-g estimate but not computed here; the referee's brute-force value at a = 4, n = 2, 3 is 8e-6 α⁴.
5. On a ring of N cells there is a winding term of order λ^N that depends on the bit parity of the cells it crosses (referee, section 3). It is negligible for long rings, but two-cell rings are not a valid numerical test of this bond.
6. Not refereed by a human.

## Reproducing

From this folder, with Python 3, numpy, scipy and mpmath:

    python3 bond9.py; python3 bond10.py                 # seconds
    python3 ring_checks.py 7 6                          # about 80 seconds
    python3 two_cell_checks.py                          # about 5 minutes
    python3 ring_three_cells.py                         # about 3 minutes
    python3 ed_two_cell.py                              # about 2 minutes
    python3 mobile_selector.py 4                        # about 1 minute; n = 9 takes longer
    cd certificate; python3 run9.py 9 35; python3 run9.py 10 42; python3 run9.py 11 48; python3 run9.py 12 56; cd ..   # 7 to 40 seconds each
    python3 crosscheck_certificate.py 9 certificate/full_loop_bands_n9.json 230    # about 10 minutes
    python3 member.py 9 2.0 100 1e-5 0.3 100            # about 5 minutes
    python3 member.py 10 2.0 100 1e-6 0.1 110           # about 10 minutes
    python3 scan_mp.py 8 80 2.0 100 1e-4 0,0.3          # n, digits, a, D, η list, f_γ list

## Relation to the rest of the repository

This addresses item 2 of "The conjecture is false at large d" in `STATUS.md` (smallest counterexamples) and the outlook of `contrib/small_d_note` (an encoding outside the single-token product class of C86). If it survives review, the smallest candidates become d = 9 (frozen) and d = 10 (mobile) in place of 11 and 13. The frozen chain has been through one subagent referee pass, and the mobile bond was corrected after it and re-checked in the referee's addendum. Both come with the same infrared inputs, the same kind of certificate for the Bloch data, and no finite-K budget. The exact invariance may also be useful for the larger alphabets, since it removes the finite-K budget and the coloured vacancies.
