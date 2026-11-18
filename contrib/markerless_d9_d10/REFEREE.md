# Referee report: markerless d = 9 (frozen) and d = 10 (mobile)

Referee: Claude (Fable 5.2) subagent, 12 November 2026. Adversarial brief. No human has checked this report.
Reviewed: `README.md` as of 22:24 (27229 bytes), the scripts it names, `certificate/`, `data/`.
All referee code and outputs are in `referee/`. Nothing outside `referee/` and this file was written.
I did not read or run `ring_checks.py`, `twocell9.py`, `two_cell_checks.py`, `ed_two_cell.py`, `ring_three_cells.py`.

## Summary of verdicts

| item | verdict |
| :-- | :-- |
| 1. Bonds from the README text | h9 confirmed (exact). h10 confirmed up to one README/code mismatch (J·K against K) |
| 2. Exact invariance of F, Q, heads, occupancy | confirmed, L ≤ 7 (d = 9) and L ≤ 6 (d = 10) |
| 3. Second-order logical Hamiltonian | confirmed for n = 2, 3, 4; padding is faithful (three-cell ring) |
| 4. Fermi-point count | confirmed at n = 5, 9, 10 and at the certified n = 9 member |
| 5. Physics claims, frozen d = 9 | no new hole found; same status as C84 with two caveats |
| 5. Physics claims, mobile d = 10 | **not confirmed as stated**: the selector fails in h10 once vacancies are present |
| 6. Certificate | (a), (b), (c), (d) confirmed; (e) plausible, README is honest about it |

## 1. Bonds rebuilt from the README text

`referee/ref_bonds.py` was written from the README formulas before opening `bond9.py` or `bond10.py`. `referee/compare_bonds.py` compares at 12 random points with every coupling generic (including K, t_d, the three phases, J, t, μ).

- h9: maximum entry difference 0.0. 82 nonzero off-diagonal entries, as stated.
- h10: differs in exactly the four entries ⟨L_a R_b|h|L_a R_b⟩. The README writes `h10 = J [ w + K Π_bad ]`; `bond10.py` has `emb = J*w9 + K*bad9`. With K replaced by K/J in my version the difference is 4e-16. Harmless for the physics, but the text and "the definition" disagree.

Places where the text needed a choice: `bond9.bond` returns `(w, o, bad)` and has no K (h9 has to be assembled by the reader); it also has undocumented parameters `sw` and `phiB`. "+ h.c." has to be read as applying to all of w. The on-site normalisation of h10 (half on each bond) is not written. None of these led to a mismatch.

Verdict: confirmed, with the J·K mismatch.

## 2. Exact invariance

The three-line argument is correct. I read it as: the head rule changes x L_a H y into x H R_a y, and none of the three affected ordered pairs can gain or lose an (L, R) adjacency for any x, y. For h10 the same holds on the squeezed word because material terms act only on physically adjacent material pairs.

`referee/ring_invariants.py` (own enumeration, vectorised) checks every nonzero off-diagonal bond entry on every ring word:

| bond | L | words | transitions | violations of F, Q, heads, occupied |
| :-- | --: | --: | --: | :-- |
| h9 | 2 to 7 | 4 782 969 at L = 7 | 33 894 126 at L = 7 | 0 |
| h10 | 2 to 6 | 1 000 000 at L = 6 | 6 000 000 at L = 6 | 0 |

A negative control (`ring_invariants_control.py`) confirms the checker is not vacuous: a second head rule, a C84-style contextual walker, and a charge-violating flip each produce hundreds of violations at L = 4.

One thing the README does not state, and should: the cell pattern itself is conserved. The bond operator B_i = [x_i is not an L]·[x_{i+1} is not an R] (a cell boundary) is conserved bond by bond by h9 on every word, and for h10 the cyclic boundary pattern of the squeezed word is conserved. I checked both exhaustively (0 violations, same runs). This is what makes "frozen cells" true here; in C84 it follows from on-site marker projectors, here from two-site diagonal operators. It also shows that in F = 0 every cell has exactly one head unless the whole ring is headless.

Verdict: confirmed.

## 3. Second-order logical Hamiltonian

Method (`referee/chain.py`, `second_order.py`): the many-site Hamiltonian is built literally from my 81 × 81 bond, with no assumption about cells. With H = H0 + gC + λV and P the code space, the exact coefficients are M0 = −PVRVP, M1 = +PVRCRVP, M2 = −PVRCRCRVP, with R = (H0 − E0)^−1 by sparse LU on the token sector. M2 is tested against const + Σ h_mm' c†_m c_m' with Jordan–Wigner order along the chain.

Padded two cells [R_0 | cell | cell | L_0], relative to the largest hopping entry:

| n | a, D, t_d | first order | bit dependence of M0 / M0 | M1 | Gaussian residual | T against `scalar9` | μ^L, μ^R against `scalar9` |
| --: | :-- | --: | --: | --: | --: | --: | --: |
| 2 | 2.6, 3, 0.07 | 0 | 2e-16 | 0 | 1e-11 | 1.4e-15 | 7e-13, 5e-12 |
| 2 | 2.2, 0.4, 0.3 | 0 | 2e-16 | 0 | 1e-13 | 2.6e-15 | 2e-14, 2e-14 |
| 3 | 2.6, 3, 0.07 | 0 | 4e-16 | 0 | 8e-12 | 2.0e-15 | 4e-12, 4e-12 |
| 3 | 2.2, 0.4, 0.3 | 0 | 4e-16 | 0 | 2e-13 | 1.5e-15 | 5e-14, 6e-14 |
| 4 | 2.6, 3, 0.07 | 0 | 6e-16 | 0 | 2e-11 | 3e-16 | 6e-12, 5e-12 |
| 4 | 2.2, 0.4, 0.3 | 0 | 5e-16 | 0 | 2e-13 | 6e-16 | 7e-14, 7e-14 |

n = 4 (108 600 states) goes one step beyond the author's brute-force check. With t_d = 0 the blocks μ^L and μ^R vanish separately (1e-14 to 2e-10 of T, the double-precision floor set by the large constant in M2).

I also implemented the README's scalar formula from the text (`readme_scalar.py`, dense Kronecker-sum solves, not a Sylvester eigen-solver). It agrees with the brute force at the same level. "μ^R likewise on the right factor" cannot be read literally: the species labels have to be exchanged (start from R_BA, pass through R_BB), otherwise the middle resolvent has the wrong species. The orientation of the complex edges is also not stated; it has to be taken from the bond.

Is the padding faithful? Yes. Three cells on a genuine ring (`item3_rings.py`, n = 1 and n = 2, 82 080 states): the one-particle block is block-circulant with μ^L + μ^R on the diagonal and T between neighbours to between 1e-14 and 7e-12 of T; M0 is three times the padded value; and the full 64-dimensional M2 is Gaussian to 2e-11 once the seam hop carries the sign (−1)^(N−1). Changing the pad letters to R_1 and L_1 changes nothing (difference exactly 0). A head is not a frozen letter and cannot be used as a pad.

Two cells on a ring are not faithful: a pair created at one boundary can be annihilated at the other, which gives a bit-dependent λ² g⁰ term (relative size 3e-5 to 3e-3) and non-Gaussian λ² g² terms. On a ring of N cells the analogous winding term is of order λ^N and depends on the bit parity of the cells it crosses. It is negligible for large N but it rules out two-cell rings as numerical tests.

Verdict: (a) to (e) confirmed.

## 4. Fermi points

`referee/item4_fermi.py`: d(k) = det h(k) is real; I count its sign changes on a uniform grid of 3000 to 5000 points and bisect. Since z^n d is a polynomial of degree 2n in z, 2n sign changes prove exactly 2n simple zeros on the circle. Chirality comes from the inertia of h(k) (LDL† pivot signs) at the midpoints between zeros. A winding-number count on |z| = 0.999 and 1.001 is an independent second count.

| member | sign changes | right, left | root separation | pair-sum separation | zeros inside 0.999, 1.001 |
| :-- | --: | :-- | --: | --: | :-- |
| n = 5, a = 2.02, D = 20, η = 1e-2, f_γ = 0.3 | 10 | 5, 5 | 0.151223 | 8.1186e-4 (55 pairs) | 0, 10 |
| n = 5, a = 2.00, same | 10 | 5, 5 | 0.141869 | 1.3378e-3 | 0, 10 |
| n = 9, a = 2, D = 100, η = 1e-5, f_γ = 0.3 | 18 | 9, 9 | 0.0180182283002 | 2.23995624279e-4 (171) | 0, 18 |
| n = 10, a = 2, D = 100, η = 1e-6, f_γ = 0.1 | 20 | 10, 10 | 0.0181039521295 | 1.47281654213e-4 (210) | 0, 20 |
| n = 9 certified (a = 4, D = 1000, e = 35, F = 83) | 18 | 9, 9 | 0.0168084783571 | 8.83303929e-5 (171) | 0, 18 |

All agree with the README and with `data/member_*.json` to the digits printed; n = 9 is identical at 100 and 130 digits. The README's n = 5 row is the a = 2.02 point (it says so now).

Independent source of T and μ: my own implementation of the README formula at n = 9 and 100 digits (622 s) agrees with `scalar9.blocks` to 1.6e-83 (T), 4e-84 (μ^L) and 3.5e-65 (μ^R) of the largest hopping entry; at n = 5 to 2e-54, 2e-55, 1e-47.

Two remarks. Equal numbers of right- and left-movers are automatic for 2n simple zeros of a periodic Hermitian family, so that line carries no information. The number of negative eigenvalues runs 0, 1, …, n, …, 1, 0 around the circle: all left-movers sit in one cluster and all right-movers in the other, about π apart, which is why the pair sums are nearly degenerate.

`scalar9.blocks` takes about 11 s at n = 9 and 100 digits, not minutes.

Verdict: confirmed.

## 5. Physics

### Frozen chain, d = 9

I looked for unintended processes and found none at second order in the port. Two adjacent heads always belong to different cells, so the port and its conjugate act only at a cell boundary with the left head at the right end and the right head at the left end. The switch and the link use N_L and N_R with the same amplitude for every allowed neighbour, so a cell's token dynamics does not depend on the neighbouring cell's state. A token cannot enter a neighbouring cell (boundary bonds are conserved). The three-cell ring computation above includes every one of these contexts literally.

Competitors that exist only without markers: the fully headless rings (all L or all R, energy 0, 2^L-fold degenerate) and, in F ≥ 1, headless cells L…LR…R. At zero port the energy of any configuration is a sum over cells of (ν − 2cos(π/(b'+1))) or 0, so everything loses to e_b < 0 for any K ≥ 0. I agree with the README here.

Differences from C84 that the README should flag:

1. The frozen pattern is protected by conserved bond operators, not on-site projectors. True, checked, but not stated.
2. The certified members use φ_A = π/(2n) and a = 4, D = 1000; the tabulated a = 2, D = 100 family uses φ_A = π/(2(n + 1)) and has a token gap of 0.08, below the clock gap. The inherited finite-g bound needs a pair gap of 3 and does not cover that family. The README says the family is uncertified, which is enough, but the two families should not be read as one.
3. Three-cell processes at order λ⁴ exist here as in C84 (a blind token that crossed its cell can annihilate with the neighbour's token at the other boundary). The README lists them as not computed.

Verdict: no new gap. Same conditional status as C84, with a weaker (inherited) port budget and no finite-K budget.

### Mobile chain, d = 10

The structural claim is right: no term creates an (L, R) junction, vacancies only block, and the legal space is (legal material word) × (any vacancy placement). I checked the squeezed boundary pattern exhaustively for L ≤ 5.

The energetic claim is not right as written. In C88 the on-site energies are contact terms, so a fast vacancy fluid multiplies the whole material Hamiltonian by the contact probability c and the selector is unchanged. In h10 the on-site energies are genuine and only the two-site terms need contact. To leading order the effective material Hamiltonian is J(c·w + o), not c·J(w + o). The cell energy becomes J(ν − 2c·cos(π/(b + 1))), which is positive as soon as c < ν/(2cos(π/(b + 1))) = 0.9225 for n = 9 (0.9348 for n = 10). At half occupancy c = 1/2 − 2/π² = 0.297 and the cell energy is about +1.2 J. Then every head costs energy and the fully headless ring, which is also in F = 0, lies below the patterned sector.

Literal ED of my h10 (`item5_mobile_selector.py`), one n = 9 cell on a ring plus N_v vacancies against the headless ring R_0^10 with the same N_v, at the README's ν = 1.7702928:

| N_v (density) | ε_cell, t/J = 1 | ε_cell, t/J = 5 | ε_cell, t/J = 0.2 |
| :-- | --: | --: | --: |
| 0 (0) | −0.1487 | −0.1487 | −0.1487 |
| 1 (0.09) | −0.047 | +0.009 | −0.098 |
| 2 (0.17) | +0.101 | +0.189 | −0.014 |
| 3 (0.23) | +0.256 | +0.364 | +0.081 |
| 6 (0.375) | +0.665 | +0.798 | +0.350 |

ε_cell = E_0(one cell) − E_0(headless); the README formula gives −0.1487 at every density. The same sign change occurs at n = 4. This is a single-cell finite-size computation, but the effect is of order one and matches the leading-order estimate.

Consequences. (i) At the tabulated ν the patterned sector is not the lowest sector of h10 at any appreciable vacancy density. (ii) The same rescaling changes a and D seen by the tokens (a_eff = Δ/c + 2cos, D_eff = D/c), so the certified Bloch data are not those of the mobile chain at the stated couplings. (iii) The ν window that selects b = 10 is only ±0.020 wide, so the walker's kinetic renormalisation would have to be known to about 1 %.

This looks repairable: embed J·o in the material–material block as C88 does, or scale ν, Δ, D by c. The first option gives up "a vacancy next to a head changes no diagonal energy". Either way claim 2 needs new text and at least a small vacancy computation. One point in the bond's favour that survives: at zero port the cell Hamiltonian is bit-independent for every J/t, so there is no analogue of C88's t²/(JK) bit field.

Verdict: structure confirmed; the c = n + 1 claim at the stated parameters is not confirmed.

### "Why seven letters fail, and why nine"

- The seven-letter result is established at zero field only (`seven_check.py` has `f=0.0`). Scanning the uniform field with the author's own T and μ (`item5_seven_field.py`) gives 2 Fermi points at n = 2 for fields outside the dimerisation gap (up to about 1.5 ‖T‖), where the README's "at most n mod 2 bands" predicts none; n = 3 stays at 2 of 6. The design still fails at the tested point, but the stated bound is a zero-field statement, tested at one parameter point and n ≤ 3. It is not a no-go: with nearly flavour-independent S (tokens close to critical) a field could sit inside several dimerised bands at once.
- Reason 1 (four data letters): plausible, not demonstrated. "Energetics alone does not do it" is asserted.
- Reason 2 (a walker letter distinct from data): convincing. My control with a contextual walker breaks F and the boundary pattern immediately.
- Reason 3 (blind carrier needs its own letters; a second head rule lets tokens cross): the second half is convincing (my second-head-rule control breaks both invariants). The first half is an argument about this design, not about alternatives.
- Reason 4 (blind letters of both charges): plausible and consistent with the measured cancellation of μ at t_d = 0, but not tested with one blind species.

## 6. Certificate

(a) Diff against `codex/work_C84`: `certify_tree.py`, `rational_bands.py`, `full_loop_integer.py` differ only in `m=n+2;S=2*n+3` → `m=n+1;S=2*n+2`; `certify_roots.py` differs in one metadata line; `loop_scan.py` and `structured.py` are identical. Confirmed. The header comments still describe C84's vertex counts.

(b) `run9.py 9 35` in `referee/cert_rerun/` (10 s): `certificate_n9.json`, `rational_bands_n9.json`, `full_loop_bands_n9.json`, `root_certificate_n9.json` and `root_certificate_exact_n9.json` are identical to the stored files apart from timing fields. Confirmed. I did not rerun n ≥ 10. The README tables agree with the stored JSON for all seven n, except that "n = 14 … above 6.5e-7" rounds a lower bound of 6.4947e-7 upwards.

(c) The graph. Reading `full_loop_integer.py` and `certify_tree.py`: path of 2n + 2 vertices, the first n + 1 with diagonal D + 4 (X_n … X_0), the rest with diagonal 4 (D_0 … D_n), path entries −1, closing entry H[X_n, D_n] = −i t_d with t_d = β 2^−e and β = 1/d_{n+1}, source at (vertex 0, vertex m) = (X_n, D_0), conversions at D_0 … D_{n−1}. That is the README graph. Gauge: with u(X_j) = e^{−ijφ_A} on the A graph and u = 1 on the data rail, the express edges become −1 and the link becomes −t_d e^{i(ψ_A − nφ_A)}; with φ_A = π/(2n), ψ_A = π this is −i t_d, equal to species B's −t_d e^{iπ/2}. So the two graphs are gauge-equivalent, T picks up one overall phase and μ none. Numerical test against the literal letter model with matched flux (`item6_gauge.py`, n = 2, 3, three parameter points each): T equals the C84-form −2 Re C with no phase to 1e-15, μ equals −2i Im A with the same sign to between 6e-14 and 3e-10 of T, the imaginary part of T and the real part of μ vanish, and the Bloch spectra agree in the same k orientation. Confirmed.

(d) `crosscheck_certificate.py` compares the certificate's floating roots with roots of the `scalar9` Bloch matrix at the certified parameters, allowing a shift of 0 or π. That does support the statement, given that `scalar9` is tied to the letter model at small n. My own version (table in item 4, last row) uses a different root finder and needs no shift: all 18 roots agree with `full_loop_bands_n9.json` to 15 digits, and the separations lie above the certified lower bounds 0.0168080 and 8.8310e-5. Confirmed for n = 9; n = 10 not repeated.

(e) Inherited estimates. C84 (10) is a three-line Neumann tail that uses only the source norm √2 α², the pair conversion norm 2g and a pair gap of 3. All three hold here at a = 4 (Gershgorin gives 4 − 2t_d), so the constant 1/4 transfers; it is in effect re-derived. The largest entry of my brute-force λ² g⁴ coefficient at a = 4 is 8e-6 α⁴ (n = 2, 3; matrix dimension 16 and 64), far below α⁴/4. The 2^50 λ⁴ remainder rests on a contour argument needing an isolated vacuum band (gap above 1/32; the clock gap is 0.236 down to 0.061 at n = 20), ‖V‖ = √2 and λ ≤ 2^−13. Those hold and the band is exactly degenerate here, so the transfer is plausible, but I did not re-derive it. The README's tag "[inherited, hypotheses checked, constants not re-derived]" is accurate. It is a two-cell statement only.

## Errors, overclaims and ambiguities

1. **Major.** "With (I1)–(I4), c = n + 1. The n = 9 member gives c = 10 > 9, and the n = 10 member gives c = 11." Together with "the same midpoint ν = 2 − (δ_b + δ_{b−1})/2 selects b" and "the vacancy dictionary and (I4) are taken over from C88 on the strength of the identical legal-space structure." The dictionary does not carry over because h10's on-site energies are genuine while C88's are contact terms (section 5).
2. "h10 = J [ w + K Π_bad ] (material-material pairs only)". The code has K, not J·K.
3. "and μ^R likewise on the right factor." Needs the species labels exchanged.
4. "The edges are −e^{iφ_s} along the express rail … and −t_d e^{iψ_s} between D_n and X_n." Orientation not given.
5. "n = 14 (e = 68, pair-sum separation above 6.5e-7)". The certified bound is 6.4947e-7.
6. "2n roots on the unit circle to 1e-68 or better". The n = 9 value is 8.2e-68.
7. "its hopping is a reflection times a positive matrix, the flavour chains dimerize, and at most n mod 2 bands cross the Fermi level." Zero field only; with a field n = 2 has one crossing band.
8. "Members with 2n simple Fermi points were found at every n tried: 2, 5, 8, 9, 10, 11, 12, 14, 16 and 20, certified at 9, 10, 11, 12, 14, 16 and 20." I found no stored output for an n = 2 member in `data/`, and did not test n = 2 myself.
9. Omission: conservation of the cell-boundary bonds is used ("fixed cell pattern") but neither stated nor proved.
10. "Sectors with headless stretches or no heads have zero or positive energy on the dead sites and lose to e_b < 0." True for h9. In F = 0 the only headless configuration is the fully headless ring; in h10 it does not lose.

## Overall assessment

The frozen nine-letter construction did what the README says in every check I could run: the bond is as written, the invariance is exact, the second-order logical Hamiltonian is a free-fermion bilinear given by the scalar formula (now checked literally to n = 4 and on a genuine ring), the scalar formula is reproduced by independent code at n = 9, and the Fermi-point counts and separations are reproduced by a different method, including at the certified n = 9 member. The certificate is C84's pipeline with the stated one-line change and the graph it assumes is the graph of this bond. So d = 9 has the same conditional status as C84's d = 11, with no finite-K budget and with a port budget that is inherited.

The mobile ten-letter claim is not at that level. The legal-space structure is right, but the bond as defined does not select the cell pattern when vacancies are present, and the certified Bloch data are not those of the dressed mobile chain. I would not quote d = 10 until the on-site terms are redefined or rescaled and a vacancy computation is shown.

What I would most want checked next:

1. The mobile chain with contact on-site energies (or rescaled ν, Δ, D): the selector against the headless and wrong-length sectors at half occupancy, and whether a certified member survives the dressing.
2. A re-derivation of the two-cell remainder constant for this bond, and the three-cell λ⁴ terms.
3. An independent rerun of the certificates at n ≥ 10 and an audit of `certify_roots.py` and `rational_bands.py`, which I did not read.

## Not done

- No rerun of `ed_two_cell.py`, `two_cell_checks.py`, `ring_checks.py`, `ring_three_cells.py`, `member.py`.
- Certificates for n ≥ 10 not rerun; `certify_roots.py`, `rational_bands.py`, `loop_scan.py`, `structured.py` not audited beyond the diff.
- The n = 10 cross-check and the 150/160-digit member runs not repeated.
- No brute force beyond n = 4, and none at member parameters (double precision is not enough there).
- Mobile ED only for one cell, bits frozen, no port and no tokens.
- The seven-letter design only at the author's parameter point.
- The infrared inputs (I1)–(I4) were out of scope.

## Files in `referee/`

`ref_bonds.py`, `compare_bonds.py` (item 1); `ring_invariants.py`, `ring_invariants_control.py` (item 2); `chain.py`, `second_order.py`, `readme_scalar.py`, `item3_run.py`, `item3_rings.py`, `item3_n4.py` (item 3); `item4_blocks.py`, `item4_fermi.py`, `item4_compare_blocks.py` (item 4); `item5_mobile_selector.py`, `item5_seven_field.py` (item 5); `item6_gauge.py`, `cert_rerun/` (item 6); JSON and log outputs beside each script. Block files without `a` in the name are a = 2.0.

## Addendum: the corrected mobile bond (README as of 22:45, 33893 bytes)

The author replaced h10 by C88's form, h10 = J embed(b9) − t (vacancy hops) − (μ/2)(n_occ ⊗ 1 + 1 ⊗ n_occ), so that the material on-site energies are contact terms. Checks are in `referee/addendum_checks.py` (output `addendum_checks.json`, `.log`); the earlier sections of this report are unchanged and describe the first version.

1. **Bond from text.** My rebuild of the corrected h10 agrees with `bond10.bond(p)` exactly (maximum entry difference 0.0 at 12 random points, K and J generic, 100 off-diagonal entries). `bond10.bond(p, onsite='genuine')` equals my first-version rebuild exactly, now with J·K as the README always wrote it.
2. **Invariance.** F, Q, head number, occupied number and the squeezed cell-boundary pattern: 0 violations on all 10^L ring words, L = 2 to 5 (500 000 transitions at L = 5).
3. **Selector.** One cell plus N_v vacancies against the headless ring, ν at the frozen midpoint. ε_cell is negative in every run: n = 4 with N_v ≤ 6 and n = 9 with N_v ≤ 4, at t/J = 25, 5, 1 and 0.2. I also measured the contact probability c in the headless ground state:

| n = 9, N_v | c | ε_cell/frozen, t/J = 25 | t/J = 5 | t/J = 1 | t/J = 0.2 |
| --: | --: | --: | --: | --: | --: |
| 1 | 0.900 | 0.902 | 0.909 | 0.940 | 1.040 |
| 2 | 0.802 | 0.805 | 0.817 | 0.865 | 1.002 |
| 4 | 0.626 | 0.630 | 0.646 | 0.708 | 0.856 |

   At t/J = 25 the ratio equals c to 0.2–0.7 % (n = 9) and 0.4–1.3 % (n = 4, including half occupancy, 0.294 against c = 0.291). So the rescaling is uniform at large t/J, with a deviation that grows like J/t. My numbers reproduce the README's table and `data/mobile_selector_n*.json` where they overlap.
4. **Is the major objection resolved?** Yes, at the level at which I raised it: to first order in J/t the corrected bond is the frozen bond times cJ, and a head no longer costs energy. What remains is what C88 also leaves open, and the README now says so (claim 2, "What is not done" item 3). The paragraph "What this does and does not show" is accurate about scope. Two things in it still need changing:
   - "At zero port the cell Hamiltonian is bit-independent for every J/t" is true only at f = 0. The field is now a contact term, no longer f times the conserved charge. For one n = 3 cell with three vacancies and f = 1e-3, the energy of a single 1-bit is (0.3905, 0.3916, 0.3905) f by position at t/J = 25 and (0.376, 0.408, 0.376) f at t/J = 1, against c = 0.391; it is exactly bit-independent at f = 0. That is a flavour-dependent field of about 0.03 (J/t) f. There is indeed no K-dependent bit field. Taking the charge term out of embed() (it commutes with everything) would restore exact bit independence.
   - "the schedule needs J/t small" should be quantified. The selector window at n = 9 is ±0.020 in ν, which allows roughly a 13 % non-uniform shift of ε_cell; the measured shift is under 1 % at t/J = 25 but 4 to 13 % at t/J = 1 for N_v ≤ 4, so t/J = 1 is already marginal below half occupancy. Separately, the certificate tolerates a normalised perturbation of 2^−35, and the first non-uniform corrections to T, μ and the field are of relative order J/t. No J/t budget exists; the mobile member is an implicitly bounded family, as in C88. The one-cell test also cannot see the quantity the selector actually discriminates, the 1.8e-4 density gap between neighbouring cell lengths.
5. **Earlier errors 2 to 10.** Fixed: 2 (J·K), 3 (species labels), 4 (edge orientation, which matches my implementation), 5 (6.4e-7), 6 (1e-67), 7 (zero field stated), 9 (conservation of B_i stated; the proof is correct), 10 (selector sentence split into h9 and h10). Not fixed: 8, n = 2 is still listed among the members with no stored output. Small point: the "8e-6 α⁴" quoted from this report is the largest matrix entry of the λ² g⁴ coefficient, not its norm.

Revised verdict for the mobile chain: the bond now has C88's conditional status at d = 10, with exact invariance and no colours, and with the same unbudgeted dependence on J/t. I have not checked wrong-length competitors, several cells, tokens or the port with vacancies, and n = 9 was not run at half occupancy.
