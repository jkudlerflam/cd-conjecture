# Margins audit: Fermi points, centres and quartic margins of the 23-letter chain, n = 1 to 23

## Claim

Setting: the frozen-cell 23-letter bond of `contrib/explicit_candidate/`, at leading order, meaning the λ²g² coefficient of the logical chain h(k) = μ + ε + T e^{ik} + T† e^{−ik}, with α_n⁴ removed.

For every n from 1 to 23, at one explicit point inside the original alphabet (schedule below), in 115- to 313-digit arithmetic:

1. Herm(−T) is positive definite and the whitened field is below 1/4, so every band crosses on both arcs. There are exactly 2n simple Fermi points, two per band.
2. The sea centres are distinct modulo π. The pair-transfer margin δ_PT = min dist(2θ_a − 2θ_b, 2πℤ) is at least 4.2 × 10⁻⁸.
3. Every number-conserving quartic momentum misses resonance. δ₄, the minimum of |p_i + p_j − p_k − p_l| mod 2π over distinct unordered pairs with repeated indices allowed (pair transfer, exchange, umklapp, particle-hole umklapp), is at least 4.2 × 10⁻⁸. For n ≥ 2, δ₂ ≥ 0.0140 and the umklapp margin δ_U ≥ 3.8 × 10⁻³.
4. No margin collapses. The two smallest δ₄ (n = 15, 21) are near-coincidences of neighbouring centres, and one dyadic step in the field or loop lifts both above 3.6 × 10⁻⁷.

## Status

Numerical, high precision, not certified. Errors are estimated by reruns at 60 more digits, not enclosed in intervals. C33's certified n = 23 point is reproduced to every printed digit. Upgrade: C33-style rational inertia and pair-phase tests at every n.

## How it was established and checked

- Author: Claude (Opus 5.5) subagent for Jonah Kudler-Flam, 2026-10-28.
- AI assistance: none beyond this subagent (no Codex call, since the work is repository-dependent code).
- Checked by: internal cross-checks only, against a literal two-cell computation from `ycell.py` (n ≤ 4), C33's and C37's float codes (n = 2 to 4), all C33/C37 positivity brackets, a brute-force k grid (n = 2 to 6), and 60-digit reruns (n = 5, 12, 17, 23). Not refereed.

## Method

Coefficient (`leading.py`). A diagonal ±1 gauge (the Jordan–Wigner string of the data left of a neutral token) makes each token Hamiltonian A_s data-independent, and a conversion at stub site j becomes c_j or c_j† times a projector on d(j). Then K = −Z†R₀CR₀CR₀Z with R₀ = (A_sL ⊗ 1 + 1 ⊗ A_sR)⁻¹, evaluated on vacuum and one-particle states in the product eigenbasis of the A_s, where every restricted source is rank one. Literal phases are kept, and nothing of size 2ⁿ is formed. K is quadratic plus a scalar to 10⁻¹³ (literal check), so T_xy = ⟨A_x|K|B_y⟩ and μ (one-particle blocks of both cells minus the vacuum value) are its Gaussian part.

Bands (`bands.py`). With Herm(−T) = LL† and W = L⁻¹, h becomes 2cos k − 2 sin k A + B with ‖A‖ < 0.36, ‖B‖ < 0.31, the same zeros and inertia. Pencil roots are refined by Newton's method to 60 digits. The inertia between roots (smallest whitened eigenvalue 0.011) rises 0 → n and falls n → 0 in every run, fixing the band assignment.

Schedule (`run_schedule.py`), the C33/C37 family. t = t_x = t_sw = 1, stub energy 4 above the vacuum (Δ = 4 − 2cos(π/(n+2))), D_x = D_M = 1000, δ = 0, φ^x_0 = π/(2n), φ^x_± = 0, port phase π/2, stub and switch phases 0, C33's loop phases (φ^d_0, φ^d_±) = (π, π/2). With β = 1/det C_express and F = Herm(−T)/(2β²), the loop is t_d = β 2^{−p}, with p set so the whitened loop ‖WμW†‖ is about 10⁻³ (C33's n = 23 value is 8.5 × 10⁻⁴), and the field is ε = 2β² 2^{−f} with f = ⌈−log₂ λ_min(F) + 2⌉, between 1/8 and 1/4 of the narrowest band (physical field λ²g²α_n⁴ε). At n = 23 this gives (p, f) = (137, 214), against C33's (136, 214). A second run uses whitened loop κ = 1/8.

## Results

v is in units of λ²g²α_n⁴ t per cell. δ_U is the smaller of the one-band (4k_F) and two-band (2k_F^a + 2k_F^b) umklapp distances. The last two columns are the κ = 1/8 run. Full lists of roots, centres, k_F and velocities are in `results/n*.json`.

| n | p | f | log10 λ_min(F) | δ₂ | δ₄ | δ_PT | δ_U | centre spread | log10 v_min, v_max | δ₄ at κ = 1/8 | δ_PT at κ = 1/8 |
| --: | --: | --: | --: | --: | --: | --: | --: | --: | :-- | --: | --: |
| 1 | – | 10 | -2.41 | 2.891e+00 | 5.0e-01 | – | 5.0e-01 | 0 | -13.8, -13.8 | 5.0e-01 | – |
| 2 | 2 | 20 | -5.20 | 7.706e-02 | 4.1e-04 | 4.1e-04 | 4.1e-03 | 2.1e-04 | -22.6, -20.8 | 1.0e-02 | 1.4e-02 |
| 3 | 8 | 29 | -7.98 | 2.302e-02 | 2.0e-06 | 2.0e-06 | 4.8e-03 | 2.1e-04 | -31.4, -29.5 | 1.0e-02 | 1.7e-02 |
| 4 | 14 | 38 | -10.77 | 2.214e-02 | 4.8e-06 | 4.8e-06 | 3.8e-03 | 2.5e-04 | -40.2, -38.1 | 1.6e-03 | 2.0e-03 |
| 5 | 21 | 48 | -13.55 | 1.802e-02 | 7.0e-06 | 7.0e-06 | 7.7e-03 | 2.5e-04 | -48.9, -47.3 | 2.6e-04 | 3.2e-03 |
| 6 | 27 | 57 | -16.33 | 1.764e-02 | 5.9e-06 | 1.6e-05 | 6.8e-03 | 2.9e-04 | -57.7, -56.0 | 8.6e-05 | 2.5e-03 |
| 7 | 33 | 66 | -19.11 | 1.651e-02 | 6.7e-06 | 6.7e-06 | 6.7e-03 | 3.4e-04 | -66.4, -64.7 | 6.5e-05 | 1.2e-03 |
| 8 | 40 | 75 | -21.90 | 1.637e-02 | 4.9e-06 | 4.9e-06 | 6.0e-03 | 2.0e-04 | -75.2, -73.4 | 2.4e-05 | 2.1e-04 |
| 9 | 46 | 84 | -24.68 | 1.576e-02 | 1.4e-06 | 2.8e-06 | 5.7e-03 | 2.5e-04 | -84.0, -82.1 | 1.5e-05 | 1.8e-04 |
| 10 | 53 | 94 | -27.46 | 1.539e-02 | 9.0e-07 | 9.0e-07 | 9.6e-03 | 1.4e-04 | -92.4, -91.2 | 9.7e-06 | 5.6e-04 |
| 11 | 59 | 103 | -30.24 | 1.515e-02 | 2.4e-06 | 2.4e-06 | 9.0e-03 | 1.8e-04 | -101.2, -99.9 | 6.5e-06 | 1.6e-04 |
| 12 | 65 | 112 | -33.02 | 1.509e-02 | 1.1e-06 | 1.1e-06 | 8.2e-03 | 2.2e-04 | -110.0, -108.6 | 1.4e-05 | 7.4e-04 |
| 13 | 72 | 121 | -35.80 | 1.491e-02 | 3.8e-07 | 3.8e-07 | 7.6e-03 | 1.4e-04 | -118.9, -117.3 | 2.8e-06 | 1.7e-05 |
| 14 | 78 | 131 | -38.58 | 1.461e-02 | 1.2e-06 | 1.2e-06 | 1.2e-02 | 2.0e-04 | -127.4, -126.3 | 9.4e-06 | 3.4e-04 |
| 15 | 85 | 140 | -41.37 | 1.456e-02 | 9.7e-08 | 9.7e-08 | 1.1e-02 | 1.4e-04 | -136.2, -135.0 | 3.5e-06 | 3.7e-05 |
| 16 | 91 | 149 | -44.15 | 1.453e-02 | 1.2e-06 | 1.2e-06 | 1.0e-02 | 1.8e-04 | -145.0, -143.8 | 1.2e-06 | 7.6e-05 |
| 17 | 98 | 158 | -46.93 | 1.447e-02 | 5.8e-07 | 5.8e-07 | 9.4e-03 | 1.3e-04 | -153.9, -152.5 | 1.3e-05 | 3.2e-05 |
| 18 | 104 | 168 | -49.71 | 1.422e-02 | 3.0e-07 | 1.8e-06 | 8.7e-03 | 1.5e-04 | -162.3, -161.4 | 9.8e-07 | 2.4e-04 |
| 19 | 111 | 177 | -52.49 | 1.423e-02 | 4.2e-07 | 4.2e-07 | 9.7e-03 | 1.0e-04 | -171.2, -170.1 | 1.6e-06 | 6.3e-05 |
| 20 | 117 | 186 | -55.27 | 1.421e-02 | 6.1e-07 | 6.1e-07 | 1.1e-02 | 1.5e-04 | -180.0, -178.9 | 1.0e-06 | 1.5e-04 |
| 21 | 124 | 195 | -58.05 | 1.420e-02 | 4.2e-08 | 4.2e-08 | 1.1e-02 | 1.0e-04 | -188.8, -187.6 | 5.9e-06 | 4.1e-05 |
| 22 | 130 | 205 | -60.83 | 1.400e-02 | 8.5e-07 | 8.5e-07 | 6.1e-03 | 1.3e-04 | -197.3, -196.5 | 4.4e-07 | 1.6e-05 |
| 23 | 137 | 214 | -63.61 | 1.402e-02 | 1.8e-07 | 2.0e-07 | 6.9e-03 | 8.9e-05 | -206.1, -205.2 | 8.1e-06 | 1.9e-05 |

- Centres lie within 2.6 × 10⁻⁴ of k = 0 (2.2 × 10⁻² at κ = 1/8). Pair transfer between neighbouring centres, usually near an extremum of the centre-versus-band curve, sets δ₄ in 18 of 22 cases.
- At κ = 1/8, δ_PT ≥ 1.6 × 10⁻⁵ everywhere, and δ₄ is set instead by crowding of the n(2n + 1) pair sums, since all Fermi points sit in two arcs (|k| from 1.23 to 1.91 at n = 23). The stronger loop is not uniformly better (n = 22 drops to 4.4 × 10⁻⁷). δ₄ of order 10⁻⁷ to 10⁻⁶ near n = 20 seems natural for one common field.
- Velocities. Every Fermi velocity lies between 2 and 10^{2.35} times λ_min(Herm(−T)), and v_max/v_min ≤ 10^{2.04}, although Herm(−T) spans 61 decades at n = 23. The symbol is graded, so every zero is a hybridized crossing with slope set by the narrowest scale (confirmed by finite differences at n = 3, 4). log10 v_min falls by about 8.7 per flavour. The infrared scale is the narrowest band, not the bandwidth.
- Tolerance. Whitened velocities are 1.98 to 2.00, so a quadratic correction E moves no Fermi point by more than δ₄/8 if ‖E‖ ≲ λ_min(Herm(−T))δ₄/4. With C33's g⁴ remainder bound this needs g ≲ 10⁻¹⁰⁸ at n = 23. C33 used 10⁻¹¹⁰, while C37 takes 10^{−2.5×10¹⁶} because it budgets against its proved δ_n, not the actual margin.

## Comparison with C33 and C37

- Code level. T and μ agree with C33's `structured.py` to 1.3 × 10⁻¹¹ of |T|, and with C37's except μ at n = 4 (1.2 × 10⁻⁶ of |T|), where the vacuum scalar is 10¹⁰|T| and my double-precision literal check hits the same floor. That is cancellation, not disagreement, and C37's modular certificate does not use it.
- Positivity. log10 λ_min(F) and λ_max(F) lie inside every C33 and C37 certified bracket, n = 1 to 23 (`compare.json`).
- C33 at n = 23 (t_d = β2⁻¹³⁶, ε = 2β²2⁻²¹⁴). All 46 roots agree to 1.1 × 10⁻¹⁵ and the centres to 6 × 10⁻¹⁶. δ₂ = 0.0140196381522 and δ₄ = 3.99638535 × 10⁻⁷ agree in every digit C33 prints, with the same witness (pair transfer between bands 8 and 9), above C33's certified 0.0140193 and 3.1558 × 10⁻⁷. The minimum speed, 10^{−62.357} in units of 2β², is the same.
- C37 (loop phases (π, 0, 0), t_d = 10⁻²⁰⁰⁰, field 10⁻¹⁰⁰⁰). With C37's phases and a small loop and field in C37's order (whitened 5 × 10⁻²⁶ and 10⁻⁵), every pair sum is distinct, as C37's certificate says. To first order the pair-transfer splitting scales with t_d and the two-band umklapp splitting with the field, so at C37's literal values they are 10⁻¹⁹⁹⁴ to 10⁻¹⁸⁹³ and 10⁻⁹⁷⁷ to 10⁻⁷⁹⁴. They are positive and far above C37's proved δ_n (10^{−1.07×10¹⁵} at n = 23), but of no quantitative use. At n = 22 C37's phases give a linear pair-transfer coefficient a thousand times below its neighbours (3 × 10⁻⁷), the only near-collapse I found.

## Flags and fixes

The smallest δ₄ at the main schedule are n = 21 (4.2 × 10⁻⁸) and n = 15 (9.7 × 10⁻⁸), both pair transfer between neighbouring centres.

- n = 21 (`fixscan_n21.json`). f + 1 gives 3.7 × 10⁻⁷, p − 3 gives 4.5 × 10⁻⁷, κ = 1/8 gives 5.9 × 10⁻⁶. Other single steps give 2.1 × 10⁻⁸ to 1.2 × 10⁻⁷.
- n = 15 (`fixscan_n15.json`). f + 1 gives 3.9 × 10⁻⁷, f + 3 gives 8.3 × 10⁻⁷, κ = 1/8 gives 3.5 × 10⁻⁶.
- For quantitative use, scale t_d with β and the field with λ_min(F) as here, not C37's literal values.

## Verdict

The finite-dimensional gap is closed at the high-precision numerical level for every n from 1 to 23, at one explicit schedule inside the 23 letters, in agreement with C33 at n = 23 and with C37's nonvanishing certificate. The margins are small, δ₄ between 4 × 10⁻⁸ and 2.4 × 10⁻⁶ for n ≥ 9, but none vanishes, and the two worst cases are repaired by one dyadic step. Open are certification at n ≠ 23, the finite-g transfer at every n (C33 did it at n = 23), and the infrared and frozen-cell questions of `contrib/explicit_candidate/` §6.

## Reproducing

numpy, scipy and the standard library, run here with `python3 -B`. `ycell.py` and C33's and C37's `structured.py` are imported read-only.

- `literal_check.py 4` (1 min) → `literal_check.json`.
- `run_schedule.py` (1 h on 9 cores) → `results/n*.json`, `run_schedule_stdout.txt`.
- `band_check.py`, `compare.py`, `summarize.py`, `make_table.py` (a minute each) → `band_check.json`, `compare.json`.
- `precision_check.py n …` (n = 23 about 12 min) → `precision_check.json`.
- `fixscan.py n` (n = 21 about 25 min) → `fixscan_n*.json`.

## Relation to the rest of the repository

This addresses item 1 of `contrib/explicit_candidate/` §6 at leading order in g for n ≤ 23. STATUS could add that distinct centres and non-resonant margins now have an independent numerical check through n = 23 (δ₄ of order 10⁻⁷), and that C37's literal point gives margins of order 10⁻¹⁹⁰⁰.
