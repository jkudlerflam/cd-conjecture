# C192 checkpoint: direct attempt on d = 2 (coordinator, Claude, 2026-11-09)

Preliminary. Nothing here is a verdict yet. All numbers come from C176's converged runs (codex/work_C176/typhon_out/out/, last line of each jsonl), symmetric line φ = π/2, |w| = 1. Scripts: sq2.py and pairfit.py in this folder, run from codex/work_C176 with run tags as arguments.

## Question

In C173's effective model, is the pseudospin (valley) sector gapless (c = 2) or does it have a small gap (c = 1)? Standard two-component Luttinger liquid language, valleys A and B, pseudospin backscattering g₁ cos(√8 φ_s) of dimension 2K_s.

## Calibrated findings so far

1. Charge parameter. Small-q slope of the density structure factor, fit c + a|q| + bq² with a cos² taper to R = 90, 120.
   - Controls: hard rods read π·slope = 0.27 to 0.28 (exact 0.281), two XX legs read 0.96 to 0.99 (exact 1).
   - Candidates at L = 256 (n = 1/8, 5/32, 3/16): 1.65 to 1.84. So K_c ≈ 1.8, with 1/K_c ≈ 0.55 to 0.6. The R = 60 window is noisier (1.2 to 1.5).
   - An earlier version of the fit had no constant term and read the controls 30% low. That was the bug.
2. Inter-valley pair correlator. D(r) e^{−iπr} at pair separation 2 and 3 is real and positive at every r from 3 to 123 (L = 256, n = 3/16). It decays as r^{−0.69} between r = 12 and r = 83 (0.0114 to 0.0030 at separation 2). Three-component fits give 0.62 to 0.69 for the momentum-π piece on every window, and a noisy 0.7 to 2.3 for the same-valley pieces at ±2k₁.
   - In a two-component liquid the momentum-π pair has exponent 1/K_c + K_s. An exponent of 0.69 is incompatible with K_s ≥ 0.69 at any K_c. With K_c = 1.8 the pseudospin part is at most about 0.14.
   - That is the Luther–Emery form (exponent 1/K_c, pseudospin pinned), or a liquid with K_s^eff far below 1.
3. Valley structure factor. P_j = c†_{j+1}c_j + h.c. is 2cos k in momentum space, about +1.94 in valley A and −1.94 in valley B, and it is also the nearest-neighbour part of the charge current. Free two-valley fermions give π·slope = 2[cos²(k₁+k_F) + cos²(k₁−k_F)] = 3.46 at n = 3/16. The data give 4.2 (R = 120) to 5.4 (R = 60), a ratio of 1.2 to 1.55, drifting down with R. Read as K_s it says K_s^eff > 1 and gapless at r of 12 to 40. No control validates this estimator yet.
4. The lowest-order g-ology estimate from C175's threshold scattering lengths gives a gapless pseudospin sector for n above about 0.08 at |w| = 1, with K_s ≈ 1.7 to 1.9. It is not reliable at k_F ≈ 0.3, where effective-range corrections are of order one.

## The tension

- Items 1 and 2 together look like a pair liquid (c = 1): pairs at momentum π with quasi-long-range order, exponent near 1/K_c.
- The entropy reads c ≈ 1.9 from L = 128 to 256, and L·Δ_P is flat near 2.4. A pair liquid with a short pseudospin length would give c = 1 and L·Δ_P growing with L.
- Item 3 says the valley fluctuations are not suppressed.
- No self-consistent (K_c, K_s) of the standard two-component liquid fits all of these. KT monotonicity still holds: K_s^eff < 1 at any scale with g₁ ≠ 0 means an eventual gap and c = 1.

## Next steps

- Odd-even test from the energies: compare Δ_P = E(N+1) + E(N−1) − 2E(N) with the smooth second difference [E(N+8) + E(N−8) − 2E(N)]/64 at L = 256 (N = 32, 40, 48 and their ±1 partners). The excess is the cost of one unpaired particle.
- Check the same-valley pair exponent and the single-particle envelope with r-window averaging.
- Read C189 (K_s estimators with controls, infinite chain) and C191 (finite-momentum phase shifts, KT flow) when they report. Don't duplicate them.
- Then write the REPORT section and tell Jonah the verdict.

## Added findings (second checkpoint)

5. Odd-even energy test, symmetric line, |w| = 1. Δ_P = E(N+1) + E(N−1) − 2E(N), and E'' is the smooth curvature from ΔN = 4, 6, 8.

| L | L·Δ_P (three densities) | L·E'' (middle N) | L·(Δ_P − E'') |
|---|---|---|---|
| 128 | 2.283, 2.337, 2.383 | 0.816 | 1.521 |
| 192 | 2.280 (N = 30) | 0.799 | 1.480 |
| 256 | 2.181, 2.239, 2.381 | 0.791 | 1.448 |

A fit excess = a/L + b gives a = 1.593, b = −0.00057. So there is no pair-breaking gap above about 3e-4. Also v_c/K_c = (2/π)·0.79 ≈ 0.50, and with v_F = 1.84 at n = 5/32 the Galilean estimate is K_c ≈ 1.91, consistent with the structure factor.

6. Valley imbalance p against the valley offset δ (C181-type runs off the line, n = 3/16).

| run | p | δ | p/δ |
|---|---|---|---|
| p154_L256 | −0.1064 | 0.1193 | 0.892 |
| p152_L256 | −0.1848 | 0.1967 | 0.939 |
| p150_L256 | −0.2808 | 0.2742 | 1.024 |
| p154_L128 | −0.1286 | 0.1193 | 1.079 |
| p152_L128 | −0.1968 | 0.1967 | 1.000 |
| p150_L128 | −0.2643 | 0.2742 | 0.964 |

The free two-valley value is 0.769. The response is finite and linear down to the smallest offset, 1.16 to 1.33 times free at L = 256. A pseudospin gap larger than about δ_min·(small factor) would show as a threshold, and none is seen.

7. Scattering channels (from C175). Unlike-valley even channel a_+ = 20.71 (weakly attractive, carries the shallow bound state, E_b = ε''/a_+² = 0.0175 against the exact 0.0195). Unlike-valley odd channel a_− = 0.289 (near hard core). Like-valley odd-wave a_s = −3.17 (strong attraction, k_F a_s ≈ −0.93). Pseudospin SU(2) is strongly broken in the direction that favours K_s > 1.

8. Species-basis estimate. Each valley alone has K_0 ≈ 1.6 from the odd-wave attraction. Fixing the inter-valley density coupling to reproduce K_c = 1.85 gives K_s ≈ 1.4, so the backscattering cos(√8 φ_s) has dimension 2K_s ≈ 2.8 and is irrelevant.

## The tension, sharpened

Items 1, 3, 5, 6, 8 all read as a gapless pseudospin sector with K_s > 1 (c = 2 stable). The one outlier is the momentum-π pair correlator, which decays as r^(−0.69). In the standard two-component liquid its exponent is 1/K_c + K_s ≥ 0.53 + K_s, so 0.69 would need K_s ≈ 0.15. Still open.

## Third checkpoint: the outlier is resolved (frame), K_s > 1

9. Frame. C191 found that the two valleys are weakly coupled only as mutually bosonic fields. In that frame the lattice boson is e^{iθ_A} (plus valley B), the inter-valley pair is e^{i√2 θ_c}, and the Jordan–Wigner fermion carries a string over the other valley. So
   - boson one-body exponent η_b = 1/(4K_c) + 1/(4K_s),
   - pair exponent at momentum π, η_pair = 1/K_c,
   - fermion exponent η_b + min(K_c, K_s).
   In the fermion frame the boson exponent is at least 1, so a measured η_b below 1 settles the frame by itself.

10. Test, predicted before looking (η_b ≈ 0.32). onebody.py, pairs of sites symmetric about the chain centre, fit to tan(πr/2L)^(−η).
   - c = 2 control (two XX legs): η_b = 0.500 to 0.502 at L = 192, 256 (exact 1/2).
   - Candidate, largest window, L = 128, 192, 256: n = 3/16: 0.325, 0.324, 0.326. n = 5/32: 0.335, 0.333, 0.335. n = 1/8: 0.347, 0.349, 0.351. Residual rms below 0.005 in log.
   - The boson momentum peak at the valley momenta grows as W^0.70 (2.59, 3.46, 4.22 at W = 64, 96, 128).

11. Pair exponent with the open-chain conformal form (pairconf.py). Control c2 (d = 3): 0.995 to 1.007 (exact 1). Candidate at L = 256, largest windows: 0.619 (3/16), 0.618 (5/32), 0.612 (1/8), with 0.635 and 0.621 at L = 128 and 192 for 3/16.

12. K_s. 1/K_s = 4η_b − 1/K_c, with 1/K_c = 0.62 (pairs) or 0.55 (structure factor).

| n | 4η_b | K_s (pairs) | K_s (structure factor) |
|---|---|---|---|
| 3/16 | 1.30 | 1.47 | 1.34 |
| 5/32 | 1.34 | 1.39 | 1.27 |
| 1/8 | 1.40 | 1.27 | 1.17 |

   K_s > 1 needs only K_c < 3.3 (n = 3/16) or K_c < 2.5 (n = 1/8), against a measured 1.6 to 1.8. The backscattering cos(√8 φ_s) has dimension 2K_s ≈ 2.3 to 2.9, irrelevant. No drift toward 1 between L = 128 and 256 (η_b moves by 0.001 to 0.004).

13. Reading. On the symmetric line the effective model is a c = 2 two-component liquid, K_c ≈ 1.6 to 1.8, K_s ≈ 1.2 to 1.5. The hidden pseudospin gap of the Yang–Gaudin template is not there. K_s falls with density, consistent with C191's separatrix near n ≈ 0.07.
