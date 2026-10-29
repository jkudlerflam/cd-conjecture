# Red team of the d = 13 mobile chain (C88)

## Claim under attack

C88 (`codex/work_C88/`) builds a mobile chain with eleven material letters and two vacancy colours, d = 13, whose pin-free length-15 sector carries C85's certified n = 13 material member. The conditional infrared count is c = 13 + 1 = 14 > 12 = d - 1. The member (`MEMBER.json`) fixes g = 2^-272, λ = 2^-608, K ≥ 2^2560, half occupancy, and states the mobile K condition "with fixed positive J, t". This report tries to show that the chain does not have c = 14 at the theoretical-physics standard (spec `specs/C97_redteam_d13.md`).

## Status

No fatal failure was found. Two of the attacks give precise failures of the member as written, and both are fixable inside the construction because J/t and K are free parameters with no lower bound on J/t and no upper bound on K.

- F1 (serious as stated, fixable). The vacancy fluid moves the cell-length selector at second order in J/t. At the scheduled μ_sel the globally selected length leaves 15 for 14 at J/t ≈ 0.0056. A fixed positive J/t of order one does not select length 15.
- F2 (serious as stated, fixable). The dictionary "mobile payload band = α_contact × frozen band", on which the use of the frozen certificate rests, holds only when D J/t is small (D = 1000 is the rail energy). The per-step deviation of the virtual rail amplitude is about -0.16 D J/t, and it is of order one at D J/t ≳ 1. Staying inside the certificate's allowance needs J/t ≲ 2^-47 (normalized 2^-35 allowance) or J/t ≲ 2^-175 (the raw tolerance of PROOF (30)), and then K ≳ 2^2593 through C93's colour condition.

Everything else I attacked is harmless at J/t this small: the scale ledger, the vacancy–payload couplings and umklapp terms at half occupancy, phase separation, and the fixed-marker-number defect pairs. The repaired schedule (J/t below both thresholds, K raised accordingly) is self-consistent. The verdict on c = 14 is unchanged. It remains a conditional candidate resting on the same inherited infrared inputs as before.

## How it was established and checked

- Author: Claude (Opus 5.5) red-team subagent for Jonah Kudler-Flam, 2026-10-29. No Codex call. Python 3 with NumPy and SciPy only, all local, small sizes.
- Read: STATUS.md (smallest-candidate section), REPORT §9.28 and §9.41 to §9.59, C88 RESULT and PROOF, C91 and C93 RESULT, C94 RESULT, and the referee reports `../referee_C84_C85/` and `../referee_C88/`. Nothing outside this folder was modified.
- Tags: [proved] means exact arithmetic or an algebraic argument, [numerical] means floating-point evidence at the stated sizes, [heuristic] means a physics argument not made rigorous, [assumed] means an input from elsewhere.
- Scripts (each writes the JSON of the same name, all run in under a minute except `selector_shift.py`, which takes about a minute for the N = 12 gap ring):
  - `gapfluid.py`: gap-space vacancy fluid (C91's H0), ground state, static contact susceptibility.
  - `selector_shift.py` → `selector_shift.json`: attacks (iii) and (iv), with a validation against exact diagonalization.
  - `carrier_toy.py` → `carrier_toy.json`: attacks (i) and (iv), the virtual rail traversal with mobile vacancies.
  - `hierarchy.py` → `hierarchy.json`: attacks (i), (ii) and (iii), the scale ledger, lock-in dimensions and compressibility.

## Model used for the mobile computations

[proved, C88 eq. (8) and C91] On the legal subspace at K = ∞ and zero port, the chain is H = t H_gap + J Σ_j c_j b_j. H_gap moves a vacancy from gap j to gap j+1 with unit amplitude, c_j = 1(g_j = 0) is the contact on material bond j, and b_j is the literal material bond operator. In C88 every material term is a contact term, the D-scale and μ_sel diagonals included (PROOF after (1)). For the vacuum clock, b_j is the unit clock hop on bond j, and the marker's μ_sel is split half to each of its two bonds.

[numerical] The periodic gap ring reproduces the free-fermion contact probability on the physical ring to 1e-15 at N = 6, 8, 10, 12 (α = 0.29266, 0.29473, 0.29568, 0.29620, tending to 1/2 - 2/π² = 0.29736). This validates the mapping.

[numerical] The static contact susceptibility S_ij = Σ_n ⟨0|c_i|n⟩⟨n|c_j|0⟩/(E_n - E_0) (in units 1/t) has a local part s0 = 0.09852, s1 = 0.00238 and a uniform part -A/N with A = 0.07818 at N = 12. Beyond nearest neighbours the residual spread is below 6e-5. The values at N = 8, 10, 12 agree to about 1e-3. A agrees with C95's density-relaxation coefficient 0.0786/t.

## Findings

### (i) Scale hierarchy (`hierarchy.json`)

[proved, exact integers for d_15] log2 of the scales in units of J: β = -149.57, E_unit = -858.58, W_ref = J α_c λ² E_unit = -2076.33, λ² = -1216, zero-port splitting 2^84/K = -2476, vacuum-tail J/K = -2560. Relative to W_ref, the zero-port splitting is 2^-399.7, the A and Z grammar terms are 2^-483.7 or smaller, and the finite-g tail is 2^-244.9. The λ² scalar part is 2^860 times W_ref but payload-blind, since a g⁰ round trip on the path-graph carrier crosses every bit an even number of times. Its sign factors z_s(a) therefore square to one, whatever the vacancies do in between [proved]. Harmless.

[proved] Colour field t² w_r/(JK), w_r ≤ 0.686. At K = 2^2560 it stays below PROOF (30)'s raw material tolerance (2^-162.6 W_ref) only for t/J < 2^160.8, and below the normalized 2^-35 allowance for t/J < 2^224.6. C93's K ≥ max(2^2560, 2^2243 (t/J)²) covers both. The colour field against the λ² E_unit band is harmless as long as K grows with (t/J)². Once F2 forces t/J up to about 2^175, K = 2^2560 is no longer enough and K ≳ 2^2593 is needed. This is fixable because nothing bounds K from above (err_K and the F4 remainder both decrease with K).

### F2. The payload dictionary needs D J/t small (`carrier_toy.json`)

[numerical] Toy: one rail excitation on material letters 0..4 of an 8-letter, 8-vacancy gap ring, with hops -J c_k and diagonal (D/2)(c_(k-1) + c_k) J, as in C88. I compare the exact traversal amplitude ⟨vac, r|(H - E_vac)^-1|vac, 0⟩ with C88/C93's fast-vacancy dictionary c_k → α.

| D | t/J | lowest level of H - E_vac (J) | per-step deviation | deviation / (D J/t) |
| --: | --: | --: | --: | --: |
| 10 | 1 | 1.02 | -0.574 | -0.057 |
| 10 | 100 | 2.41 | -0.0157 | -0.157 |
| 10 | 10^5 | 2.44 | -1.6e-5 | -0.159 |
| 30 | 10 | 5.82 | -0.315 | -0.105 |
| 30 | 10^4 | 8.33 | -4.9e-4 | -0.162 |
| 100 | 10 | 10.68 | -0.581 | -0.058 |
| 100 | 10^5 | 28.96 | -1.6e-4 | -0.163 |

The dictionary is exact as J/t → 0, with a per-step correction of about -0.16 D J/t. For D J/t ≳ 1 it fails at order one. The amplitude ratio over four steps falls to 0.003 at D = 100, t = J, and the rail gap is set by t instead of αDJ (1.44 J against αDJ ≈ 29 J). Mechanism: the diagonal is itself a contact term, so a rail letter fenced by vacancies costs nothing, and slow vacancies decide which bonds the traversal can use.

[heuristic, application to C88] The real carrier has 29 vertices with diagonals 1004 and 4. The certified material member therefore describes the mobile chain only if 29 × 0.16 × 1004 × (J/t) is below the allowance. That gives J/t ≲ 2^-47 with the normalized 2^-35 allowance, or J/t ≲ 2^-175 with PROOF (30)'s raw tolerance τ/W_ref = 2^-162.6. The two rails have different diagonals, so the correction is not a pure congruence of the Bloch matrix, and Sylvester inertia does not protect the Fermi-point count. I did not compute its congruence-invariant part. At "fixed positive J, t" of order one, the mobile payload band is not the certified band, and the 26-Fermi-point count does not transfer. Classification: serious as stated, fixable by adding J/t ≤ 2^-175 (or the sharper bound a whitened estimate would give) to the member. C91 already used J/t = 2^-5000 for its ultraviolet budget, so the repair is in the campaign's spirit, but MEMBER.json does not state it.

### (ii) Vacancy–payload couplings and umklapp (`hierarchy.json`)

[proved, kinematics] Material density ρ_m, b = 15. A cell harmonic of squeezed momentum 2πℓ/15 is uniform when ℓρ_m/15 is an integer, with charge dimension (ℓ/15)² K_f. At ρ_m = 1/2 the first one is ℓ = 30, dimension 4K_f. This coincides with the ordinary half-filling lattice umklapp, and it is the "at least four" of C88 §10. Relevant lock-in (dimension below 2 at K_f = 1) occurs only at ρ_m ∈ {1, 15/16, 15/17, 5/6, 15/19, 3/4, 5/7}, that is, six or fewer vacancies per cell. Those densities should be avoided, and half occupancy does avoid them.

[heuristic] K_f at small J/t. The mean contact energy per bond is e_b = E_15/15 = -0.00483 J (attractive and small, because μ_sel nearly cancels the clock energy). An XXZ estimate gives K_f - 1 ≈ 1.5e-3 (J/t). The umklapp dimension stays 4.00 for every J/t ≤ 1. Harmless.

[heuristic] Payload operators. The h = 0 non-forward quadratic and quartic terms are excluded by the certified root and pair-sum margins (inherited). The forward vacancy–payload density coupling from dα/du (C93) has strength of order W_ref. A multicomponent Luttinger liquid stays stable when g² < v_v v_p, and here g²/(v_v v_p) ~ W_ref/(30 t) ~ 2^-2081 J/t. Vacancy-mediated payload–payload interactions are of relative size W_ref/t. All harmless.

### (iii) Phase separation and crystallization (`hierarchy.json`)

[heuristic, mean field in the contact energy] f(ρ) = -(2t/π) sin πρ + J e_b [ρ² - sin²(πρ)/π²] is convex for J|e_b| sin πρ < πt/2, so the attractive instability needs J/t > π/(2|e_b|) ≈ 326. The d0 = 4464 wrapper failed because its length-fluid hopping was the smallest scale and its attraction was 10^4 to 10^5 times larger. Here the vacancy hopping t is the largest scale, and the contact attraction is 0.005 J per bond. A crystal locked to the cells needs one of the low-vacancy commensurate densities above, not half occupancy. The strong local feature is the marker contact, +μ_sel/2 = 0.94 J per marker bond, which attracts vacancies to markers. At J ≪ t it is a weak periodic potential with no gap opening at ρ_m = 1/2 (its only uniform harmonic is the dimension-4 term). It does drive F1. Harmless.

### F1. The selector moves at J/t ≈ 0.0056 (`selector_shift.json`)

[proved, standard second-order perturbation theory for J ≪ t with the material clock as spectator] E = J α Σ⟨b_j⟩ - (J²/t) Σ_ij S_ij ⟨b_i b_j⟩ + O(J³/t²), with the material in the product of cell clock ground states.

[numerical, validation] Exact diagonalization of the coupled gap-fluid and clock model (two cells of length 3 on a 6-gap ring, and two of length 4 on an 8-gap ring) gives second-order coefficients -0.27200 and -0.42550 at J/t = 0.02. The formula with the exact finite-ring S gives -0.27176 and -0.42489. The residual scales as J³ (residual/J³ ≈ -0.012 and -0.030, constant across J/t = 0.02, 0.05, 0.1).

[numerical, scheduled μ_sel = 1.889188] Per material site, e(l) = J α c1(l) - (J²/t) c2(l). c2 decreases with l, so vacancies favour shorter cells: a marker bond releases μ_sel/2 when a vacancy sits there, and shorter cells have more marker bonds. Against l = 15:

| l | α_c [c1(l) - c1(15)] (J) | c2(l) - c2(15) (J²/t) | J/t at which e(l) = e(15) |
| --: | --: | --: | --: |
| 14 | 9.56e-6 | 1.719e-3 | 0.0056 |
| 13 | 4.72e-5 | 3.679e-3 | 0.0128 |
| 12 | 1.29e-4 | 5.930e-3 | 0.0218 |
| 10 | 5.49e-4 | 1.157e-2 | 0.047 |
| 16 to 21 | positive | negative | never |

So for J/t above about 0.0056, the ground state over all marker numbers at fixed N prefers length-14 cells. Those carry twelve active bits, which the certificate does not cover. This agrees in scale with C94's b = 15 numerics, where ε fell from 3.16e-5 at J/t = 0.04 to 8.2e-7 at 0.0025 against α γ ≈ 8.4e-6. The second-order shift acts as an O(J/t) renormalization of μ_sel, so it can also be removed by retuning μ_sel at fixed small J/t. Classification: serious for the member as written ("fixed positive J, t"), fixable by J/t ≪ 0.0056 (implied by F2's J/t ≤ 2^-175) or by an O(J/t) retuning of μ_sel. The higher orders are not bounded here, and the all-order version is C92/C94/C95's open charge comparison.

### (iv) Transfer of the frozen-sector results (`selector_shift.json`)

[numerical] At fixed N and fixed marker number M (both conserved), wrong lengths come in pairs. The second-order term raises their cost: for (14, 16) the cost is J α (9.0e-4) + (J²/t)(2.2e-4), and similarly for (13, 17) and (12, 18). Inside a fixed-M sector the vacancies do not destabilize length 15. F1 is purely a grand-canonical (marker-number) instability.

[assumed, inherited] Pins: the referee's replacement bound (R1) plus C92's contact bound ⟨A_j⟩ ≥ 1/4096 give a pin cost at least J/4096, which is far above every J/t-suppressed correction once J/t ≪ 0.0056. Colour walls cost JK. Marker-free and all-vacancy sectors are excluded (C92). I found nothing new here beyond F1 and F2. What remains is the all-order, volume-uniform dressed-defect comparison already listed as open.

### (v) Other observations

- [numerical] In the mobile chain the gap of a virtual rail excitation is min(αDJ, O(t)), not D J (toy table above). Any estimate that uses D J as a mobile energy denominator, for example in a dressed-gap threshold for C48, needs t ≫ D J. That is the same condition as F2.
- [heuristic] A cosmetic point: C88 §10 quotes the charge dimension (Q_h/30π)² at K_f = 1. At finite J/t it should carry K_f. This changes nothing at the J/t that F1 and F2 require.
- [assumed] I did not attack the inherited inputs (C48 massive matching, Lemmas N and TP, Gaussian completeness, conformal-observable control). They are the same as for the d = 14 and d = 24 candidates.

## Gap table

| # | Item | Size | Class |
| :-- | :-- | :-- | :-- |
| F1 | Second-order vacancy shift of the length selector | crossing to l = 14 at J/t = 0.0056 | serious as stated, fixable (J/t ≪ 0.0056 or retune μ_sel) |
| F2 | Payload dictionary α × frozen band | per-step deviation -0.16 D J/t, order one for D J/t ≳ 1 | serious as stated, fixable (J/t ≲ 2^-175 raw, 2^-47 normalized) |
| F2' | K needed once F2 is imposed | K ≳ 2^2593 via C93 | fixable (K is unbounded above) |
| H1 | K-suppressed payload terms (splitting, A, Z, colour) | ≤ 2^-399.7 W_ref at t ≈ J | harmless |
| H2 | Lock-in and umklapp at ρ_m = 1/2 | dimension 4K_f ≈ 4 | harmless (relevant only for ρ_m ≥ 5/7) |
| H3 | Forward vacancy–payload coupling | g²/(v_v v_p) ~ 2^-2081 J/t | harmless |
| H4 | Phase separation | needs J/t > 326 | harmless |
| H5 | Fixed-M wrong-length pairs | second order raises cost | harmless |
| H6 | Rail gap set by t at slow vacancies | min(αDJ, O(t)) | part of F2 |

## Verdict

The d = 13 chain survives this red team as a conditional c = 14 candidate, but not with the member as written. "Fixed positive J, t" is false in two precise ways. The vacancy fluid moves the selected cell length to 14 at J/t ≈ 0.0056 (second order, checked against exact diagonalization). The frozen certificate describes the mobile payload band only when D J/t is small, with per-step error about 0.16 D J/t. Both are cured by making J/t extremely small (J/t ≲ 2^-175 suffices for the raw tolerance) and raising K to about 2^2593. Nothing I found bounds J/t from below or K from above, so the repaired schedule is consistent. No relevant coupling, phase separation or lock-in appears at half occupancy. MEMBER.json should state the J/t condition explicitly. The inherited infrared inputs are untouched.

## Reproducing

From this folder: `python3 hierarchy.py`, `python3 carrier_toy.py`, `python3 selector_shift.py` (the last writes `selector_shift.json` and its log). `gapfluid.py` is a library imported by the other two.
