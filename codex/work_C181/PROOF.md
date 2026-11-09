# C181 derivations: the partially polarized pair route off φ = π/2

This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

Standard: theoretical physics. Tags [proved], [derived], [numerical], [heuristic], [assumed], [literature] open the paragraph they apply to. Every numerical claim names its file in this folder.

Model (3.2) of C173 in units |t₂| = 1: blockade, NN hop w = u e^{iφ}, NNN hop −1, V₂ = −2, V₃ = −2/3. Band ε(k) = 2u cos(k − φ) − 2 cos 2k. Valley a is the lower one off π/2, valley b the upper one, δ = ε_b − ε_a ≥ 0.

## 1. The valley offset is a Zeeman field, and the valley masses stay equal to first order

[proved, algebra] Put φ = π/2 − ϑ. Then cos(k − φ) = sin(k + ϑ) = sin k + ϑ cos k + O(ϑ²), so

    ε(k) = ε₀(k) + 2uϑ cos k + O(ϑ²),   ε₀(k) = 2u sin k − 2 cos 2k.      (1.1)

The symmetric valleys sit at sin k = −u/4, cos k = ±c, c = √(1 − u²/16). The first-order shift moves each minimum by ±2ucϑ without moving it to first order in energy curvature, so

    δ = 4ucϑ + O(ϑ³),   1/m_{a,b} = 8 − u²/2 ± 2ucϑ + O(ϑ²),      (1.2)

and the relative mass imbalance is η ≡ (1/m_a − 1/m_b)/(1/m_a + 1/m_b) = δ/(16 − u²) + O(ϑ²). At u = 1 the slope is dδ/dϑ = 3.873, at u = 2 it is 6.928 [numerical check, fewbody_w2_R100.out and fewbody_w1_R400.out reproduce both slopes to 4 digits]. At the offsets that matter here (δ up to about 0.6), η is at most 0.05. So off π/2 the leading effect is a pseudospin field coupling to N_a − N_b, the analogue of the Zeeman field of the attractive Gaudin–Yang gas, with a valley-mass imbalance one order smaller than in C164's accidental two-valley bands [derived].

The symmetry Σ (C180) is lost off π/2. The antiunitary Θ = I·K (inversion times conjugation) survives in the parity class and keeps each valley (C180 claim 1), so the antiunitary constraints on each valley's own data stay.

## 2. Few-body data as functions of φ

Code: fewbody.py (pair by banded relative-coordinate diagonalization over K; trimer by nbody_box.py, gap coordinates with hard walls at gap R, parabolic minimum over K near K_pair + k_a), ad_scat.py (atom–dimer levels in the same box).

### 2.1 The pair

[numerical, fewbody_w2_R100.out, fewbody_w1_R400.out] The ab pair binding relative to ε_a + ε_b, E_b(φ) = ε_a + ε_b − min_K E₂(K), is constant to four digits over the whole range:

| u | φ | δ | E_b | 1/M_pair |
| :-- | :-- | :-- | :-- | :-- |
| 2 | π/2 | 0 | 0.12255 | 2.251 |
| 2 | 1.555 | 0.109 | 0.12260 | 2.248 |
| 2 | 1.550 | 0.144 | 0.12264 | 2.246 |
| 2 | 1.540 | 0.213 | 0.12275 | 2.241 |
| 2 | 1.52 | 0.352 | 0.12310 | 2.225 |
| 2 | 1.48 | 0.629 | 0.12432 | 2.167 |
| 1 | π/2 | 0 | 0.019527 | 3.374 |
| 1 | 1.565 | 0.0225 | 0.019526 | 3.374 |

So the pair is a robust composite, and what changes with φ is only how it competes with two a atoms: E₂ − 2ε_a = δ − E_b. Pairs lie below two a atoms iff δ < E_b, i.e. ϑ < 0.0050 at u = 1 and ϑ < 0.0177 (φ > 1.5531) at u = 2. This reproduces C175 section 5 and extends it to φ = 1.48 [numerical].

### 2.2 The trimer aab

[numerical, same files] B3 = (E_pair,min + ε_a) − E₃,min:

| u | φ | δ | B3 | 1/M₃ |
| :-- | :-- | :-- | :-- | :-- |
| 2 | π/2 | 0 | 0.01888 | 1.464 |
| 2 | 1.565 | 0.040 | 0.01806 | 1.476 |
| 2 | 1.560 | 0.075 | 0.01738 | 1.485 |
| 2 | 1.555 | 0.109 | 0.01672 | 1.493 |
| 2 | 1.550 | 0.144 | 0.01609 | 1.501 |
| 2 | 1.540 | 0.213 | 0.01490 | 1.516 |
| 2 | 1.520 | 0.352 | 0.01276 | 1.538 |
| 2 | 1.500 | 0.490 | 0.01089 | 1.553 |
| 2 | 1.480 | 0.629 | 0.00925 | 1.559 |
| 1 | π/2 | 0 | 0.000597 | 2.314 |
| 1 | 1.568 | 0.011 | 0.000587 | 2.315 |
| 1 | 1.565 | 0.022 | 0.000576 | 2.316 |

(R = 100 at u = 2, converged per C175; R = 400 at u = 1, C175's converged value 6.0e-4 reproduced.) The trimer stays bound everywhere, and B3 falls with the offset, roughly as B3(0)(1 − 0.8δ) at u = 2. [heuristic] The atom binds to the pair partly through virtual b-components that the offset makes more expensive.

### 2.3 Atom–dimer threshold lengths

[numerical, ad_scat_w2.out; the wall convention adds a common offset of order half the dimer size, about 4 sites] At u = 2 the two lowest atom–dimer scattering levels at K = K_pair + k_a give effective lengths of about 9 and 17 to 19 sites at R = 120 and 180, at φ = π/2 and 1.545. At φ = π/2 the two channels (atom left, atom right of the pair) are exchanged by Θ, so the splitting of the two lengths, about 8 sites, is 2|a_o| with a_o the off-diagonal (transmission) length. Atoms are therefore not impenetrable to pairs at threshold in the pseudospin sense: the exchange coupling J_⊥ ∝ n³|a_o|/M is nonzero. The trimer's shallow binding at u = 1 means an atom–dimer length of order 1/√(2μB3) ≈ 95 sites (1/μ = 1/m_a + 1/M_pair = 10.9), so at u = 1 the atom–dimer channel is nearly free for n ≳ 0.01.

## 3. Dilute phase diagram in (δ, μ)

### 3.1 Which species appears first

[derived, from the data of section 2] Measure the chemical potential from the a-valley bottom, μ̃ = μ − ε_a. A species X with N_X particles and energy E_X (relative to N_X ε_a) appears from the vacuum at μ̃ = E_X/N_X. With E_atom = 0, E_pair = δ − E_b and E_aab = δ − E_b − B3, the first species to appear is

- pairs if δ < E_b − 2B3,
- aab trimers if E_b − 2B3 < δ < E_b + B3,
- a atoms if δ > E_b + B3.

(Pair against trimer per particle: (δ − E_b)/2 < (δ − E_b − B3)/3 iff δ < E_b − 2B3. Trimer against atom: (δ − E_b − B3)/3 < 0 iff δ < E_b + B3.) At u = 2 the trimer window is δ in (0.086, 0.139), φ in (1.5507, 1.5584). At u = 1 it is δ in (0.0183, 0.0201), a sliver of width 2e-3 in δ, 5e-4 in φ. At u = 2 the bound tetramer (C175) makes the pair gas self-bound, which moves the pair–trimer boundary up by about twice the extra binding per particle of the pair liquid (at least 0.002 per particle from the tetramer, larger in the bulk). Each of these single-species dilute liquids has c = 1.

### 3.2 Mixtures at small density are not two-component

[derived, with C164's lattice coefficient as input] When two species X and Y are both dilute they are impenetrable at threshold (generic 1D threshold law), and the mixture maps to a pseudospin chain on the particle order with an exchange J_⊥ of order n³|a_o|/M and an Ising part. C164 showed that unequal masses add a ferromagnetic (segregating) Ising coupling of order 0.3 η_XY² n²/M per X–Y interface, which dominates the n³ terms as n → 0, so the mixture phase separates and the grand-canonical transitions between the species of 3.1 are first order. The relevant mass imbalances are large: pair against atom η = 0.38 (u = 1) and 0.45 (u = 2), trimer against pair 0.19 and 0.21. This confirms C175 section 5, but it bounds the obstruction to a small density: C164's criterion n Δa ≳ 0.1 η² (Δa the length difference that drives mixing) gives n of a few times 1e-3 for Δa of order the pair size. Above that the dilute expansion says nothing about miscibility.

[derived] What the dilute expansion does fix: at small n there is no c = 2 phase in any direction of the (δ, μ) plane. This is C175's conclusion with the windows made explicit, and with the caveat that it covers only n ≲ 0.01 (u = 2) and, because the u = 1 pairs and trimers are 21 and about 95 sites across, n ≲ 1/a_+ ≈ 0.05 and n ≲ 0.01 (u = 1).

## 4. Finite density: the partially polarized liquid

### 4.1 Bosonization with imbalance

[derived, from C180's operator table] Take a homogeneous state with both valleys occupied, ρ_a ≠ ρ_b, incommensurate. Its low-energy theory has two bosons (φ_a, θ_a), (φ_b, θ_b) with a general 2 × 2 velocity and Luttinger structure (Σ is broken off π/2, so no a ↔ b relation, and the drift terms ∂φ_+∂θ_− are allowed). Every nonlinear vertex carries momentum P = n_a k_a + n_b k_b + 2π(p_aρ_a + p_bρ_b) (C180 claim 3). In particular:
- inter-valley backscattering cos 2√2φ_− carries 2π(ρ_a − ρ_b) ≠ 0, so the term that opens the pairing gap at π/2 oscillates and is cut off. This is the 1D FFLO mechanism (as in the spin-imbalanced attractive Hubbard and Gaudin–Yang gases).
- the two-pair vertex O₄ of C179 carries 2(k_a + k_b) + (density harmonics), which off π/2 is 2K_pair ≠ 2π at generic densities.
- composite locking: an operator with p_aρ_a + p_bρ_b ∈ ℤ, for example ρ_a = 2ρ_b (the trimer, p = (1, −2)), is resonant on a line of fixed valley ratio and, if relevant, opens a plateau of finite width in δ (as in Burovski, Orso and Jolicoeur's composite liquids of mass-imbalanced Fermi mixtures, cited from memory, not checked here).

So a homogeneous imbalanced two-valley state is c = 2 on an open set of (δ, μ) away from these plateaus and resonance tongues, provided it is thermodynamically stable (positive-definite compressibility matrix) [derived]. The magnetization curve m(δ) = ρ_a − ρ_b at fixed n then has plateaus at m = 0 (the pair liquid, width set by the spin gap Δ_s), possibly at ρ_a = 2ρ_b (trimers), and at m = n (saturation), with c = 2 between them, unless a first-order jump replaces part of it.

### 4.2 Where the window should be

[heuristic] At u = 1, n ≈ 0.17 the pair gap is tiny (C175's Yang–Gaudin estimate 1e-4), the trimer is about 95 sites across, and so both lower plateaus should be invisible at L ≤ 256. The c = 2 window should then run from ϑ of order 1e-5 up to the saturation offset. C175's ring ED (polarized sectors 0.46 to 0.66 above the balanced one for N = 4 to 8 at L = 24, 30) puts δ_sat at roughly 0.1 to 0.25, so φ from about 1.51 up to π/2. This is consistent with C173's open-chain readings, c ≈ 2 at φ = 1.56 and 1.54 and one valley by φ = 1.5. At u = 2 the pair gap is of order E_b = 0.12, so the window should start at φ ≈ 1.553 or below, and its existence is the clean test.

### 4.3 Commensurate total density

[derived] At density n = 1/q the total-density vertex cos(2qφ_+) (dimension of order q²K₊/2 in the convention where K₊ → 1/2 for the dilute two-valley gas) is resonant. For q = 4 it is relevant when K₊ < 1/4, plausible for the compact pairs at u = 2. The u = 2 runs are therefore at n = 9/32 and 3/8 rather than 1/4.

### 4.4 FFLO signature

[derived, Luttinger counting] In the partially polarized phase the pair correlation oscillates at K_pair ± q with q = π(ρ_a − ρ_b) (the difference of the two valley Fermi wavevectors, k_{F,σ} = πρ_σ for a single band with two pockets), the single-particle function is a power law with two valley wavevectors, and the parity gap closes as 1/L. In the pair liquid the pair correlation sits at K_pair and G_f decays exponentially.

## 5. Numerics

### 5.1 Ring magnetization curves (ring_mag.py with C178's kring.py; files ringmag_*.out, summary by ringsum.py)

For each φ, the lowest level in every momentum sector at fixed N, and the ground state's valley imbalance M (fermion n(k) weight on cos k > 0 minus cos k < 0; off π/2 the lower valley has cos k < 0, so M < 0 means excess a). M from n(k) is not quantized and under-reads in paired states (one conversion reads |M| ≈ 1.8, two read ≈ 2.6 to 3.7).

[numerical, ringmag_w2_L28.out, ringmag_w2_L32.out] |w| = 2, n ≈ 0.28. L = 28, N = 8: balanced (|M| < 0.01) from π/2 down to φ = 1.53 (δ = 0.28), then one conversion (|M| = 1.78) from 1.52 to 1.48 (δ = 0.35 to 0.63). L = 32, N = 9: minimal (|M| ≈ 0.8) down to 1.53, then |M| = 1.04 at 1.52 and 2.6 at 1.50 and 1.48.
[numerical, ringmag_w2_L36.out, partial at write-up, φ = π/2 to 1.555] L = 36, N = 10: balanced down to 1.56, one conversion (|M| = 1.78) already at 1.555 (δ = 0.11), so the first conversion moves from δ ≈ 0.3 (L = 28, 32) to δ ≈ 0.1 (L = 36). The ring file continues to be written by a local run after write-up.
[numerical, ringmag_w2V0_L32.out] Control V₂ = V₃ = 0 (two-component at π/2 by C175), L = 32, N = 9: minimal |M| ≈ 0.9 down to 1.53, then |M| = 2.7 at 1.52 and 3.5 at 1.50 and 1.48, with several imbalance sectors within 0.004 of each other at 1.52.
[numerical, ringmag_w1_L30.out, ringmag_w1_L36.out] |w| = 1, n = 1/6. L = 36, N = 6: balanced down to φ = 1.54 (δ = 0.12), one conversion (|M| = 1.85) from 1.53 to 1.50 (δ = 0.16 to 0.27), and the next sector (|M| = 3.7) 0.057 above at 1.50. L = 30, N = 5: minimal imbalance down to 1.50, with the next sector 0.018 above there.

Reading [numerical, against R1 to R4]. R1 passes. R2's staircase criterion is met formally at L = 28 and 32 (an intermediate |M| ground state over φ = 1.48 to 1.52, with no jump to full polarization), and so is R4 for the control. But the control shows that the test does not discriminate at these sizes: the first conversion sits at δ ≈ 0.3 both with and without the attractions, so it is set by the finite-size cost of a conversion (of order the level spacing of 9 particles), not by a pair gap. A pair gap of order E_b = 0.12 cannot be resolved this way at L ≤ 32. R3 fails as registered: at |w| = 1 the first conversion comes only at φ ≈ 1.535 (L = 36), and the ring is not saturated at 1.50 (two of six particles converted). At L = 36 the steps come at δ ≈ 0.12 and 0.3, with no jump to full polarization. [heuristic] The rings rule out a metamagnetic jump from balanced to fully polarized at both |w|, which is one of the ways the partially polarized phase could have failed.

### 5.2 Local open-chain previews (local/dmrg/, summarize.py)

[numerical, local/dmrg/loc_w2_*.jsonl] |w| = 2, N/L = 9/32. The open chains are not uniform. The bulk density is 0.345, against N/L = 0.281, with depleted edges (0.13 to 0.2 in the outer eighths), at φ = π/2 and at 1.50 alike, so the liquid is self-bound near n ≈ 0.35 and the fixed-N chain squeezes it into the middle. In the bulk the entanglement profile is flat (S = 1.58 to 1.61 over the middle 40 sites at L = 64), and the fermion n(k) is a single broad hump between the valleys (centred at −1.6 at π/2, at −1.9 at φ = 1.50), the signature of tightly bound composites rather than two valley Fermi seas. Two-size estimates between L = 32, 64, 96 at φ = 1.50 are dominated by the edge regions (values near 4, with c_prof between 0.9 and 1.8) and are not meaningful. The |w| = 2 point at this density is a dense clustered liquid, and the Typhon runs at n = 9/32 measure that liquid at its own density with depleted edges, not a dilute pair–atom mixture.

[numerical, same files] Correlations at |w| = 2, L = 64, N = 18 (bulk n = 0.345). The boson and fermion one-body functions fall to the noise floor (about 5e-3) within 8 sites, and the pair functions D₂, D₃ are of order 1e-4 beyond 2 sites, at φ = π/2 and at 1.50 alike. So neither single particles nor pairs are coherent in the bulk, which points to larger clusters. On the ring at L = 32, N = 9 the ground state for φ ≤ 1.50 has the trimer composition (6, 3) (section 5.1, |M| ≈ 2.6 of a possible 3), and at L = 28, N = 8 it has (5, 3), the composition closest to 2:1. [heuristic] At |w| = 2 the off-symmetric state at this density is a cluster liquid, a trimer-rich one at φ ≤ 1.52, so the composite locking of section 4.1 (alternative M4(c)) pre-empts the partially polarized two-component liquid at the sizes reached. This is indicative only, since the open-chain p from n(k) (0.11) does not resolve the composition when every one-body function is short-ranged.

[numerical, local/dmrg/loc_w1_p1.5550_*.jsonl] |w| = 1, φ = 1.555 (δ = 0.061), n = 3/16, L = 64 (N = 11, 12, 13, truncation 4e-8 at χ = 256) and L = 96 (N = 18, χ = 128 only, truncation 2.5e-6, not converged). The bulk density is 0.223 (edges 0.10). The boson one-body function stays at 0.02 to 0.08 out to r = 30 with a two-valley beat, the boson and fermion n(k) have two equal humps (p = 0.009 from fermions, 0.006 from bosons), so the chain is still valley-balanced at this offset and size, consistent with the ring, where the first conversion at L = 36 comes at δ ≈ 0.12. The fermion function decays faster (1e-3 by r = 18), which [heuristic] for hard-core particles in two valleys reflects the Jordan–Wigner string over the total density, not a single-particle gap. The parity gap E(13) + E(11) − 2E(12) = −0.0035, so L·Δ_P = −0.22: no parity gap at L = 64. The pair function D₂ decays slowly (0.04 at r = 0, 0.005 to 0.01 out to r = 26) at momentum π. The preview c_two(64 → 96) is 1.20, 1.19, 1.37 for W = 1, 6, 12, with L = 96 unconverged. Grand-canonical μ_eff ≈ [E(13) − E(11)]/2 = −2.143.

[numerical, local/dmrg/loc_w1_p1.5300_*, loc_w1_p1.5000_*, locctl*] |w| = 1, n = 3/16 off the balanced point.
- Valley imbalance (boson n(k) weights, lower valley share W_a): L = 64: 0.503, 0.557, 0.617 at φ = 1.555, 1.53, 1.50. L = 96: 0.537 (p = 0.075) at 1.555 and 0.578 (p = 0.157) at 1.53. The imbalance grows continuously with the offset and with L at fixed offset (the finite-size conversion cost falls), and no chain jumps to one valley. The density stays homogeneous in the bulk (variation below 3 percent).
- Correlations at φ = 1.53, L = 96. The boson one-body function stays at 0.02 to 0.08 out to r = 44, beating between the two valley wavevectors. D₂ decays from 0.04 to about 0.005 at r = 40. The dominant pair momentum is within the Fourier resolution (2π/48) of π, so the FFLO shift q = π(ρ_a − ρ_b) ≈ 0.08 is not resolved (M6 untestable at these sizes).
- Entanglement. Fixed-window two-size estimates at φ = 1.53 are c_two(64 → 96) = 0.81, 0.99, 1.15 (W = 1, 6, 12). At φ = 1.555, where L = 64 is still balanced and L = 96 is not, they are 1.56, 1.54, 1.71. The controls at the same sizes: c = 1 control 1.07, 1.01, 1.06 (c_prof 0.93, 0.95), c = 2 control 3.46, 1.54, 2.33 (c_prof 1.87, 1.99). The fixed-window estimator is thus calibrated for c = 1 but not for c = 2 at these sizes. A window proportional to L (cuts with |l/L − 1/2| ≤ 1/4, relwin.py, chosen after seeing the controls) reads 1.00 (c = 1 control), 1.79 (c = 2 control), 1.56 (φ = 1.555) and 1.06 (φ = 1.53). The profile fits of the cusp chains (2.1 to 3.3) are inflated by the depleted edges and are not used.
- Reading [numerical, against the preview rule registered in RESULT.md]: φ = 1.53 fails the preview (c_two 0.99 and 1.15 for W = 6 and 12), and the preview is uncalibrated in any case. The calibrated relative-window estimator also reads c ≈ 1 at φ = 1.53.
- [heuristic] Caveat that limits this reading. In the partially polarized liquid the second mode is the gas of excess a particles (pseudospin solitons), at density ρ_a − ρ_b ≈ 0.03 here. The chains hold only one to three of them (N_a − N_b = pN ≈ 1.4 at L = 64 and 2.8 at L = 96 by the n(k) weights), so their entanglement contribution is not yet in the logarithmic regime, and a c = 2 phase would read close to 1 at these sizes. The test needs L(ρ_a − ρ_b) ≫ 1, which favours larger offsets (φ = 1.50 to 1.52, where p ≈ 0.2 to 0.25) and L ≥ 128 to 256.

## 7. Status of the registered test

Typhon job 17604471 (resubmitted from 17604449 with a smaller memory request) had not started at write-up, with an estimated start of 2026-11-13. RESULT.md lists its outputs and the collection command. Its |w| = 1 points at φ = 1.50 and 1.54 (p ≈ 0.2) at L = 128 are the most informative given the caveat above.

## 6. Qubit bond and the C179 vertex list off J_x + J_y = 0

### 6.1 The bond

[derived, bond.py] With w = J_x + J_y + 2iD, Δ = J_x − J_y, V = 4J_z, μ = 4J_z + 2h_z = |t₂|(μ_eff + 6) and |t₂| = Δ²/(2V) (C173 PROOF 0.1, 3.1; the collaborator's validation note uses the same μ map), the effective point (u, φ, μ_eff) at V/Δ = r and Δ = 1 is

    J_x + J_y = (u/2r) cos φ,  D = (u/4r) sin φ,  J_x − J_y = 1,  J_z = r/4,  h_z = [(μ_eff + 6)/(2r) − r]/2.

bond.py reproduces C173's candidate at φ = π/2 to all digits. Off the stratum the only change is J_x + J_y = (u/2r) cos φ ≠ 0.

### 6.2 Vertices off the stratum

[derived] Off J_x + J_y = 0 the symmetry U of C179 (Σ of C180) is lost. Going through C179's list:
- Inter-valley backscattering cos 2√2φ_−: momentum 2π(ρ_a − ρ_b), nonzero in any imbalanced state. Cut off.
- The two-pair vertex O₄ = (ψ_a†ψ_b†)² (order Δ⁴/V³ in the qubit chain, relevant on the stratum by C179's K_ρ): momentum 2(k_a + k_b) − 2π. From (1.1) both minima shift by −(u²/2)m ϑ to first order, so k_a + k_b − π = −u²ϑ/(8 − u²/2) and the detuning is q₄ = 2u²ϑ/(8 − u²/2). This is small: 0.0042 at u = 1, φ = 1.555 and 0.0083 at φ = 1.54 (checked numerically, the pair minimum K_p of fewbody.py gives 2K_p − 2π = −0.0050 and −0.0097), and 0.048 at u = 2, φ = 1.535. O₄ is cut off only if its locking length ξ₄ satisfies q₄ξ₄ ≫ 1. With C179's estimate of a few thousand sites at V/Δ = 40 this holds by a factor of 10 to 30 at u = 1, φ ≤ 1.555. At V/Δ = 8 (ξ₄ of order 100) it fails at u = 1. So off-stratum protection at u = 1 needs large V/Δ or larger ϑ.
- ΔN = 2 pair creation at separation 2, amplitude 2Δ(J_x + J_y)/V, which is nonzero off the stratum, 7.9e-4|t₂| at u = 1, φ = 1.555, V/Δ = 40, and at separation 3, 2Δ³/(3V²) = |t₂|/30 (collaborator, cusp_validation). Their undressed momentum is k_a + k_b ≈ π (or 2k_a, 2k_b, all far from 0). Dressed by density harmonics they are resonant only on lines p_aρ_a + p_bρ_b = 1/2 of the (ρ_a, ρ_b) plane, which cut tongues.
- The valley offset itself is the field δ, and the mass imbalance η = δ/(16 − u²) is marginal and harmless.

So at generic (δ, n) in the partially polarized region no vertex is resonant, and a c = 2 phase there, if it exists, is open in B₂, subject to q₄ξ₄ ≫ 1 [derived, conditional on the phase].
