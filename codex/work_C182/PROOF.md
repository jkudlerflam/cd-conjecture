# C182 PROOF: cluster bands and two-component cluster liquids of the Fibonacci-cusp model

This run was done by a Claude Opus 5.5 agent (2026-11-09), not by Codex (Astra out of credits).

Model and conventions. C173 PROOF eq. 3.2 in the c frame, units |t₂| = 1: hard-core particles, no two adjacent, H = Σ_j [w c†_{j+1}c_j + h.c. − (c†_{j+2}c_j + h.c.) − 2 n_j n_{j+2} − (2/3) n_j n_{j+3}], w = |w|e^{iφ}, ε(k) = 2|w| cos(k − φ) − 2 cos 2k. Particles never pass, so statistics do not matter. For 0 < |w| < 4 and φ near π/2 the band has two minima, valley a (the lower one, energy ε_a) and valley b (ε_b = ε_a + δ, δ ≥ 0 the valley offset). At φ = π/2, k_a = −arcsin(|w|/4), k_b = π − k_a, ε_min = −2 − |w|²/4.

## 1. Cluster bands and the valley symmetry

[proved] An N-particle cluster state is ψ(x₁,…,x_N) = e^{iKx₁} f(g₁,…,g_{N−1}) with gaps g_i = x_{i+1} − x_i ≥ 2. Its band E_N(K) is the bottom of the relative operator at total momentum K. At φ = π/2 the valley exchange Σ (C180 claim 1, C173 PROOF 3.3) acts on single-particle momenta as k → π − k, so on N-particle momenta as K → Nπ − K. Hence E_N(K) = E_N(Nπ − K).
- N even, balanced valley content (n_a = n_b): the cluster sits at the Σ-invariant momentum K = Nπ/2 (mod 2π), and its band minimum there is unique unless its curvature is negative.
- N odd, or unbalanced content: the clusters (n_a, n_b) and (n_b, n_a) sit at K and Nπ − K and are degenerate. Every trimer band (aab, abb) and every pentamer band therefore has two degenerate minima at φ = π/2. This is the "two-valley question one level up" in its unavoidable form.

Dirichlet bound [proved, min-max]. Truncating every gap at R is a restriction to a subspace, so the box eigenvalues are upper bounds on the infinite-line spectrum bottom at that K. A box energy below the scattering threshold therefore certifies a bound cluster, and per-particle energies from boxes are upper bounds.

## 2. The pair band

### 2.1 Relative problem

[proved] For two particles, r = x₂ − x₁ ≥ 2, (h_K f)(r) = Σ_s c_s f(r + s) + V(r) f(r), c_d = a_d e^{−iKd} + a_{−d}, a_{±1} = w, w*, a_{±2} = −1, V(2) = −2, V(3) = −2/3. The truncated operator is a Hermitian pentadiagonal matrix, diagonalized by banded LAPACK (pairband.py). At φ = π/2 and K = π it reproduces the collaborator's closed-form E_pair = −2[p + 1/p + u²p/(1 + p)²] to all printed digits at u = 1, 2 (PROGRESS).

### 2.2 The ab pair band at φ = π/2 has one minimum, at K = π, for every 0 < |w| < 4

[proved algebra, numerical evaluation on a dense grid] The collaborator's eq. (8), 1/m_pair = −(u/2)C₁ − 2C₂ with the exact wavefunction (DERIVATION.md eqs 5 to 7), evaluated at 4000 points of u ∈ (0, 4), is positive and monotonically decreasing, from 4 at u → 0 to 0.643 at u → 4 (pairmass_exact.py, pairmass_exact.out). Values: 3.374 (u = 1), 2.251 (2), 1.296 (3), 0.931 (3.5), 0.695 (3.9). The banded numerics give the same curvature (0.93138 at u = 3.5, 0.69469 at u = 3.9).
- Consistency check. As u → 0 the bound pair becomes two free particles, and 1/m_pair → (8 − u²/2)/2 = 4. Near u = 4 the single-particle valleys merge into a quartic minimum, so the free-pair value (8 − u²/2)/2 → 0, but the bound pair keeps a finite mass. Binding lowers the pair mass at small u and raises it near u = 4, crossing at u ≈ 3.5.

[numerical] The full K scan (361 points plus golden refinement, R = 2000 to 4000, R → 2R changes below 2e-15 at the bound minimum; sym_pairs.out) finds, at every |w| ∈ {0.3, 0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.25, 3.5, 3.6, 3.7, 3.8, 3.9}, a single bound minimum at K = π. The other local minima of the two-body lower edge are the like-valley thresholds K = 2k_a, 2k_b at 2ε_min, above it by E_b.

So the pair band never has two minima at the symmetric point. A liquid of ab pairs is a single-species liquid, c = 1 if it is gapless.

### 2.3 Like pairs

[numerical] Same-valley pairs (K = 2k_a) are unbound for |w| ≤ 3.6 (box energy above threshold by O(1/R²)), and bound at 3.7 (by 0.0022) and 3.8 (by 0.0125) (sym_pairs.out), in line with C175's like-pair resonance near |w| ≈ 3.63. There the aa and bb pairs form a pair band with two degenerate minima, but it lies above the ab pair by E_b − B_aa = 0.131 (3.7) and 0.100 (3.8). The dilute pair liquid is still the ab one.

## 3. Onset species and the pre-emption lemma

[derived, standard] In the grand-canonical ensemble the density first becomes nonzero at μ₀ = inf_N E_N^min/N, and the dilute phase is a gas of the cluster that attains it. If the infimum is approached only as N → ∞, the onset is first order, into a liquid at finite density.

[proved, elementary] Pre-emption lemma. Let clusters X and Y (N_X, N_Y particles) have equal per-particle energy e at some parameter point. If X and Y bind, E_XY = E_X + E_Y − B with B > 0, then E_XY/(N_X + N_Y) = e − B/(N_X + N_Y) < e. So the point where X and Y coexist at the dilute onset belongs to the window of XY or of a larger cluster. A dilute two-species liquid at the onset needs two clusters on the lower envelope that do not bind.

Valley offset as a field. [derived, input from C175 phi_scan and C181 fewbody] Write E(m, n) = mε_a + nε_b − B(m, n) for a cluster with m a-particles and n b-particles. The binding B changes slowly with φ (the ab binding by 0.1 percent from φ = π/2 to 1.54 at |w| = 2, the trimer binding by about 15 percent over δ = 0 to 0.14 there). Then the per-particle energy relative to ε_a is e(m, n) = [nδ − B(m, n)]/(m + n), a straight line in δ with slope n/(m + n). The onset species is the lower envelope of these lines, which is the "Zeeman field" picture of C179 claim 10 seen at the onset. Pure a atoms (e = 0) win for δ > δ* = sup_{m,n≥1} B(m, n)/n.

Consequence for pairs and atoms [derived]. Pairs and atoms have equal onset where δ = E_b. The aab trimer has B(2, 1) = E_b + B₃ with B₃ > 0 (C175, C181 F2), so by the lemma it pre-empts that point. In the dilute limit the trimer is the onset species for E_b − 2B₃ < δ < E_b + B₃, and the would-be pair–atom coexistence is replaced by two lines, pair/trimer at δ = E_b − 2B₃ and trimer/atom at δ = E_b + B₃, each subject to the lemma again (pentamer aaabb and tetramer aaab).

## 4. The two-valley trimer band, one level up

[proved, §1] At φ = π/2 the trimer band always has two degenerate minima (aab at K₃ and abb at 3π − K₃). So a dilute trimer gas would pose C175's question again, with the trimer as the particle.

[derived] It never arises at the onset, and its unlike channel is reactive.
- Onset. Trimers beat pairs per particle iff E₃/3 < E₂/2, i.e. iff B₃ > E_b/2, with B₃ = E₂ + ε_min − E₃ the trimer binding against pair + atom. The computed B₃ is far below E_b/2 wherever it is resolved (§5), so trimers are never the onset species at φ = π/2.
- Reactive unlike channel. An aab and an abb trimer have total momentum K₃ + (3π − K₃) = 3π ≡ π, the same as three ab pairs, and the same valley content. Three far-separated pairs at K = π have energy 3E₂ (an infimum over the six-body spectrum, so E₆(π) ≤ 3E₂). Then 3E₂ < 2E₃ iff B₃ < E_b/2, which is the same inequality as above. So whenever trimers are not the onset species, two unlike trimers can rearrange into three pairs with energy release E_b − 2B₃ > 0. The C175 pseudospin construction for a trimer gas (impenetrable unlike scattering with a 2 × 2 length matrix) has no stable threshold in the unlike channel. The like channel (aab + aab, two ab bonds and two spare a's) has no such decay.
- Same structure for every odd cluster. An unlike pair of odd clusters at φ = π/2 is valley-balanced and has the momentum of a set of ab pairs, so it is open to rearrangement into pairs whenever pairs are lower per particle.

## 5. The cluster ladder and the dilute onset

Per-particle energies e_N = E_N^min/N − ε_a (box upper bounds, R in parentheses; local/ladder_w3_summary.out, local/cl_w3_*_small.json, local/l6_small.jsonl; |w| = 1, 2 numbers from C175 and C181 as cited).

| |w|, δ | e₂ | e₃ | e₄ | e₅ | e₆ | onset |
| :-- | :-- | :-- | :-- | :-- | :-- | :-- |
| 1, 0 (C175) | −0.0098 | −0.0067 | ≤ −0.0088 (R 80, unresolved) | | | pairs (tetramer open) |
| 2, 0 (C175) | −0.0613 | −0.0469 | −0.0635 (R 60) | −0.0570 (R 36) | ≤ −0.0600 (R 14, box too small) | tetramers or larger |
| 3, 0 | −0.0975 (exact) | −0.0877 (R 40 = R 60) | −0.1107 (R 26) | −0.1087 (R 18) | −0.1194 (R 14) | collapse |
| 3, 0.75 E_b | −0.0247 | −0.0324 | −0.0388 | −0.0450 | | collapse |

[numerical] At |w| = 3 the hexamer upper bound at K = π already lies below the tetramer, and at δ = 0.75 E_b the ladder falls monotonically through N = 5. [heuristic, from the trend and the lemma of §3] The dilute gas at |w| ≳ 2 is self-bound: the onset is first order, into a liquid at finite density, and there is no dilute regime in which a two-species question could be posed. The trimer pair binding B₃ at |w| = 3 is 0.068 against E_b/2 = 0.0975, so trimers are not the onset species at φ = π/2 (§4).

Consequence [derived from §3 and the table]. The would-be pair–atom coexistence at δ = E_b is pre-empted by trimers wherever B₃ > 0 (all |w| checked), and at |w| ≥ 3 by ever larger clusters. No line in (|w|, δ) found so far carries two coexisting dilute species that do not bind. The Typhon ladder (job 17604443) completes this table for |w| = 1 to 3.5 and δ/E_b = 0 to 2 with N ≤ 5 (N ≤ 4 at |w| ≤ 1.5) and the hexamers at R = 16, 18.

## 6. Many-body phase diagram, local runs (C176's fibdmrg.py, open chains at fixed N)

All local runs are small (L ≤ 120, χ ≤ 192) and serve as cheap diagnostics. The full screen is the Typhon job (§7).

### 6.1 Equation of state at φ = π/2 (local/eos/eos.out, L = 40, N = 1 to 20)

[numerical] The chemical potential on even N, μ₂(N) = [E(N) − E(N − 2)]/2, and the energy per particle e/N:
- |w| = 1: μ₂ rises monotonically from −2.221 (N = 4) to −1.914 (N = 16). No coexistence up to n = 0.5, as C176's L = 48 scan found. The C173/C176 candidate region is a uniform phase.
- |w| = 0.5: μ₂ is flat to 0.003 between n = 0.2 and 0.3 (−1.953, −1.950, −1.952, −1.942) and then rises. At N/L = 0.35 the ground state is a dense droplet with mid-chain density 0.47 (local/b4, the same energy as the EoS run to 1e-9). So the small-|w| model has a first-order (coexistence) region between a dilute pair gas near n ≈ 0.2 and a dense phase near n ≈ 0.47.
- |w| = 2: e/N − ε_min is lowest at N = 4 to 6 (−0.0587) and rises beyond, with μ₂ nearly flat (−3.061, −3.059, −3.049, −3.033). Clusters of a few particles are self-bound, and open chains at N/L = 1/6 to 0.3 form droplets of mid-chain density 0.23 to 0.35 (local/mini, local/b2, local/b3).
- |w| = 3: e/N falls monotonically through N = 8 (−0.095, −0.110, −0.119, −0.123 at N = 2, 4, 6, 8), μ₂ falls. The dilute gas collapses, and open chains at N/L = 1/4 settle at n_mid = 0.29.

[heuristic] Mechanism. At w = 0 each leg of the zigzag ladder is the ferromagnetic Heisenberg chain (C173 claim 4), the boundary of phase separation. Any further attraction tips it. The inter-leg V₃ does so when |w| is small, the kinetic energy of w stabilizes a uniform gas near |w| = 1, and the w-assisted ab binding drives clustering again at |w| ≳ 2. This is consistent with the cluster ladder of §5.

### 6.2 Uniform liquids above the coexistence region (local/b2, local/b4)

[numerical] Where the chain is homogeneous (n_mid within 0.005 of N/L), the paired or clustered liquids read one mode:
| point | n | c(40 → 80), W = 1 / W = 5 | c_prof(80) | valley ratio |
| :-- | :-- | :-- | :-- | :-- |
| |w| = 3, π/2 | 0.3 | 0.87 / 0.85 | 0.60 | 0.98 |
| |w| = 3, π/2 | 0.4 | 0.87 / 0.96 | 0.92 | 0.98 |
| |w| = 3, δ = E_b | 0.3 | 0.87 / 0.85 | 0.61 | 0.92 |
| |w| = 2, π/2 | 0.4 | 0.85 / 0.95 | 0.89 | 0.98 |
| |w| = 2, δ = E_b | 0.4 | 0.85 / 0.95 | 0.89 | 0.99 |
| c1 control | 0.2 | 1.10 / 1.01 | 0.93 | |
| c2 control | 0.2 | 2.17 / 1.80 | 1.87 | |

The truncation is below 1e-10 in all five model rows. Both valleys are equally occupied (valley ratio near 1, even at δ = E_b), so these are balanced paired liquids with the valley sector gapped, c = 1.

### 6.3 Inhomogeneous readings (not usable for c)

[numerical] At |w| ≤ 1 and at |w| = 2 inside the coexistence or droplet regime, the two-size estimates are high (2.4 to 3.9) but the mid-chain density differs from N/L by 0.03 to 0.12 and drifts with L, so they fail the homogeneity condition registered in Q10 and Q11. The |w| = 2, δ = E_b lead of local/mini (Q10) falls in this class: at n = 1/4, n_mid = 0.328, 0.313, 0.297 at L = 40, 80, 120, inside the |w| = 2 droplet regime of 6.1. By the registered rule it is a droplet effect, not a lead. At |w| = 1, n = 0.35 the chain is uniform in the EoS sense (μ monotonic) but its edges are depleted over 3 to 4 sites, so the bulk density differs between L = 40 and 80 (0.414, 0.388). Its profile fit c_prof(80) = 2.0 sits in the BCS-crossover regime of the C173 candidate (pair size about 21 sites, n·a₊ ≈ 7), where a pair gap below the finite-size spacing reads as c ≈ 2.

## 7. The Typhon screen (pending at the time of writing)

Job 17604443 runs (i) the cluster ladder N ≤ 5 (N ≤ 4 at |w| ≤ 1.5) at |w| = 1 to 3.5 and δ/E_b = 0 to 2, with the trimer box checked at 1.5R, which completes the table of §5 and the lemma's bookkeeping off φ = π/2 (pair/trimer, trimer/atom, and pentamer, tetramer channels; summarize.py prints B₃, B(2+3) = E₂ + E₃ − E₅ and the aaab binding); (ii) the balanced hexamer at |w| = 2 to 3.5 with R = 16, 18; (iii) a DMRG screen at L = 48, 96 on the grid of Q7, with the parity gap at L = 96; (iv) stage 2 at the best six points by the registered rule. Given §6.1, screen points at N/L below the coexistence or droplet density will be inhomogeneous, and their c readings should be discarded by the homogeneity condition used in §6.3 (n_mid within 0.02 of N/L).
