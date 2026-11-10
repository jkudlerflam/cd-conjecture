# C188 derivations: the SU(2)_1/ordered boundary of the cyclic stratum

Claude Opus 5.5, agent C188. Tags: [proved], [derived], [heuristic], [numerical], [assumed].

## 1. Setting and exact facts

h = α XZ + β YX + γ ZY on a ring or open chain. C154 Prop. 2 [proved there]: the on-site rotations σ_j = P^{-j} τ_j map it to the period-3 XYZ chain

    H = Σ_j Σ_a J^a_j τ^a_j τ^a_{j+1},   (J^x, J^y, J^z)_j = (α, β, γ), (γ, α, β), (β, γ, α) for j ≡ 0, 1, 2 (mod 3),   (1.1)

with a twist P at one bond if 3 ∤ L. Symmetry: A_4 (twisted Klein group and P), acting on (τ^x, τ^y, τ^z) by its real 3-dimensional irrep, time reversal Θ, and σ-translation, which is τ-translation composed with P.

Lemma 1 [proved; checked numerically at L = 9, 10 to 7e-14, see PROGRESS]. A uniform π rotation about σ^x maps (α, β, γ) to (−α, −β, γ). Hence

    spec H(1, β, −γ) = − spec H(1, β, γ).   (1.2)

So the "ferromagnetic" half of the fundamental triangle (γ < 0) is the negative of the antiferromagnetic half, and its ground state is the top state of the AF chain. In the τ frame of (1, β, −|γ|) every bond can be gauged to the purely ferromagnetic form −(J^x τ^xτ^x + J^y τ^yτ^y + J^z τ^zτ^z) with all J^a ≥ 0 and the strong axis rotating from bond to bond (x, y, z with weights (1, β, |γ|) cyclically permuted). The ferromagnetic Heisenberg point is (1, 1, −1).

Parametrization used throughout: β ∈ (0, 1], t = −γ/β ∈ [0, 1]. The edges of this half-triangle are t = 0 (the YX line (1, β, 0), which is self-dual under (1.2)), t = 1 (the XZ–FM edge β = −γ) and β = 1 (the AF–FM edge through (1, 1, 0)).

Classical energy [derived]. In the ferromagnetic gauge a uniform product state of direction n has energy per bond −(1 + β + tβ)|n|²/3, independent of the direction, because Σ over one period of the cyclically permuted diagonal couplings is (1 + β + tβ) times the identity. So the classical ground manifold is the whole sphere, and the direction (cube axes, giving 6 states, or body diagonals, giving 8) is selected by fluctuations. The A_4 invariants of a vector n are |n|², n_x n_y n_z (odd under Θ) and Σ n_a⁴, so the selection is by the quartic cubic anisotropy, which is a fluctuation effect of order (1 − β)² and (β − tβ)² at small anisotropy.

## 2. Field theory of the boundary

### 2.1 What can sit next to SU(2)_1 in this stratum

[derived, input C154 §3.3] At an SU(2)_1 point of the stratum the infrared fields are g (j = 1/2, dimension 1/2) and the currents J^a, J̄^a. σ-translation acts as the centre g → −g times P. The allowed perturbations of dimension ≤ 2 with spin 0 reduce to the isotropic marginal coupling λ J·J̄:
- tr g (dimerization) is odd under the centre;
- N^a = tr(g σ^a) (Néel) is a 3 of A_4, which has no invariant;
- the symmetric traceless part of J^a J̄^b is the 5 of SO(3), which restricts to 1' + 1'' + 3 of A_4 and has no invariant;
- chiral currents and spin-2 operators have no invariant.
The first A_4-invariant anisotropy has dimension 4 (for example Σ_a (J^a J̄^a)², or the cubic part of N⁴), irrelevant.

So inside the SU(2)_1 basin the only flowing coupling is λ. For λ > 0 (the Heisenberg sign) it is marginally irrelevant, for λ < 0 the chain dimerizes. The isotropic line λ_∥ = λ_⊥ is the separatrix of the Kosterlitz–Thouless flow of an XXZ-type chain, between the critical fan and the Néel-massive region. With anisotropic marginal couplings forbidden by A_4, the flow cannot leave the separatrix, so Néel order is not reached perturbatively from SU(2)_1 however large λ > 0 is [derived]. On finite rings a large bare λ > 0 shows up as strong logarithmic corrections, x_t = 1/2 − λ/4 + …, x_s = 1/2 + 3λ/4 + … (Affleck, Gepner, Schulz, Ziman 1989), that is r_s = x_s/x_t large and the triplet low. That mimics a soft triplet at small L without any order, and must be separated from true order by the L dependence (x_t rising back toward 1/2·2πv in the marginal case, falling to zero faster than 1/L in the ordered case).

### 2.2 The four ways to leave the SU(2)_1 phase

[derived] Moving inside the 2-dimensional stratum, the SU(2)_1 phase can end only in one of four ways.
1. λ changes sign. The transition is the c = 1 line λ = 0 (SU(2)_1 with no marginal coupling, the analogue of J_2c of the J_1–J_2 chain), and the other side is dimerized (c = 0). C154 found this near the XZ corner.
2. The velocity vanishes. The uniform susceptibility of the τ magnetization is χ = 1/(2πv) at the SU(2) fixed point in spin-1/2 units (Eggert, Affleck, Takahashi: χ = 1/(π²J) with v = πJ/2), so a ferromagnetic instability with the emergent SU(2) intact is v → 0, a point with z ≠ 1 (as at the XXZ ferromagnetic point or the J_1–J_2 ferro/singlet point). Not a CFT, no c.
3. A level crossing with a state outside the basin. First order, with a jump in ∂e_0/∂(parameter). The ordered state (ferromagnetic in the ferromagnetic gauge, or Néel-like) is a distinct, gapped, A_4- and Θ-breaking state, so a crossing is the natural mechanism.
4. A continuous transition with z = 1. Then the transition line is a CFT C* with exactly one relevant symmetric scalar (one tuning in a 2-dimensional stratum), and SU(2)_1 is the infrared end of C* perturbed by it. This is a nontrivial massless flow, so Zamolodchikov's theorem gives c* > 1 strictly. Conversely, any c > 1 point on the boundary that is not an isolated junction must be of this kind.

So the question "does the boundary carry c > 1" is the same as "is the boundary continuous with z = 1 anywhere along a segment".

### 2.3 Candidates for C*, and why the generic transition is first order

[derived, with one heuristic step marked] C* must (a) carry the LSM anomaly of the projective Klein group with translation, (b) have one relevant A_4 × Θ × translation-invariant scalar, (c) flow to SU(2)_1 for one sign and to the ordered phase (which breaks A_4 and Θ) for the other.
- SU(2)_k, k odd (k even fails (a): the mixed centre/SO(3) anomaly is k mod 2). Scalars tr_j g with dimension 2j(j+1)/(k+2). Centre-even, A_4-invariant ones are the singlet parts of integer j: j = 1 (dimension 4/(k+2), relevant), j = 2 (12/(k+2), relevant for k ≥ 5). So only k = 3 satisfies (b), with c = 9/5. Its relevant singlet is Φ^{aa}_{j=1}. One sign gives the massless flow SU(2)_3 → SU(2)_1 (Affleck–Haldane, spin-3/2 chains), the other a massive phase. That massive flow is SO(3)-diagonal-symmetric, and its vacua break only the centre (dimerization, required by the anomaly). A gapped SO(3)-symmetric phase is stable against the dimension-4 cubic anisotropy, so the massive side is dimerized, not ordered. SU(2)_3 fails (c) and can only sit on an SU(2)_1/dimer boundary.
- Products SU(2)_1 × X with X critical. If X's order parameter σ_X is an A_4 singlet, X ordering breaks nothing that the ordered phase breaks, so X is not the transition to the ordered phase. If σ_X carries the 3 of A_4 at the momenta of N, the coupling σ_X·N is invariant with dimension Δ_σ + 1/2 < 2 for every unitary minimal X, so it is a second relevant invariant: the decoupled point has codimension 2 (a point, not a line), and the generic flow from it locks σ_X and N together, which is direct N order or a first-order jump.
- [heuristic] I know of no other unitary CFT with this symmetry and anomaly content and exactly one relevant invariant whose two sides are SU(2)_1 and an A_4-breaking gapped phase. If none exists, the boundary is first order (or z ≠ 1) along its whole length, apart from isolated junctions with the dimer line, where codimension-2 points such as SU(2)_3 × (tuning of λ) could in principle sit.

So the field-theoretic expectation is case 3 (or case 2 where the velocity vanishes), with c > 1 possible only at isolated points where a dimerized region touches the ordered region. The numerics below test this along the boundary.

## 3. Numerics on rings

Tools [numerical]. ringk.c is a matrix-free Lanczos in a σ-momentum sector with C155's conventions, validated against RingTable at L = 12 to 1e-10 (PROGRESS). Ghosts are removed by the Cullum–Willoughby test (ringk_run.py). L = 12, 18 with all momenta (scan_L12_18*.jsonl, scan_table.out); L = 24 in sectors 0, 1, 4, 8, 12 (local24.jsonl, local24_table.out); L = 30 on Typhon. Momentum classes (absolute σ-momenta, which is what analyze_cross.py uses): S = the SU(2)_1 singlet, at momentum 0 for L ≡ 0 mod 12 and π for L ≡ 6 mod 12 (the τ-frame Heisenberg singlet carries π L/2); N = an A_4 triplet at {π, π ± 2π/3}; F = an A_4 triplet at {0, ± 2π/3}. An ordered state has a fixed real-space pattern, so its absolute momenta do not depend on L: F is the τ'-uniform ferromagnet of the ferromagnetic gauge (period 3 in σ), N a τ'-staggered pattern. SU(2)_1 excitations sit at fixed momenta relative to S instead: the Néel triplet tr(g σ^a) is at relative {π, π ± 2π/3}, which is N for L ≡ 0 mod 12 and F for L ≡ 6 mod 12. So at L = 12, 24 the soft SU(2)_1 triplet shows up in the N sectors and at L = 18, 30 in the F sectors.

### 3.1 Where the boundary is

[numerical, scan_table.out, L = 18] The SU(2)_1 ground state persists up to t_c(β) = 0.56 to 0.57 (β = 1, 0.9), 0.62 (0.8), 0.66 (0.7), 0.72 (0.5), 0.79 (0.4), 0.82 (0.3), 0.87 (0.2), 0.92 (0.1). The curve runs from the XZ corner, where it approaches the XZ–FM edge t = 1, to the AF–FM edge β = 1 at t ≈ 0.565. This is Q1, with t_c(1) slightly above the registered window (0.5 to 0.6, inside it).

### 3.2 The velocity vanishes linearly at the boundary

[numerical] v = x_c/2π from the current level at L = 18 is linear in t on every cut up to the crossing. β = 1: 1.66, 1.37, 1.08, 0.78, 0.48, 0.18 at t = 0.30 to 0.55 (steps 0.05), constant decrement 0.29 to 0.30, extrapolating to zero at t = 0.581. β = 0.7: decrement 0.20 to 0.21, zero at t ≈ 0.68. β = 0.1: 0.357 at t = 0.3 down to 0.006 at t = 0.9. The L = 24 values agree (β = 1, t = 0.5: v = 0.457 against 0.480 at L = 18).

### 3.3 The crossing is not first order, and the gaps close as 1/L²

[numerical, fine_cross_table.out, analyze_cross.py on local24.jsonl] Per-site slope jump |∂(E_S − E_other)/∂t|/L at the crossing t_x:

| β | t_x at L = 12, 18, 24 | jump at L = 12, 18, 24 |
| :-- | :-- | :-- |
| 1.0 | 0.5655, 0.5628, 0.5664 | 0.114, 0.048, 0.015 to 0.024 |
| 0.7 | 0.6605, 0.6593, 0.6635 | 0.082, 0.033, 0.018 |
| 0.4 | 0.7768, 0.7846, 0.7926 | 0.054, 0.017, 0.0035 |

The jump falls at least as fast as L^-2, so the energy density has a continuous first derivative at the boundary in the thermodynamic limit. Q2 registered a jump that does not shrink with L. That fails: the transition is not a first-order level crossing.

[numerical, gap_at_tx.jsonl, local24.jsonl] L² times the current-level gap at the crossing: 14.5, 12.2, 11.5 (β = 1), 12.1, 10.3, 8.3 (β = 0.7), 7.5, 7.3, 5.4 (β = 0.4) at L = 12, 18, 24. A z = 1 critical point would give L × gap constant. These numbers show L² × gap constant or slowly falling, so the dynamic exponent is z ≥ 2 [numerical], matching the linear vanishing of v in 3.2.

Reading [derived plus numerical]. The SU(2)_1/ordered boundary is case 2 of §2.2: the ferromagnetic instability of the SU(2)_1 liquid, where χ ∝ 1/v diverges and the current (magnon) branch becomes quadratic. The current multiplet and the Néel triplet descend together, which is why r_s and r_c grow near the boundary (r_c = 4.2 at β = 1, t = 0.55, L = 24) while c_cas falls. The transition is the analogue of the XXZ chain at Δ = −1 or the J_1–J_2 ferro/singlet point: continuous in the energy, with z = 2 and no conformal invariance. There is no c to assign to it, and the SU(2)_1 side has c = 1 up to the boundary.

### 3.4 No c > 1 tower

[numerical] On the SU(2)_1 side near the boundary the Casimir c falls below 1: c_cas(18, 24) = 0.92, 0.83, 0.43 at β = 1, t = 0.45, 0.50, 0.55; 0.80 at β = 0.7, t = 0.6; 0.81 at β = 0.4, t = 0.7 (Heisenberg control 1.004). At L = 12, 18, c_cas decreases toward the boundary on every cut with β ≥ 0.3 (scan_table.out). The SU(2)_3 signature (r_c near 10/3 with c_cas rising toward 9/5) appears nowhere: where r_c passes through 10/3 (β = 1, between t = 0.53, r_c = 3.0, and t = 0.55, r_c = 4.2, at L = 24) c_cas(18, 24) lies between 0.83 and 0.43 and is falling. This is Q3 [numerical].

[numerical, ring_cS.out, xxz_control.out] The half-ring entanglement slope c_S(12, 18) rises above 1 within about 0.05 of the crossing (1.14 at β = 1, t = 0.55; 1.17 and 1.37 at β = 0.4, t = 0.70, 0.75). The XXZ control, which has c = 1 exactly, does the same as it approaches its ferromagnetic point: c_S(12, 18) = 1.04, 1.11 and c_cas = 0.94, 0.84 at Δ = −0.9, −0.95. So the excess is the known finite-size signature of a c = 1 liquid near a z = 2 ferromagnetic point, with a crossover length that grows as v → 0, and not a second soft sector. The open-chain pairs in §4 test this at larger L.

### 3.5 The XZ corner

[numerical, scan_table.out] For β ≤ 0.2, c_cas(12, 18) exceeds 1 (β = 0.2: 1.25 at t = 0.3 up to 7 at t = 0.85). Along each such cut (c_cas − 1) v is nearly constant (0.17 to 0.20 at β = 0.2 for t = 0.3 to 0.7), so the excess is an additive, t-independent Casimir energy from another sector divided by the vanishing B velocity. This is C154 §3.4's gapped-but-slow A Majorana (velocity 4) of the XZ point, whose Casimir contribution survives until L exceeds its correlation length. At β = 0.1 the L = 18 singlet lies below the triplet (r_s < 1) for t between 0.35 and 0.8, the dimerized pocket of C154 extended into the ferromagnetic half. The L = 24 check of the corner is in §5.

### 3.6 The ordered side

[numerical] With the labels above, every size tells the same story at the boundary: the state that crosses below S is the soft SU(2)_1 Néel triplet (N at L = 12 and 24, F at L = 18), continuously connected to the level whose x_t collapses in 3.2. Deeper in, the ground state becomes the F triplet, the τ'-uniform ferromagnet, at every L (β = 1, t ≥ 0.65 to 0.7: absolute momenta {0, ± 2π/3} at L = 12 and 18). In between there is a narrow window whose ground-state momentum depends on L (L = 12: N for 0.566 < t < 0.617; L = 18: N for 0.585 < t < 0.633; L = 24: N only for 0.566 < t < 0.573, then F, then a momentum-π singlet at t = 0.60). Deep in the ferromagnetic region the F triplet and a second F triplet are split by x ≈ 0.08 to 0.14 at L = 18, consistent with a ferromagnet with cubic anisotropy whose tunnelling splitting falls with L. The intermediate window may be an incommensurate (spiral) sliver next to the z = 2 point, as in frustrated ferromagnetic chains, or a finite-size effect of the anisotropy selection. I did not resolve it, and it does not bear on c at the boundary, which is fixed from the SU(2)_1 side.

### 3.7 L = 30

[numerical, typhon_out2/ring30.jsonl, job 17605308, 300 Lanczos steps, converged to 1e-14] At β = 1, t = 0.5: E_S = −28.418643 (sector 15), the Néel triplet at x_t = 1.130 (sector 10), the current at x_c = 2.795 (sector 14, v = 0.445). Across L = 12, 18, 24, 30: x_t = 1.150, 1.113, 1.122, 1.130 (stable) and c_cas = 0.84, 0.83, 0.90 for (12, 18), (18, 24), (24, 30), creeping up toward 1 as expected for SU(2)_1 with a slow marginal coupling. At t = 0.45: E_S = −29.074229, x_t = 2.051 (2.039, 2.045 at L = 18, 24), c_cas(24, 30) = 0.93 with v(24) (0.92 at (18, 24)). A first attempt with 130 steps (job 17604564, typhon_out/ring30_job1_unconverged.jsonl) missed the singlet entirely. With a random start of overlap ~ dim^{-1/2} = 1.7e-4 and a relative gap ~ 3e-3, 130 steps are below what the Chebyshev bound requires. Those numbers are not used.

## 4. Open chains

[numerical, dmrg_tau.py (real τ-frame MPO, validated against exact ED at L = 10), typhon_out2/dmrg_L*.jsonl, dmrg_local_L48.jsonl, dmrg_local_L192.jsonl; summarize_dmrg.py] c_pair(L1, L2) = 6[S̄(L2) − S̄(L1)]/ln(L2/L1), with S̄ the mean of the four central cuts. Truncation below 1e-8 at every SU(2)_1-side point.

| point | t | c_pair(48, 96) | c_pair(96, 192) |
| :-- | :-- | :-- | :-- |
| Heisenberg (1, 1, 1), control | | 0.916 | 0.934 |
| β = 1 | 0.50 | 0.928 | 0.931 |
| β = 1 | 0.54 | 0.827 | 0.862 |
| β = 1 | 0.56 | 0.52 | |
| β = 0.7 | 0.55 | 0.956 | |
| β = 0.7 | 0.62 | 0.862 | |
| β = 0.4 | 0.65 | 1.002 | |
| β = 0.4 | 0.75 | 0.761 | 0.793 |
| β = 0.7 | 0.62 (again, L = 192) | 0.862 | 0.884 |
| β = 0.2 (corner) | 0.60 | 1.144 | 1.057 |
| β = 0.2 (corner) | 0.80 | 1.322 | 0.939 |
| β = 0.2 (corner) | 0.85 | 1.666 (L = 96 not converged, dE = 5e-6) | |

Away from the corner the SU(2)_1-side points read the same as the Heisenberg control (SU(2)_1 with its logarithmic corrections) or below it as the boundary is approached, never above 1.01, and the (96, 192) values do not exceed the (48, 96) values by more than the control does. At β = 1, t = 0.56, within 0.006 of the ring crossing, c_pair falls to 0.52: the open chain at L = 96 already sees the ordered state. On the ordered side the bulk magnetization in the τ frame is large (|⟨τ⟩| = 0.92 at β = 1, t = 0.7 and 0.84 at β = 0.7, t = 0.8, L = 96) and S does not grow logarithmically. The corner rows (β = 0.2) start above 1.2 at (48, 96) and fall at (96, 192): 1.144 → 1.057 at t = 0.6 and 1.322 → 0.939 at t = 0.8, with the within-chain CC fits falling too (1.26, 1.15, 1.07 and 1.44, 0.97, 0.81 at L = 48, 96, 192). That is the slow A-Majorana crossover of §3.5 and §5, now seen on open chains, and it fails the candidate criterion (it falls with L). The t = 0.85 corner point lies 0.02 from the boundary; its L = 96 run is not converged and its within-chain fit is 0.51, so its pair value carries no information. At the crossings themselves (β = 0.7, t = 0.66; β = 0.4, t = 0.79) the open chains are partly ordered (|⟨τ⟩| = 0.05 to 0.15) and the CC fits are not meaningful. With these caveats this is Q7, passed. The ring c_S excess of §3.4 does not survive on open chains of 48 to 192 sites.

## 5. The XZ corner at L = 24

[numerical, local24_table.out] c_cas falls with L at every corner point: β = 0.2: 1.25 → 1.11 (t = 0.3), 1.52 → 1.17 (t = 0.6), 3.32 → 0.97 (t = 0.8), for (12, 18) → (18, 24); β = 0.1: 1.97 → 1.38 (t = 0.3), 4.13 → 2.13 (t = 0.6). At β = 0.1, t = 0.6, (c_cas − 1)v drops from 0.58 to 0.23, as the A-Majorana Casimir term freezes out. r_s(24) = 0.92 and 0.80 at β = 0.1, t = 0.3, 0.6, falling with L, so the dimerized pocket of C154 (marginally relevant isotropic coupling, c = 0) extends into the ferromagnetic half near the corner. Near the corner the ordered region therefore borders the dimerized pocket, and the λ = 0 line (c = 1) meets the v = 0 line at a junction point close to XZ. The junction sits in the region where v_B → 0 and the A sector is slow. Its c is not defined (v = 0 there), so it is not a c > 1 CFT either, but I did not resolve its structure.

## 6. Summary

[derived plus numerical] The SU(2)_1/ordered boundary of the cyclic stratum is a line of z = 2 ferromagnetic instabilities of the SU(2)_1 liquid. It is neither first order nor a z = 1 multicritical line. The SU(2)_1 side has c = 1 up to the boundary (ring Casimir c, ring spectroscopy and open-chain entanglement all read at or below the Heisenberg control). No SU(2)_3 or other c > 1 tower appears at any size from 12 to 30. The field-theoretic reason is §2: A_4 leaves the isotropic marginal coupling as the only flowing coupling of SU(2)_1, so the liquid can end only by dimerizing (c = 1 line), by v → 0 (this boundary), or by a level crossing. SU(2)_3, the one WZW candidate with a single relevant invariant, borders a dimerized phase and so cannot sit on a boundary with the ferromagnet.
