# Route A assembled: unbounded c at fixed d0, as it stands

## Claim

This is an integration document, not a new counterexample. Route A asks for a one-site translation-invariant, nearest-neighbour, Hermitian bond h_n on ℂ^{d0}, with d0 fixed, such that h_n flows to n + 1 compact bosons (Virasoro c = n + 1). The document writes that chain as one argument, with every component, interface and energy scale, and numbers at n = 2 and 3. No located incompatibility kills route A. The decisive missing piece is still a fixed-alphabet port interface and compiler for the second-order spectral data (G1 and G2 in §7). Two constraints found here restrict how the pieces may be combined (§8).

## Status

Heuristic as a whole. Every finite-dimensional ingredient it cites is proved or numerically checked in its source, the infrared step is standard bosonization input, and the spatial compiler does not exist. New in this folder are the third-order operator-class decomposition at n = 1, 2, 3 [numerical], the filter-conditioning sensitivities [numerical], the observation that every traceless effective term of C19's cell carries ρ² at every order [proved, §4], a leading-order argument that the existing helix compensator also equalizes the second-order bond [heuristic, §4], and the worked hierarchy [arithmetic on stated assumptions].

## How it was established and checked

- Author: Claude (Opus 5.5) agent acting as the C23 integrator, 2026-10-28. No Codex run was used. Astra C21, C22 and C23 ran in parallel. Their PROGRESS files and Astra's own `codex/work_C23/RESULT.md` (snapshot of 2026-10-28 00:07) were read and are cited where they agree or differ. The parallel Claude folders `contrib/routeA_compile/` and `contrib/routeA_spatial/` did not exist at the final read; `contrib/routeA_demonstrator/` had code but no README.
- Checks: `scales.py` recomputes every scale quoted in §3 from the source formulas and asserts agreement with C18's table, C19's ρ² and C20's vertex counts. `third_order.py` rebuilds C19's final cell from `codex/work_C19/uniform_cell.py` (imported read-only) and reproduces C19's non-Gaussian third-order residuals at n = 1, 2 to six digits before extending to n = 3.
- Not independently reviewed.

## 1. The chain in one paragraph

A selector makes the chain break spontaneously into cells of lengths b and b + 1 (C7 walker, `codex/work_C7/C7_result.md`), with the boundaries mobile so that the length pattern is a gapless c = 1 fluid (`codex/work_C8/selector_wrapper.md`, `length_fluid.md`). Each cell holds an initialized history-state core whose exactly degenerate code is the Fock space of n complex fermions, with two qubit ports at its ends whose compressions onto the code vanish. The bond between ports of neighbouring cells, λΣ_a F_a^R ⊗ F_a^L, then acts at second order, and C19's spectral weights make that second-order term an n-flavour free-fermion hopping with independent phases and chemical potentials and exact fermion parity (`codex/work_C19/RESULT.md`). The Fermi points are displaced and incommensurate, so the payload is a stable c = n Gaussian theory (`contrib/ir_control/`), and together with the length fluid the chain is an (n + 1)-boson phase (`codex/work_C13/RESULT.md`, blind check `codex/work_C20/RESULT.md`). Since d0 never depends on n, c is unbounded at fixed d0.

## 2. Parts list and interface ledger

Table 1. Components for the n-flavour target.

| Component | Alphabet | Bond terms | Source | Status |
|:--|:--|:--|:--|:--|
| Logical target H_* | n complex fermions per cell | hopping t_a e^{iθ_a}, fields μ_a (C19 eq. 4) | C19 §2, ir_control | proved kinematics; IR heuristic (verdict b) |
| Fermi-point design | none | Sidon points (δ₂ ≥ π/6n, δ₄ ≥ π/576n²), C13 mixed margin, C20 existence margin 1/(4M(n)) | ir_control, C13 §4.3, C20 §2 | proved margins |
| Cell spectral data | abstract cell, dim 16(n+2)2^n (one-port-at-a-time variant 16(n+2)²2^n) | reflection R H_diag R with k = n+1 poles and Jacobi weights; exact parity Q | C19 §3 to §7; `contrib/graded_route_a/`; `contrib/separated_port_cells/` | proved; ED at n = 1, 2 |
| Port interface | two qubit port letters | λΣ_a F_a^R F_a^L; on-cell checks | C18 §7, §8; minimal_d0 §2 | bare Bell checks fail (proved); diagonal escape proved; hopping escape open |
| History clock | 9 or 15 letters × work space w; two idle port letters: 15w+4 | hopping rules with work gates; conserved-defect penalties; input pins | C8 RUN7 §36 to §37; C18 §2 to §6; `contrib/graded_route_a/sweeping_clock_interface.md` | proved |
| Length selection, frozen | none extra (clock-selected) or C7 {#,L,X,R} | chemical potential on the clock | C18 §5; C7 thm on period b | proved; frozen only |
| Mobile wrapper | 4 × (6q+3) × 4 = 96q + 48 for a q-letter core | C7 walker, arrow wrapper with vacancies, length filter, two helix tracks, endpoint compensators M_L, M_R, marker hop C_new | C8 `selector_wrapper.md`, `arrow_wrapper.md`, `helix_equalizer.md`, `joint_fluid.md` | proved for a positive frustration-free core, first order |
| Infrared | none | none | C13, C20, ir_control, `contrib/research_continuation/IR_AUDIT.md` | heuristic, standard inputs stated |
| Spatial compiler | unknown | unknown | C18 §9 obligations; C21, C22 in progress; `contrib/routeA_demonstrator/` | missing |

Table 2. Interfaces between components.

| Interface | State | Source |
|:--|:--|:--|
| Spectral weights → second-order target with fields via exchange, no 3-cell term at λ² | checked (proved; ED n = 1, 2) | C19 §7, §10 |
| Spectral weights → exact parity of the whole cell chain | checked (proved) | C19 §4, §6; graded_route_a §3 |
| Two ports → separate one-port gadgets on shared memory | checked (proved) | separated_port_cells |
| Cell data → fixed-alphabet port terms | open. Bare positive Bell checks give scalar resolvents (C18 §7); unitary-decorated paths and trees give scalar return measures (C21, C22 PROGRESS); syndrome-conditioned work energies give diagonal ZZ only (C18 §8) | C18, C21, C22 |
| Energy program → history clock (sine weights, Schur error) | checked for timewise-annihilating codes (C18 §6, programmable_cell P1); open for the hopping program | C18 §6 |
| Clock → length selection, frozen | checked | C18 §5 |
| Clock-selected length → mobile cells | incompatible as built (§8) | C8 RUN7 verdict; C18 §11 |
| Core → mobile wrapper | checked at first order for positive frustration-free (FF) cores; second order assumed (heuristic argument in §4); non-FF cores not covered | selector_wrapper eq. 2 to 10; helix_equalizer §2 |
| Wrapper → joint ring spectrum | checked at first order | joint_fluid; C13 §6 |
| Mobile chain → (n+1)-boson IR | assumed (standard bosonization; gluing data assumed) | C13 §3 to §4; C20 §4 |
| Second-order logical chain → c = n at finite λ | assumed; no finite-coupling numerics for a parity-only incommensurate route-A cell | routeA_incommensurate has scripts but no data |

## 3. The scale hierarchy

### 3.1 Scales as functions of n and ℓ

Notation. k = n + 1 poles at E_e = eΔ_src, λ the port coupling, u the marker hopping, J the clock hopping, t the C7 head hopping, ℓ the cell length, Δ_min the smallest excluded local gap. For the cell length no compiler exists, so the table uses the register floor ℓ_min(n) = n + ⌈log₂ 16(n+2)⌉ + 2 (the C19 register at one qubit per site plus two ports), which gives ℓ = 10 and 12. Any working compiler will need a longer cell, which only lowers every ℓ-dependent entry.

Table 3. Values from `scales.py` (J = t = 1).

| Scale | Formula | n = 2 | n = 3 | Source |
|:--|:--|:--|:--|:--|
| Filter visibility ρ² | 1/(16k² r_max²), r_max² = 24, 1200, 41160 | 5.79e-6 | 9.49e-8 | C19 eq. 16 [proved] |
| Clock length C | (2m+4)2^m + m − 4, m = ℓ − 2 | 5124 | 24582 | C8 RUN7 §37 [proved] |
| Clock gap Γ | ≈ 3π²J/C² | 1.13e-6 | 4.90e-8 | C18 §5 [proved] |
| Clock length margin γ_len (frozen) | ≈ 9π²J/(8C²) | 5.89e-7 | 2.43e-8 | C18 §5 [proved] |
| Initialization gap g_in | between Jπ²/C³ (endpoint pin) and the proved bound ~ J C⁻⁵ | 7.3e-11 to 1.7e-16 | 6.6e-13 to 6.5e-20 | C18 §4 [proved] |
| C7 neighbouring-length margin η_b (b = ℓ) | → 2π²t/b⁴ | 1.34e-3 | 6.87e-4 | selector_wrapper eq. 6 [proved] |
| C7 head gap | ≈ 3π²t/b² | 0.284 | 0.200 | C7 [proved] |
| Wrapper scale U | ≤ U_max(b) → π²t/b⁴ | 6.1e-4 | 3.2e-4 | selector_wrapper eq. 7 [proved sufficient] |
| Helix overlap, gap, compensator (φ = 1/b, r_p = 1/2) | cos(φ/2)^{2b}; 1/((b+1)(2b+1)); (1/r₁+1/r₂)/sin²(φ/2) | 0.975, 4.3e-3, 1.6e3 | 0.979, 3.1e-3, 2.3e3 | helix_equalizer §3 [proved] |
| Payload margins, Sidon design | δ₂, δ₄, δ_mix; v_min | 0.698, 0.0233, 0.349; 1.44 | 0.433, 0.0120, 0.205; 1.38 | ir_control recipe, C13 eq. 5 [numerical, bounds proved] |
| Payload margins, C19 geometric design | δ₂; v_min; pair momentum min 2θ_a | 0.099; 0.101; 0.214 | 0.0134; 0.0134; 0.029 | C19 eq. 7 [numerical] |
| C20 existence margin (no flavour U(1), all mixed vertices) | 2π/(4M(n)), M = 49, 217 | 0.0321 | 0.00724 | C20 eq. 16 to 17 [proved] |
| Second-order logical bandwidth W | λ²ρ² max t_a /Δ_src | | | C19 eq. 23 [proved] |
| Length-fluid bandwidth J_len | u A_b² β_wrap β_helix | | | helix_equalizer eq. 11 [proved] |
| Third-order traceless coefficient | c₃ ρ² λ³/Δ_src², c₃ = 2.26, 4.47, 6.52 at n = 1, 2, 3 | 4.47 | 6.52 | this folder [numerical] |
| Third-order non-Gaussian part | (in units of ρ²λ³/Δ_src²) | 2.7e-4 | 5.9e-3 | this folder, matches C19 §10 at n = 1, 2 [numerical] |
| Weight-error sensitivity S(n) | cross-flavour hopping per unit relative weight error, worst (typical) | 30 (3.4) | 160 (23) | this folder [numerical] |

Two lessons follow. ρ² falls by a factor between 110 (n = 1 to 2) and about 35 (large n) per added flavour, and the Cauchy matrix C_ef = 1/(e+f) behind the filters has condition number 1.4 × 10³ and 4.6 × 10⁴ at n = 2, 3, although the cross-flavour sensitivity S(n) grows only about fivefold per flavour. And C19's geometric Fermi-point design is unusable beyond n = 2, with v_min δ₂ = 1.8 × 10⁻⁴ at n = 3 against 0.60 for the Sidon design. Switching costs nothing, since the hopping bound, not the field bound, sets ρ in every case computed (n ≤ 12).

### 3.2 The ordering and a worked chain

For the frozen, clock-selected version the sufficient ordering is

J ≫ Γ ≳ γ_len ≫ g_in ≫ M ≳ Δ_src ≫ λ ≫ W,

with M the norm of the energy program. Each "≫" is implemented in `scales.py` with a safety factor 0.1 and the following conditions. (i) The Schur error of C18 §6, about M/g_in for M ≈ Δ_src, must be below the weight accuracy ε_w = 0.1 v_min δ₂/S(n), so Δ_src = ε_w g_in. (ii) Third order must sit below the margins, c₃ λ/Δ_src ≤ 0.1 v_min δ₂. (iii) For an exactly annihilated code there is no first-order leakage and no splitting; otherwise both must be below W v_min δ₂. (iv) The program frustration must be below the length margin, which (i) implies. With the Sidon design:

| | n = 2 | n = 3 |
|:--|:--|:--|
| ε_w | 3.3e-3 | 3.7e-4 |
| λ/Δ_src | 2.2e-2 | 9.1e-3 |
| Δ_src (pin-lift g_in / proved g_in) | 2.5e-13 / 5.5e-19 | 2.5e-16 / 2.4e-23 |
| λ | 5.5e-15 / 1.2e-20 | 2.3e-18 / 2.2e-25 |
| W in units of J | 7.2e-22 / 1.6e-27 | 2.0e-27 / 1.9e-34 |

If the energy events sit at clock times of small sine weight f_t (endpoint f_t ≈ 2π²/C³), M must exceed Δ_src by 1/f_t and Δ_src falls by that factor. Astra C23 (RESULT §7 to §8), using C20's existence margins, its rigorous cubic bound and a stricter RG budget, reaches W ≈ 1.4 × 10⁻⁵³ (n = 2) and 1.7 × 10⁻⁶² (n = 3). Both tables are feasibility arithmetic on assumed resources, not parameters of a constructed Hamiltonian, and C23 rightly stresses that the resources stop being independent once a spatial parent exists.

The mobile version replaces the top of the chain by t ≫ U ≤ U_max(b) ≫ Δ_wrap, where Δ_wrap collects the wrapper's own vacancy, length-filter and helix gaps (selector_wrapper eq. 10). The core's clock hopping J_c is free there, because the wrapper bound only needs positivity, and the core's source poles are then placed below Δ_wrap, Δ_src ≪ Δ_wrap, after which the chain continues as above. The length-fluid bandwidth is set comparable to W by u ≈ W/(A_b² β).

### 3.3 Is there an ordering for every n?

Yes, at the physics standard and conditional on the compiler. Every condition is an upper bound on a smaller scale in terms of larger ones, and none bounds a small scale from below by a larger one, so choosing scales from the top down succeeds at each fixed n, where all margins, ρ², c₃ and S are positive and finite. Nothing needs to be uniform in n, because each h_n is a separate chain whose c is set at its own fixed couplings. The one quantity that does not shrink with λ or u is a second-order mismatch between the two cell lengths, which scales like W itself. §4 argues that in the arrow-wrapped architecture it is controlled by other tunable ratios. In a design whose two lengths have different source spectra it would have to be removed by construction.

## 4. Corrections the pieces generate on each other

A structural fact controls most of this. At ρ = 0 C19's cell is H′ ⊗ I_logical, the ports act trivially on the logical register, and the cell is analytic in ρ. Exact parity forbids an O(ρ) effective term, since a single insertion of D_e is odd on one cell. Hence, for C19's cell, every traceless effective term at every order in λ carries ρ² [proved, from C19 §4 to §6]. The k-th order then has relative size (λ/Δ_src)^{k−2} times a coefficient, rather than the (λ/Δ_src)C_n/ρ² of C19's conservative sufficient condition. At third order the coefficient is c₃ = 2.3, 4.5, 6.5, far below the naive bound r_max² = 24, 1200, 41160 [numerical], and a factor 2.6 × 10⁵ (n = 2) and 1.1 × 10⁷ (n = 3) below C23's rigorous operator-norm ratio 27λ/(4Δρ²) (C23 eq. 24). The coupling window is therefore λ/Δ_src ≈ 10⁻², not C23's sufficient 10⁻¹⁶ to 10⁻²⁰. Higher-order coefficients are assumed similar. A compiler must preserve the spectator property at ρ = 0, or its violation enters at relative size 1/ρ².

Table 4. Corrections, relative to the target hopping coefficient.

| Correction | Size | Effect | Verdict |
|:--|:--|:--|:--|
| Third order, same flavour, number conserving | 1.40, 1.33 × λ/Δ_src (n = 2, 3) | shifts t_a, θ_a, μ_a and the Fermi points | renormalizes; keep shifts below δ/16 (C20) or retune |
| Third order, same-flavour pairing | 0.18, 0.20 × λ/Δ_src | momentum 2θ_a | renormalizes only with displaced centres; this is the routeA_numerics pairing mechanism |
| Third order, cross-flavour conversion and pairing | 0.38, 0.39 and 0.10, 0.15 × λ/Δ_src | momenta p_b − p_a, p_a + p_b | renormalizes if below v_min δ₂ |
| Third order, non-Gaussian | 2.7e-4, 5.9e-3 × λ/Δ_src (operator norm in units of ρ²) | quartic; forward parts deform the Gaussian form | renormalizes if ε log(1/δ₄) ≪ 1 (ir_control H4) |
| Finite clock coupling (C18 §6 Schur bound) | weight errors ε_w → cross-flavour hopping S(n) ε_w | Gaussian | renormalizes if S ε_w ≪ v_min δ₂; zero first order for timewise-annihilated codes |
| Code splitting from a non-FF program or initialization | ε_split | on-cell fields, flavour mixing, anomalous bilinears | must be ≪ W v_min δ₂; structural unless the code is exact |
| Residual length–payload coupling | see below | density-current cross stiffness plus mixed vertices at ±2πν ± p_i ± p_j | renormalizes if small and G > 0 (C13 §3, C20 §3) |
| Marker hop × port insertion | ~ uλ/Δ_min | length-dependent hopping | tunable: relative size (λ/Δ_min)/(A_b²β) when J_len ≈ W |
| Two marker hops | ~ u²/Δ_min | length interactions, payload-dependent through compensators | tunable |
| Virtual excursions (leakage) | O(λ/Δ_min) beyond the target | quasi-local retarded kernels | renormalizes at the physics standard (IR_AUDIT §9); mathematically open (ir_control H6) |

Every row renormalizes couplings in the sufficient regime of §3.2 [heuristic], given a Sidon-type design, exact parity and a second-order length equalization. No row was found that destabilizes the target there.

The length–payload row needs care. First-order helix equalization does not by itself equalize the second-order bond (C19 §11), and endpoint rescaling cannot equalize cells whose source gaps differ (C23 eq. 27). In C8's arrow wrapper, however, both lengths carry the same core and differ by one vacancy. An end vacancy only hides the port. An interior vacancy switches off one core bond (`tagged_vacancy.md` eq. 5) and changes the source spectrum, but it costs an extra 2U against hopping w < U/(k+1), so its weight is O((w/U)²). Place the source poles far below the wrapper's own excitations, Δ_src ≪ Δ_wrap, which the top-down ordering permits. Then vacancy and helix excitations are far off shell during an excursion, each port insertion picks up the endpoint amplitude ⟨χ_p|P_occ|χ_p⟩ = r_p and each compensator insertion its expectation m_p, and the two-cell bond is (r_p m_p)²(r_p′ m_p′)² K₂. The existing choice m_p = 1/r_p equalizes it, with residuals O(Δ_src/Δ_wrap) and O((w/U)²). The logical part of the first is further suppressed, because Σ_e X_ea = 0 removes the constant and linear terms of the shifted Cauchy kernel [heuristic, derived here; G3 is the check]. So for this architecture the mismatch appears tunable, not structural. It would be structural (order W) whenever the two lengths have different source spectra, which includes the clock-selected cell. C23 §6.3 offers a length fluid faster than W as a fallback that averages a residual mismatch into a retuned target.

## 5. The large-n accumulation point

Normalize each h_n by its largest coefficient, which in the mobile construction is either the C7 head hopping t (with K_sel = 4t + 9U) or, if they are chosen larger, the core's positive clock and legality terms J_c. The wrapper requires U ≤ U_max(b) → π²t/b⁴, and every scale that carries the logical physics (endpoint rewards, vacancy hopping, filters, energy program, λ, u) sits below U, so all of them vanish relative to the norm as b(n) → ∞, while μ/t → 2. If t dominates, a convergent subsequence tends to the C7 selector at its dilute threshold μ = 2t (C23 eq. 33), with head band bottom tk² (z = 2) and cells of diverging length costing π²t/l², plus whatever core terms were kept of order t. If J_c dominates, the limit keeps positive clock constraints whose zero-energy histories exist at every core length and whose tower gaps J_c/C² become arbitrarily small in long cells. Either way h_∞ is gapless with a ground-state degeneracy growing exponentially in L, hence non-conformal with extensive zero-temperature entropy [heuristic]. C13 §9 and C23 §9 leave the limit undetermined among such options; the point here is that each natural normalization gives a non-conformal limit. Physical logical velocities, about bW, vanish faster than any selector scale and the conformal onset length diverges, as the c-theorem argument in STATUS requires. The exact h_∞ depends on compiler coefficients that do not yet exist [assumed].

## 6. What the numerics already taught

- `contrib/routeA_numerics/`. A G1-type cell reproduces the SU(3) exchange coefficient 27/8192 but breaks SU(3) at relative order λ, and at λ = 0.2 the spectrum drifts away from SU(3)₁, so symmetric targets are dropped. An XX target with Fermi points at ±k_F is gapped (ξ ≈ 1.25/λ) by a generated pairing γ ≈ 0.4λ. Table 4 shows the same mechanism in C19's cell (same-flavour pairing 0.14 to 0.20 λ/Δ_src at n = 1 to 3). Nothing in a parity-only cell suppresses it, so the defence is kinematic: Sidon pair momenta 2θ_a are 2.40 and 2.01 at n = 2, 3, against 0.21 and 0.029 for C19's geometric design.
- `contrib/encoded_chain_numerics/`. Route B gives c = 1.00 for r = 2 up to λ/Δ = 0.8, and c ≈ 2 for r = 4 at accessible scales. Generated O(λ²) chemical potentials can dope a sector away from umklapp, and degenerate Fermi points leave marginal non-abelian couplings whose sign matters. Route A's distinct Fermi points remove the latter issue.
- `contrib/mobile_cells/`. A d = 6 mobile toy has a unique ground state on every ring tested and c = 2 from three estimators; c_p = 2 is not yet established.
- `contrib/routeA_incommensurate/` contains a graded cell with a displaced incommensurate target and its ED and VUMPS drivers, but no data yet.

## 7. Remaining gaps, ranked, with the next calculation

1. G1. A fixed-alphabet port interface giving parity-odd, logical-operator-valued second-order exchange with an exact frustration-free code. Known no-gos: bare positive Bell checks (C18 §7), unitary-decorated paths and trees (C21, C22), boundary moments local in the register, which C19's original first moments violate for a product layout (C21), at least ⌈n/2⌉ + 1 poles (C19 §8). C21 reports delayed Jacobi filters with k = n + s poles whose first s − 1 moments are scalar, which eases the moment constraint. C22's merged 17-state bus with uniform flux fails at n = 2 (nonnormal hopping matrix, gapped half-filled band); non-uniform fluxes are open. The live candidate is a syndrome-mixing positive check such as `contrib/routeA_demonstrator/cell.py`'s H_pair, which annihilates Φ⁺ and couples Φ⁻ to |01⟩, |10⟩ through s^±_c (n = 1, d = 10, no results yet). Next calculation: put a syndrome-mixing check on C18's 15-letter clock with a two-mode swap program at n = 2, compute the source spectral measure in the invariant sector as in C18 §8, and test for independent η ⊗ γ hopping with a nonzero phase.
2. G2. The whole-cell compiler. Initialize n modes, derive n from ℓ, and produce n + 1 poles and weights, or any nonresonant Gaussian target, with relative accuracy ε_w ≤ 0.1 v_min δ₂/S(n) after the sine-weight correction. Next calculation: the complete rule table at n = 2, with the Schur error measured against ε_w = 3.3 × 10⁻³.
3. G3. Second-order composition with the wrapper. Next calculation: wrap C18 §8's small cell (exact code, nonzero second-order ZZ, legal space 1920), with its clock in frustration-free Laplacian form, in the tagged-vacancy wrapper at lengths b and b + 1. Compute both second-order bonds in the exact invariant sector, and the mismatch with m_p = 1/r_p as a function of Δ_src/Δ_wrap and w/U. The same model gives the uλ mixed term. Confirm that G1's port terms keep the core frustration-free, which the wrapper theorem needs.
4. G4. Finite-coupling numerics on a parity-only incommensurate route-A chain. Next calculation: complete routeA_incommensurate with the Sidon design at n = 2, logical chain h₂ + h₃ from cluster ED, two calibrated estimators at λ/Δ_src = 0.02 to 0.1. Prediction from Table 4: c = 2.
5. G5. Switch C19's target to the Sidon design and rerun its exact checks. Short.
6. G6. Second-order infrared transfer. Restate IR_AUDIT eqs. (5) to (10) with W = λ²ρ²/Δ_src and the ρ² structure of §4, and compute the fourth-order coefficient at n = 1, 2 to test that it stays of order one.
7. G7. Gluing data for C20 §4 (translation character, winding lattice) for a parity-only payload with one Jordan–Wigner grading inside the C8 wrapper. C13 §6.2 has the free reference.
8. G8. The alphabet. The mobile version costs d0 = 96q_core + 48 (helix_equalizer eq. 1). With the nine- or fifteen-letter clock, two port letters and a single work qubit (q_core = 9w + 4 or 15w + 4, w = 2) that is already 2160 or 3312. The estimates d0 ≈ 45 and 130 in minimal_d0 apply to frozen cells. Count after G1 and G2.
9. G9. Scope. Irrational density gives moving zero-mode sectors (C13 §7). This is a question of definition, not of construction.

## 8. Does anything kill route A?

No. The located obstructions restrict implementations, and two constraints follow from combining the pieces. C23 reached (a) independently.

(a) The clock-selected length (the d0 ≈ 45 architecture) is a frozen selector. Cells of lengths b and b + 1 would carry clocks of different length and unrelated histories (C23 §2.2 makes the same point), so a local marker hop has no useful overlap and the second-order source spectra differ, which is C23's rescaling obstruction. The mobile version must use the C7 walker with the arrow wrapper, which keeps one fixed-length core and moves only vacancies, at a factor of about 96 in d0.

(b) C8's wrapper theorem uses positivity of every non-reward term (selector_wrapper eq. 2) and a frustration-free core (`tagged_vacancy.md` eq. 1). So the core's clock must be in Laplacian rather than hopping form, and its energy program must annihilate the code termwise, or the theorem must be extended with the core's frustration energy below tη_b − U. A non-frustration-free energy program, one of C19 §9's two suggested escapes, therefore carries an extra obligation in the mobile setting. Syndrome-mixing positive checks satisfy it automatically.

Neither constraint nor any other result found here excludes the remaining construction. The fixed-alphabet interface and compiler (G1, G2) are a finite construction problem with located no-gos and live candidates.

## Reproducing

numpy and scipy only.

- `python3 -B scales.py` (seconds) writes `scales_results.json` and prints a summary (saved as `scales_stdout.txt`). It asserts C19's ρ² at n = 1, 2, C18's table entries and C20's counts M(1..8).
- `python3 -B third_order.py 1 2 3` (about 20 s at n = 3, roughly 2 GB peak) writes `third_order_results.json` (stdout saved as `third_order_stdout.txt`). The `-B` flag keeps Python from writing bytecode into `codex/work_C19/`.

## Relation to the rest of the repository

This is the C23 assembly requested in `specs/C23_routeA_assembly.md`, written alongside Astra's closed-book `codex/work_C23/RESULT.md`. The two agree on the verdict (no kill; compiler decisive; FF wrapper hypothesis; frozen clock-selected cells). They differ in three places: the third-order window (numerical c₃ and the ρ² structure here, against C23's rigorous but loose bound), the second-order equalization (tunable in the arrow-wrapped geometry here, open there), and the accumulation point (non-conformal in every natural normalization here, undetermined there). It does not change "Where things stand". It suggests three amendments to STATUS step 5 and the route-A paragraph: the d0 ≈ 45 estimate is for frozen cells only; the second-order coupling window is λ/Δ_src ≲ 10⁻² at n ≤ 3, not λ/Δ ≪ ρ²; and C19's geometric Fermi points should be replaced by the Sidon design.
