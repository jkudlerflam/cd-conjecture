# A folded-bus spatial cell for route A

## Claim

Setting: one-site translation-invariant, nearest-neighbour, Hermitian bond on C^23, frozen cells, route-A coupling at second order through the marker-bulk bond. One logical complex fermion is stored in each physical data qubit of a cell, so a cell of n+1 sites carries n flavours. There is no history state and no controller.

1. [proved] An explicit 23-letter bond has an exactly degenerate 2^n code per cell, exact on-site U(1) charge and parity, zero first-order port compressions, no three-cell term at second order, and no odd orders in λ on open chains. The vacuum token is C22's C7-type selector, so d0 = 23 includes length selection.
2. [proved in a limit, numerical at finite parameters] Astra C22 showed that a linear bus fails: the two ports see the flavours in opposite orders, the hopping matrix is non-normal and the bands dimerize. Here the right port reaches the data only through a flavour-blind express track ending at the left port's junction. The hopping matrix then tends to a phase times a positive Gram matrix, and all n bands cross a small common field (checked through T for n ≤ 4, through the bands for n ≤ 3).
3. [numerical] A weak link that closes a loop separates the n Fermi-sea centres, while particle-hole symmetry keeps the real part of the one-cell matrix at zero.
4. [numerical] ED of two and three cells at n = 1, 2 (Hilbert spaces up to 373248) reproduces the exact second-order bond with error proportional to λ⁴.
5. [heuristic] At each fixed n, with n-dependent small parameters, the logical chain is n free bands with exact U(1), distinct centres and non-resonant Fermi points, inside the stable region of `contrib/ir_control/`, so c = n at d0 = 23 on compatible rings. The scales shrink exponentially (or faster) in n, and the cells are frozen.

## Status

Numerical and heuristic. The finite-dimensional statements in 1 are proved. Upgrades would be (i) a proof of the Gram limit with explicit control at finite D_x, (ii) a band check for n ≥ 4 in extended precision, and (iii) a mobile, equal-payload version.

## How it was established and checked

- Author: Claude (Opus 5.5) subagent for Jonah Kudler-Flam, 2026-10-28.
- AI assistance: all design, code and runs by Claude. No Codex call. It builds on Astra's finished C22 (`codex/work_C22/RESULT.md`) and uses its 17-letter alphabet, gauge, selector and exchange formula. Interim C21 notes were relayed by the coordinator.
- Checked by: internal cross-checks only (ED against the exact resolvent, two extraction methods, real-rail and particle-hole controls). Not independently rerun.

## Argument

### 1. Starting point

[inherited, C22 §1-5] C22's 17 letters (L_a, R_a, X_{s,a} with s = v, 0, +, -, and M_∅, M_s) carry one token per cell. The vacuum token selects the length. Excited species have charge 0, ±1, neutral hopping carries Z on the departure data qubit, and a conversion g turns 0 into ± while flipping the data bit. C22 proved that uniform rail phases give T equal to a phase times a real matrix, hence a common Fermi-sea centre. It also proved that at n = 2 the leading T = [[a,b],[c,0]] with |c| > |b| gaps both bands, because the two ports see data site i at distances n-i and i.

### 2. The bond

[definition] Add six letters Xx_{s,a} for s = 0, +, -, a token on an express track. The vacuum species never uses it. Letters keep C22's arrows (Xx = (L,R)), so the legal-cell space is still one token per cell, with dimension (7n+4)2^n. Terms, each with its Hermitian conjugate:

- vacuum hop on data sites and into the marker, amplitude -t (C22);
- excited hop on data sites (the stub), -t e^{iφ_s}(-1)^{a[s=0]}; excited species may not enter the marker from the stub;
- express hop Xx R → L Xx and Xx M_∅ → L M_s, -t_x e^{iφ^x_s}(-1)^{a[s=0]};
- junction switch (M_*, Xd_{s,a}) ↔ (M_*, Xx_{s,a}) on every marker-bulk pair, -t_sw, so the switch acts only on the first site of a cell;
- weak loop link Xd_{s,a} M_∅ → L_a M_s, -t_d e^{iφ^d_s}(-1)^{a[s=0]};
- on-site energy Δ for excited Xd, Δ + D_x for excited Xx and M_s, C22's conversion g, marker potential μ_sel, a uniform data field ε, and penalty K on forbidden pairs;
- port λ[F0^M F0^X + e^{iφ}F+^M F-^X + e^{-iφ}F-^M F+^X] on (marker, Xd) pairs.

Working point: t = t_x = t_sw = 1, Δ = 0, D_x = D_M = 4, g = 0.05, φ^x_0 = 0.4, φ^d_0 = 0.9, φ = π/2, other phases zero. Particle-hole asymmetry is zero.

[proved] Every move preserves the outer arrows and the species charge, so C22's legal-sector argument, the charge Q = data + token charge and the parity Π(-1)^{Q_site} (logical parity on the code) carry over. Each F changes vacuum to excited. Hence G†FG = 0, no three-cell term arises at second order, and (-1)^{N_exc} on alternate cells sends V to -V. The vacuum sector is C22's path of b = n+1 vertices, so its selector window and margin Θ(t b⁻⁴) are inherited. [numerical] The code is the full ground space of every tested cell, with the gaps in §5.

### 3. Why folding repairs the exchange

At t_d = 0 the excited token graph is a path: M, x(n-1), ..., x(0), d(0), d(1), ..., d(n-1). The left source sits at the junction d(0) and the right source at M. Both see data site j at graph distance j beyond d(0).

[proved, structure] On a tree, the resolvent between M and d(j) factorizes as G_{M,d0}G_{d0,dj}/G_{d0,d0}. The transfer factor G_{M,d0}/G_{d0,d0} is flavour-independent. As D_x → ∞ it becomes energy-independent across the window of stub energies, with relative variation of order nW/D_x [heuristic estimate]. In that limit the right-port filters are a constant times the left-port filters f_j(τ), so T → e^{iΦ} c Γ with Γ_ij = ∫_0^∞ f_i f_j dτ. Γ is positive definite, because f_j(τ) starts at order τ^{2j+1} at short times, so the f_j are linearly independent. Then h(k) = μ + Te^{ik} + T†e^{-ik} equals μ ± 2R_s at k = -Φ and π - Φ, where R_s is the Hermitian part of e^{-iΦ}T. If R_s > 0 and ‖μ‖ < 2λ_min(R_s), every band changes sign on both arcs, which gives at least 2n Fermi points.

[numerical] Eigenvalues of R_s/|T| from one-particle matrix elements (`hermpart.json`):

| n | D_x | asymmetry ‖R-R^T‖/‖R‖ | eig(R_s)/\|T\| |
| --: | --: | --: | :-- |
| 2 | 4 | 0.10 | 0.0176, 1.22 |
| 3 | 4 | 0.10 | 2.1e-4, 0.064, 2.02 |
| 4 | 4 | 0.07 | -6.1e-6, 9.6e-4, 0.075, 2.14 |
| 4 | 10 | 0.05 | 1.3e-5, 3.9e-3, 0.127, 2.13 |
| 4 | 20 | 0.03 | 4.5e-5, 5.0e-3, 0.155, 2.12 |

At n = 4, D_x = 4 is too small and D_x = 10 suffices, which fits D_x growing with n. [numerical] Real express rails (φ^x_0 = 0, t_d = 0) give |T| ≈ 10⁻¹⁹, C22's sine cancellation. The express phase is what makes the hopping nonzero.

### 4. Distinct Fermi-sea centres

[numerical] At t_d = 0, T is real to 10⁻⁹ after removing its phase and μ vanishes to roundoff, so the centres coincide. Breaking particle-hole symmetry fails. Charged energies Δ ± 0.3 give |μ|/|T| ≈ 2.4 × 10³ at n = 2, because [heuristic] one-cell terms at the directly coupled left port carry no express suppression. The loop link works. The cycle M → express → junction → stub → M carries a species-dependent flux, μ becomes iA with A real antisymmetric and ∝ t_d (|Re μ| ≤ 10⁻¹⁷|T|), h(k) loses its k → -k antiunitary symmetry, and the centres separate linearly in t_d. The link reopens the reversed-order route, so it must stay weak. At n = 2, D_x = 10, Δ = 0.5, t_d = 10⁻² already makes R_s indefinite.

### 5. Logical chain at the working point

[numerical] `logical_n*.json`. The quadratic part is the Hilbert-Schmidt projection, referred to half filling, and ε = ½λ_min(R_s). The non-Gaussian entry is the norm of the remainder over |T|.

| n | t_d | cell gap | poles k | \|T\| (λ=1) | γ_min/γ_max | Fermi points | centres | δ₂ | δ₄ | v_min/\|T\| | non-Gaussian/\|T\| |
| --: | --: | --: | --: | --: | --: | --: | :-- | --: | --: | --: | --: |
| 1 | 10⁻³ | 0.70 | 9 | 3.6e-7 | 1 | 2 | 0 | 2.64 | 1.01 | 1.94 | 2.1e-3 |
| 2 | 10⁻³ | 0.27 | 40 | 1.1e-8 | 0.014 | 2+2 | ±0.0017 | 0.54 | 6.8e-3 | 0.047 | 5.9e-3 |
| 3 | 10⁻⁵ | 0.12 | 78 | 1.3e-10 | 1.0e-4 | 2+2+2 | 0.0013, 0, -0.0014 | 0.68 | 2.6e-3 | 9.5e-4 | 2.2e-2 |

Here δ₂ is the smallest separation of Fermi points, and δ₄ is the smallest |p₁+p₂-p₃-p₄| mod 2π over distinct pairs. That includes inter-flavour Cooper transfer, which is 2(centre_a - centre_b). [numerical] The non-Gaussian part scales as g². At g = 0.01 it is 2.3 × 10⁻⁴ (n = 2) and 8.7 × 10⁻⁴ (n = 3), with the band data unchanged. [numerical] k ≥ n+1 by a wide margin, so the pole bound 2n ≤ 4(k-1) is not the binding constraint. The binding constraint is the conditioning of Γ.

### 6. Exact diagonalization

[numerical] `ed_check.py` builds open chains of M legal cells from the cell matrices. It diagonalizes each exact charge sector (dense, or shift-invert Lanczos), forms the canonical effective Hamiltonian by polar decomposition against the unperturbed code, and compares it with M E0 + λ² Σ_bonds H2. The ED point uses Δ = 0.5, D_x = 2, g = 0.3, t_d = 0.05 for a larger T. The IR point is the working point with t_d = 10⁻³.

| point | n | M | dim | λ | ‖error‖ | error/λ⁴ | T relative error |
| :-- | --: | --: | --: | :-- | :-- | --: | :-- |
| ED | 1 | 2 | 484 | .2, .1, .05 | 4.3e-5, 2.7e-6, 1.7e-7 | 0.027 | 3.5e-4, 1.1e-4, 2.9e-5 |
| ED | 2 | 2 | 5184 | .2, .1, .05 | 1.2e-5, 7.4e-7, 4.6e-8 | 0.0074 | 1.3e-2, 3.1e-3, 7.9e-4 |
| ED | 1 | 3 | 10648 | .2, .1, .05 | 2.3e-5, 1.4e-6, 9.0e-8 | 0.014 | |
| ED | 2 | 3 | 373248 | .2, .1 | 1.9e-5, 1.2e-6 | 0.012 | |
| IR | 1 | 2 | 484 | .1, .05 | 1.8e-6, 1.1e-7 | 0.018 | 2.8e-4, 7.0e-5 |
| IR | 2 | 2 | 5184 | .1, .05 | 5.0e-7, 3.2e-8 | 0.0050 | 2.0e-3, 4.9e-4 |
| IR | 1 | 3 | 10648 | .1, .05 | 1.3e-6, 8.2e-8 | 0.013 | |

Eigenpair residuals are below 6 × 10⁻¹⁴, and the retained band is at least 0.27 below the next level. A constant error/λ⁴ confirms the absence of odd orders and of second-order three-cell terms. The ED hopping matrix converges to T as λ². These are full legal-cell spaces, with no excitation truncation. No central charge is measured.

### 7. Scale hierarchy and infrared conditions

[heuristic, with the numbers above] In decreasing size:

- the cell gap, 0.70, 0.27, 0.12, 0.07 for n = 1 to 4. At Δ = 0 it is of order t n⁻³ minus t_sw²/D_x, so Δ ≳ t_sw²/D_x + |g| gives an n-independent bound;
- the selector margin Θ(t b⁻⁴) (C22);
- J_n = λ²|T|, with |T| ∝ g²α_n⁴τ_x², α_n = √(2/(n+2)) sin(π/(n+2)) (C22) and express transfer τ_x ~ (t_x/D_x)^{n+1};
- the narrowest band, 10⁻² to 10⁻⁵ of the widest for n = 2 to 4. That is a factor 50 to 100 per flavour, consistent with C21 (far flavours enter only through high boundary moments);
- the field window |ε| < 2γ_min, and the loop splitting ∝ t_d.

For `ir_control` (verdict b) and C20: exact U(1) forbids pairing, which answers `contrib/routeA_numerics/`. Distinct centres remove the common-centre Cooper channel, and δ₂, δ₄ > 0 at fixed n. ε_NG/(v_min δ₂) is 0.009 at n = 2, g = 0.01. Higher orders need λ ≲ 10⁻² at n = 2. All conditions are met at n = 2 and 3 by explicit choices. For general n, take D_x large, then t_d, g, ε and λ small, in that order [heuristic, fixed n]. J_n then falls faster than exponentially, which is allowed for existence but puts n ≥ 4 beyond double-precision band checks here.

### 8. d0, and what remains

[proved count] d0 = 23 for the complete frozen-cell bond, selector included (C22's failing design has 17, C18's one-qubit cell 64). If the heuristic conditions hold for n ≥ 23, this beats c ≤ d0 - 1 on compatible rings. Open: (i) quantitative control of the Gram limit at finite D_x; (ii) frozen cells, whose payload grows with length (n = b - 1), so C20's equal-payload mobile pair does not exist as is; (iii) the very narrow window, with ir_control's leakage step assumed. No obstruction specific to spatial designs was found. The two that came up, reversed port orders and common centres from tree-like phases, are removed by the fold and by one loop.

## Reproducing

numpy and scipy only, in this folder. `ycell.py` holds the literal bond and the cell, code, exact second-order bond, fermion fit and bands. `python3 logical_chain.py n [t_d ...]` writes `logical_n{n}.json` (seconds for n ≤ 2, minutes for n = 3). `python3 hermpart.py` writes `hermpart.json` (about 10 minutes). `python3 ed_check.py [ir] n M λ ...` writes `ed_n{n}_M{M}[_ir].json`. M = 3, n = 2 takes about 30 minutes per λ on 6 threads, with stdout in `ed_n2_M3_stdout.txt`.

## Relation to the rest of the repository

This supplies the repair that C22 left open, a mechanism that breaks the common-phase relation and repairs the non-normal exchange. It addresses the spatial route-A item in STATUS "Where help is most useful", for frozen cells. It escapes C22's tree theorem through conversions and the loop, and C18's static-Bell-check obstruction because the port creates a mobile token. STATUS should not change until the mobile and IR steps are done.
