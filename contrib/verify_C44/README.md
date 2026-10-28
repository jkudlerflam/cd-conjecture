# Rule-table rebuild of C44 and the n = 1 central charge

## Claim

The C44/C45 referee (`contrib/referee_C44_C45/`) asked for (ii) an independent rebuild of C44's source formula (11), transport overlap (9) and vanishing first-order port coupling, and (i) a calibrated c = 2 for the n = 1 member from two estimators. Rebuilt from the 46-letter rule text, all three agree with C44 at n = 1, 2 and extend to n = 3. The projected n = 1 chain gives c = 2 from finite-entanglement scaling and from open-chain interval entanglement, each calibrated on exactly solvable references.

## Status

Numerical. The rebuild is an independent check at C44's diagnostic parameters. The central charge is that of the projected model, not of the 4464-state chain, so it tests part (b) of the referee's input S (Gaussian stability of the retained theory) at n = 1 and not part (a), the elimination of the auxiliary sectors.

## How it was established and checked

- Author: Claude (Opus 5.5) subagent for Jonah Kudler-Flam, 2026-10-28.
- AI assistance: all code and runs by Claude, partly on Typhon. No Codex, since every step needs the repository's rule text.
- Checked by: internal cross-checks only (independent evaluation routes, a basis-rotation roundoff test, exact free-fermion references). Nothing from C32, C39, C44 or C45 was imported or rerun, and their saved outputs were read only for comparison.

## Method and results

### Rebuild of (9), (11) and the first-order coupling

`core46.py` builds the letters, arrows and moves from C32 §1 and C39 §1, §2 alone (Z on rightward neutral departures, e^{iπ/2n} on neutral express hops, conversion g, rank-one vacuum edges in q, Q and 2Q allocations, onsite terms split over bonds). The kernel is 2^n at n = 1, 2, 3, and the bare bond C39 (6) equals C44's `bare_n1`, `bare_n2` to 8e-17 (C44 orders data bits little-endian). `wrap.py` adds the ordered vacancy fluid with D_j gating whole core bonds, both helix tracks, the compensator C8 (8a) and the selector head.

The second-order bond is computed (a) from exact eigenbases of the charge ⊗ core block and the full 2^N helix track, (b) from C44 (11) as written, and (c) as (b) with the charge ground energy removed exactly and each charge scale (ground, tunnelling partner, fast levels) diagonalized at its own norm. The full helix confirms the reached space of C44 (10) (weight outside it below 6e-15). Route (c) is invariant to 2.5e-16 under random rotations of the core basis and tends to the s/δ → 0 bond linearly in t = s/δ. A plain double-precision eigendecomposition fails at small s, with errors up to 4e-9 at n = 1 and 9e-7 at n = 2, against payload hoppings of 2e-8 and 4e-10.

At C44's own s, its saved bonds agree with (c) to between 2e-13 and 5e-11 at n = 1, and between 2e-11 and 3.4e-9 at n = 2. The n = 2 residual grows with t. Imposing C44's tunnelling gap, which differs from the reflection-split value by 1.3e-7 relative (the doublet is split by only 2.5e-9), removes part of it, leaving at most 1.8e-9. That is still above the n = 2 payload hopping at ρ = 0.1 (4.3e-10). C44's scalar density interactions and its four Kbar − K_bare values are reproduced to every quoted digit. Its traceless length-payload numbers hold at n = 1 but are roundoff at n = 2 (traceless E_LR 5.7e-12 and 3.3e-10 in C44, 3e-14 and 5e-15 here). The occupation filter factors out exactly (residual 1.5e-13 at J_f = 434), the buses never mix, and a Sylvester solve with the full helix agrees with (c) to 7.5e-11. At n = 3 (4032 states per bus) the payload hopping is 9.0e-12 and the traceless mixed terms stay below 1.4e-13.

The first-order port compression is identically zero at all four length pairs for n = 1, 2, 3, since every port term turns a vacuum token into an excited one. `transport.py` builds the marker move literally (selector letters, typed vacancies with H_0 → H_2, both helix tracks, the biased code). Its overlap equals A_b² β_v cos(1/2b)^{2b} times the identity to 1e-16 at n = 1, 2, 3, matches C44's β values, and vanishes without the conversion. A plain eigensolve breaks the reflection symmetry of χ_2 at tiny tunnelling splittings (by 9e-5 at n = 3, ρ = 0.1), so the scripts diagonalize the two reflection sectors separately.

### Central charge of the projected n = 1 chain

The approximation is the second-order projected Hamiltonian on the tiling manifold, in cell coordinates. Each cell carries a length bit and the payload qubit, with bond Σ_{pq} Π_p ⊗ Π_q ⊗ K_pq − j_l(|10⟩⟨01| + h.c.) and a marker chemical potential. K_pq is the rebuilt bond at C44's diagnostic core parameters in the s/δ → 0 limit, exact in ρ, g and D. Dropped are the O(t) remainder, the λ⁴, λ²u_mov and u_mov² terms, readout leakage, the ring gluing and the physical coordinate. The uniform payload hopping phase is removed by the exact gauge c_j → e^{ijφ}c_j. Without it, uniform U(1)-broken MPS appear unable to follow the phase winding, and VUMPS stalls.

The hierarchy matters. The length fluid's density interaction V is attractive and, as s/δ → 0, scales as ρ⁴, with V/|t| = −49, −6.3 and −0.32 at ρ = 0.16, 0.1 and 0.05 (|t| = 1.96e-8 is the payload hopping). At C44's finite diagnostic t the fast term gives V/|t| between −3e4 and −3e5 instead. The length fluid survives only if |V| < 2j_l.

Units are |t|, with j_l = 1, length density ν between 0.29 and 0.38, the payload half filled, and the resonant densities 1/4, 1/2, 3/4 avoided. M0 is the decoupled free reference (c = 2) and C0 the payload alone (c = 1). M1 is the scheduled point ρ = 0.05 (V = −0.32), M2 is M1 with every traceless length-payload term multiplied by 300 (to about 0.3), and C1 is ρ = 0.1 (V = −6.3). FES fits S = (c/6) ln ξ with one-site VUMPS over χ = 16 to 96 (16 to 24 for one species, where larger χ did not converge from random starts), and FES2 takes ξ from the two leading transfer-matrix scales, one per species. The profile fits open-chain S(ℓ) with Friedel terms for ℓ ≥ 8, at L = 32 and 48 with χ = 192 for two species (where M0 matches its exact profile to 4e-4 at L = 32) and at L = 48 and 96 with χ = 128 and 160 for one species. Calibrated values (in parentheses) are 2 × model/M0, or C1/C0, at identical χ or L.

| Model | FES | FES2 | Profile, smaller L | Profile, larger L |
| :-- | :-- | :-- | :-- | :-- |
| M0, exact c = 2 | 1.97 | 1.99 | 1.84 (exact 1.85) | 1.88 (exact 1.90) |
| M1, scheduled | 1.95 (1.97) | 1.99 (2.00) | 1.87 (2.03) | 1.88 (2.00) |
| M2, stress test | 1.96 (1.98) | 1.99 (2.00) | 1.87 (2.03) | 1.87 (1.99) |
| C0, exact c = 1 | 0.95 | n/a | 0.948 (exact 0.948) | 0.967 (exact 0.967) |
| C1, ρ = 0.1 | 0.95 (1.00) | n/a | 0.947 (0.999) | 0.967 (1.000) |

With χ ≥ 32 only, the calibrated FES values for M1 and M2 lie between 1.98 and 2.00, and across fit windows ℓ ≥ 4 to ℓ ≥ 12 the calibrated profile values stay between 1.95 and 2.04. C1 locks every cell to the long length (ν = 1), as V < −2j_l predicts, and both estimators return c = 1 for it.

## Reproducing

numpy and scipy, `python3 -B` from this folder. Item (ii) is `run_source.py 1 2` and `run_source.py 3` (minutes, about 20 at n = 3), plus `compare_C44.py`, `gap_attribution.py`, `roundoff_check.py` and `transport.py` (seconds to minutes). Item (i) starts with `run_c.py prepare`, which writes `kpq.npz`. Then `vumps_fes.py MODEL χ...` and `run_c.py MODEL finite χ L` produce `results*/` (minutes to hours, with `jobs_*.sh` the Typhon arrays), and `fitc.py` writes `c_results.json`.

## Relation to the rest of the repository

This answers the referee's items (i) and (ii) at the level stated. It does not test S(a), the microscopic chain, or the lacunary schedule at large n. "Where things stand" should change only to note that C44's core-dependent formulas have an independent rebuild and that the projected n = 1 chain is calibrated at c = 2.
