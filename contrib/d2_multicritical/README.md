# Where Ising lines meet c = 1 regions in qubit chains

## Claim

Class: periodic, one-site translation-invariant, nearest-neighbour Hermitian qubit chains in the C15 convention plus a transverse field,

H = Σ_j [Jx XX + Jy YY + Jz ZZ + D (XY − YX) − h Z_j + hx X_j],

or in Kitaev variables τ = −(Jx + Jy), δ = Jx − Jy, V = 4Jz (μ = V − 2h). The quantity is the total Virasoro c.

Families: XYZ in a longitudinal field (catalog A7, containing A6 and the dual ANNNI chain A8) in the slices δ/τ = 0.2, 0.5, 1, 2; its extension into the one-unitary-Z₂ stratum by the axial DM term D; and XXZ at Δ = 2 in a tilted field, an antiunitary-only family that contains A5 as the XY coupling goes to zero.

In these families every point where an Ising line meets a c = 1 region is non-conformal (a band edge or Lifshitz point with z > 1, or a classical multiphase point) or has c = 1 (an eight-vertex point, or the Kosterlitz–Thouless end of A6). No Ising line runs inside a c = 1 region, so there is no Ising × Luttinger (c = 3/2) or two-boson (c = 2) point. Every ring value above 1.3 occurs where a free-fermion control with c = 1 gives the same overshoot.

## Status

Numerical, plus a standard bosonization argument for why the lines do not meet. Rings to L = 20, VUMPS to χ = 64. Tracking the tilted-field Néel line at L ≥ 24, and scanning the TR-only and no-symmetry strata, would upgrade it.

## How it was established and checked

- Author: Claude (Anthropic, Opus 5.5), a subagent of the coordinating session, 2026-10-27/28. No Codex/GPT was used.
- Code: `ring.py` (momentum- and parity-resolved ED, orbit builder copied from `../d2_classification/ed.py`; its self-test matches dense diagonalization in every (k, P) sector to 10⁻¹³), `vumps.py` (copied unchanged from `../d2_followup_220_1282/`), and the scripts below.
- Checked by: nobody else. Every estimator was run on the controls below.

## Calibration

c_S is the chord fit of the ring ground-state entropy over l ≥ L/4 at the largest L, c_half the slope of S(L/2) against (1/3) ln L, and c_fit six times the slope of S against ln ξ over the last four bond dimensions of VUMPS.

| control (exact c) | c_S | c_half | VUMPS c_fit |
|---|---|---|---|
| Ising (1/2), L ≤ 20 | 0.501 | 0.499 | 0.498 (d2_followup) |
| XX (1) | 1.004 | 0.992 | 1.000 (d2_followup) |
| cluster-Ising, three-site (3/2), L ≤ 18 | 1.56 | 1.46 | 1.500 (d2_followup) |
| two XX chains as one d = 4 site (2) | 2.008 | 1.983 | 1.94 (d2_followup) |
| Ising ⊗ XXZ(Δ = 1/2), d = 4, two velocities (3/2) | 1.506 | 1.492 | 1.42 (ξ ≤ 41) |
| free chiral Kitaev, V = 0, δ = 0.5, D = 0.35 (1) | 1.13 to 2.37, L = 20 | | |

The d = 4 rings are built exactly from two qubit rings. In the Ising ⊗ XXZ spectrum the extra Majorana tower is visible directly (x_σ = 0.126 in Ising-velocity units), while a single-velocity Casimir fit gives 1.81, so Casimir values are not used for verdicts. The last row is the main warning: in a chiral Fermi sea, where c = 1 exactly, ring entropies at L ≤ 20 overshoot to 2.4.

## Lines tracked and what they meet

**1. A7, the band-edge Ising line.** In every slice the polarized phase is bounded by an Ising line where the Bogoliubov gap at k = 0 closes, located by crossings of L·|E(0,−) − E(0,+)| at L = 12/16/20. At δ = 0.5, h_I = 0.9995 (V = 0, exact 1), 1.97 (V = 2), 2.96 (V = 4), 3.973 (V = 6) and 4.974 (V = 8). Near h_I at V = 6 and L = 20 the tower gives σ/ε ≈ 0.16 and L₋₁σ/ε ≈ 1.14 (Ising 1/8 and 9/8), with v_I ≈ 0.45. Floating lobes rise from the eight-vertex points at h = 0 (V = 3 for δ = 0.5). Across the four slices c_S(16) never exceeds 1.22 (`data/maps/kit_d*.png`).

On rings the lobe seems to reach the Ising line for V ≳ 3. It does not. Along V = 6 the density estimate K_L(L = 20) rises from 0.65 at h = 2.8 to 1.35 at h = 3.95, and VUMPS separates the phases.

| point (δ = 0.5) | ring c_S, L = 20 | VUMPS |
|---|---|---|
| V = 6, h = 2.0 (lobe) | 1.0 | c_fit = 0.997, ξ = 43 to 164 |
| V = 6, h = 3.7 | 0.87 | gapped, S = 0.263, ξ = 6.7 to 6.9 for χ = 32 to 64 |
| V = 8, h = 4.8 | 0.95 | gapped, S = 0.256, ξ = 7.8 |

At V = 6, cells of 2, 3, 4, 5 and 7 sites all converge to the same gapped state (energy equal to 10⁻¹⁰), and the ring energies approach it from below as e^(−L/ξ) with ξ ≈ 7. So a gapped topological strip (here 3.7 ≤ h < 3.97 at least) separates the Ising line from the lobe. In the dual ANNNI slice the two approach only at the classical multiphase point (V/h → 2, h → ∞). As δ → 0 the Ising line becomes the saturation edge of XXZ in a field, z = 2.

The reason is general for band-edge lines [heuristic, standard bosonization]. The line is the onset of a dilute liquid with K → 1, where the pairing source cos 2θ (dimension 1/K) is relevant, while a floating phase needs K < 1/2. An Ising mode inside a floating phase would need a second, independent low-energy degree of freedom, and in a one-orbital chain the onset fermions are the particles that form the liquid.

**2. Z₂ stratum, adding D.** The k = 0 sectors barely feel D. On the trivial side the quasiparticle gap at k = 0 is 0.170 (V = −1, h = 0.62) and 0.145 (V = 2, 4) for every D from 0 to 0.45 at L = 16, and at V = −1 the k = 0 crossing stays at h_I = 0.5543 with σ/ε = 0.127. At V = 0 this is exact, since D couples to a conserved current. The line ends at D*, where the quasiparticle branch at k → 0⁺ goes soft and the ground state moves to nonzero momentum (a chiral Fermi sea). From E_qp(2π/16) = 0, D* = 0.25 exactly at V = 0, and D* ≈ 0.40, 0.25, 0.23 (±0.03) at V = −1, 2, 4.

The endpoint is a Lifshitz point. At V = 0 the σ level stays at L·ε₁ = 0.78 (L = 10 to 18) while higher levels collapse with z_eff = 1.8 to 2.8, toward the cubic branch z = 3 of C15/C16. In every interacting case the first soft momentum is k → 0⁺, as in the free chain (a dip near k ≈ 3π/4 is a three-quasiparticle state, present at V = 0 too). Inside the chiral region on the continuation h = h_I (V = −1, D = 0.5), c_S(20) = 1.037 with no feature at h_I. The largest chiral value in the D maps (1.86 at V = −2.25, h = 0.25, D = 0.6, L = 16) falls as 2.06, 1.86, 1.75, 1.67 for L = 14 to 20, inside the range of the free chiral control and next to the ferromagnetic point.

**3. Antiunitary-only stratum.** For XXZ at Δ = 2 in a tilted field (hx ≤ 2, h ≤ 8, L = 12/16), c_S(16) ≤ 1.21 everywhere. Above the Néel phase every hx in the window shows a c ≈ 1 region with ground-state momentum alternating between 0 and π up to saturation. The A5 Néel Ising line is not resolved at these sizes. Its expected junction with the c = 1 region is the commensurate–incommensurate point hx = 0, h = h_c1 (z = 2), because a uniform transverse field is non-resonant with the antiferromagnetic transverse correlations (momenta π and ±2πm), and its square locks the boson only where K > 1/2. In A7 itself hx couples linearly to the k = 0 σ and turns the band-edge line into a crossover (standard, not scanned).

**4. Earlier junction.** The A6 Ising line ends on the Δ = −1 segment at h ≈ 0.50 with c_S = 1.03 and c_E = 0.96 (`../d2_classification`), the c = 1 Kosterlitz–Thouless end found by Dmitriev et al.

| meeting point | nature | evidence |
|---|---|---|
| band-edge Ising and U(1) line (δ → 0) | saturation edge, z = 2 | exact |
| band-edge Ising and floating lobe (A7) | do not meet, gapped strip | VUMPS at V = 6, 8; K_L > 1/2 |
| dual ANNNI | classical multiphase point | exact limit |
| eight-vertex points | c = 1, where lobes start | exact (zero-field XYZ) |
| D endpoint | Lifshitz, z > 1 | z_eff → 3 at V = 0, same soft branch at V ≠ 0 |
| A6 junction, Δ = −1 | c = 1 (KT) | c_S 1.03, c_E 0.96 |
| tilted-field Néel line | not resolved, expected z = 2 | map, L ≤ 16 |

## Verdict

No point with c > 1. Every meeting point located is non-conformal or has c = 1. The one route to c = 3/2 these families could offer, an Ising line inside a floating phase, is blocked by a gapped strip that VUMPS resolves where rings cannot. The argument covers band-edge Ising lines in one-orbital chains. It says nothing about an Ising line of an independent degree of freedom, or about the TR-only and no-symmetry strata, which remain the places to look.

## Reproducing

numpy, scipy and numba (as in `../d2_classification`); matplotlib for pictures only. Maps are gzipped JSON lines in `data/maps/`.

- `ring.py` (self-test). Maps: `scan.py kit --fixed tau=1 delta=0.5 [D=...] --x V -4 8 49 --y h 0 6 25 --L 12 16 --out ...` and `scan.py xyz --fixed Jx=1 Jy=1 Jz=2 --x hx 0 2 21 --y h 0 8 33 ...`, then `mapplot.py`, `maplines.py` (about 1 core-hour per slice).
- Lines: `track.py kit --fixed tau=1 delta=0.5 --scan V 0 8 5 --solve h 0.5 5.5 --sig 0,-1 --L 12 16 20` (largest-h root); D version in `data/tracks/`.
- Cuts: `cut.py kit --fixed tau=1 delta=0.5 V=6 --x h 2.8 4.4 33 --L 12 16 20`, `cutplot.py`.
- Dispersions and endpoints: `qpline.py 16 'tau=1,delta=0.5,V=-1,h=0.62' 0 0.2 0.3 0.35 0.38 0.4`, `endpoint.py kit 'tau=1,delta=0.5,V=0,h=1,D=0.25'`.
- Calibration: `ringc.py ising|xx|cluster_ising|xx_pair|ising_xxz --L 10 ... 20`; VUMPS `fes.py MODEL --n 2 --chis 8 12 16 24 32 48 64` (5 to 10 min each).

Typhon job 17513250 (apollo, 32 cpus, about 10 min) computed the A7 slices. Its second array task failed, so the D slices cover every other grid point (`*.half.*`). Jobs 17513223, 17513234, 17513235, 17513398, 17513478, 17513485, 17513524 and 17513562 waited about two hours without starting and were cancelled. Everything else ran on the laptop.

## Relation to the rest of the repository

This is the numerical side of the d = 2 residual (STATUS.md, small d; C15 §4, C16 item 4, C17 item 1). It finds C17's top-ranked open mechanism absent in A7 and its D and tilted-field extensions. "Where things stand" is unchanged.
