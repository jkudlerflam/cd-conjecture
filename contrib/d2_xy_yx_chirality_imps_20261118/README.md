# A staggered-chirality Ising transition inside an incommensurate liquid, at d = 2

**Correction, 2026-11-19 (REPORT §9.209).** The c = 3/2 claim below is withdrawn. Cuts at δ = 3/16 and 7/32 show the same structure. Below the chirality line there is a region that reads as the c = 2 control does, not as a disordered Ising side:
- fitted c = 1.973 ± 0.010 at δ = 3/16, h = 1.783;
- fitted c = 1.919 ± 0.009 at δ = 1/8, h = 1.7405;
- the order parameter decays as a stable power of ξ, with an exponent that varies continuously and reaches about 1/8 at the line.

Read as an extended c = 2 phase, it forces c ≥ 2 on the line, so the 1.4 to 1.6 readings at h_I are crossovers. The phase is a wedge that closes between δ = 7/32 and 1/4. At δ = 1/4 there is no chirality phase, only a c = 1 liquid up to the gap at h = 1.836.

Table of the cuts:

| δ | c = 1 liquid to | c ≈ 2 region | κ locks at | gapped from |
|---|---|---|---|---|
| 1/8 | 1.738 | 1.739 to h_I | 1.7467 | 1.761 |
| 3/16 | 1.776 | 1.780 to h_I | 1.7875 | 1.792 |
| 7/32 | 1.805 | 1.808 to h_I | about 1.8105 | 1.812 |
| 1/4 | 1.834 | none | none | 1.836 |

The original 2026-11-18 text follows, unchanged.


## Claim

Take the translation-invariant nearest-neighbour qubit chain (d = 2)

    H = Σ_j [ a X_jY_{j+1} + b Y_jX_{j+1} + Z_jZ_{j+1} ] − h Σ_j Z_j,   a = (1+δ)/2, b = (1−δ)/2.

It is Hermitian and has only Z2 parity ΠZ, together with translation, an antiunitary A = K·Π_{j odd} Z_j and an antiunitary inversion. Fix δ = 1/8.

- The staggered vector chirality κ_j = X_jY_{j+2} − Y_jX_{j+2} orders at h_I = 1.7466 ± 0.0006, inside a gapless incommensurate liquid. The spiral pitch is about 34 sites and K ≈ 2.4 from earlier work.
- At h_I the finite-entanglement numerics read as compact boson × Ising, Virasoro c = 3/2:
  - fitted c = 1.53 to 1.59;
  - order-parameter exponent β/ν = 0.13 ± 0.01;
  - no hysteresis.
- Decoupled c = 3/2 controls run with the same code read 1.42 to 1.52, with β/ν between 0.10 and 0.14.
- The point h = 7/4, where the chain was first posed, is 0.0034 inside the ordered phase and flows to c = 1.

If c = 3/2 at h_I holds, the chain violates c ≤ d − 1 at d = 2. It is the "staggered Ising transition inside an incommensurate liquid" listed as an open mechanism in [STATUS.md](../../STATUS.md).

## Status

[numerical]. Infinite-MPS finite-entanglement scaling at χ ≤ 256 (correlation lengths up to about 700 sites), calibrated on decoupled controls with c = 1, 3/2 and 2. Three things would upgrade it:
- an independent measurement of η = 1/4;
- the Ising velocity;
- a second method (for example L ≥ 512 open chains) agreeing on c at h_I.

## How it was established and checked

- Authors: J. Kudler-Flam, with AI assistance.
- AI assistance: Claude Opus 5.5 (code, runs, analysis), 2026-11-18. Earlier exact structure and open-chain DMRG are from a Claude Fable 5.2 agent (2026-11-17/18), not included here.
- Checked by: no second person or second code.
  - The rotated real frame was checked against the original complex frame.
  - The infinite-MPS energies are consistent with ED at L = 12 and 16 from the earlier work.
  - The h_I estimate agrees with an independent open-chain correlation-ratio estimate (1.7465 ± 0.0015).

## Method

- **Solver.** Two-site iDMRG (TeNPy 1.1.1), 2-site cell, parity conserved, in the frame where site j is rotated by jπ/2 about z, where H is real. At each χ, 20 mixer sweeps are followed by mixer-free 20-sweep blocks until E, S and κ change by less than 1e-11, 5e-6 and 5e-6.
  - TeNPy VUMPS diverged on this model for χ ≥ 48 and was not used.
- **Branches.** κ is pinned at χ = 16 and released (also run unpinned, and with the bond order pinned instead).
- **Estimators.**
  - c from S against (1/6) ln ξ, with ξ the largest finite transfer-matrix length.
  - β/ν from d ln κ / d ln ξ.
  - h_I is where that exponent crosses −1/8.

| χ | 128 | 160 | 200 | 256 |
|---|---|---|---|---|
| crossing in h | 1.7482 | 1.7472 | 1.7468 | 1.7467 |

Main readings (c fitted over χ ≥ 96):

| h | 1.738 | 1.7405 | 1.744 | 1.7465 | 1.747 | 1.749 | 1.750 | 1.755 |
|---|---|---|---|---|---|---|---|---|
| c_fit | 1.07 | 1.92 | 1.73 | 1.59 | 1.53 | 1.35 | 1.25 | 1.09 |
| d ln κ / d ln ξ at χ = 200 | symmetric | symmetric | −0.28 | −0.136 | −0.116 | −0.053 | −0.032 | −0.008 |

Calibration (c fitted over χ ≥ 96):

| control | true c | c_fit | d ln m / d ln ξ |
|---|---|---|---|
| XX | 1 | 0.998 | n/a |
| XX × ordered TFIM | 1 | 1.005 to 1.10 | ≈ 0 |
| XX × critical TFIM (velocity ratios 0.35, 1, 2.5, and an XXZ boson with K ≈ 2.3) | 3/2 | 1.42 to 1.52 | −0.10 to −0.14 |
| XX × XX | 2 | 1.96 | n/a |

Full χ ladders are in [tables/](tables/).

## Caveats

- **Below h_I the readings do not fall to 1 at this χ.** Near h ≈ 1.739 to 1.741 they read c ≈ 1.9 to 2.2, stable to χ = 320, with quasi-long-range bond order (XX − YY on next-nearest neighbours, staggered). This is either a separate c ≈ 2 point or narrow regime, with crossovers up to h_I (the reading consistent with the c-theorem and with c = 3/2 at h_I), or an extended c = 2 phase. The latter would conflict with the h_I readings.
- **η not measured independently.** It was not extracted, because the connected correlators carry boson-modulated oscillations.
- **Weak first order not excluded.** One beyond ξ ≈ 700 is not excluded. Alberton, Ruhman, Berg and Altman (PRB 95, 075132) discuss this for an Ising transition coupled marginally to a boson.

## Reproducing

The scripts are in [scripts/](scripts/). They need TeNPy 1.1.1, numpy, scipy and h5py.

- `submit_all.sh` lists every SLURM job; `run3.sbatch` makes them resumable and self-resubmitting.
- `ims.py` is the driver. Example:

  ```
  python ims.py model=target h=1.7465 chis=32,48,64,96,128,160,200,256 branch=pin
  ```

- Controls use `model=ctrl tA=.. JzA=.. J=.. g=.. [tB=..]`.
- `an2.py` prints the χ ladders.
- `pw.py` computes plane-wave dispersions (rough boson velocity v_b ≈ 0.16).
- Cost: about 1 to 6 CPU-hours per (h, χ ≤ 256) ladder on 4 cores.
- The states (275 MB) are not committed. They regenerate from the scripts.

## Relation to the rest of the repository

This addresses mechanism 2 under "Elsewhere at d = 2" in [STATUS.md](../../STATUS.md). It is a numerical candidate for c = 3/2 at d = 2, separate from the c = 2 corner candidate. See also [d2_floating_ising](../d2_floating_ising/).
