# A c ≈ 2 wedge in a parity-only qubit chain, at d = 2

## Result

The chain is

    H = Σ_j [ a X_jY_{j+1} + b Y_jX_{j+1} + Z_jZ_{j+1} ] − h Σ_j Z_j,   a = (1+δ)/2, b = (1−δ)/2,

with Z2 parity ΠZ, translation, and two antiunitaries. Just below the line where the staggered chirality κ_j = X_jY_{j+2} − Y_jX_{j+2} orders, there is a wedge that reads c ≈ 2 by two methods.

| point | iDMRG (χ ≤ 320) | open chains (L = 96 to 384, χ ≤ 1024) |
|---|---|---|
| δ = 3/16, h = 1.783 | 1.973 ± 0.010 | 1.96 ± 0.08 |
| δ = 1/8, h = 1.7405 | 1.919 ± 0.009 | 1.90 ± 0.12 |
| δ = 3/16, h = 1.772 (c = 1 neighbour) | 1.00 | 1.00 ± 0.01 |
| controls XX × gapped TFIM, XX × critical TFIM, XX × XX | 1.00 to 1.10, 1.42 to 1.52, 1.96 | 0.99, 1.48, 1.96 |

| δ | c = 1 liquid to | c ≈ 2 wedge | κ locks at | gapped from |
|---|---|---|---|---|
| 1/16 | gapped to 1.725, then c ≈ 1 | none seen | — | — |
| 1/8 | 1.738 | 1.739 to h_I | 1.7467 | 1.761 |
| 3/16 | 1.776 | 1.780 to h_I | 1.7875 | 1.792 |
| 7/32 | 1.805 | 1.808 to h_I | about 1.8105 | 1.812 |
| 1/4 | 1.834 | none | none | 1.836 |

A π/4 rotation maps the chain onto the family of the STATUS c = 2 candidate at V/Δ = 4, with δ = 1/8 on its tuning condition (REPORT §9.213). Inside the wedge at δ = 3/16 the C193 estimators give K_s ≈ 1.8 to 2.0, the candidate's band (REPORT §9.214).

Status: numerical, not asymptotic. Open: a crossover beyond about 400 sites, the identity of the second mode, whether the lower edge is pair binding, and where the wedge closes between δ = 1/16 and 1/8. The first version of this folder (2026-11-18) claimed c = 3/2 at h_I. That claim is withdrawn (REPORT §9.209b), because the 1.4 to 1.6 readings there are crossovers out of the c ≈ 2 wedge.

## Methods and checks

- iDMRG: TeNPy 1.1.1, 2-site cell, parity conserved, in the real frame where site j is rotated by jπ/2 about z. c from S against (1/6) ln ξ over χ ≥ 96. κ pinned at χ = 16 and released.
- Open chains: finite two-site DMRG. c is the slope across L of the central-quarter average of S(l) against the same average of (1/6) ln[(2L/π) sin(πl/L)]. Checked against a hand-built ED at L = 12 to 12 digits.
- Correlators: `meas_obc.py`, checked against dense vectors to 2e-14, written in the C193 format. The data are in [meas/](meas/).
- Authors: J. Kudler-Flam, with AI assistance (Claude Opus 5.5, 2026-11-18 to 24). No second person has checked it.

## Files

- [scripts/](scripts/):
  - `ims.py` (iDMRG driver), `an2.py` (iDMRG ladders), `run3.sbatch`, `submit_all.sh`;
  - `ofd.py` (open-chain driver), `otest.py` (ED and free-fermion checks), `ana_obc.py` (open-chain c), `oe.py` (odd–even energies), `runo.sbatch`, `submit_obc.sh`;
  - `meas_obc.py` and `envelope.py` (correlators).
- [tables/](tables/): χ ladders for the target, the controls and the open chains; exponents and odd–even energies; the cuts below δ = 1/8.
- The states (several GB) are not committed. They regenerate from the scripts.
