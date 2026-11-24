Correlators on the χ = 1024 open-chain states (χ = 384/512 for T3), written by `../scripts/meas_obc.py` in the C193 ckpt_meas format, so `codex/work_C193/ckpt_meas/ana193.py` reads them after gunzip.

- T1: δ = 3/16, h = 1.783 (wedge). T2: δ = 1/8, h = 1.7405 (wedge). T3: δ = 3/16, h = 1.772 (c = 1, below the lower edge).
- Frame: the rotated frame of `ims.py` (site j rotated by jπ/2 about z). Particles are down spins, b† = S⁻, n = (1 − Z)/2.
- `meas.Gb_rows[i]` = [Re, Im] of ⟨S⁻_i S⁺_j⟩ for j = i+1 … L−1, i ≤ L/2 + 5.
- `meas.Dpair["d_i0"]` = [Re, Im] of the connected ⟨S⁻_i0 S⁻_i0+d S⁺_j S⁺_j+d⟩ for j = i0+d+1 … L−d−1, d = 2, 3, i0 = L/2−3 … L/2 (already connected, so `mm` and `pp` are zero).
- `meas.Cnn[i0]` = connected ⟨n_i0 n_j⟩, j > i0. `meas.dens` = ⟨n_j⟩.
- Checked against dense vectors at L = 10 to 2.4e-14.
