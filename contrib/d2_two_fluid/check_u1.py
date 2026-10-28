#!/usr/bin/env python3
"""U(1) two-up-spin check: exact relative-coordinate levels (twobody.isolated_levels with
gL = gR = 0) against E(K) = V + c^2/V - 2mu, c = 2t cos(K/2 + phi); pair minimum Q = -2 phi
= 2 k0; and, for V < -2t, the pair onset mu_p = (V + 4t^2/V)/2 against the all-up energy per
site V - mu_p = (V^2 - 4t^2)/(2V) < 0."""
import numpy as np
from twobody import isolated_levels, pair_band_u1
rng = np.random.default_rng(7)
K = np.linspace(-np.pi, np.pi, 2001, endpoint=False)
err = 0; qerr = 0; worst_up = -np.inf; n = 0
for _ in range(200):
    phi = rng.uniform(-np.pi, np.pi); V = rng.uniform(-6, 6); mu = 0.0
    E1, E2, ed = isolated_levels(K, 1.0, phi, V, mu, 0.0, 0.0)
    # S level (free single) is also an eigenvalue of the truncated problem; the pair is the
    # isolated level that is not the single
    single = -2 * np.cos(K + phi) - mu
    pair = pair_band_u1(K, 1.0, phi, V, mu)
    cand = np.where(np.abs(E1 - single) < 1e-9, E2, E1)
    m = ~np.isnan(pair) & (np.abs(pair - ed) > 1e-6)
    err = max(err, np.nanmax(np.abs(cand[m] - pair[m])) if m.any() else 0)
    if V < -2:
        Q = K[np.nanargmin(pair)]
        qerr = max(qerr, abs((Q + 2 * phi + np.pi) % (2 * np.pi) - np.pi))
        mup = (V + 4 / V) / 2
        worst_up = max(worst_up, V - mup); n += 1
print(f"max |E_numeric - E_formula| = {err:.2e}; max |Q - (-2 phi)| = {qerr:.2e} (grid 2pi/2001 = {2*np.pi/2001:.1e})")
print(f"V < -2t points: {n}; largest all-up energy per site at the pair onset, V - mu_p = {worst_up:.3f} (< 0 means the up-polarized state is already below the vacuum)")
