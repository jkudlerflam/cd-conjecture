#!/usr/bin/env python3
"""Low-lying Sz = 0 ring levels exactly at u = r (frustration-free line), scaled as L^2 (E - E0) with
one-site momentum labels from G = T o prod X. Compared with the exact one-defect branch over the Neel
state, E_-(k) = 4 r - 4 sqrt(r^2 - sin^2(k/2)) (cell momentum k), whose L^2 E at the lowest k is
L^2 * E_-(4 pi / L). A second soft branch would show up as extra levels with L^2 E = O(1).
Written for codex/work_C153 (Claude Opus agent, 2026-11-07)."""
import numpy as np
from spec_ed import spectrum
for b in [0.0, 1.0, 2.0, 4.0]:
    r = np.sqrt(1 + b * b / 4)
    for L in [12, 16, 20]:
        w, ps = spectrum(L, r, b, k=10)
        k1 = 4 * np.pi / L
        e1 = 4 * r - 4 * np.sqrt(r * r - np.sin(k1 / 2) ** 2)
        print(f'b={b} L={L} E0={w[0]:.2e} one-defect L^2E(k1)={L*L*e1:.3f} | ' +
              ' '.join(f'{L*L*(x-w[0]):.2f}[{p/(2*np.pi/L):+.0f}]' for x, p in zip(w[1:], ps[1:])), flush=True)
