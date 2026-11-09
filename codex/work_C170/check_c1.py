#!/usr/bin/env python3
"""C170: corrected second-order N^4 coefficient at the switch, tested on random harmonics against lr2.dE2.
Claim: u^2 coefficient of dE2 = c2 - c1 with
  c2 = (g^2 + gxy^2) cot^2 q / 4,
  c1 = sin^2(q/2) |tan(q/2) Gam_c - i D_c|^2,  Gam_c = Gam_xz - i Gam_yz,  D_c = D_y + i D_x,
  Gam_xz = (J_xz + J_zx)/2, D_y = (J_zx - J_xz)/2, Gam_yz = (J_yz + J_zy)/2, D_x = (J_yz - J_zy)/2, q = -a.
Fields contribute only to the linear term."""
import numpy as np
from lr2 import dE2
rng = np.random.default_rng(1)
L = 100
for a in (2 * np.pi * 13 / 100, 2 * np.pi * 31 / 100):
    q = -a
    for t in range(4):
        v = 0.1 * rng.standard_normal(9)
        g, gxy, Jxz, Jzx, Jyz, Jzy, hx, hy, hz = v
        dJ = np.array([[g, gxy, Jxz], [gxy, -g, Jyz], [Jzx, Jzy, 0.]]); dB = np.array([hx, hy, hz])
        us = np.linspace(0, 1, 7)
        e = [dE2([np.sqrt(1 - u) * np.cos(0.3), np.sqrt(1 - u) * np.sin(0.3), np.sqrt(u)], dJ, dB, L=L, a=a)[0] for u in us]
        quad = np.polyfit(us, e, 2)[0]
        Gxz, Dy, Gyz, Dx = (Jxz + Jzx) / 2, (Jzx - Jxz) / 2, (Jyz + Jzy) / 2, (Jyz - Jzy) / 2
        Gc, Dc = Gxz - 1j * Gyz, Dy + 1j * Dx
        c2 = (g * g + gxy * gxy) / np.tan(q) ** 2 / 4
        c1 = np.sin(q / 2) ** 2 * abs(np.tan(q / 2) * Gc - 1j * Dc) ** 2
        c1_C165 = (Gxz ** 2 + Gyz ** 2) * (1 - np.cos(q)) ** 2 / (2 * (1 + np.cos(q)))
        print(f'a={a:.3f} trial {t}: numeric u^2 coeff {quad:+.6e}  corrected c2-c1 {c2-c1:+.6e}  C165 c2-c1(Gamma only) {c2-c1_C165:+.6e}')
