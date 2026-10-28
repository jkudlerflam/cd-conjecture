#!/usr/bin/env python3
"""Marginal couplings of the q = 3 logical chain with the nearest-neighbour bond to ALL orders
(exact two-cell diagonalization, effective.core_bond) and the three-cell term at second order.
Prints, per lam, the chemical-potential bilinear, spin and charge current couplings (same
normalization as marginal.py: Hubbard U = 1 -> +0.5 charge, -0.5 spin; positive = relevant), and
saves the bonds for VUMPS (data/core_bonds.npz).
Written by Claude (Opus 5.5) for contrib/encoded_chain_numerics, 2026-10-27."""
import sys, numpy as np
import cells as C, effective as EF, marginal as MG
out = {}
for g in (0.5, 0.2):
    c = C.cell(3, 1, 2, 1.0, g)
    so = EF.second_order(c)
    m3 = EF.monomials(so['T3'], c, 3)
    N3 = MG.couplings(MG.project(m3))
    for lam in ((0.1, 0.2, 0.4, 0.8) if g == 0.5 else (0.2, 0.4)):
        hb, gap, smin = EF.core_bond(c, lam, so)
        m2 = EF.monomials(hb, c, 2)
        N = (MG.couplings(MG.project(m2)) + lam ** 2 * N3) / lam ** 2
        hop = [cf for mm, cf in m2 if len(mm) == 2 and mm[0][1] != mm[1][1]]
        mu = [cf for mm, cf in m2 if len(mm) == 2 and mm[0][1] == mm[1][1]]
        t_eff = np.mean([abs(x) for x in hop]) if hop else 0
        mu_eff = np.mean([x.imag for x in mu]) if mu else 0
        print(f'g={g} lam={lam}: 2-cell gap={gap:.3f} smin={smin:.3f} hop/lam={t_eff/lam:.5f} '
              f'(alpha^2={so["alpha"]**2:.5f}) mu/lam^2={mu_eff/lam**2:.5f} '
              f'N_spin/lam^2={np.diag(N)[3]:.4f} N_par/lam^2={np.diag(N)[0]:.4f} N_perp/lam^2={np.diag(N)[1]:.4f} '
              f'max quad offdiag {np.abs(N - np.diag(np.diag(N))).max():.1e}', flush=True)
        out[f'g{g}_l{lam}'] = hb
np.savez('data/core_bonds.npz', **out)
