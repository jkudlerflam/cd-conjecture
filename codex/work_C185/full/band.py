#!/usr/bin/env python3
"""C177 (Claude Opus 5.5 agent, 2026-11-09): exact one-particle band of the full qubit chain at the Fibonacci cusp.

Bond (C152 convention) h = Jx XX + Jy YY + D (XY - YX) + Jz ZZ + (hz/2)(Z1 + 1Z), plus optional generic terms
Jxz (XZ + ZX) and a transverse field hx (hx/2)(X1 + 1X).  Delta = Jx - Jy = 1, V = 4 Jz.
A particle is a down spin (bit 1).  On a ring of L sites (momentum sectors from est.RingTable) the dressed vacuum
is the lowest even-popcount state at k = 0 and the dressed one-particle band is
    eps(k) = E_odd(k) - E_vac,
exact up to corrections that fall exponentially in L (the vacuum is gapped, correlation length of order one site).
Fit (k -> pi - k symmetric harmonics, the symmetry of the J_x + J_y = 0 stratum):
    eps(k) = c0 + a1 sin k + b2 cos 2k + a3 sin 3k + b4 cos 4k,
compared with the effective model (C173 eq. 3.2, |t2| units) eps = 2|w| sin k - 2|t2| cos 2k - mu.
So |w| = a1/2, |t2| = -b2/2, mu_eff = -c0/|t2|.  The antisymmetric harmonics (cos k, sin 2k, ...) are also fitted,
to check the symmetry numerically.
"""
import sys, json
import numpy as np
import scipy.sparse.linalg as spla
import est as E

SX, SY, SZ, I2 = E.SX, E.SY, E.SZ, E.I2


def bond(Jx, Jy, D, Jz, hz, Jxz=0.0, hx=0.0):
    k = np.kron
    return (Jx * k(SX, SX) + Jy * k(SY, SY) + D * (k(SX, SY) - k(SY, SX)) + Jz * k(SZ, SZ)
            + 0.5 * hz * (k(SZ, I2) + k(I2, SZ)) + Jxz * (k(SX, SZ) + k(SZ, SX)) + 0.5 * hx * (k(SX, I2) + k(I2, SX)))


def couplings(r, wa=1.0, mu_eff=-2.2, Dscale=1.0, Jxpy=0.0):
    """Leading-order map (C173 full_dmrg.couplings) at phi = pi/2, Delta = 1, V = r.
    Dscale multiplies the leading-order D; Jxpy adds J_x + J_y."""
    V = float(r); t2 = 1.0 / (2 * V)
    Jz = V / 4; D = Dscale * wa * t2 / 2
    Jx = (Jxpy + 1.0) / 2; Jy = (Jxpy - 1.0) / 2
    mu = t2 * (mu_eff + 6); hz = (mu - 4 * Jz) / 2
    return dict(Jx=Jx, Jy=Jy, D=D, Jz=Jz, hz=hz)


def band(h, L, parity_split=True):
    T = E.table(L)
    Ev = None; eps = np.zeros(L)
    for m in range(L):
        s = T.sector(m)
        H = T.sector_matrix(h, m)
        pc = np.array([bin(int(x)).count('1') for x in s['R']])
        for par in ((0, 1) if m == 0 else (1,)):
            idx = np.nonzero(pc % 2 == par)[0]
            Hs = H[idx][:, idx]
            if Hs.shape[0] <= 400:
                w = np.linalg.eigvalsh(Hs.toarray())[0]
            else:
                w = spla.eigsh(Hs, k=1, which='SA', tol=1e-13)[0][0]
            if par == 0:
                Ev = w
            else:
                eps[m] = w
    return Ev, eps - Ev


def fit(eps, L):
    k = 2 * np.pi * np.arange(L) / L
    cols = [np.ones(L), np.sin(k), np.cos(2 * k), np.sin(3 * k), np.cos(4 * k),
            np.cos(k), np.sin(2 * k), np.cos(3 * k), np.sin(4 * k)]
    A = np.stack(cols, 1)
    c, *_ = np.linalg.lstsq(A, eps, rcond=None)
    res = float(np.max(np.abs(A @ c - eps)))
    names = ['c0', 'a1', 'b2', 'a3', 'b4', 'anti_cos1', 'anti_sin2', 'anti_cos3', 'anti_sin4']
    return dict(zip(names, map(float, c))), res


if __name__ == '__main__':
    out = []
    for r in [float(x) for x in sys.argv[1].split(',')]:
        for L in [int(x) for x in sys.argv[2].split(',')]:
            cp = couplings(r, mu_eff=-12.0)
            h = bond(**cp)
            Ev, eps = band(h, L)
            f, res = fit(eps, L)
            t2 = -f['b2'] / 2; w = f['a1'] / 2
            rec = dict(r=r, L=L, couplings=cp, Evac_per_site=Ev / L, eps=eps.tolist(), fit=f, fit_resid=res,
                       w=w, t2=t2, w_over_t2=w / t2, t2_over_lo=t2 * 2 * r, mu_eff=-f['c0'] / t2,
                       a3_over_t2=f['a3'] / t2, b4_over_t2=f['b4'] / t2)
            out.append(rec)
            print(json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in rec.items() if k not in ('eps', 'couplings', 'fit')}),
                  {k: round(v / t2, 5) for k, v in f.items()}, flush=True)
    if len(sys.argv) > 3:
        json.dump(out, open(sys.argv[3], 'w'))
