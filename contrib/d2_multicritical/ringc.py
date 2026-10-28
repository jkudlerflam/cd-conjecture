#!/usr/bin/env python3
"""Ring estimators of c at one point (and on the calibration controls).

  c_S     chord fit of the ring ground-state entropy, l >= L/4, at the largest L.
  c_half  slope of S(L/2) against (1/3) ln L over the sizes given.
  c_Cas   -6A/(pi v) from E0/L = e + A/L^2 + B/L^4, with v the velocity from the
          lowest level at relative momentum 4 pi/L in the ground-state parity
          sector, read as the stress tensor (x = 2) of a single-velocity CFT.  For
          several velocities the same fit gives sum_a c_a v_a / v.
  towers  the lowest levels per (relative momentum, parity), in units 2 pi v/L.

usage: ringc.py MODEL --L 10 12 14 16 18 20 [--out data/ringc_X.json]
MODEL is a control (ising, xx, cluster_ising, xx_pair, ising_xxz) or FAM:k=v,...
For the decoupled d = 4 controls the ring of d = 4 sites is the product of two
qubit rings, so its spectrum and entropies are built exactly from the two qubit
chains (energies and momenta add, entropies add).
"""
import os
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '1')
import sys, json, argparse
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ring
from scan import bond_of

X, Y, Z, I2 = ring.X, ring.Y, ring.Z, ring.I2
ISING = ring.bond(np.diag([-1., 0, 0]), [0, 0, -1.])
XX = ring.bond(np.diag([1., 1., 0]), [0, 0, 0])
XXZ05 = ring.bond(np.diag([0.5, 0.5, 0.25]), [0, 0, 0])
CLUSTER = -np.kron(np.kron(X, Z), X) + np.kron(np.kron(Y, Y), I2)
PAIRS = dict(xx_pair=(XX, XX), ising_xxz=(ISING, XXZ05))


def qubit_data(L, h, k=4):
    sp, key, psi = ring.spectrum(L, h, k=k, want_gs=True)
    return sp, key, ring.entropies(psi, L)


def combine(L, A, B):
    """Product of two qubit rings (sectors labelled (m, P) with P the product parity)."""
    spA, kA, SA = A; spB, kB, SB = B
    out = {}
    for (ma, pa), ea in spA.items():
        for (mb, pb), eb in spB.items():
            key = ((ma + mb) % L, (pa or 1) * (pb or 1))
            vals = (np.asarray(ea)[:, None] + np.asarray(eb)[None, :]).ravel()
            out.setdefault(key, []).extend(vals[np.isfinite(vals)].tolist())
    out = {k: np.sort(v)[:6] for k, v in out.items()}
    gs = ((kA[0] + kB[0]) % L, (kA[1] or 1) * (kB[1] or 1))
    return out, gs, SA + SB


def data(model, L):
    if model in PAIRS:
        hA, hB = PAIRS[model]
        return combine(L, qubit_data(L, hA), qubit_data(L, hB))
    if model == 'cluster_ising':
        return qubit_data(L, CLUSTER)
    if model == 'ising':
        return qubit_data(L, ISING)
    if model == 'xx':
        return qubit_data(L, XX)
    fam, rest = model.split(':', 1)
    kw = {s.split('=')[0]: float(s.split('=')[1]) for s in rest.split(',') if s}
    return qubit_data(L, bond_of(fam, **kw))


def analyse(recs):
    Ls = np.array(sorted(recs), float)
    E0 = np.array([min(np.nanmin(v) for v in recs[int(L)][0].values()) for L in Ls])
    out = dict(L=Ls.tolist(), E0=E0.tolist())
    vs = []
    for L in Ls.astype(int):
        sp, gs, S = recs[L]
        e0 = min(np.nanmin(v) for v in sp.values())
        m0, P0 = gs
        cand = []
        for dm in (2, -2):
            key = ((m0 + dm) % L, P0)
            if key in sp:
                cand.append(np.nanmin(sp[key]) - e0)
        vs.append(L * min(cand) / (4 * np.pi))   # stress tensor: x = 2 at relative momentum 2
    out['v'] = vs
    A = np.vstack([np.ones_like(Ls), Ls ** -2, Ls ** -4]).T
    coef = np.linalg.lstsq(A, E0 / Ls, rcond=None)[0]
    out['A'] = float(coef[1])
    out['c_Cas'] = float(-6 * coef[1] / (np.pi * vs[-1]))
    Lmax = int(Ls[-1])
    out['c_S'] = ring.c_chord(recs[Lmax][2], Lmax)
    half = np.array([recs[int(L)][2][int(L) // 2 - 1] for L in Ls])
    ev = Ls % 2 == 0
    if ev.sum() >= 2:
        out['c_half'] = float(3 * np.polyfit(np.log(Ls[ev][-4:]), half[ev][-4:], 1)[0])
    # tower table at Lmax: levels in units 2 pi v / L, relative momentum and parity
    sp, gs, S = recs[Lmax]
    e0 = min(np.nanmin(v) for v in sp.values())
    tab = []
    for (m, P), v in sp.items():
        dm = (m - gs[0]) % Lmax
        dm = dm - Lmax if dm > Lmax // 2 else dm
        for e in v:
            if np.isfinite(e):
                tab.append((float((e - e0) * Lmax / (2 * np.pi * vs[-1])), int(dm), int(P)))
    tab.sort()
    out['tower'] = tab[:40]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('model')
    ap.add_argument('--L', type=int, nargs='+', default=[10, 12, 14, 16, 18, 20])
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    recs = {L: data(a.model, L) for L in a.L}
    res = analyse(recs)
    res['model'] = a.model
    print(a.model, 'c_Cas=%.4f c_S=%.4f c_half=%s v=%s' % (res['c_Cas'], res['c_S'], res.get('c_half'),
                                                           np.round(res['v'], 4).tolist()))
    print('tower (x, dk, P):', [(round(x, 3), d, p) for x, d, p in res['tower'][:16]])
    if a.out:
        json.dump(res, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
