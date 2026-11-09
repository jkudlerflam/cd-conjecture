#!/usr/bin/env python3
"""C177 (Claude Opus 5.5 agent, 2026-11-09): the two-magnon problem at total momentum K, in the full qubit chain
(ring ED) and in a two-body effective model, for the coordinator's check of C175's bound pair at K = pi.

Effective two-body model on a ring of L sites: two hard-core bosons that never pass (blockade d >= 2 on both sides),
one-particle hops t_m (m = 1..4) taken from the exact full-chain band (band.py convention
eps(k) = c0 + sum_m [t_m e^{-ikm} + c.c.]: t_1 = i a1/2, t_2 = b2/2, t_3 = i a3/2, t_4 = b4/2), interactions V(d) at
d = 2, 3 (and L - 2, L - 3).  Relative coordinate d = x2 - x1 in [2, L-2], psi(x1, x1 + d) = e^{iK x1} f(d):
  (H f)(d) = sum_m [t_m e^{-iKm} f(d+m) + t_m^* e^{iKm} f(d-m) + t_m f(d-m) + t_m^* f(d+m)] + (2 c0 + V(d)) f(d),
with f = 0 outside [2, L-2]; bosonic states satisfy f(d) = e^{iKd} f(L-d) (the two labellings of one configuration).
Hops by m >= 2 that would pass the other particle land outside [2, L-2] and are dropped.

Full chain: E_2(K) - E_vac from momentum-resolved ring ED (est.RingTable), even parity, at the band-fit couplings
(mu_eff = -12, where four-particle states are far above two-particle ones).
"""
import sys, json
import numpy as np
import scipy.sparse.linalg as spla
import est as E
import band as B


def eff2(L, K, tm, c0=0.0, V=None, nev=4):
    V = V or {}
    ds = np.arange(2, L - 1)
    n = len(ds)
    H = np.zeros((n, n), complex)
    idx = {d: i for i, d in enumerate(ds)}
    for i, d in enumerate(ds):
        H[i, i] += 2 * c0 + V.get(d, 0.0) + V.get(L - d, 0.0)
        for m, t in tm.items():
            for dd, amp in ((d + m, t * np.exp(-1j * K * m) + np.conj(t)), (d - m, np.conj(t) * np.exp(1j * K * m) + t)):
                if dd in idx:
                    H[i, idx[dd]] += amp
    assert np.abs(H - H.conj().T).max() < 1e-12
    # bosonic subspace: P f(d) = e^{iKd} f(L-d)
    P = np.zeros((n, n), complex)
    for i, d in enumerate(ds):
        P[i, idx[L - d]] = np.exp(1j * K * d)
    w, U = np.linalg.eig(P)
    Bs = U[:, np.abs(w - 1) < 1e-8]
    Q, _ = np.linalg.qr(Bs)
    Hb = Q.conj().T @ H @ Q
    return np.sort(np.linalg.eigvalsh((Hb + Hb.conj().T) / 2))[:nev]


def tm_from_fit(f):
    return {1: 0.5j * f['a1'], 2: 0.5 * f['b2'], 3: 0.5j * f['a3'], 4: 0.5 * f['b4']}


def band_min(tm, c0=0.0, nk=200001):
    k = np.linspace(-np.pi, np.pi, nk)
    e = c0 + sum((t * np.exp(-1j * k * m) + np.conj(t) * np.exp(1j * k * m)).real for m, t in tm.items())
    return float(e.min())


def full2(h, L, nev=3):
    T = E.table(L)
    out = {}
    for m in (0, L // 2):
        s = T.sector(m)
        H = T.sector_matrix(h, m)
        pc = np.array([bin(int(x)).count('1') for x in s['R']])
        idx = np.nonzero(pc % 2 == 0)[0]
        Hs = H[idx][:, idx]
        w = np.sort(spla.eigsh(Hs, k=nev, which='SA', tol=1e-12)[0])
        out[m] = w
    return out[0][0], out[L // 2]


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'check':  # C175: E - 2 eps_min = -0.01953 at L = 240, |w| = 1 (bare model)
        tm = {1: 1j, 2: -1.0}
        em = band_min(tm)
        for L in (120, 240, 400):
            e = eff2(L, np.pi, tm, V={2: -2.0, 3: -2 / 3})
            print('bare |w|=1 L', L, 'E - 2 eps_min =', np.round(e - 2 * em, 5), flush=True)
    elif mode == 'full':
        T = {d['r']: d for d in json.load(open('tune.json'))}
        bands = {}
        out = []
        for r in [float(x) for x in sys.argv[2].split(',')]:
            for L in [int(x) for x in sys.argv[3].split(',')]:
                cp = B.couplings(r, mu_eff=-12.0, Dscale=T[r]['Dscale'])
                h = B.bond(**cp)
                Ev, eps = B.band(h, 12)
                f, res = B.fit(eps, 12)
                t2 = -f['b2'] / 2
                ev, e2 = full2(h, L)
                tm = tm_from_fit(f)
                eb = eff2(L, np.pi, tm, c0=f['c0'], V={2: -2.0 * t2, 3: -2 / 3 * t2})
                em = band_min(tm, f['c0'])
                rec = dict(r=r, L=L, t2=t2, fit=f, E2_full=(e2 - ev).tolist(), E2_effbare=eb.tolist(), two_eps_min=2 * em)
                out.append(rec)
                print(f"r={r} L={L}: (E2 - 2eps_min)/t2 full {np.round((e2 - ev - 2 * em) / t2, 5).tolist()} "
                      f"eff(bare V) {np.round((eb[:3] - 2 * em) / t2, 5).tolist()}", flush=True)
        json.dump(out, open(sys.argv[4], 'w'), indent=1)
