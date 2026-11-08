"""C160 Typhon task (Claude Opus agent, 2026-11-08): classical ground textures of random lattice bonds in the
doubly rescaled wedge around the pure-DM direction near h_F (PROOF section 3.3). Longer chains and more starts
than the local pilot texture_periodic.py.
Usage: python3 texture_task.py SEED NBONDS N STARTS OUTFILE
For each bond: J = -1 + eps A(D) + sc eps^2 S', B = sc eps^2 b (eps = 0.2), open chain of N sites, L-BFGS from
STARTS initial states (random and helices of several pitches about random axes). On the bulk window
[N/4, 3N/4] we list the spectral peaks of n_j and test whether every peak with more than 1e-3 of the total
weight is a harmonic m q1 (|m| <= 30, modulo 2 pi) of the strongest one. Output: one JSON line per bond.
"""
import sys, json
import numpy as np
from scipy.optimize import minimize

seed, nb, N, starts, outf = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
rng = np.random.default_rng(1000 + seed)


def A(D):
    return np.array([[0, D[2], -D[1]], [-D[2], 0, D[0]], [D[1], -D[0], 0]])


def eg(p, J, B):
    n_ = len(p) // 2
    th, ph = p[:n_], p[n_:]
    n = np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], 1)
    E = np.sum(np.einsum('ia,ab,ib->i', n[:-1], J, n[1:])) + np.sum(n @ B)
    G = np.zeros_like(n); G[:-1] += n[1:] @ J.T; G[1:] += n[:-1] @ J; G += B
    dth = np.stack([np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), -np.sin(th)], 1)
    dph = np.stack([-np.sin(th) * np.sin(ph), np.sin(th) * np.cos(ph), 0 * th], 1)
    return E, np.r_[np.sum(G * dth, 1), np.sum(G * dph, 1)]


def ground(J, B, eps):
    best = None
    for s in range(starts):
        if s < starts // 4:
            p0 = np.r_[np.arccos(rng.uniform(-1, 1, N)), rng.uniform(0, 2 * np.pi, N)]
        else:
            q = eps * (0.2 + 2.0 * s / starts)
            ax = rng.normal(size=3); ax /= np.linalg.norm(ax)
            e1 = np.cross(ax, [1, 0, 0] if abs(ax[0]) < 0.9 else [0, 1, 0]); e1 /= np.linalg.norm(e1)
            e2 = np.cross(ax, e1); j = np.arange(N)
            tilt = rng.uniform(0, 0.8)
            n = np.cos(tilt) * (np.outer(np.cos(q * j), e1) + np.outer(np.sin(q * j), e2)) + np.sin(tilt) * ax
            p0 = np.r_[np.arccos(np.clip(n[:, 2], -1, 1)), np.arctan2(n[:, 1], n[:, 0])]
        r = minimize(eg, p0, args=(J, B), jac=True, method='L-BFGS-B', options=dict(maxiter=200000, gtol=1e-13, maxcor=40))
        if best is None or r.fun < best[0]:
            best = (r.fun, r.x)
    p = best[1]; th, ph = p[:N], p[N:]
    return best[0], np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], 1)


def peaks(n):
    b = n[N // 4: 3 * N // 4]; M = len(b)
    F = np.fft.rfft((b - b.mean(0)) * np.hanning(M)[:, None], n=16 * M, axis=0)
    P = np.sum(np.abs(F) ** 2, 1); k = np.fft.rfftfreq(16 * M) * 2 * np.pi
    tot = P.sum()
    if tot < 1e-12 * M:
        return dict(kind='uniform')
    loc = np.where((P[1:-1] > P[:-2]) & (P[1:-1] >= P[2:]))[0] + 1
    # weight of a peak: power within +- 2 resolution cells
    res = 2 * np.pi / M
    pk = []
    for i in loc:
        w = P[np.abs(k - k[i]) < 2 * res].sum() / tot
        if w > 1e-3:
            pk.append((float(k[i]), float(w)))
    pk.sort(key=lambda a: -a[1])
    q1 = pk[0][0]
    unexplained = 0.0; worst = None
    for kk, w in pk[1:]:
        ok = any(abs(kk - abs(np.mod(m * q1 + np.pi, 2 * np.pi) - np.pi)) < 2 * res for m in range(2, 31))
        if not ok:
            unexplained += w; worst = kk if worst is None else worst
    return dict(kind='textured', q1=q1, npeaks=len(pk), unexplained_weight=unexplained, first_unexplained=worst, peaks=pk[:8])


eps = 0.2
with open(outf, 'a') as fo:
    for t in range(nb):
        D = rng.normal(size=3); D /= np.linalg.norm(D)
        S = rng.normal(size=(3, 3)); S = (S + S.T) / 2; S -= np.trace(S) / 3 * np.eye(3)
        b = rng.normal(size=3); sc = rng.uniform(0.1, 0.8)
        J = -np.eye(3) + eps * A(D) + sc * eps ** 2 * S; B = sc * eps ** 2 * b
        E, n = ground(J, B, eps)
        rec = dict(seed=seed, t=t, eps=eps, sc=sc, D=D.tolist(), S=S.tolist(), b=b.tolist(), N=N, E_per_site=E / N, **peaks(n))
        fo.write(json.dumps(rec) + '\n'); fo.flush()
