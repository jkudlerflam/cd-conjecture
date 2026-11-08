"""Re-minimization of the 48 C160 texture bonds (coordinator, Claude Opus 5.5, 2026-11-08).

Why: the C160 Typhon texture test (texture_task.py) failed its own criterion (46/48 bonds with
non-harmonic peak weight > 1e-2). Two suspected causes: (i) the peak test counts Hanning sidelobes at
q1 +- 2.5 bins as separate peaks and double-counts main-lobe power in their +-2-bin weights, so a clean
helix scores about 0.5; (ii) L-BFGS from random or badly pitched helix starts on an open chain of 2000
sites can stop in multi-domain metastable states. This script separates the two.

For one bond (J = -1 + eps A(D) + sc eps^2 S, B = sc eps^2 b, read from the C160 jsonl):
 1. best uniform cone helix n_j = c a + s (cos(qj) e1 + sin(qj) e2) from the exact per-site energy;
 2. ring minimizations on all cells (P, w) with P <= PMAX and q = 2 pi w/P within +-35% of the helix q
    (plus P = 1, 2), started from the optimal helix: best single-period energy e_ring;
 3. open chain of N sites, minimized from (a) the tiled best ring texture and (b) STARTS old-style starts;
    bulk energy density e_bulk on the middle half for each;
 4. harmonic test on the best open-chain bulk: least-squares fit by harmonics m q (|m| <= 40), with q
    optimized near the ring value; residual fraction r = |n - fit|^2 / |n - <n>|^2.
Consistency with C160 M1 on the lattice: r small (single period), and e_bulk(best) >= e_ring - tol.
Usage: python3 remin.py JSONL LINE PMAX N STARTS OUT.json
"""
import sys, json
import numpy as np
from math import gcd
from scipy.optimize import minimize

jf, line, PMAX, N, STARTS, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
rec = json.loads(open(jf).read().splitlines()[line])
eps, sc = rec['eps'], rec['sc']
D, S, b = np.array(rec['D']), np.array(rec['S']), np.array(rec['b'])
A = np.array([[0, D[2], -D[1]], [-D[2], 0, D[0]], [D[1], -D[0], 0]])
J = -np.eye(3) + eps * A + sc * eps ** 2 * S
B = sc * eps ** 2 * b
rng = np.random.default_rng(7 + line)


def ang2n(p):
    m = len(p) // 2
    th, ph = p[:m], p[m:]
    return np.stack([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)], 1), th, ph


def n2ang(n):
    return np.r_[np.arccos(np.clip(n[:, 2], -1, 1)), np.arctan2(n[:, 1], n[:, 0])]


def eg(p, ring):
    n, th, ph = ang2n(p)
    nn = np.roll(n, -1, 0) if ring else n[1:]
    nl = n if ring else n[:-1]
    E = np.sum(np.einsum('ia,ab,ib->i', nl, J, nn)) + np.sum(n @ B)
    G = np.zeros_like(n)
    if ring:
        G += np.roll(n, -1, 0) @ J.T + np.roll(n, 1, 0) @ J
    else:
        G[:-1] += n[1:] @ J.T
        G[1:] += n[:-1] @ J
    G += B
    dth = np.stack([np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), -np.sin(th)], 1)
    dph = np.stack([-np.sin(th) * np.sin(ph), np.sin(th) * np.cos(ph), 0 * th], 1)
    return E, np.r_[np.sum(G * dth, 1), np.sum(G * dph, 1)]


def relax(n0, ring):
    # small random kick off the poles keeps the angle chart regular
    r = minimize(eg, n2ang(n0), args=(ring,), jac=True, method='L-BFGS-B',
                 options=dict(maxiter=400000, gtol=1e-12, maxcor=50))
    return r.fun, ang2n(r.x)[0]


def frame(x):
    a = np.array([np.sin(x[0]) * np.cos(x[1]), np.sin(x[0]) * np.sin(x[1]), np.cos(x[0])])
    t = np.array([1.0, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(a, t); e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    c, s = np.cos(x[2]), np.sin(x[2])
    e1, e2 = c * e1 + s * e2, -s * e1 + c * e2
    return a, e1, e2


def helix_e(x):
    a, e1, e2 = frame(x)
    al, q = x[3], x[4]
    c, s = np.cos(al), np.sin(al)
    sym = 0.5 * (e1 @ J @ e1 + e2 @ J @ e2)
    asy = 0.5 * (e1 @ J @ e2 - e2 @ J @ e1)
    return c * c * (a @ J @ a) + s * s * (sym * np.cos(q) + asy * np.sin(q)) + c * (a @ B)


best = None
for _ in range(60):
    x0 = np.r_[np.arccos(rng.uniform(-1, 1)), rng.uniform(0, 2 * np.pi), rng.uniform(0, np.pi),
               rng.uniform(0, np.pi / 2), rng.uniform(0, 1.0)]
    r = minimize(helix_e, x0, method='Nelder-Mead', options=dict(maxiter=20000, xatol=1e-12, fatol=1e-15))
    if best is None or r.fun < best.fun:
        best = r
xh = best.x
a, e1, e2 = frame(xh)
al, qh = xh[3], abs(np.mod(xh[4] + np.pi, 2 * np.pi) - np.pi)
if xh[4] < 0:
    e2 = -e2
eh = best.fun


def helix(Nn, q):
    j = np.arange(Nn)
    return np.cos(al) * a + np.sin(al) * (np.outer(np.cos(q * j), e1) + np.outer(np.sin(q * j), e2))


# 2. ring cells
cells = [(1, 0), (2, 1)]
dq = min(0.35 * qh, 0.08)
lo, hi = qh - dq, qh + dq
for P in range(3, PMAX + 1):
    for w in range(1, P // 2 + 1):
        if gcd(w, P) == 1 and lo <= 2 * np.pi * w / P <= hi:
            cells.append((P, w))
ring = []
for P, w in cells:
    q = 2 * np.pi * w / P
    n0 = helix(P, q) if P > 2 else np.tile(a, (P, 1))
    n0 = n0 + 1e-3 * rng.normal(size=n0.shape)
    n0 /= np.linalg.norm(n0, axis=1)[:, None]
    E, n = relax(n0, True)
    ring.append((E / P, P, w, n))
ring.sort(key=lambda t: t[0])
e_ring, Pb, wb, nb = ring[0]


def bulk_e(n):
    i0, i1 = N // 4, 3 * N // 4
    e = np.einsum('ia,ab,ib->i', n[i0:i1], J, n[i0 + 1:i1 + 1]) + n[i0:i1] @ B
    return e.mean()


# 3. open chain
trials = []
tiled = np.tile(nb, (N // Pb + 1, 1))[:N]
E, n = relax(tiled + 1e-4 * rng.normal(size=tiled.shape), False)
trials.append(('tiled_ring', E, bulk_e(n), n))
for s in range(STARTS):
    if s % 2 == 0:
        n0 = rng.normal(size=(N, 3))
    else:
        q = qh * rng.uniform(0.5, 1.5)
        n0 = helix(N, q) + 0.05 * rng.normal(size=(N, 3))
    n0 /= np.linalg.norm(n0, axis=1)[:, None]
    E, n = relax(n0, False)
    trials.append((f'start{s}', E, bulk_e(n), n))
trials.sort(key=lambda t: t[1])
lab, Ebest, ebest, nbest = trials[0]


# 4. harmonic residual
def resid(n, q, M=40):
    x = n[N // 4: 3 * N // 4]
    j = np.arange(N // 4, 3 * N // 4)
    cols = [np.ones_like(j, dtype=float)]
    if q > 1e-9:
        for m in range(1, M + 1):
            cols += [np.cos(m * q * j), np.sin(m * q * j)]
    X = np.stack(cols, 1)
    coef, *_ = np.linalg.lstsq(X, x, rcond=None)
    var = np.sum((x - x.mean(0)) ** 2)
    return float(np.sum((x - X @ coef) ** 2) / max(var, 1e-300)), float(var / len(x))


def spec_peaks(n):
    x = n[N // 4: 3 * N // 4]
    M = len(x)
    w = np.blackman(M)          # sidelobes below -58 dB, unlike Hanning's -31 dB
    F = np.fft.rfft((x - x.mean(0)) * w[:, None], n=8 * M, axis=0)
    P = np.sum(np.abs(F) ** 2, 1)
    k = np.fft.rfftfreq(8 * M) * 2 * np.pi
    loc = np.where((P[1:-1] > P[:-2]) & (P[1:-1] >= P[2:]))[0] + 1
    loc = [i for i in loc if P[i] > 1e-4 * P.max()]
    return sorted([float(k[i]) for i in loc], key=lambda kk: -P[np.argmin(abs(k - kk))])[:6]


rr0, var0 = resid(nbest, 0.0)
best_r = (rr0 if var0 > 1e-14 else 0.0, 0.0)
if best_r[0] > 0:
    cands = {2 * np.pi * wb / Pb} if Pb > 1 else set()
    for kk in spec_peaks(nbest):
        for m in range(1, 9):
            cands.add(kk / m)
    res_bin = 2 * np.pi / (N // 2)
    for qc in cands:
        if qc < 0.5 * res_bin:
            continue
        loc_best = (np.inf, qc)
        for q in np.linspace(qc - 2 * res_bin, qc + 2 * res_bin, 81):
            rr, _ = resid(nbest, q)
            if rr < loc_best[0]:
                loc_best = (rr, q)
        qq = loc_best[1]
        for q in np.linspace(qq - res_bin / 20, qq + res_bin / 20, 81):
            rr, _ = resid(nbest, q)
            if rr < loc_best[0]:
                loc_best = (rr, q)
        if loc_best[0] < best_r[0]:
            best_r = loc_best

# 5. single-period cells matched to the open-chain texture (long periods near onset)
e_uni = min(t[0] for t in ring if t[1] == 1)
xb = nbest[N // 4: 3 * N // 4]
var_bulk = float(np.mean(np.sum((xb - xb.mean(0)) ** 2, 1)))
e_ring2, cell2 = None, None
qf = best_r[1]
if qf and qf > 0:
    res_bin = 2 * np.pi / (N // 2)
    cand = []
    for P in range(3, 1201):
        w = int(round(qf * P / (2 * np.pi)))
        if w >= 1 and gcd(w, P) == 1 and abs(2 * np.pi * w / P - qf) < 0.5 * res_bin:
            cand.append((abs(2 * np.pi * w / P - qf), P, w))
    cand.sort()
    for _, P, w in cand[:12]:
        seg = nbest[N // 4: N // 4 + P]
        for n0 in (seg, helix(P, 2 * np.pi * w / P)):
            n0 = n0 + 1e-4 * rng.normal(size=n0.shape)
            n0 = n0 / np.linalg.norm(n0, axis=1)[:, None]
            E, _n = relax(n0, True)
            if e_ring2 is None or E / P < e_ring2:
                e_ring2, cell2 = E / P, (P, w)

res = dict(e_uniform=float(e_uni), var_bulk=var_bulk, e_ring_matched=None if e_ring2 is None else float(e_ring2),
           cell_matched=cell2, seed=rec['seed'], t=rec['t'], sc=sc, q_helix=float(qh), cone=float(al), e_helix=float(eh),
           e_ring=float(e_ring), P_ring=int(Pb), w_ring=int(wb), ncells=len(cells),
           e_bulk_best=float(ebest), best_label=lab, E_open_best=float(Ebest),
           trials=[(t[0], float(t[1]), float(t[2])) for t in trials],
           resid=float(best_r[0]), q_fit=float(best_r[1]), old_E_per_site=rec['E_per_site'],
           old_unexplained=rec['unexplained_weight'])
json.dump(res, open(out, 'w'))
np.save(out.replace('.json', '.npy'), nbest.astype(np.float32))
print(json.dumps({k: v for k, v in res.items() if k != 'trials'}))
