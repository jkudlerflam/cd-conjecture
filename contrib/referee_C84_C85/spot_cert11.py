"""(6) Spot checks of the C84 n=11 and C85 n=13 certificates (inactive endpoint).
Adapted from contrib/referee_C80_C82/spot_cert.py: the carrier has p = n+1 data
vertices and p+1 express/marker vertices (2p+1 in all); only D_0..D_(n-1) receive
conversion sources; D_(p-1) is inactive but closes the loop.
(a) The integer full-loop matrices F, M against the referee's floating-point
evaluation of the folded-graph formula (C80 PROOF eq. (5), C75 eq. (3)) at
D = 1000, t_d = beta 2^-eta, up to one common factor.
(b) Bloch scan of the whitened matrices on a fine grid: inertia changes,
roots and margins compared with the reported ones.
Usage: python3 spot_cert.py <scratch work folder> <n> <eta_power> <tag>"""
import sys, json
import numpy as np
from scipy.linalg import solve

W, n, eta, tag = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
D = 1000; p = n + 1; S = 2 * p + 1
d0, d1 = 1, D + 4
for _ in range(2, p + 2):
    d0, d1 = d1, (D + 4) * d1 - d0
beta = 1 / d1
td = beta * 2.0 ** -eta
Hg = np.zeros((S, S), complex)
for i in range(S):
    Hg[i, i] = 4 + (D if i <= p else 0)
for i in range(S - 1):
    Hg[i, i + 1] = Hg[i + 1, i] = -1
Dj = lambda j: p + 1 + j
Hg[0, Dj(p - 1)] = -1j * td; Hg[Dj(p - 1), 0] = 1j * td
I = np.eye(S); A = np.kron(Hg, I) + np.kron(I, Hg)
e = np.zeros(S * S, complex); e[Dj(0)] = 1          # |M, D_0>
x = solve(A, e)
U = np.zeros((S * S, n), complex); V = np.zeros_like(U)
for j in range(n):
    m = np.zeros((S, S)); m[Dj(j), Dj(j)] = 1
    U[:, j] = np.kron(m, I) @ x; V[:, j] = np.kron(I, m) @ x
RU = solve(A, U); RV = solve(A, V)
C = U.conj().T @ RV
Am = U.conj().T @ RU + V.conj().T @ RV
T = -C.real
mu = -1j * Am.imag
a = np.load(f'{W}/full_loop_n{n}.npz', allow_pickle=True)
bits = int(a['bits']); sc = 2 ** (2 * bits)
Fi = np.array([[int(v) / sc for v in row] for row in a['F']])
Mi = np.array([[int(v) / sc for v in row] for row in a['M']])
ratioT = T / (-Fi)
mask = abs(Mi) > 0
ratioM = (mu[mask] / (1j * Mi[mask])).real
out = dict(n=n, eta_power=eta,
           T_ratio_over_beta2=float(ratioT.mean() / beta ** 2),
           T_ratio_rel_spread=float(np.ptp(ratioT) / abs(ratioT.mean())),
           M_ratio_over_beta2=float(ratioM.mean() / beta ** 2),
           M_ratio_rel_spread=float(np.ptp(ratioM) / abs(ratioM.mean())),
           T_entry_dynamic_range=float(abs(T).max() / abs(T).min()),
           onsite_to_hopping_norm=float(np.linalg.norm(mu, 2) / np.linalg.norm(T, 2)))
print(out, flush=True)
w = np.load(f'{W}/full_loop_whitened_n{n}.npz', allow_pickle=True)
Tw, muw, fw = w['T'].astype(complex), w['mu'].astype(complex), w['field'].astype(complex)
B = muw + fw
ks = np.linspace(-np.pi, np.pi, 200001)
assert np.allclose(B, B.conj().T)
ins = []
for c0 in range(0, len(ks), 5000):
    kk = ks[c0:c0 + 5000]
    h = B[None] + Tw[None] * np.exp(1j * kk)[:, None, None] + Tw.conj().T[None] * np.exp(-1j * kk)[:, None, None]
    ins.append(np.sum(np.linalg.eigvalsh(h) < 0, axis=1))
ins = np.concatenate(ins)
crossings = [(float(ks[i]), int(ins[i - 1]), int(ins[i])) for i in range(1, len(ks)) if ins[i] != ins[i - 1]]
roots = np.array([c[0] for c in crossings])
rep = json.load(open(f'{W}/full_loop_bands_n{n}.json'))
rr = np.array(rep['roots'])


def margins(r):
    r = np.sort(np.mod(r, 2 * np.pi))
    d2 = min(min(np.diff(r)), 2 * np.pi - r[-1] + r[0])
    s = []
    for i in range(len(r)):
        for j in range(i, len(r)):
            s.append(np.mod(r[i] + r[j], 2 * np.pi))
    s = np.sort(np.array(s)); diffs = np.diff(s)
    d4 = min(min(diffs[diffs > 1e-12]), 2 * np.pi - s[-1] + s[0])
    return float(d2), float(d4), int(np.sum(diffs <= 1e-12)), len(s)


out.update(grid_crossings=len(crossings), reported_roots=len(rr),
           max_root_mismatch=float(np.max(abs(np.sort(roots) - np.sort(rr)))) if len(roots) == len(rr) else None,
           grid_step=float(ks[1] - ks[0]),
           inertia_up_then_down=[c[2] for c in crossings],
           recomputed_margins_from_reported_roots=margins(rr),
           reported_delta2=rep.get('delta2'), reported_delta4=rep.get('delta4'))
print(out, flush=True)
import os
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'spot_{tag}.json'), 'w'), indent=1)
