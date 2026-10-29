"""(4) Spot checks of the certified n=15 numbers.
(a) The integer full-loop matrices F, M (C75 full_loop_n15.npz) against the
referee's floating-point evaluation of PROOF eq. (3) at D=1000,
t_d = beta 2^-84, Delta - E_v = 4.
(b) Band crossings recomputed on a fine momentum grid from C75's whitened
matrices, compared with the reported roots and margins.
Argument: path to a scratch copy of codex/work_C75."""
import sys, json
import numpy as np
from scipy.linalg import solve

W = sys.argv[1]
n = 15; D = 1000; S = 2 * n + 1
d0, d1 = 1, D + 4
for _ in range(2, n + 2): d0, d1 = d1, (D + 4) * d1 - d0
beta = 1 / d1  # float, ~1e-48
td = beta * 2.0 ** -84
Hg = np.zeros((S, S), complex)
for i in range(S): Hg[i, i] = 4 + (D if i <= n else 0)
for i in range(S - 1): Hg[i, i + 1] = Hg[i + 1, i] = -1
Dj = lambda j: n + 1 + j
Hg[0, Dj(n - 1)] = -1j * td; Hg[Dj(n - 1), 0] = 1j * td
I = np.eye(S); A = np.kron(Hg, I) + np.kron(I, Hg)
e = np.zeros(S * S, complex); e[Dj(0)] = 1
x = solve(A, e)
U = np.zeros((S * S, n), complex); V = np.zeros_like(U)
for j in range(n):
    m = np.zeros((S, S)); m[Dj(j), Dj(j)] = 1
    U[:, j] = np.kron(m, I) @ x; V[:, j] = np.kron(I, m) @ x
RU = solve(A, U); RV = solve(A, V)
C = U.conj().T @ RV
Am = U.conj().T @ RU + V.conj().T @ RV
T = -C.real            # units 2 alpha^4 (g^2 coefficient), up to the stated constant
mu = -1j * Am.imag
a = np.load(f'{W}/full_loop_n15.npz', allow_pickle=True)
bits = int(a['bits']); sc = 2 ** (2 * bits)
Fi = np.array([[int(v) / sc for v in row] for row in a['F']])
Mi = np.array([[int(v) / sc for v in row] for row in a['M']])
ratioT = T / (-Fi)
ratioM = (mu / (1j * Mi)).real
mask = abs(Mi) > 0
out = dict(T_ratio_mean=float(ratioT.mean()), T_ratio_rel_spread=float(np.ptp(ratioT) / abs(ratioT.mean())),
           beta2=beta ** 2, T_ratio_over_beta2=float(ratioT.mean() / beta ** 2),
           M_ratio_mean=float(ratioM[mask].mean()), M_ratio_rel_spread=float(np.ptp(ratioM[mask]) / abs(ratioM[mask].mean())),
           M_ratio_over_beta2=float(ratioM[mask].mean() / beta ** 2),
           T_entry_dynamic_range=float(abs(T).max() / abs(T).min()))
print(out)
# (b) bands from C75's whitened matrices
w = np.load(f'{W}/full_loop_whitened_n15.npz', allow_pickle=True)
Tw, muw, fw = w['T'].astype(complex), w['mu'].astype(complex), w['field'].astype(complex)
B = muw + fw
ks = np.linspace(-np.pi, np.pi, 400001)
ev_prev = None; crossings = []
for k in ks:
    h = B + Tw * np.exp(1j * k) + Tw.T * np.exp(-1j * k)
    s = int(np.sum(np.linalg.eigvalsh(h) < 0))
    if ev_prev is not None and s != ev_prev: crossings.append((float(k), ev_prev, s))
    ev_prev = s
roots = np.array([c[0] for c in crossings])
rep = json.load(open(f'{W}/full_loop_bands_n15.json'))
rr = np.array(rep['roots'])
def margins(r):
    r = np.sort(np.mod(r, 2 * np.pi))
    d2 = min(min(np.diff(r)), 2 * np.pi - r[-1] + r[0])
    s = []
    for i in range(len(r)):
        for j in range(i, len(r)): s.append(np.mod(r[i] + r[j], 2 * np.pi))
    s = np.sort(np.array(s)); diffs = np.diff(s)
    d4 = min(min(diffs[diffs > 1e-12]), 2 * np.pi - s[-1] + s[0])
    return float(d2), float(d4), int(np.sum(diffs <= 1e-12))
out.update(grid_crossings=len(crossings), max_root_mismatch=float(np.max(abs(np.sort(roots) - np.sort(rr)))) if len(roots) == len(rr) else None,
           grid_step=float(ks[1] - ks[0]), reported_delta2=rep['delta2'], reported_delta4=rep['delta4'],
           recomputed_margins_from_reported_roots=margins(rr),
           inertia_sequence=[c[1:] for c in crossings][:3])
print(out)
json.dump(out, open('spot_n15.json', 'w'), indent=1)
