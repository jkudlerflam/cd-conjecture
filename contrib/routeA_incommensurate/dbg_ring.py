import numpy as np, itertools
import cell as CE
from ed_check import eigbasis
cell = CE.build(CE.INC); cb = eigbasis(cell)
D, N, L = cb['D'], cb['N'], 3
e = cb['e']
def V(x):
    x = x.reshape([D] * L + [-1]); y = np.zeros_like(x)
    for j in range(L):
        j1 = (j + 1) % L
        for c in range(3):
            t = np.moveaxis(np.tensordot(cb['FR'][c], x, ([1], [j])), 0, j)
            t = np.moveaxis(np.tensordot(cb['FL'][c], t, ([1], [j1])), 0, j1)
            y += t
    return y.reshape(D ** L, -1)
diag = sum(np.meshgrid(*([e] * L), indexing='ij')).ravel()
idx = list(itertools.product(range(N), repeat=L))
flat = np.array([np.ravel_multi_index(i, [D] * L) for i in idx])
X = np.zeros((D ** L, len(flat)), complex); X[flat, np.arange(len(flat))] = 1
VX = V(X)
den = np.where(diag > 1e-9, diag, np.inf)
H2 = -(VX / den[:, None]).conj().T @ VX            # -P V Q H0^-1 Q V P   (lam = 1)
H2 = 0.5 * (H2 + H2.conj().T)
# ring from the predicted bond
h2p = CE.expected_h2(cell); hb = h2p.reshape(N, N, N, N)
Il = np.eye(N ** L).reshape([N] * L + [N ** L]); Hl = np.zeros((N ** L, N ** L), complex)
for j in range(L):
    j1 = (j + 1) % L
    t = np.tensordot(hb, Il, ([2, 3], [j, j1]))
    Hl += np.moveaxis(t, [0, 1], [j, j1]).reshape(N ** L, N ** L)
J0 = cell['pred']['J0']
print('||H2 - Hl||/J0 =', np.abs(H2 - Hl).max() / J0)
w1 = np.linalg.eigvalsh(H2); w2 = np.linalg.eigvalsh(Hl)
print('spectra diff/J0', np.abs(w1 - w2).max() / J0, 'rel', np.abs((w1 - w1[0]) - (w2 - w2[0])).max() / J0)
# open chain check (no wrap)
# --- full ring: Rayleigh quotient of the first-order dressed second-order ground state, lam = 0.05
lam = 0.05
w2, U2 = np.linalg.eigh(Hl)
phi = np.zeros(D ** L, complex); phi[flat] = U2[:, 0]
v = V(phi[:, None])[:, 0]
psi = phi - lam * v / np.where(diag > 1e-9, diag, np.inf)
Hpsi = diag * psi + lam * V(psi[:, None])[:, 0]
print('dressed Rayleigh quotient', np.vdot(psi, Hpsi).real / np.vdot(psi, psi).real, ' 2nd-order E0', lam ** 2 * w2[0])
import json
d = json.load(open('data/bigring_INC_0.05.json'))
print('run_ring lowest', min(min(s['E']) for s in d['results'][0]['sectors']))
# --- replicate run_ring's Davidson with this V, sector Q = +1
from scipy import linalg
udiag = np.real(np.diag(cb['U']))
par = np.ones([D] * L)
for j in range(L):
    sh = [1] * L; sh[j] = D
    par = par * udiag.reshape(sh)
par = par.ravel()
pl = np.ones([N] * L)
for j in range(L):
    sh = [1] * L; sh[j] = N
    pl = pl * udiag[:N].reshape(sh)
pl = pl.ravel()
Q = 1; k = 32
sl = np.where(pl == Q)[0]
wl, Ul = linalg.eigh((Hl * lam ** 2)[np.ix_(sl, sl)])
sec = np.where(par == Q)[0]
Qm = np.zeros((D ** L, k), complex)
idx = np.array(np.unravel_index(sl, [N] * L))
fl = np.ravel_multi_index(tuple(idx), [D] * L)
Qm[fl, :] = Ul[:, :k]
Qm = Qm[sec]
dsec = diag[sec]
def mvs(v):
    full = np.zeros((D ** L, v.shape[1]), complex); full[sec] = v
    return (diag[:, None] * full + lam * V(full))[sec]
for it in range(30):
    HQ = mvs(Qm); Hs = Qm.conj().T @ HQ
    w, s = linalg.eigh(0.5 * (Hs + Hs.conj().T)); Qm = Qm @ s; HQ = HQ @ s
    R = HQ - Qm * w[None, :]; res = np.linalg.norm(R, axis=0).max()
    print(it, 'res %.2e' % res, 'E0 %.8f' % w[0])
    if res < 1e-10: break
    T = R / np.where(dsec[:, None] > 0.5, dsec[:, None] - w[None, :], 1.0 + lam)
    for _ in range(2): T = T - Qm @ (Qm.conj().T @ T)
    Qb, _ = np.linalg.qr(np.hstack([Qm, T])); HQb = mvs(Qb); Hs = Qb.conj().T @ HQb
    w2, s2 = linalg.eigh(0.5 * (Hs + Hs.conj().T)); Qm = Qb @ s2[:, :k]
print('sector 2nd order E0', wl[0])
