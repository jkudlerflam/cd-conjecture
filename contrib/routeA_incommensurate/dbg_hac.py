import sys, numpy as np, time
sys.argv = ['x', 'fes', '--out', '/tmp/x.json']
import bigchain as B, cell as CE, vumps as V, dressing as DR
from ed_check import eigbasis
from scipy.sparse.linalg import LinearOperator, eigsh
lam = 0.1
cell = CE.build(CE.INC); cb = eigbasis(cell); B.EVALS = cb['e'].copy()
D, N = cb['D'], cb['N']; J = cell['pred']['J0'] * lam ** 2; B.DELTA = 10 * J
W = [B.mpo(cb, lam)] * 2
heff = CE.expected_h2(cell); heff = (heff - np.trace(heff) / 16 * np.eye(16)) / cell['pred']['J0']
p0, _ = V.run([V.mpo_from_term(heff, 2, d=N)] * 2, 2, [3], d=N, tol=1e-10, maxit=400, mode='damped', log=lambda r: None)
emb = [np.concatenate([x, np.zeros((x.shape[0], D - N, x.shape[2]))], 1) for x in p0.AL]
As = [DR.apply_mpo(x, DR.dressing_mpo(cb, lam)) for x in emb]
AL, AR, C = DR.canonical(As)
psi = V.UMPS(AL, AR, C)
print('E', V.energy_AL(psi.AL, W), 'chi', AL[0].shape[0])
Ls, eL, _ = V.left_envs(psi.AL, psi.C, W, 1e-13)
Rs, eR, _ = V.right_envs(psi.AR, psi.C, W, 1e-13)
print('eL eR per cell', eL / 2, eR / 2)
k = 0
f = lambda v: B.hac(Ls[k], W[k], Rs[1], v)
x = psi.AC[k] / np.linalg.norm(psi.AC[k])
th0 = np.real(np.vdot(x, f(x)))
print('<AC|H|AC>', th0)
t0 = time.time()
e_d, v_d = B.lowest_patched(f, x, 1e-10)
print('davidson', e_d, 'overlap', abs(np.vdot(v_d, x)), time.time() - t0, B.STATS)
n = x.size
op = LinearOperator((n, n), matvec=lambda v: f(v.reshape(x.shape)).ravel(), dtype=complex)
t0 = time.time()
w, U = eigsh(op, k=4, which='SA', v0=x.ravel(), tol=1e-12, ncv=60, maxiter=20000)
o = np.argsort(w)
print('eigsh', w[o], 'overlap', abs(np.vdot(U[:, o[0]], x.ravel())), time.time() - t0)
# C problem
fc = lambda c: V.hc(Ls[k], Rs[k], c)
c0 = psi.C[k] / np.linalg.norm(psi.C[k])
print('<C|H|C>', np.real(np.vdot(c0, fc(c0))))
e_c, v_c = B.lowest_patched(fc, c0, 1e-10)
print('C lowest', e_c, 'overlap', abs(np.vdot(v_c, c0)))
