import sys, numpy as np, time
sys.argv = ['x', 'fes', '--out', '/tmp/x.json']
import bigchain as B, cell as CE, vumps as V, dressing as DR
from ed_check import eigbasis
lam = 0.1
cell = CE.build(CE.INC); cb = eigbasis(cell); B.EVALS = cb['e'].copy()
D, N = cb['D'], cb['N']
W = [B.mpo(cb, lam)] * 2
heff = CE.expected_h2(cell); heff = (heff - np.trace(heff) / 16 * np.eye(16)) / cell['pred']['J0']
p0, _ = V.run([V.mpo_from_term(heff, 2, d=N)] * 2, 2, [3], d=N, tol=1e-10, maxit=400, mode='damped', log=lambda r: None)
emb = [np.concatenate([x, np.zeros((x.shape[0], D - N, x.shape[2]))], 1) for x in p0.AL]
print('E code-only', V.energy_AL(emb, W))
Wd = DR.dressing_mpo(cb, lam)
for scale in (1.0,):
    Wd2 = Wd.copy(); Wd2[0, 1:] *= scale
    As = [DR.apply_mpo(x, Wd2) for x in emb]
    t0 = time.time()
    AL, AR, C = DR.canonical(As)
    print('scale', scale, 'chi', AL[0].shape[0], 'E(AL)', V.energy_AL(AL, W), 'sv', np.round(np.diag(C[0]).real[:12], 4), 'nonzero', int(np.sum(np.diag(C[0]).real > 1e-12)), 't', round(time.time() - t0, 1), flush=True)
    # consistency AL C = C AR
    print('   |AL C - C AR|', max(np.linalg.norm(np.tensordot(AL[k], C[(k + 1) % 2], ([2], [0])) - np.tensordot(C[k], AR[k], ([1], [0]))) for k in range(2)))
print('pred 2nd order', -0.7 * lam ** 2)
