"""C190 diagnostic: spectrum and Davidson convergence of one two-site problem at the chain centre.
Usage: python3 diag_local.py <config.json> <ckpt.npz> [site]"""
import sys
import json
import time
import numpy as np
import scipy.sparse.linalg as spla
import z2dmrg as Z
import run190 as R

cfg = json.load(open(sys.argv[1]))
psi = Z.MPS.load(sys.argv[2])
mpo, gph, shift = R.build(cfg)
L = psi.L
j0 = int(sys.argv[3]) if len(sys.argv) > 3 else L // 2 - 1
dt = np.result_type(mpo.dtype, psi.M[0].dtype)
psi.M = [m.astype(dt) for m in psi.M]
Z.right_canonicalize(psi)
w = mpo.w
Le = np.zeros((1, w, 1), dtype=dt); Le[0, 0, 0] = 1
for j in range(j0):
    Z._qr_move_right(psi, j)
    Le = Z.env_left(Le, psi.M[j], mpo.W[j], psi.n0[j], psi.n0[j + 1], mpo)
Re = np.zeros((1, w, 1), dtype=dt); Re[0, 2, 0] = 1
for j in range(L - 1, j0 + 1, -1):
    Re = Z.env_right(Re, psi.M[j], mpo.W[j], psi.n0[j], psi.n0[j + 1], mpo)
D = psi.dims
ts = Z.TwoSite(Le, mpo.W[j0], mpo.W[j0 + 1], Re, psi.n0[j0], D[j0], psi.n0[j0 + 2], D[j0 + 2], mpo)
th = np.tensordot(psi.M[j0], psi.M[j0 + 1], axes=([2], [0]))
x0 = ts.pack_dense(th)
print('size', ts.size, 'D', D[j0], D[j0 + 2])
op = spla.LinearOperator((ts.size, ts.size), matvec=ts.matvec, dtype=dt)
t0 = time.time()
ev = spla.eigsh(op, k=8, which='SA', tol=1e-12, return_eigenvectors=False, ncv=60, maxiter=20000)
ev = np.sort(ev)
print('lowest local eigenvalues (relative to the lowest):', ev[0], (ev - ev[0]).tolist(), f'{time.time() - t0:.1f}s')
print('t2 units:', ((ev - ev[0]) / 0.0125).tolist())
print('diag range', ts.diag.min(), ts.diag.max())
for mmax, keep in ((32, 4), (64, 8), (128, 16)):
    for start in ('ckpt', 'rand'):
        x = x0 if start == 'ckpt' else np.random.default_rng(1).standard_normal(ts.size).astype(dt)
        t0 = time.time()
        e, v, nmv, rn = Z.davidson(ts.matvec, x, ts.diag, res_tol=1e-10, maxit=2000, mmax=mmax, keep=keep)
        print(f'davidson mmax={mmax} keep={keep} start={start}: nmv={nmv} res={rn:.2e} e-e0={e - ev[0]:.2e} '
              f'{time.time() - t0:.1f}s')
