"""C190 diagnostic: Davidson preconditioner floor sigma on one centre two-site problem."""
import sys, json, time
import numpy as np
import z2dmrg as Z
import run190 as R
cfg = json.load(open(sys.argv[1])); psi = Z.MPS.load(sys.argv[2])
mpo, gph, shift = R.build(cfg)
L = psi.L; j0 = L // 2 - 1
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
x0 = x0 + 1e-3 * np.linalg.norm(x0) * np.random.default_rng(3).standard_normal(x0.size) / np.sqrt(x0.size)
print('size', ts.size)
t0 = time.time(); ts.build_kron4(); print(f'kron4 build {time.time()-t0:.2f}s', 'e0', ts.k4e0)
for sigma in [float(a) for a in sys.argv[3:]]:
    for mmax in (32,):
        t0 = time.time()
        e, v, nmv, rn = Z.davidson(ts.matvec, x0, ts.diag, res_tol=1e-9, maxit=600, mmax=mmax, keep=6, sigma=sigma,
                                   prec=lambda r, th: ts.precond4(r, th, sigma))
        print(f'KRON sigma={sigma:g} mmax={mmax}: nmv={nmv} res={rn:.2e} e={e:.13f} {time.time()-t0:.1f}s', flush=True)
        x1 = np.random.default_rng(5).standard_normal(ts.size).astype(x0.dtype)
        e, v, nmv, rn = Z.davidson(ts.matvec, x1, ts.diag, res_tol=1e-9, maxit=600, mmax=mmax, keep=6, sigma=sigma,
                                   prec=lambda r, th: ts.precond4(r, th, sigma))
        print(f'KRON rand-start sigma={sigma:g} mmax={mmax}: nmv={nmv} res={rn:.2e} e={e:.13f} {time.time()-t0:.1f}s', flush=True)
