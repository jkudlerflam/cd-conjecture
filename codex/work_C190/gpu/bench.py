"""C190g profile and benchmark (Claude Opus 5.5 agent, 2026-11-09).

Usage: python3 bench.py <numpy|cupy> --L 256 --chis 256,384,512,768 [--nb 3] [--nmv 80] [--cpu 1] [--svd-dev]

Builds the p153 MPO at the collaborator's field, a random Z2-block MPS with every bond at min(chi, 2^dist) in
mixed canonical form around the centre bond, and the exact environments for it.  Then performs nb consecutive
right-moving two-site updates at the centre, exactly as dmrg() does them (TwoSite, kron4 preconditioner, Davidson,
split, env_left), with a fixed number nmv of matvecs per Davidson solve (res_tol = 0), and reports the time of
each part.  The matvec count does not depend on the state, so a random state is fine for timing (production
solves use 50 to 120 matvecs).  With --cpu k the first k bonds are repeated with the original z2dmrg.py on the
same environments and tensors (host, OMP threads as set), the matvecs of the two codes are compared on the same
vector, and the SVD of a production-size block is compared between cusolver and LAPACK.
"""
import os
import sys
import time
import json
import argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(os.path.join(HERE, 'z2dmrg.py')):
    sys.path.insert(0, os.path.dirname(HERE))
import z2dmrg as Z0
import z2dmrg_xp as ZX
import models as Mo

ap = argparse.ArgumentParser()
ap.add_argument('backend')
ap.add_argument('--L', type=int, default=256)
ap.add_argument('--chis', default='256,384,512,768')
ap.add_argument('--nb', type=int, default=3)
ap.add_argument('--nmv', type=int, default=80)
ap.add_argument('--cpu', type=int, default=0)
ap.add_argument('--svd-dev', action='store_true')
ap.add_argument('--out', default=None)
a = ap.parse_args()

ZX.USE_GPU = a.backend == 'cupy'
ZX.SVD_HOST = not a.svd_dev
xm = ZX._cp if ZX.USE_GPU else np
sync = ZX.sync


def now():
    sync()
    return time.time()


def mem():
    if not ZX.USE_GPU:
        return None
    pool = ZX._cp.get_default_memory_pool()
    free, total = ZX._cp.cuda.runtime.memGetInfo()
    return dict(pool_used_GB=pool.used_bytes() / 1e9, pool_total_GB=pool.total_bytes() / 1e9,
                dev_used_GB=(total - free) / 1e9, dev_total_GB=total / 1e9)


def build_state(L, chi, j0, mpo, seed=5):
    rng = np.random.default_rng(seed)
    D = [1] + [int(min(chi, 2 ** min(j, L - j))) for j in range(1, L)] + [1]
    n0 = [1] + [(d + 1) // 2 for d in D[1:L]] + [1]
    for j in range(1, L):
        if D[j] == 1:
            n0[j] = 1
    M = []
    for j in range(L):
        T = np.zeros((D[j], 2, D[j + 1]), dtype=complex)
        for qa in (0, 1):
            for s in (0, 1):
                qb = qa ^ s
                sa, sb = ZX.bsl(n0[j], D[j], qa), ZX.bsl(n0[j + 1], D[j + 1], qb)
                shp = T[sa, s, sb].shape
                T[sa, s, sb] = rng.standard_normal(shp) + 1j * rng.standard_normal(shp)
        M.append(ZX.to_mod(T, xm))
    psi = ZX.MPS(M, n0, 0)
    t0 = now()
    for j in range(0, j0):                     # rescale the pushed factor each step (no overflow at L = 256)
        ZX._qr_move_right(psi, j)
        psi.M[j + 1] /= xm.linalg.norm(psi.M[j + 1])
    for j in range(L - 1, j0 + 1, -1):
        ZX._rq_move_left(psi, j)
        psi.M[j - 1] /= xm.linalg.norm(psi.M[j - 1])
    nrm = xm.linalg.norm(xm.tensordot(psi.M[j0], psi.M[j0 + 1], axes=([2], [0])))
    psi.M[j0] /= nrm
    t_can = now() - t0
    W = [ZX.to_mod(x, xm) for x in mpo.W]
    w = mpo.w
    Lenv = [None] * (L + 1)
    Renv = [None] * (L + 1)
    Lenv[0] = xm.zeros((1, w, 1), dtype=complex); Lenv[0][0, 0, 0] = 1.0
    Renv[L] = xm.zeros((1, w, 1), dtype=complex); Renv[L][0, 2, 0] = 1.0
    t0 = now()
    for j in range(j0):
        Lenv[j + 1] = ZX.env_left(Lenv[j], psi.M[j], W[j], psi.n0[j], psi.n0[j + 1], mpo)
    for j in range(L - 1, j0 + 1, -1):
        Renv[j] = ZX.env_right(Renv[j + 1], psi.M[j], W[j], psi.n0[j], psi.n0[j + 1], mpo)
    t_env = now() - t0
    return psi, W, Lenv, Renv, dict(t_canon=t_can, t_envs=t_env, n_env=L)


def one_bond(mod, psi_M, n0, dims, j, W, Le, Re, mpo, chi, nmv, timers):
    """One right-moving update at bond j with module mod (ZX or Z0); returns A, B, k0, Lnew, v."""
    t0 = now()
    ts = mod.TwoSite(Le, W[j], W[j + 1], Re, n0[j], dims[j], n0[j + 2], dims[j + 2], mpo)
    th = (xm if mod is ZX else np).tensordot(psi_M[j], psi_M[j + 1], axes=([2], [0]))
    x0 = ts.pack_dense(th)
    t1 = now(); timers['setup'] += t1 - t0
    ts.build_kron4()
    t2 = now(); timers['kron_eigh'] += t2 - t1
    tm = {'mv': 0.0, 'pc': 0.0}

    def mvf(x):
        s = now(); out = ts.matvec(x); tm['mv'] += now() - s
        return out

    def pcf(r, thv):
        s = now(); out = ts.precond4(r, thv, 1e-2); tm['pc'] += now() - s
        return out
    e, v, k, rn = mod.davidson(mvf, x0, None if mod is ZX else ts.diag, res_tol=0.0, maxit=nmv - 1, mmax=32,
                               keep=6, sigma=1e-3, prec=pcf)
    t3 = now()
    timers['matvec'] += tm['mv']; timers['precond'] += tm['pc']
    timers['davidson_other'] += (t3 - t2) - tm['mv'] - tm['pc']
    timers['nmv'] += k
    thv = ts.unpack_dense(v)
    # cutoff -1 keeps chi states (a random environment gives a nearly product local ground state, and the
    # timing must be done at the full bond dimension)
    A, B, k0, disc = mod.split(thv, n0[j], dims[j], n0[j + 2], dims[j + 2], chi, -1.0, 'right')
    t4 = now(); timers['split_svd'] += t4 - t3
    Ln = mod.env_left(Le, A, W[j], n0[j], k0, mpo)
    t5 = now(); timers['env'] += t5 - t4
    pr = mod.schmidt_probs_center(A, n0[j], dims[j], k0, A.shape[2])
    timers['entropy_svd'] += now() - t5
    timers['total'] += now() - t0
    return A, B, k0, Ln, ts, x0


def newt():
    return dict(setup=0.0, kron_eigh=0.0, matvec=0.0, precond=0.0, davidson_other=0.0, split_svd=0.0, env=0.0,
                entropy_svd=0.0, total=0.0, nmv=0)


def svd_compare(B):
    out = {}
    Bh = ZX.asnp(B)
    t0 = time.time(); U, s, Vh = np.linalg.svd(Bh, full_matrices=False); out['t_host'] = time.time() - t0
    out['rec_host'] = float(np.linalg.norm((U * s) @ Vh - Bh) / np.linalg.norm(Bh))
    if ZX.USE_GPU:
        cp = ZX._cp
        Ud, sd, Vd = cp.linalg.svd(B, full_matrices=False); sync()
        t0 = now(); Ud, sd, Vd = cp.linalg.svd(B, full_matrices=False); out['t_dev'] = now() - t0
        out['rec_dev'] = float(cp.linalg.norm((Ud * sd) @ Vd - B) / cp.linalg.norm(B))
        out['orth_dev'] = float(cp.max(cp.abs(Ud.conj().T @ Ud - cp.eye(Ud.shape[1]))))
        out['ds_max_rel'] = float(np.max(np.abs(cp.asnumpy(sd) - s)) / s[0])
        e1 = cp.asnumpy(cp.linalg.eigvalsh(B @ B.conj().T)); sync()
        t0 = now(); ev, Ue = cp.linalg.eigh(B @ B.conj().T); out['t_eigh_dev'] = now() - t0
        t0 = time.time(); np.linalg.eigh(Bh @ Bh.conj().T); out['t_eigh_host'] = time.time() - t0
    out['shape'] = list(B.shape)
    return out


def main():
    L = a.L
    mpo, gph, shift, info = Mo.candidate('p153', L, -19.97589125)
    res = dict(backend=a.backend, L=L, nmv=a.nmv, nb=a.nb, svd_host=ZX.SVD_HOST,
               threads=os.environ.get('OMP_NUM_THREADS'), rows=[])
    if ZX.USE_GPU:
        cp = ZX._cp
        props = cp.cuda.runtime.getDeviceProperties(0)
        res['gpu'] = props['name'].decode() if isinstance(props['name'], bytes) else str(props['name'])
    for chi in [int(c) for c in a.chis.split(',')]:
        j0 = L // 2 - a.nb // 2 - 1
        t0 = time.time()
        psi, W, Lenv, Renv, info_b = build_state(L, chi, j0, mpo)
        info_b['t_build'] = time.time() - t0
        row = dict(chi=chi, j0=j0, build=info_b, mem_after_build=mem())
        tm = newt()
        n0, dims = list(psi.n0), psi.dims
        # warm-up (cupy kernel compilation and plan caches) on a copy, not timed
        if ZX.USE_GPU:
            one_bond(ZX, psi.M, n0, dims, j0, W, Lenv[j0], Renv[j0 + 2], mpo, chi, 4, newt())
        per = []
        cpu_inputs = []
        for k in range(a.nb):
            j = j0 + k
            if k < a.cpu:
                cpu_inputs.append((j, [ZX.asnp(psi.M[j]).copy(), ZX.asnp(psi.M[j + 1]).copy()],
                                   list(psi.n0), list(psi.dims), ZX.asnp(Lenv[j]), ZX.asnp(Renv[j + 2])))
            t1 = newt()
            A, B, k0, Ln, ts, x0 = one_bond(ZX, psi.M, psi.n0, psi.dims, j, W, Lenv[j], Renv[j + 2], mpo, chi,
                                           a.nmv, t1)
            if k == 0:
                row['svd'] = svd_compare(ZX._theta_block(ts.unpack_dense(x0), psi.n0[j], psi.dims[j],
                                                          psi.n0[j + 2], psi.dims[j + 2], 0))
                if a.cpu:
                    xv = x0 / xm.linalg.norm(x0)
                    row['mv_dev_vs_host'] = None
                    row['_xv'] = ZX.asnp(xv)
                    row['_yv'] = ZX.asnp(ts.matvec(xv))
            per.append(dict(j=j, Da=psi.dims[j], Dc=psi.dims[j + 2], size=ts.size, **t1))
            psi.M[j], psi.M[j + 1] = A, B
            psi.n0[j + 1] = k0
            Lenv[j + 1] = Ln
            del ts
        row['mem_peak'] = mem()
        if ZX.USE_GPU:
            row['mem_pool_peak_GB'] = ZX._cp.get_default_memory_pool().total_bytes() / 1e9
        row['bonds'] = per
        for key in ('total', 'matvec', 'precond', 'davidson_other', 'kron_eigh', 'split_svd', 'env', 'entropy_svd',
                    'setup'):
            row[key + '_mean'] = float(np.mean([p[key] for p in per]))
        row['t_per_matvec'] = float(np.mean([p['matvec'] / p['nmv'] for p in per]))
        # CPU comparison with the original code on the same inputs
        if a.cpu:
            Wh = mpo.W
            cper = []
            for (j, (Mj, Mj1), n0c, dimsc, Le, Re) in cpu_inputs:
                Mloc = {j: Mj, j + 1: Mj1}
                t1 = newt()
                s0 = time.time()
                one_bond(Z0, Mloc, n0c, dimsc, j, Wh, Le, Re, mpo, chi, a.nmv, t1)
                cper.append(dict(j=j, **t1))
                if j == j0:
                    ts0 = Z0.TwoSite(Le, Wh[j], Wh[j + 1], Re, n0c[j], dimsc[j], n0c[j + 2], dimsc[j + 2], mpo)
                    yh = ts0.matvec(row['_xv'])
                    row['mv_dev_vs_host_rel'] = float(np.linalg.norm(yh - row['_yv']) / np.linalg.norm(yh))
            row['cpu_bonds'] = cper
            for key in ('total', 'matvec', 'precond', 'davidson_other', 'kron_eigh', 'split_svd', 'env',
                        'entropy_svd', 'setup'):
                row['cpu_' + key + '_mean'] = float(np.mean([p[key] for p in cper]))
            row['cpu_t_per_matvec'] = float(np.mean([p['matvec'] / p['nmv'] for p in cper]))
            row['speedup_bond'] = row['cpu_total_mean'] / row['total_mean']
            row['speedup_matvec'] = row['cpu_t_per_matvec'] / row['t_per_matvec']
        row.pop('_xv', None); row.pop('_yv', None)
        res['rows'].append(row)
        print(json.dumps({k: v for k, v in row.items() if k not in ('bonds', 'cpu_bonds')}), flush=True)
        for p in per:
            print('  dev bond', json.dumps(p), flush=True)
        for p in row.get('cpu_bonds', []):
            print('  cpu bond', json.dumps(p), flush=True)
        del psi, W, Lenv, Renv
        if ZX.USE_GPU:
            ZX._cp.get_default_memory_pool().free_all_blocks()
        if a.out:
            json.dump(res, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
