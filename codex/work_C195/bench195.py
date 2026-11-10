"""C195 benchmark at L = 512 with the offloaded storage (Claude Opus 5.5 agent, 2026-11-10).

Usage: python3 bench195.py --L 512 --chis 768,1024,1536 [--nb 3] [--nmv 100] [--out f.json]

For each chi: a random Z2-block MPS of p153c with every bond at min(chi, 2^dist), brought to mixed canonical
form around the centre bond site by site on the GPU, with every environment built on the GPU and stored on the
host, as in an offloaded sweep (so the host memory after the build is that of a production sweep at this chi).
Then nb consecutive right-moving mid-chain bond updates, each exactly the step of z2dmrg_xp.dmrg's offloaded
sweep (upload of one environment and one tensor, TwoSite, kron4 preconditioner, Davidson with a fixed number
nmv of matvecs, split with the device SVD, env_left, download of one environment and one tensor), plus the
entropy SVD of a left-moving step.  Reports seconds per bond and per part, device memory (cupy pool and
memGetInfo) and host memory (VmRSS, VmHWM).
"""
import os
import sys
import gc
import time
import json
import argparse
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import z2dmrg_xp as ZX
import models as Mo

ap = argparse.ArgumentParser()
ap.add_argument('--L', type=int, default=512)
ap.add_argument('--chis', default='768,1024,1536')
ap.add_argument('--nb', type=int, default=3)
ap.add_argument('--nmv', type=int, default=100)
ap.add_argument('--out', default=None)
a = ap.parse_args()

cp = ZX._cp
ZX.USE_GPU = True
ZX.SVD_HOST = False
ZX.OFFLOAD = True


def now():
    cp.cuda.Device().synchronize()
    return time.time()


def host_mem():
    out = {}
    for line in open('/proc/self/status'):
        if line.startswith(('VmRSS', 'VmHWM')):
            out[line.split(':')[0]] = int(line.split()[1]) / 1e6
    return out


def dev_mem():
    pool = cp.get_default_memory_pool()
    free, total = cp.cuda.runtime.memGetInfo()
    return dict(pool_used_GB=pool.used_bytes() / 1e9, pool_total_GB=pool.total_bytes() / 1e9,
                dev_used_GB=(total - free) / 1e9, dev_total_GB=total / 1e9)


def rand_tensor(Dl, Dr, n0l, n0r, rng):
    T = cp.zeros((Dl, 2, Dr), dtype=complex)
    for qa in (0, 1):
        for s in (0, 1):
            sa, sb = ZX.bsl(n0l, Dl, qa), ZX.bsl(n0r, Dr, qa ^ s)
            shp = T[sa, s, sb].shape
            if shp[0] and shp[1]:
                T[sa, s, sb] = rng.standard_normal(shp) + 1j * rng.standard_normal(shp)
    return T


def build(L, chi, j0, mpo, W):
    rng = cp.random.default_rng(5)
    D = [1] + [int(min(chi, 2 ** min(j, L - j))) for j in range(1, L)] + [1]
    n0 = [1] + [(d + 1) // 2 for d in D[1:L]] + [1]
    M = [None] * L
    psi = ZX.MPS(M, n0, 0)
    psi.M = [np.zeros((D[j], 2, D[j + 1]), dtype=complex) if False else None for j in range(L)]
    w = mpo.w
    Lenv = [None] * (L + 1)
    Renv = [None] * (L + 1)
    Lenv[0] = np.zeros((1, w, 1), dtype=complex); Lenv[0][0, 0, 0] = 1.0
    Renv[L] = np.zeros((1, w, 1), dtype=complex); Renv[L][0, 2, 0] = 1.0
    # left part: QR to the right, rescaling the pushed factor, environments on the GPU, stored on the host
    carry = rand_tensor(D[0], D[1], n0[0], n0[1], rng)
    dLe = cp.asarray(Lenv[0])
    for j in range(j0):
        psi.M[j] = carry
        psi.M[j + 1] = rand_tensor(D[j + 1], D[j + 2], n0[j + 1], n0[j + 2], rng)
        ZX._qr_move_right(psi, j)
        psi.M[j + 1] /= cp.linalg.norm(psi.M[j + 1])
        dLe = ZX.env_left(dLe, psi.M[j], W[j], psi.n0[j], psi.n0[j + 1], mpo)
        Lenv[j + 1] = cp.asnumpy(dLe)
        psi.M[j] = cp.asnumpy(psi.M[j])
        carry = psi.M[j + 1]
    left_c = carry
    carry = rand_tensor(D[L - 1], D[L], n0[L - 1], n0[L], rng)
    dRe = cp.asarray(Renv[L])
    for j in range(L - 1, j0 + 1, -1):
        psi.M[j] = carry
        psi.M[j - 1] = rand_tensor(D[j - 1], D[j], n0[j - 1], n0[j], rng)
        ZX._rq_move_left(psi, j)
        psi.M[j - 1] /= cp.linalg.norm(psi.M[j - 1])
        dRe = ZX.env_right(dRe, psi.M[j], W[j], psi.n0[j], psi.n0[j + 1], mpo)
        Renv[j] = cp.asnumpy(dRe)
        psi.M[j] = cp.asnumpy(psi.M[j])
        carry = psi.M[j - 1]
    psi.M[j0], psi.M[j0 + 1] = left_c, carry
    nrm = cp.linalg.norm(cp.tensordot(psi.M[j0], psi.M[j0 + 1], axes=([2], [0])))
    psi.M[j0] /= nrm
    # the centre pair goes to the host too, as in a sweep
    psi.M[j0] = cp.asnumpy(psi.M[j0]); psi.M[j0 + 1] = cp.asnumpy(psi.M[j0 + 1])
    del dLe, dRe, carry, left_c
    return psi, Lenv, Renv


def main():
    L = a.L
    mpo, gph, shift, info = Mo.candidate('p153', L, -19.97589125)
    W = [cp.asarray(x) for x in mpo.W]
    pr = cp.cuda.runtime.getDeviceProperties(0)
    res = dict(L=L, nmv=a.nmv, nb=a.nb, gpu=pr['name'].decode() if isinstance(pr['name'], bytes) else str(pr['name']),
               host=os.uname()[1], threads=os.environ.get('OMP_NUM_THREADS'), rows=[])
    for chi in [int(c) for c in a.chis.split(',')]:
        j0 = L // 2 - a.nb // 2 - 1
        t0 = time.time()
        psi, Lenv, Renv = build(L, chi, j0, mpo, W)
        row = dict(chi=chi, j0=j0, t_build=time.time() - t0, host_after_build=host_mem(), dev_after_build=dev_mem())
        env_bytes = sum(x.nbytes for x in Lenv + Renv if x is not None)
        mps_bytes = sum(x.nbytes for x in psi.M)
        row['host_env_GB'], row['host_mps_GB'] = env_bytes / 1e9, mps_bytes / 1e9
        print(json.dumps(row), flush=True)
        per = []
        dLe, dMj = cp.asarray(Lenv[j0]), cp.asarray(psi.M[j0])
        for k in range(a.nb + 1):                        # the first bond is a warm-up (kernel compilation)
            j = j0 + k
            tm = dict(up=0.0, setup=0.0, kron=0.0, matvec=0.0, precond=0.0, dav_other=0.0, split=0.0, env=0.0,
                      down=0.0, entropy=0.0)
            nmv = a.nmv if k else 4
            s0 = now()
            dRe, dMj1 = cp.asarray(Renv[j + 2]), cp.asarray(psi.M[j + 1])
            s1 = now(); tm['up'] = s1 - s0
            D = psi.dims
            ts = ZX.TwoSite(dLe, W[j], W[j + 1], dRe, psi.n0[j], D[j], psi.n0[j + 2], D[j + 2], mpo)
            th = cp.tensordot(dMj, dMj1, axes=([2], [0]))
            x0 = ts.pack_dense(th)
            s2 = now(); tm['setup'] = s2 - s1
            ts.build_kron4()
            s3 = now(); tm['kron'] = s3 - s2
            acc = {'mv': 0.0, 'pc': 0.0}

            def mvf(x):
                q = now(); out = ts.matvec(x); acc['mv'] += now() - q
                return out

            def pcf(r, thv):
                q = now(); out = ts.precond4(r, thv, 1e-2); acc['pc'] += now() - q
                return out
            e, v, kmv, rn = ZX.davidson(mvf, x0, None, res_tol=0.0, maxit=nmv - 1, mmax=32, keep=6, sigma=1e-3,
                                        prec=pcf)
            s4 = now()
            tm['matvec'], tm['precond'] = acc['mv'], acc['pc']
            tm['dav_other'] = (s4 - s3) - acc['mv'] - acc['pc']
            thv = ts.unpack_dense(v)
            A, B, k0, disc = ZX.split(thv, psi.n0[j], D[j], psi.n0[j + 2], D[j + 2], chi, -1.0, 'right')
            del ts, th, thv, v, x0
            s5 = now(); tm['split'] = s5 - s4
            psi.n0[j + 1] = k0
            dLe = ZX.env_left(dLe, A, W[j], psi.n0[j], psi.n0[j + 1], mpo)
            s6 = now(); tm['env'] = s6 - s5
            psi.M[j] = cp.asnumpy(A)
            Lenv[j + 1] = cp.asnumpy(dLe)
            Renv[j + 2] = None
            s7 = now(); tm['down'] = s7 - s6
            ZX.schmidt_probs_center(A, psi.n0[j], D[j], k0, A.shape[2])
            s8 = now(); tm['entropy'] = s8 - s7
            tm['total'] = s8 - s0
            tm['nmv'] = kmv
            tm['size'] = int(2 * D[j] * D[j + 2])
            psi.M[j + 1] = B
            dMj = B
            del A, dRe
            if k:
                per.append(tm)
                print('  bond', j, json.dumps({kk: round(vv, 4) if isinstance(vv, float) else vv for kk, vv in tm.items()}),
                      flush=True)
        row['dev_peak'] = dev_mem()
        row['host_peak'] = host_mem()
        for key in per[0]:
            row[key + '_mean'] = float(np.mean([p[key] for p in per]))
        row['t_per_matvec'] = float(np.mean([p['matvec'] / p['nmv'] for p in per]))
        row['bonds'] = per
        res['rows'].append(row)
        print(json.dumps({k: v for k, v in row.items() if k != 'bonds'}), flush=True)
        del psi, Lenv, Renv, dLe, dMj, B
        gc.collect()
        cp.get_default_memory_pool().free_all_blocks()
        if a.out:
            json.dump(res, open(a.out, 'w'), indent=1)


if __name__ == '__main__':
    main()
