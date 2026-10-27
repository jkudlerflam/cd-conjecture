#!/usr/bin/env python3
"""Lower-memory momentum-block ED for L = 16 (spin-3/2 chart), built on ced.Sector.

Differences from ced.Sector.blocks: one CSR per bond component, accumulated site by site with
int32 indices (instead of one COO list over all sites), and a matrix-free sum over components
in the Lanczos step.  Spin labels come from comparing the Sz = 0 and Sz = 1 spectra (a level of
the Sz = 0 block that is absent from the Sz = 1 block at the same momentum is a singlet), so no
Casimir vectors on the full basis are needed.  The ground-state half-chain entropy uses
ced.Sector.entropies at the single cut L/2.

usage: ed_lean.py --L 16 --points pts.json --momenta 0 1 4 8 --nlev 6 --sz 0 --out f.jsonl
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'counterexample_search'))
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh, LinearOperator
from ced import Sector
import models


def block_component(sec, h, k):
    L, d = sec.L, sec.d
    k = k % L
    allowed = (k * sec.rperiod.astype(np.int64)) % L == 0
    ids = np.flatnonzero(allowed)
    reps = sec.reps[ids]
    inverse = np.full(len(sec.reps), -1, dtype=np.int64)
    inverse[ids] = np.arange(len(ids))
    n = len(ids)
    real = (k == 0 or 2 * k == L)
    dt = np.float64 if real else np.complex128
    M = sparse.csr_matrix((n, n), dtype=dt)
    for site in range(L):
        j = (site + 1) % L
        pi, pj = d ** site, d ** j
        pair = ((reps // pi) % d) * d + (reps // pj) % d
        rows, cols, data = [], [], []
        for inp in range(d * d):
            ci = np.flatnonzero(pair == inp)
            if not len(ci):
                continue
            a, b = divmod(inp, d)
            for op in np.flatnonzero(np.abs(h[:, inp]) > 1e-13):
                c, e = divmod(int(op), d)
                new = reps[ci] + (c - a) * pi + (e - b) * pj
                bi = np.searchsorted(sec.basis, new)
                bi[bi >= len(sec.basis)] = 0
                assert np.all(sec.basis[bi] == new)
                ri = inverse[sec.rid[bi]]
                use = ri >= 0
                v = h[op, inp] * np.sqrt(sec.rperiod[ids[ci[use]]] / sec.period[bi[use]])
                if not real:
                    v = v * np.exp(2j * np.pi * k * sec.shift[bi[use]] / L)
                elif k:
                    v = v * np.cos(np.pi * sec.shift[bi[use]])
                rows.append(ri[use].astype(np.int32)); cols.append(ci[use].astype(np.int32)); data.append(v.astype(dt))
        if data:
            M = M + sparse.coo_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
                                      shape=(n, n)).tocsr()
        del rows, cols, data
    M.sum_duplicates()
    M.eliminate_zeros()
    return ids, M


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--L', type=int, required=True)
    p.add_argument('--points', required=True)
    p.add_argument('--momenta', type=int, nargs='+', required=True)
    p.add_argument('--nlev', type=int, default=6)
    p.add_argument('--sz', type=int, default=0, help='2 Sz of the sector (0 or 2 for even L)')
    p.add_argument('--entropy', action='store_true')
    p.add_argument('--out', required=True)
    a = p.parse_args()
    m = models.spin_s(1.5)
    P0, P1, P2, P3 = m['P']
    comps = [P1, P2, P3]
    pts = json.load(open(a.points))
    t0 = time.time()
    sec = Sector(4, a.L, comps, m['charges'], [a.sz], cache=False)
    print('sector dim %d reps %d  %.0fs' % (len(sec.basis), len(sec.reps), time.time() - t0), flush=True)
    for k in a.momenta:
        t1 = time.time()
        ids, Ms = None, []
        for h in comps:
            ids, M = block_component(sec, h, k)
            Ms.append(M)
        n = len(ids)
        print('k=%d dim %d nnz %s built %.0fs' % (k, n, [M.nnz for M in Ms], time.time() - t1), flush=True)
        for (t, eta) in pts:
            c = [np.sin(t) + eta, np.cos(t), np.sin(t) - eta]
            t2 = time.time()
            dt = Ms[0].dtype
            op = LinearOperator((n, n), matvec=lambda v: c[0] * (Ms[0] @ v) + c[1] * (Ms[1] @ v) + c[2] * (Ms[2] @ v), dtype=dt)
            rng = np.random.default_rng(1234 + k)
            v0 = rng.normal(size=n).astype(dt)
            e, v = eigsh(op, k=a.nlev, which='SA', tol=1e-10, v0=v0, ncv=max(40, 3 * a.nlev + 10), maxiter=50000)
            ix = np.argsort(e); e, v = e[ix], v[:, ix]
            res = float(np.max(np.linalg.norm(np.stack([op.matvec(v[:, i]) - e[i] * v[:, i] for i in range(len(e))], 1), axis=0)))
            rec = dict(L=a.L, sz=a.sz, k=k, theta=t, eta=eta, a=[0.0] + c, E=e.tolist(), residual=res, dim=n, seconds=time.time() - t2)
            if a.entropy:
                psi = sec.expand(k, v[:, 0].astype(complex))
                rec['S_half'] = sec.entropies(psi, cuts=[a.L // 2])[0]
                del psi
            with open(a.out, 'a') as f:
                f.write(json.dumps(rec) + '\n')
            print('  pt (%.4f, %.3f) E0=%.8f res=%.1e %.0fs' % (t / np.pi, eta, e[0], res, time.time() - t2), flush=True)
        del Ms


if __name__ == '__main__':
    main()
