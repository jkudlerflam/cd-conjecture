"""One DMRG ground state of the open vacuum-clock model. Returns energy and the contact profile."""
import os
for _k in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_k, '2')
import time
import numpy as np
from qdmrg_sd import QDMRG
from model import Model


def ground(lengths, pins, V, J, mu, G, chis=(32, 64, 128), etol=1e-12, maxsweeps=10, log=None, t=1.0, on_chi=None):
    m = Model(lengths, pins, V, J, mu, G, t)
    dm = QDMRG(m.bonds(), m.qloc, m.N, m.config(), m.allowed)
    rng = np.random.default_rng(1)
    hist = []
    t0 = time.time()
    for chi in chis:
        Eold = None
        for sw in range(maxsweeps):
            nz = 1e-4 * 0.1 ** sw if sw < 2 else 0.0
            E = dm.sweep(chi, cutoff=1e-14, tol=1e-10, noise=nz, rng=rng)
            if log:
                log(f'chi={chi} sw={sw} E={E:.13f} trunc={dm.trunc:.1e} t={time.time()-t0:.0f}s')
            if nz == 0 and Eold is not None and abs(E - Eold) < etol:
                break
            Eold = E
        pr = dm.prob.reshape(m.N, m.G + 1, 2).sum(2)
        hist.append(dict(chi=chi, E=float(E), trunc=float(dm.trunc), maxchi=int(max(dm.chis())), sweeps=sw + 1,
                         contact=pr[:m.N - 1, 0].tolist(), pgmax=float(pr[:, m.G].max()), time=time.time() - t0))
        if on_chi:
            on_chi(hist)
    prob = dm.prob.reshape(m.N, m.G + 1, 2).sum(2)        # P(g_j = g)
    contact = prob[:m.N - 1, 0]
    return dict(E=hist[-1]['E'], hist=hist, contact=contact.tolist(), pgmax=float(prob[:, m.G].max()),
                first_order=m.first_order(contact), frozen=m.frozen(), N=m.N, V=V, J=J, G=G,
                lengths=list(lengths), pins=list(pins), time=time.time() - t0)
