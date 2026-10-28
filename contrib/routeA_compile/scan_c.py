"""Random parameter scan: largest Fermi-point count (c_max) of the
exact order-lambda^2 logical one-body chain, over the chemical potential."""
import sys, json, time
import numpy as np
import scipy.sparse.linalg as sla
from onebody import logical_h, blocks, max_c, prepared
from model import assemble, DEFAULT

def carrier_gap(n, p):
    ch, P = prepared(n)
    H0, V = assemble(P, p)
    mask = np.ones(ch.dim, bool); mask[ch.code] = False
    Q = np.where(mask)[0]
    return float(sla.eigsh(H0[Q][:, Q], k=1, which='SA', return_eigenvectors=False)[0])

def sample(rng):
    p = dict(DEFAULT)
    p.update(eT=rng.uniform(0.5, 2), eH=rng.uniform(0.5, 2), eD=rng.uniform(0.5, 2),
             wT=rng.uniform(-.7, .7), wH=rng.uniform(-.7, .7), wD=rng.uniform(-.7, .7),
             tauA=rng.uniform(.05, .5), tauB=rng.uniform(.05, .5),
             kappa=rng.uniform(.2, 3), phiT=rng.uniform(0, np.pi))
    return p

if __name__ == '__main__':
    n, trials, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rng = np.random.default_rng(seed)
    prepared(n)
    hist, best = {}, None
    t0 = time.time()
    for t in range(trials):
        p = sample(rng)
        g = carrier_gap(n, p)
        if g < 0.05:
            continue
        E0, h = logical_h(n, p)
        T, M = blocks(h, n)
        c, mu = max_c(T, M, nk=3000)
        hist[c] = hist.get(c, 0) + 1
        if best is None or c > best[0]:
            best = (c, p, g)
    print(json.dumps(dict(n=n, trials=trials, hist={str(k): v for k, v in hist.items()},
                          best_c=best[0], best_params=best[1], best_gap=best[2],
                          seconds=time.time()-t0)))
