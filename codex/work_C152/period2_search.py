"""Adversarial search for c_free >= 3/2 in period-2 NN quadratic chains: maximize the number of near-zero crossings,
using as smooth objective the count plus a bonus for small |E| minima of each band (local minima of min_k |E_b(k)|)."""
import numpy as np
from period2_mf import bdg, model, count, K
rng = np.random.default_rng(11)
def obj(p):
    E = bdg(*model(p)); s = np.sign(E)
    n = np.sum(s != np.roll(s, 1, axis=0)) / 4.0
    # soft term: number of local minima of |E| over k for the two lowest |E| bands, weighted by closeness to 0
    a = np.abs(E); soft = 0.0
    for b in range(4):
        x = a[:, b]; m = (x < np.roll(x, 1)) & (x < np.roll(x, -1))
        soft += np.sum(np.exp(-x[m] / 0.05))
    return n + 0.05 * soft
best_all = 0
for start in range(60):
    p = rng.normal(size=10); f = obj(p); step = 0.3
    for it in range(250):
        q = p + step * rng.normal(size=10); g = obj(q)
        if g >= f: p, f = q, g
        if it % 50 == 49: step *= 0.6
    c = count(p); best_all = max(best_all, c)
print('max c_free found over 60 adversarial starts:', best_all)
