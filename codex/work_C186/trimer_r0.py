"""C186: three-body bound states of the number-conserving cusp model past the pair boundary (r = 0, v > v_c).
Infinite line, N = 3 exactly, total momentum K scanned; threshold = min(3 eps_min, eps_min + lowest pair) = 3 eps_min
(no pair bound for v > v_c). Hard wall at extent R; bound trimers show E - 3 eps_min < 0 stable in R."""
import numpy as np, sys
from mixed_fewbody import Model, bare_valleys
u = float(sys.argv[1]); vs = [float(x) for x in sys.argv[2].split(',')]
for v in vs:
    m = Model(u, np.pi/2, v, r=0.0)
    ka, kb = bare_valleys(m)[:2]
    m.mu = float(m.eps(ka))
    for R in (80, 140):
        best = (np.inf, None)
        for K in np.concatenate([np.linspace(-np.pi, np.pi, 25, endpoint=False), [np.pi + ka, np.pi + kb, 3*ka, 3*kb]]):
            e = m.lowest(K, 3, R, Nmin=3)[0]
            if e < best[0]: best = (e, K)
        print(f"u={u} v={v} R={R}: lowest N=3 level - 3 eps_min = {best[0]:+.6f} at K = {best[1]:.4f}", flush=True)
