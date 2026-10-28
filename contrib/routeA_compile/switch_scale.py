"""c_max of the central-switch cell versus m and tau (precise one-body bond)."""
import sys, json, time
import numpy as np, scipy.linalg as la
from onebody2 import logical_h_precise
from onebody import blocks, max_c
m = int(sys.argv[1]); out = []
taus = [float(x) for x in sys.argv[2:]] or [0.1, 0.03, 0.01]
for tau in taus:
    t = time.time()
    p = dict(wD=-0.3, tauA=tau, tauB=tau)
    d0, h = logical_h_precise('switch', m, p)
    T, M = blocks(h, m)
    Ts = 0.5 * (T + T.T)
    ev = la.eigvalsh(Ts.real)
    gen = np.sort(la.eigvals(np.linalg.solve(Ts, M)).real)
    c, mu = max_c(T, M, nk=4000)
    r = dict(m=m, tau=tau, c_max=c, T_asym=float(abs(T - T.T).max() / abs(T).max()),
             T_eig_over_max=(ev / abs(ev).max()).tolist(), eig_Tinv_M=gen.tolist(),
             seconds=time.time() - t)
    print(json.dumps(r), flush=True); out.append(r)
json.dump(out, open(f'switch_scale_m{m}_' + '_'.join(sys.argv[2:]) + '.json', 'w'), indent=1)
