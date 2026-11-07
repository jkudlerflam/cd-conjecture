"""C138 continuum diagnostics (W1)-(W3), uniformity in eta.

For eta (epsilon=10^-k) and Fermi boundaries Q, report
  ratio(u)   = r_Q(Q-u)/r0(Q-u) at distances u from the Fermi point   (W2)
  Var_v      = total variation of v_Q = eps_Q'/r_Q on [-Q,Q], over eta   (W1)
  curv       = e_bulk''(m(Q))/eta^2                                    (W3)
  ZQ2eta     = Z_Q(Q)^2 eta,  Z0eta = Z_Q(0) eta
Usage: python3 cont_diag.py -> cont_diag.json
"""
import json
from pathlib import Path
import numpy as np
from continuum import Kernel, Restricted

rows = []
for k in [2, 3, 4, 5, 6]:
    eps = 10.0**(-k)
    eta = 2*np.arcsin(np.sqrt(eps/2))
    K = Kernel(eta)
    g = np.pi-eta
    b = np.pi/g
    r0 = lambda x: 1/(2*g*np.cosh(b*x))
    for Q in [0.02, 0.1, 0.3, 1.0, 2.0, 4.0, 8.0]:
        R = Restricted(K, Q, N=900)
        y = R.y
        ratio = {}
        for u in [0.0, eta, 10*eta, 0.1, 0.5, 1.0, 2.0, 3.0, 4.0]:
            if u <= Q:
                xx = Q-u
                ratio[f'{u:.4g}'] = float(R.interp(R.r, xx)/r0(xx))
        # v = eps'/r by centered differences on the mesh
        de = np.gradient(R.epsQ, y)
        v = de/R.r
        Var_v = 2*np.sum(np.abs(np.diff(v)))   # both halves
        rows.append(dict(k=k, eta=eta, Q=Q, m=R.m, ratio=ratio,
                         min_ratio=float(np.min(R.r/r0(y))),
                         Var_v_over_eta=float(Var_v/eta),
                         sup_v_over_a=float(np.max(np.abs(v))/K.a),
                         curv=float(R.curvature()/eta**2),
                         ZQ2eta=float(R.Z[-1]**2*eta), Z0eta=float(R.Z[0]*eta)))
        print(json.dumps(rows[-1]), flush=True)
Path('cont_diag.json').write_text(json.dumps(rows, indent=1))
