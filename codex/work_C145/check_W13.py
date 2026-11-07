"""C145: the reduced continuum statements behind (W1), (W3). Floating point only.
(W1') v_Q = eps_Q'/r_Q nondecreasing on [0,Q]  =>  Var(v_Q) = 2 v_Q(Q) <= 4 a tanh(2Q) (proved implication).
(W3') eta Z_Q(Q)^2 >= c_Z tanh(2Q)            =>  e_bulk'' <= a tanh(2Q)/Z_Q(Q)^2 <= (pi/c_Z) eta^2 (proved implication).
Also reports Var(v_Q)/eta, e_bulk''/eta^2 and v_Q(Q)/(2a tanh 2Q) (proved <= 1)."""
import json
import numpy as np
from continuum import Kernel, Restricted, grading
rows = []
for eps in [0.5, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
    eta = 2*np.arcsin(np.sqrt(eps/2)); K = Kernel(eta); a = np.pi*np.sin(eta)
    for Q in [0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 2.0, 3.0]:
        N = 900
        R = Restricted(K, Q, N=N, c=grading(K, Q, N))
        y = R.y
        de = np.gradient(R.epsQ, y); de[0] = 0.0
        v = de/R.r
        dv = np.diff(v)
        vQ = float(v[-1])
        # one-sided cubic slope at Q, as in Restricted.curvature
        c = np.polyfit(y[-4:]-y[-1], R.epsQ[-4:], 3); vQ1 = c[2]/R.r[-1]
        row = dict(eps=eps, eta=eta, Q=Q, m=float(R.m),
                   min_dv_rel=float(dv.min()/max(abs(vQ1), 1e-300)),
                   VarV_over_eta=float(2*np.abs(dv).sum()/eta),
                   vQ_over_2a_tanh=float(vQ1/(2*a*np.tanh(2*Q))),
                   cZ=float(eta*R.Z[-1]**2/np.tanh(2*Q)), ZQ=float(R.Z[-1]),
                   ebulk2_over_eta2=float(R.curvature()/eta**2))
        rows.append(row); print(json.dumps(row), flush=True)
json.dump(rows, open('check_W13.json', 'w'), indent=1)
