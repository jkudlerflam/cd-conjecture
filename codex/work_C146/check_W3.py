"""C146: check of the proved edge bound (W3'') and its two ingredients. Floating point only. J=1.
Proved: eta Z_Q(Q)^2 >= pi Q / (4 (1 + 16 Q/(3 pi^2)))  (PROOF Theorem 1).
Ingredients: Z_Q(Q) nondecreasing in Q; W(Q) = int_{-Q}^{Q} Z_Q has W'(Q) = 2 Z_Q(Q)^2."""
import json
import numpy as np
from continuum import Kernel, Restricted, grading
rows = []
for eps in [0.5, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
    eta = 2*np.arcsin(np.sqrt(eps/2)); K = Kernel(eta)
    Qs = np.geomspace(0.003, 4.0, 23)
    ZQ, W = [], []
    for Q in Qs:
        N = 900
        R = Restricted(K, Q, N=N, c=grading(K, Q, N))
        ZQ.append(R.Z[-1]); W.append(2*R.integrate(R.Z))
    ZQ, W = np.array(ZQ), np.array(W)
    bound = np.pi*Qs/(4*(1+16*Qs/(3*np.pi**2)))
    ratio = eta*ZQ**2/bound
    # W'(Q) = 2 Z(Q)^2 by finite differences in log Q
    dW = np.gradient(W, Qs)
    for i, Q in enumerate(Qs):
        row = dict(eps=eps, eta=eta, Q=float(Q), Q_over_eta=float(Q/eta), ZQ=float(ZQ[i]),
                   etaZ2=float(eta*ZQ[i]**2), bound=float(bound[i]), ratio_to_bound=float(ratio[i]),
                   cZ_true=float(eta*ZQ[i]**2/np.tanh(2*Q)),
                   dW_over_2Z2=float(dW[i]/(2*ZQ[i]**2)), W_over_2QZ2=float(W[i]/(2*Q*ZQ[i]**2)))
        rows.append(row)
    print(eps, 'min ratio to bound', ratio.min(), 'Z monotone', bool(np.all(np.diff(ZQ) >= -1e-9*ZQ[1:])),
          'max |dW/2Z^2-1| (interior)', float(np.max(np.abs(dW[2:-2]/(2*ZQ[2:-2]**2)-1))),
          'max W/(2QZ^2)', float(np.max(W/(2*Qs*ZQ**2))), 'min cZ', float(np.min(eta*ZQ**2/np.tanh(2*Qs))), flush=True)
json.dump(rows, open('check_W3.json', 'w'), indent=1)
