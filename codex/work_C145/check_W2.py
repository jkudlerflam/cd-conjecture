"""C145: test of the proved (W2) bound r_Q(x) >= (1 - 4 exp(-2 b d)) r_0(x), d = Q - |x|,
and of the intermediate reflection bound r_0 - r_Q <= r_0(2Q-x) + r_0(2Q+x). Floating point only."""
import json
import numpy as np
from continuum import Kernel, Restricted, grading

rows = []
for eps in [0.5, 1e-1, 1e-2, 1e-4, 1e-6]:
    eta = 2*np.arcsin(np.sqrt(eps/2)); K = Kernel(eta)
    g = np.pi-eta; b = np.pi/g
    r0 = lambda x: 1/(2*g*np.cosh(b*x))
    for Q in [0.3, 1.0, 1.5, 2.0, 4.0, 6.0]:
        N = 900
        R = Restricted(K, Q, N=N, c=grading(K, Q, N))
        y = R.y; d = Q-y
        refl = r0(2*Q-y)+r0(2*Q+y)
        viol_refl = float(np.max((r0(y)-R.r)-refl))          # must be <= 0
        lb = 1-4*np.exp(-2*b*d)
        ratio = R.r/r0(y)
        viol_W2 = float(np.max(lb-ratio))                     # must be <= 0
        m1 = d >= 1
        rows.append(dict(eps=eps, eta=eta, Q=Q, max_refl_violation=viol_refl,
                         max_W2_violation=viol_W2,
                         min_ratio_d_ge_1=float(ratio[m1].min()) if m1.any() else None,
                         max_ratio=float(ratio.max())))
        print(json.dumps(rows[-1]), flush=True)
json.dump(rows, open('check_W2.json', 'w'), indent=1)
