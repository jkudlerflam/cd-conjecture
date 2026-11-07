import json, numpy as np
from continuum import Kernel, Restricted, grading
out = []
for eps, Q in [(1e-6, 4.0), (1e-4, 6.0)]:
    eta = 2*np.arcsin(np.sqrt(eps/2)); K = Kernel(eta); g = np.pi-eta; b = np.pi/g
    for N in [900, 1800, 2700]:
        R = Restricted(K, Q, N=N, c=grading(K, Q, N)); y = R.y
        r0 = 1/(2*g*np.cosh(b*y)); ratio = R.r/r0
        lb = 1-4*np.exp(-2*b*(Q-y))
        out.append(dict(eps=eps, Q=Q, N=N, max_W2_violation=float(np.max(lb-ratio)), max_ratio=float(ratio.max()),
                        ratio_at_0=float(ratio[0])))
        print(json.dumps(out[-1]), flush=True)
json.dump(out, open('check_W2_refine.json', 'w'), indent=1)
