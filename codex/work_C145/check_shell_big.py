"""Large-L runs of check_shell.case. Usage: python3 check_shell_big.py eps L1,L2,..."""
import json, sys, time
import numpy as np
from check_shell import case
from check_attractive import bounds
eps = float(sys.argv[1]); Ls = [int(a) for a in sys.argv[2].split(',')]
eta = 2*np.arcsin(np.sqrt(eps/2)); M1 = 4*bounds(eta)['CD']; rows = []
for L in Ls:
    for n in [L//4, 3*L//8, L//2-2, L//2]:
        t = time.time()
        try:
            c = case(eta, L, n, 2.0, M1)
        except Exception as exc:
            c = dict(error=str(exc))
        c.update(eps=eps, eta=eta, L=L, n=n, etaL=eta*L, M1=M1, secs=time.time()-t)
        rows.append(c); print(json.dumps(c), flush=True)
json.dump(rows, open('check_shell_big_%s.json' % sys.argv[1], 'w'), indent=1)
