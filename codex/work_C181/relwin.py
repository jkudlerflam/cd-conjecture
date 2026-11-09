"""C181: two-size c with a window proportional to L, S averaged over cuts l/L in [1/2 - f, 1/2 + f] after removing
the conformal shape: c = 6 [<S - (c/6)... >] is solved as c = 6 (Sb2 - Sb1) / (<x>_2 - <x>_1), x = ln[(2L/pi) sin(pi l/L)].
Usage: python3 relwin.py f tagL1 tagL2 (tags in local/dmrg or a given dir via env DIR)"""
import json, sys, os, numpy as np
f = float(sys.argv[1]); d = os.environ.get('DIR', 'local/dmrg')
vals = []
for t in sys.argv[2:4]:
    r = json.loads(open(os.path.join(d, t + '.jsonl')).read().strip().split('\n')[-1]); L = r['cfg']['L']
    S = np.array(r['S']); l = np.arange(1, L); sel = np.abs(l / L - 0.5) <= f
    x = np.log((2 * L / np.pi) * np.sin(np.pi * l / L))
    vals.append((S[sel].mean(), x[sel].mean()))
print(sys.argv[2], sys.argv[3], f, round(6 * (vals[1][0] - vals[0][0]) / (vals[1][1] - vals[0][1]), 3))
