"""C163 twisted-curve check (Claude Opus agent, 2026-11-08): J = -R_z(phi) + eps X, f = eps b (c163lib.twisted_specs).
Open-chain DMRG at L = 64, 128 (chi to 128) and ring K_S(12), K_S(16).  usage: twisted.py TAG  -> twisted/TAG.json"""
import sys, os, json, time
import numpy as np
import c163lib as C, dmrg_c as DC, tk
from kest import k_est
tag = sys.argv[1]
s = [x for x in C.twisted_specs() if x['tag'] == tag][0]
J, f = C.twisted_bond(s); v = C.J_to_v(J, f)
os.makedirs('twisted', exist_ok=True)
t0 = time.time()
out = dict(spec=s, v=v.tolist(), parent={k: (val.tolist() if hasattr(val, 'tolist') else val) for k, val in tk.spiral_parent(J, f).items()})
for L in (12, 16):
    try:
        out[f'K_S{L}'] = k_est(v, L)[0]
    except Exception as ex:
        out[f'K_S{L}'] = str(ex)
json.dump(out, open(f'twisted/{tag}.json', 'w'), default=float)
out['dmrg'] = DC.run(v, (64, 128), chis=(16, 32, 64, 96, 128, 128))
out['time_s'] = time.time() - t0
json.dump(out, open(f'twisted/{tag}.json', 'w'), default=float)
print('done', tag, f'{time.time() - t0:.0f}s', flush=True)
