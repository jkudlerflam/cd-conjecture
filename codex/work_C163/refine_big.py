"""C163 large-L refinement of single bonds (Claude Opus agent, 2026-11-08).
usage: refine_big.py NAME L [CHI]      NAME in c159_refine/ (scan vector) or 'NAME@k' = NAME displaced by 0.03 along a
random T direction (seed k).  Open-chain DMRG (dmrg_c.run) at one length L, chis up to CHI; writes big/NAME_L.json."""
import sys, os, json, time
import numpy as np
import dmrg_c as DC
name, L = sys.argv[1], int(sys.argv[2]); CHI = int(sys.argv[3]) if len(sys.argv) > 3 else 192
base = name.split('@')[0]
v = np.array(json.load(open(f'c159_refine/{base}.json'))['scan']['vec'], float); v /= np.linalg.norm(v)
if '@' in name:
    rng = np.random.default_rng(int(name.split('@')[1]))
    d = np.zeros(9); d[:6] = rng.normal(size=6); d /= np.linalg.norm(d)
    v = v + 0.03 * d; v /= np.linalg.norm(v)
os.makedirs('big', exist_ok=True)
chis = tuple(c for c in (16, 32, 64, 96, 128, 160, 192, 256) if c < CHI) + (CHI, CHI)
t0 = time.time()
r = DC.run(v, (L,), chis=chis)
r.update(name=name, v=v.tolist(), L_=L, chis=chis, time_s=time.time() - t0)
json.dump(r, open(f'big/{name}_{L}.json', 'w'), default=float)
print('done', name, L, f'{time.time() - t0:.0f}s', flush=True)
