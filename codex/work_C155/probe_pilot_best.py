"""Probe the pilot's best screening point (c_obj 2.40 at rings 12/16) at larger sizes."""
import json, time, numpy as np
import est, dmrg_c
r = json.loads(open('pilot/task_0.jsonl').readline())
v = np.array(r['best_v'])
t = time.time()
e = est.estimate(est.bond_v(v), Ls=(16, 20), nev=3)
print('ED 16/20', {k: (np.round(e[k], 3) if not isinstance(e[k], list) else e[k]) for k in ['c_S', 'c_fit', 'c_cas', 'v', 'm0', 'Lgap', 'ndeg', 'z1', 'z2', 'crit']}, f'{time.time()-t:.0f}s', flush=True)
o = dmrg_c.run(v, Ls=(32, 64), chis=(16, 32, 64, 96, 96))
print('DMRG', {L: (round(d['c_prof'], 3), round(d['c_sm'], 3), f"{d['trunc']:.1e}", d['dE_last']) for L, d in o['L'].items()}, 'c_two', round(o['c_two'], 3))
S = o['L'][64]['S']; print('S(l) L=64:', np.round(S, 3).tolist())
