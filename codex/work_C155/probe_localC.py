"""Refine the best points of the local XZ and XZB restarts (local_search/task_10, task_14): ED 16/20 and DMRG L = 32, 64."""
import json, sys, numpy as np, est, dmrg_c
t = int(sys.argv[1])
r = json.loads(open(f'local_search/task_{t}.jsonl').readline())
v = r['best_v']
e = est.estimate(est.bond_v(v), Ls=(16, 20), nev=3)
print(t, r['seed'], 'ED16/20', {k: np.round(e[k], 3).tolist() for k in ('c_S', 'c_fit', 'c_cas', 'crit', 'm0', 'ndeg', 'gap1', 'v')}, flush=True)
o = dmrg_c.run(v, Ls=(32, 64), chis=(16, 32, 64, 64, 64, 64))
print(t, 'DMRG', {L: (round(d['c_prof'], 3), round(d['c_sm'], 3), f"{d['trunc']:.1e}", f"{d['dE_last']:.1e}") for L, d in o['L'].items()}, 'c_two', round(o['c_two'], 3), flush=True)
print(t, 'S(l) L=64', np.round(o['L'][64]['S'], 2).tolist())
