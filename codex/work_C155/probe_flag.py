"""Refine the one gapless uniform-sample point with c_obj > 1.1 (local_uniform_flag.json): ED 16/20 and DMRG 32/64."""
import json, numpy as np, est, dmrg_c
from crit import labels
for v in json.load(open('local_uniform_flag.json')):
    e = est.estimate(est.bond_v(v), Ls=(16, 20), nev=3)
    print('v', np.round(v, 3).tolist(), labels(v))
    print('ED16/20', {k: np.round(e[k], 3).tolist() for k in ('c_S', 'c_fit', 'c_cas', 'crit', 'm0', 'ndeg')}, flush=True)
    o = dmrg_c.run(v, Ls=(32, 64), chis=(16, 32, 64, 96, 96))
    print('DMRG', {L: (round(d['c_prof'], 3), round(d['c_sm'], 3), f"{d['trunc']:.1e}") for L, d in o['L'].items()}, 'c_two', round(o['c_two'], 3))
