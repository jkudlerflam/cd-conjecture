"""C168 follow-up (Claude Opus 5.5, 2026-11-09): open-chain DMRG of one stage-A point at one length and seed,
chi up to 192, to test whether the L = 128 readings of job 17602943 are converged ground states.
Usage: python3 dmrg_task.py IDX L SEED"""
import json, sys, numpy as np
import dmrg_c_c155 as dc
i, L, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
v = json.load(open('pts.json'))[i]
o = dc.run(v, Ls=(L,), chis=(16, 32, 64, 128, 160, 192, 192), seed=seed)
json.dump(o, open(f'dm/p{i}_L{L}_s{seed}.json', 'w'))
x = o['L'][L]
print(i, L, seed, 'E/L', x['E'] / L, 'Smid', x['Smid'], 'c_prof', x['c_prof'], 'trunc', x['trunc'], flush=True)
