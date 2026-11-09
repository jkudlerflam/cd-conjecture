"""C168: analysis of C155 stage A (12,000 uniform points on S^8, rings 12/16). Claude Opus 5.5, 2026-11-08."""
import json, glob, numpy as np
R = []
for fn in sorted(glob.glob('c155_stageA/uniform_*.jsonl')):
    for line in open(fn):
        R.append(json.loads(line))
N = len(R)
crit = np.array([r['est']['crit'] for r in R]); cobj = np.array([r['est']['c_obj'] for r in R])
cS = np.array([r['est']['c_S'] for r in R]); ccas = np.array([r['est']['c_cas'] for r in R])
ndeg = np.array([r['est']['ndeg'][-1] for r in R])
print('N', N)
for t in (0.05, 0.1, 0.15, 0.25):
    g = crit < t
    print(f'crit<{t}: {g.sum()} ({g.mean():.4f})  c_obj>1.1: {(g & (cobj > 1.1)).sum()}  c_obj>1.3: {(g & (cobj > 1.3)).sum()}  c_S>1.3: {(g & (cS > 1.3)).sum()}')
g = crit < 0.15
q = np.quantile(cobj[g], [0.1, 0.25, 0.5, 0.75, 0.9])
print('gapless c_obj quantiles 10/25/50/75/90', np.round(q, 3))
print('gapless c_S quantiles', np.round(np.quantile(cS[g], [0.1, 0.25, 0.5, 0.75, 0.9]), 3))
flag = np.where(g & (cobj > 1.1))[0]
out = []
for i in flag:
    r = R[i]; e = r['est']
    out.append(dict(idx=int(i), v=r['v'], c_obj=e['c_obj'], c_S=e['c_S'], c_cas=e['c_cas'], crit=e['crit'], m0=e['m0']))
    print(f"{i:6d} c_obj {e['c_obj']:.3f} c_S {e['c_S']:.3f} c_cas {e['c_cas']:.3f} crit {e['crit']:.3f} m0 {e['m0']} v {np.round(r['v'],3).tolist()}")
json.dump(out, open('stageA_flags.json', 'w'), indent=1)
# fraction of gapless points with c_obj in the c=1 window
print('gapless with 0.75<c_obj<1.1:', ((cobj > 0.75) & (cobj < 1.1) & g).sum())
