"""C182 stage-2 selector (rule registered in RESULT.md before the job ran).
Reads out/s_*_L48.jsonl and out/s_*_L96.jsonl, computes c_two(48->96) at W = 6, keeps points whose L = 96
truncation is below 1e-6, ranks them, and writes typhon/tasks2.txt with converged series at L = 96, 144, 192
(chi to 512) plus N +- 1 runs at L = 144 and 192 for the parity gap, for up to 6 selected points, and the two
controls at the same sizes."""
import numpy as np, json, glob, os
from analyze import Sbar
SCH = [16]*4 + [32]*4 + [64]*4 + [128]*3 + [192]*2 + [256]*2 + [384]*2 + [512]*2
rows = []
for f in sorted(glob.glob('out/s_*_L96.jsonl')):
    ser = os.path.basename(f)[:-len('_L96.jsonl')]
    f48 = f'out/{ser}_L48.jsonl'
    if not os.path.exists(f48): continue
    r96 = [json.loads(l) for l in open(f)][-1]; r48 = [json.loads(l) for l in open(f48)][-1]
    if r96['sweep'] != len(r96['cfg']['chis'])-1 or r48['sweep'] != len(r48['cfg']['chis'])-1: continue
    c = 6*(Sbar(r96['S'], 96, 6) - Sbar(r48['S'], 48, 6))/np.log(2)
    # amendment (registered in RESULT.md Q7 before any job output existed): homogeneity, |n_mid - N/L| < 0.02 at L = 96
    m = r96.get('meas')
    if m is None or abs(np.mean(m['dens'][24:72]) - r96['cfg']['N']/96) >= 0.02: continue
    rows.append((c, r96['trunc'], ser, r96['cfg']))
rows.sort(key=lambda r: -r[0])
with open('typhon/stage2_ranking.txt', 'w') as fo:
    for c, tr, ser, cfg in rows:
        fo.write(f"{ser} c_two(48->96,W=6)={c:.3f} trunc96={tr:.1e}\n")
good = [r for r in rows if r[1] < 1e-6]
sel = [r for r in good if r[0] > 1.4][:6]
if len(sel) < 6:
    sel += [r for r in good if 1.2 < r[0] <= 1.4][:6-len(sel)]
tasks = []
def add(c):
    fn = f"typhon/cfg/{c['tag']}.json"; json.dump([c], open(fn, 'w'))
    tasks.append(f"timeout 2400m python3 run.py {fn} out2 > logs/{c['tag']}.log 2>&1")
for c, tr, ser, cfg in sel:
    n = cfg['N']/96
    for L in (192, 144, 96):
        N = int(round(n*L))
        base = dict(L=L, wa=cfg['wa'], phi=cfg['phi'], chis=SCH, meas_chi=512)
        add(dict(tag=f"{ser}_L{L}", N=N, **base))
        if L > 96:
            for dN, s in ((1, 'p'), (-1, 'm')):
                add(dict(tag=f"{ser}_L{L}_N{s}1", N=N+dN, nomeas=True, **base))
for nm, extra in (('c1', dict(wa=1.0, phi=0.0, t2=0.0, V2=0.0, V3=0.0)),
                  ('c2', dict(wa=0.0, phi=0.0, V2=0.0, V3=0.0, blockade=False))):
    for L in (192, 144, 96):
        add(dict(tag=f"{nm}_n1o4_L{L}", L=L, N=L//4, chis=SCH, meas_chi=512, **extra))
open('typhon/tasks2.txt', 'w').write('\n'.join(tasks)+'\n')
print(len(sel), 'points selected;', len(tasks), 'stage-2 tasks')
