"""C160 Typhon summary (Claude Opus agent, 2026-11-08). Run in /scratch/lustre/jkudlerflam/c160:  python3 summarize.py > summary.txt"""
import json, glob, collections
import numpy as np
recs = collections.defaultdict(dict)
for f in glob.glob('dmrg/*.json'):
    r = json.load(open(f)); recs[r['label']][r['L']] = r
print('DMRG, K-stratum CSL family near h_F (v = (-1,-1,-1+kappa eps^2/2,0,0,-eps,-h eps^2,0,0))')
print('label       eps   kappa    h     L   c_prof  c_sm   Smid    trunc     dE_last   c_two(prev L)')
for lab in sorted(recs):
    Ls = sorted(recs[lab]); prev = None
    for L in Ls:
        r = recs[lab][L]
        ctwo = '' if prev is None else f'{6 * (r["Smid"] - prev[1]) / np.log(L / prev[0]):+.3f}'
        print(f'{lab:10s} {r["eps"]:.2f} {r["kappa"]:+.3f} {r["h"]:.2f} {L:4d}  {r["c_prof"]:+.3f} {r["c_sm"]:+.3f} {r["Smid"]:.4f} '
              f'{r["hist"][-1][2]:.1e} {r["dE_last"]:+.1e}  {ctwo}')
        prev = (L, r['Smid'])
print()
rows = [json.loads(l) for f in glob.glob('texture/*.jsonl') for l in open(f)]
tex = [r for r in rows if r.get('kind') == 'textured']
print(f'textures: {len(rows)} bonds, {len(tex)} textured, {len(rows) - len(tex)} uniform')
if tex:
    u = np.array([r['unexplained_weight'] for r in tex])
    print(f'unexplained (non-harmonic) peak weight: max {u.max():.2e}, median {np.median(u):.2e}, bonds above 1e-2: {(u > 1e-2).sum()}')
    for r in sorted(tex, key=lambda r: -r['unexplained_weight'])[:5]:
        print(f'  seed {r["seed"]} t {r["t"]}: q1 = {r["q1"]:.4f}, unexplained {r["unexplained_weight"]:.2e} at k = {r["first_unexplained"]}, peaks {r["peaks"][:4]}')
