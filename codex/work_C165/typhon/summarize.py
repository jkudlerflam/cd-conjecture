#!/usr/bin/env python3
"""C165 digest of the Typhon job. Pre-fixed c = 3/2 criteria (spec C165, Task 3):
 DMRG c_prof and c_sm > 1.3 at L = 128 and not falling to L = 256, truncation < 1e-7; FES pair slope > 1.3 between
 the two largest chi with chi >= 256; staggered moment vanishing continuously at the point; incommensurate peak on
 both sides. Floating + SSB flag: see stageB.py."""
import json, glob, os, numpy as np
def fes_line(f):
    d = json.load(open(f)); R = d['results']
    if not R: return f'{os.path.basename(f)}: empty'
    S = np.array([np.mean(r['S']) for r in R]); xi = np.array([r['xi'] for r in R]); ch = [r['chi'] for r in R]
    cp = 6 * np.diff(S) / np.diff(np.log(xi))
    ls = 6 * np.polyfit(np.log(xi[-3:]), S[-3:], 1)[0] if len(S) >= 3 else np.nan
    N = ''
    if 'mag' in R[-1]:
        Ns = [np.linalg.norm(np.array(r['mag']).real.T[0] - np.array(r['mag']).real.T[1]) / 2 for r in R]
        N = ' N=' + ','.join(f'{x:.3f}' for x in Ns[-4:])
    return (f"{os.path.basename(f)}: chi_max={ch[-1]} xi={xi[-1]:.1f} c_pair(last4)={','.join(f'{x:.3f}' for x in cp[-4:])} "
            f"c_ls3={ls:.3f}{N}")
print('== FES calibrations and long VUMPS'); [print(fes_line(f)) for f in sorted(glob.glob('fes_cal/*.json')) + sorted(glob.glob('fes/nfloat*.json')) + sorted(glob.glob('fes/Q2big*.json'))]
print('== DMRG')
for f in sorted(glob.glob('dm/*.json')):
    d = json.load(open(f))
    print(f"{os.path.basename(f)}: L={d['L']} c_prof={d['c_prof']:.3f} c_sm={d['c_sm']:.3f} trunc={d['trunc']:.1e} dE={d['dE_last']:.1e} "
          + (f"stag={np.round(d['stag'],4).tolist()} Sk={ {k: [round(x,3) for x in v] for k, v in d['Sk_peak'].items()} }" if 'stag' in d else ''))
print('== stage B selection'); print(open('logs/stageB.out').read() if os.path.exists('logs/stageB.out') else 'none')
if os.path.exists('stageB_selection.json'):
    sel = json.load(open('stageB_selection.json'))
    for t in sel['selected']:
        f = f'fes/{t}.json'
        print(' ', json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in sel['rows'][t].items()}), '|', fes_line(f))
    rows = sel['rows']
    print('== all screen points with c_last > 1.1 or coexist flag')
    for t, r in sorted(rows.items()):
        if r['c_last'] > 1.1 or r['coexist']:
            print(' ', t, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
