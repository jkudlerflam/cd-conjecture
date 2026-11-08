#!/usr/bin/env python3
"""C165 stage B selection (criteria fixed before any Typhon output was seen).
For each screen file fes/A*.json (two-site VUMPS, chi up to 64):
  N(chi) = |m_0 - m_1|/2 (staggered moment of the cell), c_last = 6 dS/dln xi between chi = 48 and 64.
  ordered   : N(64) >= 0.05 and N(64)/N(32) >= 0.97
  coexist   : ordered and c_last >= 0.8 and xi(64)/xi(48) >= 1.2           (floating + SSB flag)
  highc     : c_last >= 1.25 and xi(64)/xi(48) >= 1.2
  transition: along a line (same label, Dz sorted) the last disordered and first ordered point
Selected (priority coexist, highc, transition; at most 24) are continued to chi = 96..256 from their npz state;
the four highest-priority ones also get open-chain DMRG at L = 128 and 256."""
import json, glob, os, re, numpy as np
def stats(f):
    R = json.load(open(f))['results']
    if len(R) < 5: return None
    m = lambda r: np.array(r['mag']).real.T
    N = [np.linalg.norm(m(r)[0] - m(r)[1]) / 2 for r in R]; ch = [r['chi'] for r in R]
    S = [np.mean(r['S']) for r in R]; xi = [r['xi'] for r in R]
    i64, i48, i32 = len(R) - 1, len(R) - 2, len(R) - 3
    c_last = 6 * (S[i64] - S[i48]) / np.log(xi[i64] / xi[i48]) if xi[i64] > xi[i48] * 1.001 else 0.0
    return dict(N=N[i64], rN=N[i64] / max(N[i32], 1e-12), c_last=c_last, rxi=xi[i64] / xi[i48], chi=ch[-1])
def main():
    rows = {}
    for f in sorted(glob.glob('fes/A*.json')):
        tag = os.path.basename(f)[:-5]; s = stats(f)
        if s is None or s['chi'] < 64: continue
        s['ordered'] = s['N'] >= 0.05 and s['rN'] >= 0.97
        s['coexist'] = s['ordered'] and s['c_last'] >= 0.8 and s['rxi'] >= 1.2
        s['highc'] = s['c_last'] >= 1.25 and s['rxi'] >= 1.2
        rows[tag] = s
    lines = {}
    for tag in rows:
        lab, dz = re.match(r'(.*)_Dz([0-9.]+)$', tag).groups() if tag.startswith('A1') or tag.startswith('A3') else re.match(r'(.*)_hz([0-9.]+)$', tag).groups()
        lines.setdefault(lab, []).append((float(dz), tag))
    trans = []
    for lab, pts in lines.items():
        pts.sort()
        for (x0, t0), (x1, t1) in zip(pts, pts[1:]):
            if rows[t0]['ordered'] != rows[t1]['ordered']:
                trans += [t0, t1]; break
    pri = [t for t in rows if rows[t]['coexist']] + [t for t in rows if rows[t]['highc']] + trans
    sel = []
    for t in pri:
        if t not in sel: sel.append(t)
    sel = sel[:24]
    json.dump(dict(rows=rows, selected=sel), open('stageB_selection.json', 'w'), indent=1, default=float)
    tasks = []
    big = '96 128 160 192 256'
    for i, t in enumerate(sel):
        cmd = open('tasks_A.txt').read().split('\n')
        hit = [c for c in cmd if f'fes/{t}.json' in c]
        if not hit: continue
        line = hit[0]
        line = re.sub(r'(\.json) [0-9 ]+ --n', r'\1 ' + big + ' --n', line)
        tasks.append(line)
        if i < 4:
            kw = re.search(r'"([^"]+)"', line).group(1)
            for Lc, chi in ((128, 128), (256, 160)):
                tasks.append(f"python3 dmrgq.py 'fam:{kw}' {Lc} {chi} dm/B_{t}_L{Lc}.json")
    open('tasks_B.txt', 'w').write('\n'.join(tasks) + '\n')
    print('screened', len(rows), 'ordered', sum(r['ordered'] for r in rows.values()), 'coexist', sum(r['coexist'] for r in rows.values()),
          'highc', sum(r['highc'] for r in rows.values()), 'selected', len(sel))
if __name__ == '__main__':
    main()
