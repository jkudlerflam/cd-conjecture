"""C164: digest of the Typhon job (or of local/ pilots). usage: summarize.py [root]"""
import sys, json, glob, os
import numpy as np

root = sys.argv[1] if len(sys.argv) > 1 else '.'
def J(pattern):
    out = []
    for fn in sorted(glob.glob(os.path.join(root, pattern))):
        for l in open(fn):
            try:
                out.append(json.loads(l))
            except Exception:
                pass
    return out

print('== classical screen (solscreen) ==')
C = J('cl/*.jsonl')
by = {}
for r in C:
    k = (r.get('kind'), round(r.get('eps', 0), 1) if r.get('kind') in ('N', 'T') else '')
    d = by.setdefault(k, dict(n=0, stable=0, sol=0, below=0, two=0, cand=0))
    d['n'] += 1
    if r.get('stable'):
        d['stable'] += 1
        if r.get('nsol', 0) > 0:
            d['sol'] += 1
            if r['sols'][0]['E'] > 0 and r['sols'][0]['E'] < r['om_min']:
                d['below'] += 1
        d['two'] += bool(r.get('two_valley')); d['cand'] += bool(r.get('cand'))
for k, d in sorted(by.items(), key=str):
    print(k, d)
cands = [r for r in C if r.get('cand')]
print('classical candidates (E_A > 0, below the magnon gap, second valley within 30 percent):', len(cands))
for r in sorted(cands, key=lambda r: r['dE_AB'] / max(r['A']['E'], 1e-9))[:10]:
    print('  ', r['kind'], 'A', r['A'], 'B', r['B'], 'om_min', round(r['om_min'], 4), 'v', np.round(r['v'], 3).tolist())

print('== quantum optimizer (qopt) ==')
Q = J('qopt/*.jsonl')
print('restarts', len(Q), 'by kind', {k: sum(r['kind'] == k for r in Q) for k in set(r['kind'] for r in Q)})
for thr in (0.05, 0.1, 0.2, 0.3):
    print(f'  O < {thr}:', sum(r['O'] < thr for r in Q))
print('== continuation to onset (qtune) ==')
T = J('qtune/*.jsonl')
import qtune
ok = [r for r in T if qtune.quality_ok(r['hist'][-1]['info'])]
gaps = sorted(0.5 * (r['hist'][-1]['info']['EA'] + r['hist'][-1]['info']['EB']) / r['hist'][-1]['info']['bw'] for r in ok)
print('tuned', len(T), 'quality endpoints', len(ok), 'smallest gap/bw', [round(g, 4) for g in gaps[:12]])
print('== past-onset DMRG (onset_dmrg) ==')
for fn in sorted(glob.glob(os.path.join(root, 'onset/*.json'))):
    try:
        d = json.load(open(fn))
    except Exception:
        continue
    inf = d['info']
    print(os.path.basename(fn), f"kA={inf['kA']:.3f} kB={inf['kB']:.3f} gap/bw={(inf['EA'] + inf['EB']) / 2 / inf['bw']:.3f} delta={d['delta']:.4f}")
    for p in d['points']:
        sg = d.get('sign', 1)
        w = p['w'][str(sg)] if str(sg) in p['w'] else p['w'][sg]
        print(f"   mu={p['mu']:.4f} lam={p['lam']:+.4f} L={p['L']} c={p['c_prof']:.3f}/{p['c_sm']:.3f} tr={p['trunc']:.0e} "
              f"wA={w[0]:.2f} wB={w[1]:.2f} pk={[(round(q, 2), round(x, 1)) for q, x in p['peaks'][:2]]}")
    print('   bracket', d.get('bracket'), 'lam*', d.get('lam_star'))
    best = max((p['c_prof'] for p in d['points'] if p['L'] >= 64), default=float('nan'))
    print('   max c_prof (L >= 64)', round(best, 3), '-> c = 2 candidate (needs L = 256 and FES)' if best > 1.5 else '')
print('== controls ==')
for fn in sorted(glob.glob(os.path.join(root, 'ctrl/*.out'))):
    print(os.path.basename(fn), open(fn).read().strip()[-400:])
