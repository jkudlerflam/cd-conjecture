"""C163 summary of big/, fes/, twisted/, paths/ (Claude Opus agent, 2026-11-08). usage: python3 summarize.py"""
import json, glob, os, numpy as np

print('== big: open-chain DMRG at one L per file (c_prof, c_sm, Smid, truncation)')
rows = {}
for fn in sorted(glob.glob('big/*.json')):
    d = json.load(open(fn)); L = int(d['L_']); r = d['L'][str(L)] if str(L) in d['L'] else d['L'][L]
    rows.setdefault(d['name'], {})[L] = r
    print(f"  {d['name']:14s} L={L:4d} c_prof={r['c_prof']:.3f} c_sm={r['c_sm']:.3f} Smid={r['Smid']:.4f} trunc={r['trunc']:.1e} t={d['time_s']:.0f}s")
ref = {}
for fn in glob.glob('c159_refine/*.json'):
    d = json.load(open(fn)); ref[os.path.basename(fn)[:-5]] = d
print('  two-size slopes 6 dS/ln(L2/L1) using C159 L = 32, 64, 128 where available:')
for name, R in rows.items():
    S = {}
    if name in ref and isinstance(ref[name].get('dmrg'), dict):
        for k, r in ref[name]['dmrg']['L'].items():
            S[int(k)] = r['Smid']
    for L, r in R.items():
        S[L] = r['Smid']
    Ls = sorted(S)
    sl = [f"{Ls[i]}-{Ls[i+1]}: {6 * (S[Ls[i+1]] - S[Ls[i]]) / np.log(Ls[i+1] / Ls[i]):.3f}" for i in range(len(Ls) - 1)]
    print(f"    {name:14s} " + '  '.join(sl))

print('\n== fes: two-site-cell VUMPS')
try:
    import fesfit_c158 as F
except Exception:
    F = None
for fn in sorted(glob.glob('fes/*.json')):
    d = json.load(open(fn)); R = d['results']
    print(os.path.basename(fn))
    for r in R:
        print(f"  chi={r['chi']} S={np.mean(r['S']):.4f} xi={r['xi']:.2f} err={r.get('err', float('nan')):.1e}")
    if F is not None and len(R) >= 3:
        try:
            F.analyse(fn)
        except Exception as ex:
            print('  fesfit failed', ex)

print('\n== twisted curve J = -R_z(phi) + eps X')
for fn in sorted(glob.glob('twisted/*.json')):
    d = json.load(open(fn))
    if 'dmrg' not in d:
        print('  ', d['spec']['tag'], 'incomplete'); continue
    L = d['dmrg']['L']; g = lambda k: L.get(k, L.get(int(k)))
    print(f"  {d['spec']['tag']:18s} parent Delta={d['parent']['Delta']:+.3f} K_S12={d['K_S12'] if isinstance(d['K_S12'], str) else round(d['K_S12'], 2)} "
          f"K_S16={d['K_S16'] if isinstance(d['K_S16'], str) else round(d['K_S16'], 2)} | L64 c_prof={g('64')['c_prof']:.3f} "
          f"L128 c_prof={g('128')['c_prof']:.3f} c_two={d['dmrg']['c_two']:.3f} S128={g('128')['Smid']:.3f} tr={g('128')['trunc']:.0e}")

print('\n== paths (attractive window)')
allc = []
for fn in sorted(glob.glob('paths/*.json')):
    d = json.load(open(fn)); s = d['spec']
    P = sorted(d['pts'].values(), key=lambda r: r['lam'])
    line = ' '.join(f"{r['lam']:.2f}:{r['Smid']:.2f}/{r['K_S12']:.1f}" for r in P[::2])
    print(f"  {s['tag']} Delta={s['Delta']} theta={s['theta']:.2f} | lam:Smid/K_S12 {line}")
    if 'final' in d:
        f = d['final']; L = f['L']; g = lambda k: L.get(k, L.get(int(k)))
        print(f"     best lam={d['best_lam']:.3f} excess={d['excess_top'][0][0] if d['excess_top'] else 0:.4f} "
              f"final: c_prof64={g('64')['c_prof']:.3f} c_prof128={g('128')['c_prof']:.3f} c_two={f['c_two']:.3f} K_S16={f['K_S16']}")
        allc.append((f['c_two'], g('128')['c_prof'], s['tag'], d['best_lam'], f['K_S16']))
if allc:
    allc.sort(reverse=True)
    print('  largest c_two at the selected points:', allc[:6])
    print('  candidates by the pre-fixed criterion (c_two and c_prof(128) > 1.2):', [x for x in allc if x[0] > 1.2 and x[1] > 1.2])
