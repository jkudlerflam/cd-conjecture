"""C168 job digest: proxy efficiencies of the two-stage screen, stage-A refinement, sampling bound.
Claude Opus 5.5, 2026-11-08. Usage: python3 summarize.py [outdir]"""
import json, sys, os, numpy as np
d = sys.argv[1] if len(sys.argv) > 1 else 'out'


def load(fn):
    return [json.loads(l) for l in open(os.path.join(d, fn))] if os.path.exists(os.path.join(d, fn)) else []


P = load('proxy.jsonl'); S1 = load('s1.jsonl'); S2 = load('s2.jsonl')
print('== Stage P: proxy calibration (range-3 qubit chains, c = number of Fermi seas) ==')
eff = {}
for fam in ('free', 'int'):
    for ns in (1, 2):
        R = [r for r in P if r['fam'] == fam and r['nseas'] == ns and 'e1216' in r]
        if not R:
            continue
        f1 = np.array([r['e1216']['c_S'] > 1.5 for r in R])
        f2 = np.array([('e1620' in r) and r['e1620']['c_S'] > 1.5 for r in R])
        old = np.array([r['e1216']['crit'] < 0.15 and r['e1216']['c_obj'] > 1.1 for r in R])
        eff[(fam, ns)] = (f2.mean(), old.mean())
        print(f'{fam:4s} c={ns}  N {len(R):4d}  F1 {f1.mean():.3f}  F1&F2 {f2.mean():.3f}  C155 rule (crit<.15 & c_obj>1.1) {old.mean():.3f}')
print('\n== Stage S1: stage-A points with c_S(12,16) > 1.5 or C155 flag, rings 16/20 ==')
ok = [r for r in S1 if 'e1620' in r]
cs = np.array([r['e1620']['c_S'] for r in ok])
print(f'points {len(S1)}  done {len(ok)}  c_S(16,20) > 1.5: {(cs > 1.5).sum()}  > 1.3: {(cs > 1.3).sum()}  quantiles {np.round(np.quantile(cs, [.1, .5, .9]), 2) if len(cs) else None}')
print('\n== Stage S2: DMRG L = 32, 64, 128 ==')
ncand = 0
for r in sorted(S2, key=lambda r: -r.get('dmrg', {}).get('L', {}).get('128', {}).get('c_prof', -9) if 'dmrg' in r else 9):
    if 'dmrg' not in r:
        print(r['idx'], 'error', r.get('error')); continue
    Ls = r['dmrg']['L']
    g = lambda L, k: Ls.get(str(L), Ls.get(L, {})).get(k, np.nan)
    m64, m128 = g(64, 'Smid'), g(128, 'Smid')
    c2 = 6 * (m128 - m64) / np.log(2)
    cand = (g(128, 'c_prof') > 1.2) and (c2 > 1.2) and (g(128, 'trunc') < 1e-7)
    ncand += cand
    print(f"idx {r['idx']:6d} c_S(16,20) {r['e1620']['c_S']:.2f} | c_prof 32/64/128 {g(32,'c_prof'):.2f} {g(64,'c_prof'):.2f} {g(128,'c_prof'):.2f} | c_two64-128 {c2:.2f} | trunc128 {g(128,'trunc'):.1e} {'CANDIDATE' if cand else ''}")
print(f'\ncandidates meeting the pre-fixed rule: {ncand}')
N = 12000
for fam in ('free', 'int'):
    if (fam, 2) in eff and eff[(fam, 2)][0] > 0:
        e = eff[(fam, 2)][0]
        print(f'bound with eta = F1&F2 efficiency of {fam} c=2 proxies ({e:.3f}), zero candidates: p < {3.0/(e*N):.2e} (95%), {4.61/(e*N):.2e} (99%)')
