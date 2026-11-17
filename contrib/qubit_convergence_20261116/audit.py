#!/usr/bin/env python3
"""Read-only audit of C190/C195 raw records; only Python's standard library.

The canonical final entropy and the last-sweep entropy are separate estimators.
An absent predecessor/checkpoint is never interpreted as convergence.
"""
import hashlib
import json
import math
from pathlib import Path
import statistics as st

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WS = (1, 6, 12)
DIRS = ["codex/work_C190/gpu/typhon_out/" + d for d in ("outG", "outX1024", "outX1536")]
DIRS += ["codex/work_C195/" + d for d in ("typhon_out/out512", "typhon_out_b/out512", "typhon_out_x2048")]
KEYS = ("sym40_n316", "p153c", "xx1", "xx2")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sbar(S, L, W):
    return st.fmean(s for l, s in enumerate(S, 1) if abs(l - L // 2) <= W)


def ctwo(a, b):
    return {str(w): 6 * (b['Sbar'][str(w)] - a['Sbar'][str(w)]) / math.log(b['L'] / a['L']) for w in WS}


def stage_checks(hist, L):
    x = hist[-1]
    out = dict(sweeps=len(hist), res_max=x['res_max'], trunc=x['trunc'],
               dE=None, dS=None, dn=None, pass_T1_numbers=False)
    if len(hist) >= 2:
        y = hist[-2]
        out.update(dE=abs(x['E'] - y['E']),
                   dS=max(abs(x[f'Sbar{w}'] - y[f'Sbar{w}']) for w in WS),
                   dn=abs(x['n_mid'] - y['n_mid']))
        out['pass_T1_numbers'] = (out['res_max'] <= 1e-8 and out['trunc'] <= 1e-7
            and out['dE'] <= 1e-8 * max(1, L / 256) and out['dS'] <= 1e-3 and out['dn'] <= 5e-4)
    return out


def load(path):
    r = json.loads(path.read_text())
    cfg = r.get('cfg', {})
    key = next((k for k in KEYS if cfg.get('tag', '').startswith(k + '_L')), None)
    if key is None or r.get('status') != 'done':
        return None
    L = cfg['L']
    hs = [h for h in r['hist'] if 'Sbar6' in h and h['alpha'] == 0]
    stages = {}
    for chi in sorted({h['chi'] for h in hs}):
        hh = [h for h in hs if h['chi'] == chi]
        h = hh[-1]
        stages[str(chi)] = dict(chi=chi, E=h['E_phys'],
            Sbar={str(w): h[f'Sbar{w}'] for w in WS}, n_mid=h['n_mid'],
            checks=stage_checks(hh, L), sweep_seconds=h.get('t'),
            host_rss_GB=h.get('host_rss_GB'), dev_used_GB=h.get('dev_used_GB'))
    last = hs[-1]
    finalS = {str(w): sbar(r['S'], L, w) for w in WS}
    return dict(file=str(path.relative_to(ROOT)), sha256=sha(path), key=key, L=L, P=cfg['P'],
        model=cfg['model'], gauge=cfg.get('gauge'), hz=cfg.get('hz'), mu=cfg.get('mu'),
        chi=last['chi'], final_chi=cfg['final']['chi'], E=r.get('E_phys', last['E_phys']),
        Sbar=finalS, sweep_Sbar=stages[str(last['chi'])]['Sbar'],
        canonical_minus_sweep={str(w): finalS[str(w)] - last[f'Sbar{w}'] for w in WS},
        n_mid=st.fmean(r['meas']['dens'][L//4:3*L//4]), parity=r['parity'],
        stages=stages, checks=stage_checks([h for h in hs if h['chi']==last['chi']], L),
        checkpoint_local_exists=(path.parent / 'ckpt' / (path.stem + '.npz')).is_file(),
        continues=cfg.get('continues'))


def main():
    runs = [r for d in DIRS for p in sorted((ROOT / d).glob('*.json')) if (r := load(p))]
    # Highest-chi finished run, then lowest energy among starts at that chi.
    selected = {}
    for r in runs:
        k = (r['key'], r['L'], r['P'])
        if k not in selected or (r['chi'], -r['E']) > (selected[k]['chi'], -selected[k]['E']):
            selected[k] = r
    # Merge stage provenance from the selected run first; never hide source paths.
    for r in runs:
        prior = [(int(c), q, s) for q in runs if (q['key'],q['L'],q['P'])==(r['key'],r['L'],r['P'])
                 for c,s in q['stages'].items() if int(c)<r['chi']]
        if prior:
            c, q, s = max(prior, key=lambda t:(t[0], -t[2]['E']))
            step = {str(w): r['sweep_Sbar'][str(w)]-s['Sbar'][str(w)] for w in WS}
            r['C4'] = dict(from_chi=c, predecessor_file=q['file'], increments=step,
                pass_abs_S6=abs(step['6']) <= .01, predecessor_checks=s['checks'])
        else:
            r['C4'] = dict(from_chi=None, pass_abs_S6=None)
    branches, comparisons = [], []
    for key in KEYS:
        ground = {}
        for L in sorted({r['L'] for r in runs if r['key']==key}):
            sectors = [selected[key,L,p] for p in (0,1) if (key,L,p) in selected]
            g = min(sectors, key=lambda r:r['E'])
            ground[L] = g
            b = dict(key=key,L=L,P=g['P'], sectors=[dict(P=r['P'],chi=r['chi'],E=r['E'],file=r['file']) for r in sectors],
                     both_sectors=len(sectors)==2, matched_chi=len({r['chi'] for r in sectors})==1,
                     ground_order='numerical variational comparison; no certified error bound')
            if len(sectors)==2:
                b['E_even_minus_odd'] = sectors[0]['E']-sectors[1]['E']
                common = set.intersection(*[set(c for r in runs if (r['key'],r['L'],r['P'])==(key,L,p)
                                                    for c in r['stages']) for p in (0,1)])
                c = max(common, key=int)
                es = [min(s['E'] for r in runs if (r['key'],r['L'],r['P'])==(key,L,p)
                          for ch,s in r['stages'].items() if ch==c) for p in (0,1)]
                b['last_common_chi'] = int(c)
                b['common_chi_E_even_minus_odd'] = es[0]-es[1]
            branches.append(b)
        for L,a in ground.items():
            if 2*L not in ground:
                continue
            b=ground[2*L]
            comparisons.append(dict(key=key,Lsmall=L,Lbig=2*L,parities=[a['P'],b['P']],chis=[a['chi'],b['chi']],
                files=[a['file'],b['file']],c_two=ctwo(a,b),
                c_two_last_sweep=ctwo(dict(a,Sbar=a['sweep_Sbar']),dict(b,Sbar=b['sweep_Sbar'])),
                density_difference=b['n_mid']-a['n_mid'],T4=abs(b['n_mid']-a['n_mid'])<=.006,
                legacy_T2_numbers=b['checks']['pass_T1_numbers'] and min(ctwo(a,b)[str(w)] for w in (6,12))>=1.6,
                C4_large=b['C4'], all_pairings=[dict(parities=[p,q],c_two=ctwo(selected[key,L,p],selected[key,2*L,q]))
                   for p in (0,1) for q in (0,1) if (key,L,p) in selected and (key,2*L,q) in selected]))
    out=dict(source='existing saved data only; no new production', runs=runs, branches=branches,
             comparisons=comparisons, new_chi_test='UNDECIDED: no chi > 2048 production output')
    (HERE/'audit.json').write_text(json.dumps(out,indent=2)+'\n')
    for c in comparisons:
        print(c['key'],c['Lsmall'],c['Lbig'],'P',c['parities'],'chi',c['chis'],
              'c_two', {w:round(v,8) for w,v in c['c_two'].items()},'C4',c['C4_large']['pass_abs_S6'])
    for b in branches:
        if b['key']=='sym40_n316':
            print('branch',json.dumps(b))


if __name__=='__main__':
    main()
