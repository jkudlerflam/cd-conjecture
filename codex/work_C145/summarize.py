"""Aggregate the C145 numerical checks into summary_C145.json."""
import json, glob
rows = []
for f in sorted(glob.glob('check_shell_*.json')):
    rows += [r for r in json.load(open(f)) if 'min_sig_over_r_K0' in r]
by = {}
for r in rows:
    e = by.setdefault(r['eps'], dict(cases=0, min_sig_over_r=9, max_etaVar_K0=0, max_etaVar_S=0, Lmax=0))
    e['cases'] += 1; e['min_sig_over_r'] = min(e['min_sig_over_r'], r['min_sig_over_r_K0'])
    e['max_etaVar_K0'] = max(e['max_etaVar_K0'], r['etaVar_K0']); e['max_etaVar_S'] = max(e['max_etaVar_S'], r['etaVar_S'])
    e['Lmax'] = max(e['Lmax'], r['L'])
out = dict(shell={str(k): v for k, v in sorted(by.items())})
try:
    w = json.load(open('check_W13.json'))
    out['W13'] = dict(cases=len(w), min_rel_dv=min(r['min_dv_rel'] for r in w), v_monotone_all=all(r['min_dv_rel'] > -1e-8 for r in w),
                      max_VarV_over_eta=max(r['VarV_over_eta'] for r in w), max_vQ_over_2a_tanh=max(r['vQ_over_2a_tanh'] for r in w),
                      min_cZ=min(r['cZ'] for r in w), argmin_cZ=min(w, key=lambda r: r['cZ']),
                      max_ebulk2_over_eta2=max(r['ebulk2_over_eta2'] for r in w),
                      min_cZ_by_eps={str(e): min(r['cZ'] for r in w if r['eps'] == e) for e in sorted(set(r['eps'] for r in w))})
except FileNotFoundError:
    pass
json.dump(out, open('summary_C145.json', 'w'), indent=1); print(json.dumps(out, indent=1))
