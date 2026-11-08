import json, numpy as np
rows = json.load(open('doubling_A.json'))+json.load(open('doubling_B.json'))
rows = [r for r in rows if 'error' not in r]
print('cases', len(rows))
out = {}
for eps in sorted(set(r['eps'] for r in rows), reverse=True):
    R = [r for r in rows if abs(r['eps']-eps) < 1e-12*max(eps,1e-12)+1e-15]
    def mx(k, f=max):
        v = [r[k] for r in R if r.get(k) is not None and np.isfinite(r[k])]
        return f(v) if v else None
    half = [r for r in R if 2*r['n'] == r['L']]
    d = dict(eta=R[0]['eta'], n=len(R), max_LD=mx('LD'), min_LD=mx('LD', min),
             max_LD_over_cft=max(r['LD']/r['cft_piv_over4'] for r in half),
             cft_ratio_largestL=[r['LD']/r['cft_piv_over4'] for r in half if r['L'] == max(h['L'] for h in half)][0],
             max_abs_LA_ref=max(abs(r['LA_ref']) for r in R), max_L_bound_ref=mx('L_bound_ref'),
             max_ratio_Blin=max(abs(r['LB_lin_ref']-r['LB_ref'])/max(abs(r['LB_ref']),1e-300) for r in R),
             max_w_over_2a=max(max(r['max_w_over_2a'], r['max_w_over_2a_at_ystar']) for r in R),
             max_pairwidth_over_delta=mx('pair_halfwidth_over_delta'),
             max_L_bound_mix_Q=mx('L_bound_mix_Q'), max_L_thm_actual_Q=mx('L_thm_bound_actual_Q'),
             max_etaNsp_Q=mx('etaNsp_mix_Q'), max_I1_Q=mx('I1_mix_Q'),
             min_lamratio_shell=mx('min_lamL_over_lamQ_shell', min), shell_pts=sum(r.get('n_shell_pts', 0) or 0 for r in R),
             min_lamratio_deep=mx('min_lamL_over_lamQ_deep', min),
             max_etaNsp_cell_2=mx('etaNsp_cell_2'), max_logas_rel=max(abs(r['LD']/r['LD_logas']-1) for r in R if r['etaL'] < 0.2) if any(r['etaL'] < 0.2 for r in R) else None,
             max_LD_check_err=max(abs(r['LD_check_Q']-r['LD']) for r in R if 'LD_check_Q' in r),
             cont_errors=sum('cont_error' in r for r in R))
    out[str(eps)] = d
    print(eps, json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in d.items()}))
json.dump(out, open('summary_C150.json', 'w'), indent=1)
