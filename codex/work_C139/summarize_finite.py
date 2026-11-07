import json, numpy as np
from constants import K as Kc
summ = {}
for d in ['0.1', '0.5', '0.9', '0.99', '0.999']:
    R = json.load(open(f'finite_{d}.json'))
    g = np.arccos(float(d)); K = Kc(g)
    s = dict(cases=len(R), Lmax=max(r['L'] for r in R),
             max_LsupD=max(r['L_supD'] for r in R), B_over_2=R[0]['B_over_2'],
             max_LsupD_over_Bhalf=max(r['L_supD']/r['B_over_2'] for r in R),
             max_sigdev_over_D0=max(r['L_sup_sigma_dev']/r['D0'] for r in R),
             min_sigma_global=min(r['min_sigma'] for r in R),
             min_sigma_roots=min(r['min_sigma_roots'] for r in R),
             max_abs_LdF=max(abs(r['L_F_minus_Le']) for r in R), K=K,
             max_abs_LdF_over_K=max(abs(r['L_F_minus_Le']) for r in R)/K,
             max_defect_M0=max(r.get('defect_over_M0', -1) for r in R),
             half_filled_LdF_at_Lmax=[r['L_F_minus_Le'] for r in R if r['nu'] == 0.5][-1],
             casimir=-np.pi*(np.pi*np.sin(g)/(2*g))/6)
    summ[d] = s; print(d, json.dumps(s))
json.dump(summ, open('finite_summary.json', 'w'), indent=1)
