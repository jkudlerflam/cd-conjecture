import numpy as np, solscreen as S, clat
kappa, mu = 0.05, 2.3
for e in (0.3, 0.5):
    for h in (0.0, 0.03, 0.08, 0.15, 0.3):
        v = np.array([-1, -1 + mu*e*e/2, -1 + kappa*e*e/2, 0, 0, -e, -h*e*e, 0, 0])
        J, f = clat.J_f(v)
        st = clat.stable_uniform(J, f, Lmax=16, nstart=2)
        sols = S.find_solitons(J, f, st['n'], N=40)
        sz = [int(np.sum(np.linalg.norm(s['n']-st['n'],axis=1)>0.3)) for s in sols[:3]]
        print(e, h, 'stable', st['stable'], 'tex_gain', round(st['texture_gain'],5), 'om_min', round(st['om_min'],4), 'k', round(st['k_ommin']/np.pi,3),
              [(round(s['E'],4), round(s['P']/np.pi,3)) for s in sols[:3]], 'size', sz, flush=True)
