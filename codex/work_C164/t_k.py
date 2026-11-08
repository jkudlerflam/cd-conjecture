import numpy as np, time, solscreen as S
kappa, mu = 0.05, 2.3
for e, h in ((0.2, 0.0266), (0.2, 0.03), (0.4, 0.035), (0.6, 0.05)):
    v = np.array([-1, -1 + mu*e*e/2, -1 + kappa*e*e/2, 0, 0, -e, -h*e*e, 0, 0])
    t = time.time()
    o = S.analyze(v, N=60)
    print(e, h, o['stable'], 'nsol', o.get('nsol'), [ (round(s['E'],5), round(s['P']/np.pi,3)) for s in o.get('sols',[])], o.get('two_valley'), f'{time.time()-t:.1f}s')
