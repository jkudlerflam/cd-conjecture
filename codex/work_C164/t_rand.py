import numpy as np, time, json, solscreen as S, sys
rng = np.random.default_rng(int(sys.argv[1])); eps=float(sys.argv[2])
t0=time.time(); n=0; st=0; tv=0
while time.time()-t0 < 50:
    v = S.sample_near_hF(rng, eps, 'N')
    o = S.analyze(v); n+=1; st+=o['stable']
    if o['stable']:
        print(round(o['om_min'],3), o['nsol'], [(round(s['E'],4), round(s['P']/np.pi,3), round(s['edge'],3)) for s in o['sols'][:4]], o.get('two_valley'), o.get('dE_AB'), o.get('cand'))
        tv += bool(o.get('two_valley'))
print('bonds', n, 'stable', st, 'two-valley', tv, f'{(time.time()-t0)/n:.2f}s/bond')
