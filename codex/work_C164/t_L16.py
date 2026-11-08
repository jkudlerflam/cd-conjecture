import json, numpy as np, sys, qscreen as Q, time
fs = sys.argv[1:]
for fn in fs:
    for l in open(fn):
        r = json.loads(l)
        if not r.get('flag'): continue
        v = np.array(r['v']); t=time.time()
        o = Q.assess(v, 16)
        print(fn.split('/')[-1], 'L12 gap %.4f mins %s | L16 gap %.4f gap0 %.4f m0 %d mins %s flag %s  %.0fs' % (r['gap'], r['mins'][:2], o['gap'], o['gap0'], o['m0'], o.get('mins',[])[:2], o['flag'], time.time()-t), flush=True)
