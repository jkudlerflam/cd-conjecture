import json,numpy as np,sys
def last(f): return json.loads(open(f).read().strip().split('\n')[-1])
def slope(C,r,R,qs):
    s=r<=R; r=r[s]; C=C[s]; w=np.cos(np.pi*r/(2*R))**2
    y=np.array([2*np.sum(w*C*np.cos(q*r)) for q in qs])
    A=np.stack([np.ones_like(qs),qs,qs**2],1); c,a,b=np.linalg.lstsq(A,y,rcond=None)[0]
    return a
for tag in sys.argv[1:]:
    m=last(f'typhon_out/out/{tag}.jsonl')['meas']
    out=[f'{tag:14s}']
    for name,off in (('Cnn',1),('Cpp',2)):
        res=[]
        vs=[np.array(v) for v in m[name].values()]; n0=min(len(v) for v in vs); Cs=np.mean([v[:n0] for v in vs],0)
        r=(np.arange(len(Cs))+off).astype(float)
        for R in (60,90,120):
            for (lo,hi) in ((3,8),(4,10)):
                qs=np.linspace(lo*np.pi/R,hi*np.pi/R,16)
                res.append(slope(Cs,r,R,qs)*np.pi)
        out.append(name+': '+' '.join(f'{x:.3f}' for x in res))
    print('  '.join(out))
