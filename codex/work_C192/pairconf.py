# Pair exponent with the open-chain conformal form (run from codex/work_C176).
# D(r) from reference site x0 to x0+r, demodulated at momentum Q, fit |D| ~ F(r)^(-eta),
# F = d(x1-x2) d(x1+x2) / sqrt(d(2x1) d(2x2)), d(x) = sin(pi x / 2L).
import json,numpy as np,sys
def last(f): return json.loads(open(f).read().strip().split('\n')[-1])
Q=float(sys.argv[1]); d=lambda x,L: np.sin(np.pi*x/(2*L))
for tag in sys.argv[2:]:
    rec=last(f'typhon_out/out/{tag}.jsonl'); m=rec['meas']['Dpair']; L=rec['cfg']['L']
    out=[f'{tag:13s}']
    for sep in ('2','3'):
        xs=[];ys=[]
        for key,v in m.items():
            if not key.startswith(sep+'_'): continue
            x0=int(key.split('_')[1])+0.5+int(sep)/2; a=np.array(v); D=a[:,0]+1j*a[:,1]
            r=np.arange(len(D))+int(sep)+1     # first stored separation
            h=(D*np.exp(-1j*Q*r)).real; x2=x0+r
            F=d(r,L)*d(2*x0+r,L)/np.sqrt(d(2*x0,L)*d(2*x2,L))
            xs.append((r,F,h))
        res=[]
        for (lo,hi) in ((6,L//8),(10,L//5),(16,L//3),(24,int(.4*L))):
            X=[];Y=[]
            for r,F,h in xs:
                s=(r>=lo)&(r<=hi)&(h>0); X+=list(np.log(F[s])); Y+=list(np.log(h[s]))
            p=np.polyfit(X,Y,1); res.append((-p[0],np.std(np.array(Y)-np.polyval(p,X))))
        out.append(f'd={sep}: '+' '.join(f'{e:.3f}({s:.2f})' for e,s in res))
    print('  '.join(out))
