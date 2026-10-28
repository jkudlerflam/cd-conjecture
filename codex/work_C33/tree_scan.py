import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,json,time
from structured import tree,visibility
out=[]
for n in [1,2,3,4,5,6,8,10,12,16,23,32,48,64]:
 for D in [10.,100.,1000.]:
    start=time.monotonic();F,lb=tree(n,D,normalized=True)
    scale=F.max();S=(F+F.T)/2/scale
    ev=np.linalg.eigvalsh(S)
    diag=np.sqrt(np.diag(S));E=S/diag[:,None]/diag[None,:]
    ee=np.linalg.eigvalsh(E)
    r=dict(n=n,D=D,log10_T_scale_g2=float(np.log10(2*visibility(n)*scale)+2*lb/np.log(10)),log10_beta=float(lb/np.log(10)),
           herm_min_over_scale=float(ev[0]),equilibrated_min=float(ee[0]),
           asymmetry=float(np.linalg.norm((F-F.T)/scale)/np.linalg.norm(F/scale)),seconds=time.monotonic()-start)
    out.append(r);print(json.dumps(r),flush=True)
    np.savez(f'tree_n{n}_D{int(D)}.npz',F=F,logbeta=lb)
    json.dump(out,open('tree_results.json','w'),indent=2)
