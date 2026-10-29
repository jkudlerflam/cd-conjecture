"""Full marker sectors for the small frozen F3 diagnostic, without charge cuts."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import itertools,json
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
from bond13 import MAT,marker,params,material_bond

def matrix(marker_sites,p,lam):
    choices=[[i for i,a in enumerate(MAT) if marker(a)==(j in marker_sites)] for j in range(4)]
    states=list(itertools.product(*choices));ix={w:i for i,w in enumerate(states)}
    h=material_bond(p,lam);rows=[];cols=[];vals=[]
    tr={i:[(int(j),h[j,i]) for j in np.flatnonzero(h[:,i])] for i in range(121)}
    for col,w in enumerate(states):
        for i,a in enumerate(w):
            j=(i+1)%4
            for q,v in tr[a*11+w[j]]:
                z=list(w);z[i],z[j]=divmod(q,11)
                rows.append(ix[tuple(z)]);cols.append(col);vals.append(v)
    return sp.coo_matrix((vals,(rows,cols)),shape=(len(states),len(states))).tocsc()
def main():
    aa=lambda b:2*(1-np.cos(np.pi/(b+1)))
    dd=lambda b:(b+1)*aa(b)-b*aa(b+1)
    mu=2-(dd(2)+dd(1))/2
    ee=lambda b:(mu-2*np.cos(np.pi/(b+1)))/b
    delta=min(ee(1)-ee(2),ee(3)-ee(2));K=100.
    gamma=delta-1/(K-4)
    out=[]
    for lam in (0.,.001,.003):
        energies=[];resids=[]
        for markers in ([1,3],[0,3]):
            H=matrix(markers,params(n=1,K=K,musel=mu,g=.03,td=.001),lam)
            e,v=eigsh(H,k=1,which='SA',tol=1e-11,ncv=45,v0=np.ones(H.shape[0]))
            energies.append(float(e[0]));resids.append(float(np.linalg.norm(H@v[:,0]-e[0]*v[:,0])))
        lower=2*(gamma-5*np.sqrt(2)*lam)
        assert lower>0 and energies[1]-energies[0]>=lower
        rec=dict(lambda_value=lam,K=K,dimension_each=H.shape[0],selected_lengths=[2,2],wrong_lengths=[1,3],
                 energy_difference=energies[1]-energies[0],F3_lower_bound=lower,max_residual=max(resids))
        out.append(rec);print(json.dumps(rec),flush=True)
    with open('frozen_checks.json','w') as f:json.dump(out,f,indent=2)
if __name__=='__main__':main()
