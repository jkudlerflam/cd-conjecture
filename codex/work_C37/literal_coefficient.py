"""Independent direct Taylor coefficient from the literal cell, n<=3.
No subtraction of two finite-g Hamiltonians is used.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys,json,time
import numpy as np
sys.path.insert(0,'../../contrib/routeA_spatial')
from ycell import params,build_cell,code,logical
from structured import Leading

out=[]
for n in [1,2,3]:
    start=time.monotonic()
    pp=params(Dx=2,DM=2,Delta=3-2*np.cos(np.pi/(n+2)),delta=0,gp=0,gm=0,td=.03,
       phx={'0':np.pi/(2*n),'+':0.,'-':0.},phd={'0':np.pi,'+':0.,'-':0.})
    cell=build_cell(n,pp);E,G,gap=code(cell);ex=cell['exc'];H=cell['H'][np.ix_(ex,ex)]
    pp1=dict(pp,gp=1,gm=1);other=build_cell(n,pp1)
    C=(other['H']-cell['H'])[np.ix_(ex,ex)]
    w,U=np.linalg.eigh(H);w=w-E;C=U.conj().T@C@U;r=len(w)
    W=np.zeros((r*r,4**n),complex)
    for fm,fx,c in cell['terms']:
        A=U.conj().T@(cell['FM'][fm]@G)[ex];B=U.conj().T@(cell['FX'][fx]@G)[ex]
        W+=c*np.kron(A,B)
    den=(w[:,None]+w[None,:]).ravel()
    X=(W/den[:,None]).reshape(r,r,4**n)
    Y=np.einsum('ab,bck->ack',C,X,optimize=True)+np.einsum('cb,abk->ack',C,X,optimize=True)
    Y=Y.reshape(r*r,4**n)
    K=-Y.conj().T@(Y/den[:,None])
    L=logical(K,n);T0,M0,_=Leading(n,D=2,a=3,td=.03,loopflux=np.pi).calculate();sc=np.max(abs(T0))
    result=dict(n=n,excited_cell_dimension=r,pair_dimension=r*r,code_bond_dimension=4**n,
      T_relative_error=float(np.max(abs(L['T']-T0))/sc),
      M_error_over_T=float(np.max(abs(L['muA']+L['muB']-M0))/sc),
      nonGaussian_over_T=float(L['res']/sc),seconds=time.monotonic()-start)
    out.append(result);print(json.dumps(result),flush=True)
    json.dump(out,open('literal_coefficient_results.json','w'),indent=2)
