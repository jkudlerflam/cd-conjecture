import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
from ycell_reference import params,build_cell,code,second_order,one_particle_T
from finite_g import exact_T
from structured import Leading
out=[]
for n in [1,2,3]:
    D=2.;a=3.;g=.3
    E0=-2*np.cos(np.pi/(n+2))
    P=params(Dx=D,DM=D,Delta=a+E0,delta=0.,gp=g,gm=g,phx={'0':np.pi/(2*n),'+':0.,'-':0.})
    c=build_cell(n,P);E,G,gap=code(c);H2,_=second_order(c,G,E)
    T=one_particle_T(H2,n)
    S,meta=exact_T(n,D,a,g)
    err=np.max(abs(S-T))/np.max(abs(S))
    out.append(dict(n=n,relative_error=err,cell_gap=gap))
    print(out[-1],flush=True)
    assert err<1e-6
import json
json.dump(out,open('validation.json','w'),indent=2)
