"""Check the full second-order projected operator and a third-order source."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
from checks import ham
from singular_corner import heff
L=6
# V=1, h=1/2 corresponds to mu=0. Add L/4 to remove spin constant.
H0,b=ham(L,0,0,.25,0,.5)
d=H0.diagonal()+L/4
A,b=ham(L,.5,-.5,0,0,0)
p=np.where(abs(d)<1e-12)[0];q=np.where(d>1e-12)[0]
M=A.toarray()[np.ix_(p,q)]
SW=-(M/d[q])@M.T
err=0.
for N in range(4):
    E,bb=heff(L,N)
    ids=np.searchsorted(p,bb)
    err=max(err,float(np.max(abs(SW[np.ix_(ids,ids)]-E.toarray()))))
L=8
H0,b=ham(L,0,0,.25,0,.5);d=H0.diagonal()+L/4
A,b=ham(L,.5,-.5,0,0,0);A=A.toarray()
p=np.where(abs(d)<1e-12)[0];q=np.where(d>1e-12)[0]
# P A R A R A P with R=-Q/H0; two resolvents give positive sign.
source=(A[9,q]/d[q])@A[np.ix_(q,q)]@(A[q,0]/d[q])
text=f'L=6, full second-order operator maximum error: {err}\n'
text+=f'L=8, third-order source <1001|H3|vacuum>: {source}; expected 2/3\n'
open('operator_checks.txt','w').write(text)
print(text)
assert err<1e-12 and abs(source-2/3)<1e-12
