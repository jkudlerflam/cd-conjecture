from canonical_dmrg import CanonicalDMRG
from simple_dmrg import qutrit_model
from qutrit_ed import matrix
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh
import numpy as np
L,N=8,3;delta=.9;h,o=qutrit_model(delta,mu=0)
H,ss=matrix(L,N,pbc=False,U=-14/9+delta)
# The generic bond-sum ED routine gives half an onsite term at each boundary.
n0=(ss&3).astype(int);nlast=((ss>>np.uint64(2*(L-1)))&3).astype(int)
H=H+diags((-14/9+delta)/2*((n0==2).astype(int)+(nlast==2).astype(int)))
e=eigsh(H,k=1,which='SA',tol=1e-12)[0][0]
obj=CanonicalDMRG(L,N,h,o);res=obj.run([8,16,27,27],output='canonical_validation',tol=1e-11)
print('ED / canonical DMRG',e,res[-1]['E'],'error',res[-1]['E']-e)
assert abs(res[-1]['E']-e)<1e-8
