import numpy as np,json
from scipy.sparse import kron,eye,csr_matrix
from scipy.sparse.linalg import eigsh
from simple_dmrg import DMRG,qutrit_model
L=7;h,o=qutrit_model()
d=len(o);H=csr_matrix((d**L,d**L))
for i in range(L):H+=kron(kron(eye(d**i),csr_matrix(o)),eye(d**(L-i-1)))
for i in range(L-1):H+=kron(kron(eye(d**i),csr_matrix(h)),eye(d**(L-i-2)))
e=eigsh(H,k=1,which='SA',tol=1e-12)[0][0]
obj=DMRG(L,h,o);res=obj.run([8,16,27,27],output='dmrg_validation',tol=1e-11)
print('ED',e,'DMRG',res[-1]['E'],'error',res[-1]['E']-e)
assert abs(res[-1]['E']-e)<1e-8
