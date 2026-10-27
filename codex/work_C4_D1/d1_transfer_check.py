"""Direct small-size alternating six-vertex transfer-matrix check."""
import numpy as np
from d1 import tl_e, graded_swap, local_sum, periodic_tl

def transfer(u,g,L):
    e=tl_e(g);p=graded_swap(2,0)
    mon=np.eye(2**(L+1),dtype=complex)
    for j in range(L):
        w=u-(np.pi/2 if j%2 else 0)
        r=p@(np.sin(g-w)*np.eye(4)+np.sin(w)*e)
        mon=local_sum(r,(0,j+1),2,L+1)@mon
    return np.trace(mon.reshape(2,2**L,2,2**L),axis1=0,axis2=2)

g=.4;L=4;eps=1e-6
H=np.zeros((2**L,2**L),dtype=complex)
for u in [0,np.pi/2]:
    t=transfer(u,g,L)
    dt=(transfer(u+eps,g,L)-transfer(u-eps,g,L))/(2*eps)
    H+=np.linalg.solve(t,dt)
A=periodic_tl(g,L).toarray();eye=np.eye(2**L)
c=np.linalg.lstsq(np.column_stack([A.ravel(),eye.ravel()]),H.ravel(),rcond=None)[0]
print('fit Htransfer = alpha Htl + beta I:',c)
print('relative residual',np.linalg.norm(H-c[0]*A-c[1]*eye)/np.linalg.norm(H))
print('relative antihermitian part',np.linalg.norm(H-H.conj().T)/np.linalg.norm(H))
print('expected alpha',1/(np.sin(g)*np.cos(g)))
print('expected beta',-2*L/np.tan(2*g))
