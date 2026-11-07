import numpy as np, time
from circle import *
from check_c129 import energy
t=time.time()
for d in [1.001,1.1,2.,20.]:
    M=Circ(d); R=Restricted(M,np.pi/2)
    print(d, 'm',R.m,'E',R.E,'exact',e_half(M), 'r err', np.max(np.abs(R.r-M.r0(R.X))))
    err=0
    for L in range(2,11):
        for n in range(1,L//2+1):
            E,_=solve(L,n,d); err=max(err,abs(float(E)-energy(L,n,d)))
    print('  ED err',err)
print(time.time()-t)
