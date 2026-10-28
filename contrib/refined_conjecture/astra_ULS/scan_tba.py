from tba_check import TBA
import numpy as np
import sys
N=int(sys.argv[1]); M=int(sys.argv[2]) if len(sys.argv)>2 else 24
obj=TBA(N,M,1024 if N<=8 else 2048,R=max(20,5*N))
x=None
for t in [.5,1,1.5,2,2.5,3,4,5,6,7,8]:
 out,x=obj.solve(t/N**2,x,tol=3e-9)
 print(out,flush=True)
