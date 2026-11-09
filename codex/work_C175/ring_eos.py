"""Ground energies E0(N, L) of model (3.2) on rings (all momenta), phi = pi/2, for the equation of state,
pair-binding energies and the grand-canonical density at mu."""
import numpy as np, sys, json, time
from ring_k import levels
wa=float(sys.argv[1]); Ls=[int(x) for x in sys.argv[2].split(',')]; Nmax=int(sys.argv[3])
w=1j*wa
for L in Ls:
    for N in range(1,Nmax+1):
        if 2*N>L: break
        t=time.time()
        best=(np.inf,None); spec=[]
        for m in range(L):
            e=levels(L,N,w,m,k=2)
            if len(e)==0: continue
            spec.append((m,float(e[0]),float(e[1]) if len(e)>1 else None))
            if e[0]<best[0]: best=(float(e[0]),m)
        print(json.dumps(dict(w=wa,L=L,N=N,E0=best[0],m0=best[1],spec=spec,sec=round(time.time()-t,1))),flush=True)
