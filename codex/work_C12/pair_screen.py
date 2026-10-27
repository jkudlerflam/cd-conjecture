"""Closed-book two-body screen for a number-conserving three-state NN chain."""
import json
import numpy as np
from scipy.linalg import eigh_tridiagonal

# Local states 0,1,2. Bonds have 10 <-> 01: -t,
# 20 <-> 02: -tp, and 11 <-> (20,02): -g.
# Onsite doublon U and nearest-neighbor 11 energy V.
# At total momentum K, r=0 is the onsite pair and r>=1 are separated atoms.
def energies(K, U, V, g, tp, n=64):
    c=np.cos(K/2)
    diag=np.zeros(n)
    diag[0]=U-2*tp*np.cos(K)
    diag[1]=V
    off=np.full(n-1,-2*c)
    off[0]=-2*g*c
    return eigh_tridiagonal(diag,off,select='i',select_range=(0,0),check_finite=False)[0][0]

if __name__=='__main__':
    rng=np.random.default_rng(20261027)
    ks=np.linspace(0,np.pi,65)
    best=[]
    for trial in range(3000):
        U=rng.uniform(-12,8)
        V=rng.uniform(-12,5)
        g=10**rng.uniform(-1,0.8)
        tp=rng.uniform(-6,2)
        es=np.array([energies(k,U,V,g,tp) for k in ks])
        idx=np.argmin(es)
        if 1<idx<len(ks)-2 and es[idx]<-4.001:
            depth=min(es[0],es[-1])-es[idx]
            best.append(dict(U=U,V=V,g=g,tp=tp,K=float(ks[idx]),E=float(es[idx]),depth=float(depth)))
    best.sort(key=lambda p:-p['depth'])
    with open('pair_screen_results.json','w') as f: json.dump(best[:20],f,indent=2)
    print('Samples:',3000,'interior bound-band minima below atomic threshold:',len(best))
    print(json.dumps(best[:5],indent=2))
