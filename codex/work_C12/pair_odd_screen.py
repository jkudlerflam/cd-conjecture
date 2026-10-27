"""Odd-parity pair conversion: 11 -> 20 and 11 -> 02 have opposite signs."""
import json
import numpy as np
from scipy.linalg import eigh_tridiagonal

def energy(K,U,V,g,tp,n=80):
    diag=np.zeros(n)
    diag[0]=U-2*tp*np.cos(K)
    diag[1]=V
    off=np.full(n-1,-2*np.cos(K/2))
    off[0]=-2*g*np.sin(K/2)
    return eigh_tridiagonal(diag,off,select='i',select_range=(0,0),check_finite=False)[0][0]

if __name__=='__main__':
    rng=np.random.default_rng(20261028)
    ks=np.linspace(0,np.pi,65)
    best=[]
    for trial in range(5000):
        U=rng.uniform(-10,6)
        V=rng.uniform(-10,2)
        g=10**rng.uniform(-1,0.7)
        tp=rng.uniform(-4,4)
        es=np.array([energy(k,U,V,g,tp) for k in ks])
        idx=np.argmin(es)
        if 1<idx<len(ks)-2 and es[idx]<-4.001:
            best.append(dict(U=U,V=V,g=g,tp=tp,K=float(ks[idx]),E=float(es[idx]),depth=float(min(es[0],es[-1])-es[idx])))
    best.sort(key=lambda p:-p['depth'])
    with open('pair_odd_results.json','w') as f:json.dump(best[:30],f,indent=2)
    print('Samples:',5000,'interior bound-band minima below atomic threshold:',len(best))
    print(json.dumps(best[:6],indent=2))
