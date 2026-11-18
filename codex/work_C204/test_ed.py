import numpy as np, itertools
from c204lib import *
def ed(L, n, eta):
    # H' = -sum (SxSx+SySy + cos(eta) SzSz), periodic; sector with n down spins
    states=[s for s in range(1<<L) if bin(s).count('1')==n]
    idx={s:i for i,s in enumerate(states)}
    H=np.zeros((len(states),)*2)
    c=np.cos(eta)
    for i,s in enumerate(states):
        for b in range(L):
            a,bb=b,(b+1)%L
            sa=(s>>a)&1; sb=(s>>bb)&1
            H[i,i]+= -c*(0.25 if sa==sb else -0.25)
            if sa!=sb:
                t=s^(1<<a)^(1<<bb)
                H[idx[t],i]+=-0.5
    return np.linalg.eigvalsh(H)[0]
print(gregory_weights(6))
for eta in [0.3,1.0]:
    K=Kern(eta)
    for L in [5,6,7,8,9,10]:
        for n in range(1,L//2+1):
            x,res,it=bethe_solve(eta,L,n)
            Fb=-np.cos(eta)*L/4 + np.sum(-PI*np.sin(eta)*K.a1(x))
            print(eta,L,n,'%.12f %.12f %.1e'%(Fb,ed(L,n,eta),abs(Fb-ed(L,n,eta))),it)
