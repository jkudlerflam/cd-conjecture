import numpy as np, sys, time
from classical import fam, ground, diag
L=100; a=2*np.pi*13/L
t0=time.time()
for g in (0.0, 0.3, 0.6):
    for Dz in np.round(np.arange(0.9,1.5,0.05),3):
        J,B=fam(a=a,Dz=Dz,g=g)
        e,n=ground(J,B,L,q0=-a,nrand=2)
        d=diag(n)
        print(f"g={g} Dz={Dz:.3f} e={e:.6f} |N|={d['absN']:.4f} |M|={d['absM']:.4f} kpk={d['kpk']:.4f} apk={d['apk']:.4f}",flush=True)
print('t',time.time()-t0)
