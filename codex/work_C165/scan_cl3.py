import numpy as np, sys, time
from classical import fam, ground, diag
L=100; a=2*np.pi*13/L
kw=eval('dict('+sys.argv[1]+')'); Dzs=[float(x) for x in sys.argv[2].split(',')]
for Dz in Dzs:
    J,B=fam(a=a,Dz=Dz,**kw)
    e,n=ground(J,B,L,q0=-a,nrand=1)
    d=diag(n)
    print(f"{sys.argv[1]} Dz={Dz:.4f} e={e:.6f} |N|={d['absN']:.4f} |M|={d['absM']:.4f} kpk={d['kpk']:.4f} apk={d['apk']:.4f} Nvec={np.round(d['N'],3)}",flush=True)
