import numpy as np, sys
from classical import fam, ground, diag
L=100; k=int(sys.argv[3]); a=2*np.pi*k/L
kw=eval('dict('+sys.argv[1]+')')
for Dz in [float(x) for x in sys.argv[2].split(',')]:
    J,B=fam(a=a,Dz=Dz,**kw)
    e,n=ground(J,B,L,q0=-a,nrand=1)
    d=diag(n)
    print(f"a={a:.3f} {sys.argv[1]} Dz={Dz:.3f} e={e:.5f} |N|={d['absN']:.3f} |M|={d['absM']:.3f} kpk={d['kpk']:.3f} apk={d['apk']:.3f}",flush=True)
