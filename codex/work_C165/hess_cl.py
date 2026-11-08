import numpy as np, sys
from classical import fam, ground, diag, hessian
L=100; a=2*np.pi*13/L
kw=eval('dict('+sys.argv[1]+')')
for Dz in [float(x) for x in sys.argv[2].split(',')]:
    J,B=fam(a=a,Dz=Dz,**kw)
    e,n=ground(J,B,L,q0=-a,nrand=0)
    H,T=hessian(n,J,B)
    w,v=np.linalg.eigh(H)
    k=2*np.pi*np.fft.fftfreq(L); wpi=[]
    for i in range(len(w)):
        dn=np.einsum('ik,ika->ia',v[:,i].reshape(L,2),T)
        F=np.fft.fft(dn,axis=0); amp=np.sum(np.abs(F)**2,1)
        wpi.append(amp[np.abs(np.abs(k)-np.pi)<0.07].sum()/amp.sum())
    wpi=np.array(wpi); ip=np.where(wpi>0.5)[0][0]
    d=diag(n)
    print(f"{sys.argv[1]} Dz={Dz} e={e:.6f} |N|={d['absN']:.3f} w0={w[0]:.2e} w1={w[1]:.2e} lowest-pi-mode w={w[ip]:.5f} (pi weight {wpi[ip]:.2f}, index {ip})",flush=True)
