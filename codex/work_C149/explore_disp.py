import numpy as np
from bethe_ferro import Model, residual, solve
for eps in [1e-2,1e-4,1e-6]:
    eta = 2*np.arcsin(np.sqrt(eps/2)); m = Model(eta); g=np.pi-eta; b=np.pi/g
    for L in [64,512,4096]:
        n=L//2; I=np.arange(1,n+1)-(n+1)/2
        G,x,d=solve(L,n,eta)
        Z0=lambda u: np.arctan(np.tanh(b*u/2))/np.pi
        disp=L*(Z0(x)-I/L)
        lam=eta*L/(2*g*np.cosh(b*x))
        k=np.argmax(np.abs(disp))
        sel=lam>=1
        print(f"eps={eps:g} L={L} max|disp|={np.abs(disp).max():.4f} at lam={lam[k]:.3g} x={x[k]:.3f}; max on lam>=1: {np.abs(disp[sel]).max() if sel.any() else 0:.4f}; eta*B*b/(2pi)... ={b/4:.3f}")
