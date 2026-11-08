"""C149 exploration: residual of the explicit half-filled quantile configuration.
y_j = (2/b) atanh(tan(pi I_j/L)) are the quantiles of r_0 = 1/(2 gamma cosh bx), D == 0 exactly,
so Phi_j = -F_j(y)/(2 pi L). Floating point only."""
import numpy as np, sys
from bethe_ferro import Model, residual, solve
for eps in [1e-2,1e-4,1e-6]:
    eta = 2*np.arcsin(np.sqrt(eps/2)); m = Model(eta); g=np.pi-eta; b=np.pi/g
    for L in [64,512,4096]:
        n=L//2; I=np.arange(1,n+1)-(n+1)/2
        y=(2/b)*np.arctanh(np.tan(np.pi*I/L))
        Phi=-residual(m,y,L)/(2*np.pi*L)
        r0=1/(2*g*np.cosh(b*y)); lam=eta*L*r0
        pred=-eta**2*b*np.tanh(b*y)*r0/np.expm1(2*np.pi*lam)   # eta^2 rho'/(e^{2 pi lam}-1)
        j=np.argmax(np.abs(Phi))
        print(f"eps={eps:g} L={L} max L|Phi|/eta={L*np.abs(Phi).max()/eta:.4f} at lam={lam[j]:.3g}; b/(2pi)={b/(2*np.pi):.4f}; max|Phi-pred|L/eta={L*np.abs(Phi-pred).max()/eta:.3g}; lam(0)={lam[n//2]:.3g}")
