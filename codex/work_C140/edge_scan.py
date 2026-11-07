"""Structure of sigma_L near the edge, 0<Delta<1. Diagnostics, floating point. J=1."""
import sys, json
import numpy as np
from bethe_scan import solve
def a1f(x,g): return np.sin(g)/(np.pi*(np.cosh(2*x)-np.cos(g)))
def Kf(x,g): return np.sin(2*g)/(2*np.pi*(np.sinh(x)**2+np.sin(g)**2))
def r0f(x,g): return 1/(2*g*np.cosh(np.pi*x/g))
def sig(xs,x,g,L): return a1f(xs,g)-Kf(xs[:,None]-x[None,:],g).sum(1)/L
for delta in [0.1,0.5,0.9,0.99]:
    g=np.arccos(delta); b=np.pi/g
    for L in [16,17,64,256,1024]:
        for n in sorted(set([1,2,L//8,L//4,L//2-1,L//2])):
            if n<1: continue
            E,d,x=solve(L,n,delta,return_roots=True)
            sr=sig(x,x,g,L)
            # min over a fine grid
            xs=np.linspace(-x.max()-3*g-2, x.max()+3*g+2, 20001)
            s=sig(xs,x,g,L)
            j=np.argmin(sr)
            print(f"D={delta} L={L} n={n}: min L*sig at roots {L*sr.min():.4f} (root idx {j+1}/{n}), L*sig(x_n)/b={L*sr[-1]/b:.4f}, L*r0(x_n)/b={L*r0f(x[-1],g)/b:.4f}, grid min L*sig {L*s.min():.3e}, L*sig at x_n, x_n-1,x_n-2 /b: {[round(L*v/b,4) for v in sr[-3:][::-1]]}")
