import numpy as np
from bethe_scan import solve
def a1f(x,g): return np.sin(g)/(np.pi*(np.cosh(2*x)-np.cos(g)))
def Kf(x,g): return np.sin(2*g)/(2*np.pi*(np.sinh(x)**2+np.sin(g)**2))
def r0f(x,g): return 1/(2*g*np.cosh(np.pi*x/g))
for delta in [0.5,0.9,0.99]:
  g=np.arccos(delta); b=np.pi/g
  for L,n in [(256,128),(256,127),(256,100),(255,127)]:
    E,d,x=solve(L,n,delta,return_roots=True)
    xn=x[-1]
    ys=xn+np.array([-3,-2,-1.5,-1,-0.5,0,0.25,0.5,1,1.5,2,3,4,6,8,12])/b
    s=a1f(ys,g)-Kf(ys[:,None]-x[None,:],g).sum(1)/L
    print(f"D={delta} L={L} n={n} xn*b={xn*b:.3f} gaps*b: {np.round(np.diff(x[-4:])*b,3)}")
    print("  (y-xn)b :", np.round((ys-xn)*b,2))
    print("  L sig/b :", np.round(L*s/b,4))
    print("  sig/r0  :", np.round(s/r0f(ys,g),4))
