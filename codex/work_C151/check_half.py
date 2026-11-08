import numpy as np
from bethe_ferro import solve, Model
from continuum import Kernel, half_filled_bulk
from ddv import DdV, path_integral
for eta in [np.pi/2-1e-9, 0.5, 0.1, 0.01]:
  for L in [16, 64, 256]:
    n = L//2
    G, x, d = solve(L, n, eta)
    D = DdV(eta, L, x)
    I = np.arange(1, n+1)-(n+1)/2
    err = np.max(np.abs(D.theta(x)-I))
    g = np.pi-eta; b = np.pi/g; a = np.pi*np.sin(eta)
    fp = lambda z: a*b*np.sinh(b*z)/(2*g*np.cosh(b*z)**2)   # f = -a/(2g cosh bz), f' = a b sinh/(2g cosh^2)
    K = Kernel(eta)
    F = K.Delta*L/4+G
    dev = F - L*half_filled_bulk(K)
    X = x[-1]+12
    res = []
    for ep in [eta/8, eta/4, eta/2]:
        h = min(eta/8, 0.25/(L*0.5))
        F_ = lambda z: fp(z)*D.ell(z)
        val = -np.imag(path_integral(F_, [X+1j*ep, -X+1j*ep], h))/np.pi
        res.append(val)
    print(f"eta={eta:.3g} L={L} theta-I err={err:.1e} L*dev={L*dev:.8f} DdV={[round(L*v,8) for v in res]}")
