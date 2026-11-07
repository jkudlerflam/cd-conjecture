"""C138 half-filled density diagnostics (no continuum solver needed).

At n=L/2 the comparison density is r0(x)=1/(2 gamma cosh(b x)), b=pi/gamma.
For each eta (epsilon=10^-k) and L, report
  Gcum   = L sup_x |mu(-inf,x] - 1/4 - int_0^x r0|   (C136 bound: B/2)
  Dabs   = eta L sup_roots |sigma_L - r0|             ((D*) asks this to be O(1))
  Dabs_mid = same at cell boundaries t=(j+1/2)/L      (between roots)
  D1abs  = eta^2 L sup_roots |sigma_L' - r0'|          ((D1*) asks O(1))
  C136   = D/L with D = 1+3B/(pi s): the C136 absolute bound, times eta L
  dev    = L(F_L - L e_bulk(1/2))
Usage: python3 diag_half.py -> diag_half.json
"""
import json
from pathlib import Path
import numpy as np
from bethe_ferro import solve, Model
from continuum import Kernel, half_filled_bulk

rows = []
for k in [2, 3, 4, 5, 6]:
    eps = 10.0**(-k)
    eta = 2*np.arcsin(np.sqrt(eps/2))
    mod, ker = Model(eta), Kernel(eta)
    g = np.pi-eta
    b = np.pi/g
    B = np.pi/(2*eta)
    eb = half_filled_bulk(ker)
    r0 = lambda x: 1/(2*g*np.cosh(b*x))
    dr0 = lambda x: -b*np.tanh(b*x)/(2*g*np.cosh(b*x))
    for L in [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]:
        n = L//2
        G, x, d = solve(L, n, eta)
        F = mod.Delta*L/4+G
        sig = ker.a1(x)+ker.T(x[:, None]-x[None, :]).sum(axis=1)/L
        # derivative of sigma_L: a1' + (1/L) sum T'(x-x_k)
        hh = 1e-7
        Tp = (ker.T(x[:, None]-x[None, :]+hh)-ker.T(x[:, None]-x[None, :]-hh))/(2*hh)
        a1p = (ker.a1(x+hh)-ker.a1(x-hh))/(2*hh)
        dsig = a1p+Tp.sum(axis=1)/L
        # cell boundaries: midpoints between consecutive roots in counting coordinate
        z = lambda X: (mod.p(X)-(np.sum(-np.pi*np.sign(X[:, None]-x[None, :])
                                         + mod.phi(X[:, None]-x[None, :]), axis=1))/L)/(2*np.pi)
        mids = (x[1:]+x[:-1])/2
        # refine mids to exact cell boundaries by a few Newton steps on z(X)=t
        tb = (np.arange(1, n)-n/2)/L
        X = mids.copy()
        for _ in range(30):
            s_ = ker.a1(X)+ker.T(X[:, None]-x[None, :]).sum(axis=1)/L
            X = np.clip(X-(z(X)-tb)/s_, x[:-1], x[1:])
        sigm = ker.a1(X)+ker.T(X[:, None]-x[None, :]).sum(axis=1)/L
        # cumulative error, one-sided at roots
        I0 = (np.arctan(np.tanh(b*x/2))*2/(2*g*b))   # int_0^x r0 = gd(bx)/(2 gamma b)
        I0 = np.arctan(np.sinh(b*x))/(2*g*b)
        j = np.arange(1, n+1)
        Gc = L*max(np.max(np.abs((j-1)/L-0.25-I0)), np.max(np.abs(j/L-0.25-I0)))
        dens = L*eta*r0(x)
        rows.append(dict(k=k, eta=eta, L=L, etaL=eta*L, Gcum=Gc, Gcum_over_B=Gc/B,
                         Dabs=float(eta*L*np.max(np.abs(sig-r0(x)))),
                         Dabs_mid=float(eta*L*np.max(np.abs(sigm-r0(X)))),
                         D1abs=float(eta**2*L*np.max(np.abs(dsig-dr0(x)))),
                         C136_abs=float(eta*(1+3*B/(np.pi*np.sin(eta)))),
                         rel_dense4=float(np.max(np.abs(sig/r0(x)-1)[dens >= 4])) if np.any(dens >= 4) else None,
                         dev=L*(F-L*eb), dev_over_eta=L*(F-L*eb)/eta,
                         casimir=-np.pi*(np.pi/2)*np.sin(eta)/g/6))
        print(json.dumps(rows[-1]), flush=True)
Path('diag_half.json').write_text(json.dumps(rows, indent=1))
