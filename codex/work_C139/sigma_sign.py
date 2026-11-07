"""Sign of the finite counting density sigma_L = a1 - K*mu for 0<Delta<1.
Diagnostics only. J=1."""
import json, sys
import numpy as np
from bethe_scan import solve

def a1f(x, g):
    return np.sin(g)/(np.pi*(np.cosh(2*x)-np.cos(g)))
def Kf(x, g):
    return np.sin(2*g)/(2*np.pi*(np.sinh(x)**2+np.sin(g)**2))

out = []
for delta in [0.1, 0.5, 0.9, 0.99]:
    g = np.arccos(delta)
    for L in [8, 9, 16, 17, 64, 65, 256, 257, 1024]:
        for frac in [0.5, 0.4, 0.25, 0.1]:
            n = int(L*frac) if frac < 0.5 else L//2
            if n < 1: continue
            E, d, x = solve(L, n, delta, return_roots=True)
            xs = np.concatenate([np.linspace(-40, 40, 40001), x])
            sig = a1f(xs, g) - Kf(xs[:, None]-x[None, :], g).sum(1)/L
            sr = a1f(x, g) - Kf(x[:, None]-x[None, :], g).sum(1)/L
            # tail coefficient: sigma ~ (2/pi) e^{-2x}[sin g - sin 2g S]
            S = np.mean(np.exp(2*x))*n/L
            tail = np.sin(g) - np.sin(2*g)*S
            rel = np.min(sig*np.exp(2*np.abs(xs)))
            out.append(dict(delta=delta, L=L, n=n, min_sigma_roots=float(sr.min()),
                            min_sigma_scaled=float(rel), tail_coef=float(tail),
                            xmax=float(x.max())))
            print(delta, L, n, 'min sig at roots', sr.min(), 'min sig*e^{2|x|}', rel, 'tail', tail, flush=True)
json.dump(out, open('sigma_sign.json', 'w'), indent=1)
