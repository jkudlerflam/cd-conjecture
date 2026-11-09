"""C186: density onset of the mixed second-order cusp model on rings (all particle numbers, exact).

For each mu the non-vacuum ground energy E*(mu) = min over K of the lowest level that is not the exact vacuum.
The onset mu_on(L) is where E*(mu) = 0. The object that appears first is characterised by its momentum K and
<N>. Compared with the dressed one-particle onset mu_1 (infinite line, few-body code, converged in N and box).
A first object with <N> >= 2 and mu_on < mu_1 beyond finite-size effects means bound pairs or clusters.
"""
import numpy as np, json, sys, time
from scipy.optimize import brentq
from ring_mixed import Ring, model


def nonvac_lowest(ring, m):
    L = ring.L
    best = (np.inf, None, None)
    for kk in range(L):
        K = 2 * np.pi * kk / L
        e, Nx = ring.lowest(m, K, k=2)
        for ei, ni in zip(e, Nx):
            if kk == 0 and abs(ei) < 1e-9 and ni < 1e-9:
                continue                     # the exact vacuum
            if ei < best[0]:
                best = (ei, K, ni)
            break
    return best


def onset(L, u, phi, v, r=None, g=2/3, lo=None, hi=None):
    ring = Ring(L)
    def f(mu):
        return nonvac_lowest(ring, model(u, phi, v, mu, g=g, r=r))[0]
    eps_min = -2 - u ** 2 / 4
    lo = eps_min - 1.5 if lo is None else lo
    hi = eps_min + 0.6 if hi is None else hi
    flo, fhi = f(lo), f(hi)
    if flo < 0:
        return dict(L=L, u=u, phi=phi, v=v, r=(np.sqrt(v) if r is None else r), mu_on=None, note=f"E*<0 already at mu={lo}")
    mu_on = brentq(f, lo, hi, xtol=1e-6)
    e, K, N = nonvac_lowest(ring, model(u, phi, v, mu_on, g=g, r=r))
    # also the object just above onset
    e2, K2, N2 = nonvac_lowest(ring, model(u, phi, v, mu_on + 0.02, g=g, r=r))
    return dict(L=L, u=u, phi=phi, v=v, r=(np.sqrt(v) if r is None else r), mu_on=mu_on,
                K_on=K, N_on=N, E_above=e2, K_above=K2, N_above=N2, eps_min=eps_min)


if __name__ == "__main__":
    L = int(sys.argv[1])
    u = float(sys.argv[2])
    vs = [float(x) for x in sys.argv[3].split(',')]
    rmode = sys.argv[4] if len(sys.argv) > 4 else "locked"      # locked (r = sqrt v) or zero
    phi = float(sys.argv[5]) if len(sys.argv) > 5 else np.pi / 2
    for v in vs:
        t0 = time.time()
        out = onset(L, u, phi, v, r=(0.0 if rmode == "zero" else None))
        out["rmode"] = rmode; out["sec"] = round(time.time() - t0, 1)
        print(json.dumps(out), flush=True)
