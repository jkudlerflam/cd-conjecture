"""C200: general filling, trial y^+ = quantiles of the continuum counting function Z_0 (at Q_nu) at levels
I_j/L + kappa0/L, kappa0 possibly > 1/2; nodes whose level is not attained are at +infinity.
F_j evaluated with the step-free Bethe map.  Output: min_j F_j/(2 pi eta) and the index.  Floating point only."""
import numpy as np, sys, json, time
from scipy.fft import rfft, irfft, next_fast_len
from scipy.interpolate import CubicHermiteSpline
from bethe_ferro import Model
from fftsolver import Kern, solve_Q2

def counting_function(K, m, R, X):
    U = R.U; h = U.h; N = U.N
    ne = int(np.ceil(X/h))
    u = -R.Q-ne*h+h*np.arange(N+1+2*ne); M = len(u)
    offs = np.arange(-(N+ne), M-ne)*h
    g = U.w*R.r
    Lf = next_fast_len(len(offs)+N+1)
    conv = lambda kern: irfft(rfft(kern, Lf)*rfft(g, Lf), Lf)[N:N+M]
    th = -2*np.arctan(np.tanh(offs)/np.tan(m.eta))
    Z = (m.p(u)-conv(th))/(2*np.pi); s = K.a1(u)+conv(K.T(offs))
    return u, Z, s

def Fmap(m, y, L, ninf):
    n = len(y); F = L*m.p(y)+2*m.eta*ninf
    for s in range(0, n, 400):
        d = y[s:s+400, None]-y[None, :]
        idx = np.arange(s, min(s+400, n)); d[idx-s, idx] = 1.0
        ph = m.phi(d); ph[idx-s, idx] = 0.0
        F[s:s+400] -= ph.sum(axis=1)
    return F

def case(eps, L, n, kappas, X=25.0):
    eta = np.arccos(1-eps); m = Model(eta); K = Kern(eta); nu = n/L
    I = np.arange(1, n+1)-(n+1)/2
    R = solve_Q2(K, nu, hfac=10.0, withZ=False)
    u, Z, s = counting_function(K, m, R, X)
    sp = CubicHermiteSpline(u, Z, s)
    out = dict(eps=eps, L=L, n=n, Q=R.Q, lamQ=eta*L*R.r[-1], lam0=eta*L*R.r[R.U.N//2], res={})
    for k0 in kappas:
        lev = I/L+k0/L
        fin = lev < Z[-1]
        tg = lev[fin]
        i = np.clip(np.searchsorted(Z, tg), 1, len(u)-1)
        y = u[i-1]+(tg-Z[i-1])/np.maximum(s[i-1], 1e-300)
        for _ in range(40):
            y = np.clip(y, u[0], u[-1]); y = y-(sp(y)-tg)/sp(y, 1)
        ninf = n-fin.sum()
        F = Fmap(m, y, L, ninf)/(2*np.pi*eta)
        jm = int(np.argmin(F))
        out['res'][str(k0)] = [float(F.min()), jm+1, int(ninf), float(F[-1]), float((y[jm]-R.Q))]
    return out

if __name__ == '__main__':
    eps = float(sys.argv[1]); Ls = [int(a) for a in sys.argv[2:]]
    ks = [0.4, 1.0, 2.0, 5.0]
    for L in Ls:
        for n in sorted(set([max(2, L//32), L//8, L//4, 3*L//8, L//2-L//16, L//2-8, L//2-2, L//2-1])):
            if n < 8: continue
            r = case(eps, L, n, ks)
            print(f"eps={eps:g} L={L} n={n} Q={r['Q']:.3f} lamQ={r['lamQ']:.3f} lam0={r['lam0']:.2f} | " +
                  ' | '.join(f"k0={k}: {v[0]:.3f}@{v[1]} ninf={v[2]} top={v[3]:.3f} y-Q={v[4]:.3f}" for k, v in r['res'].items()), flush=True)
