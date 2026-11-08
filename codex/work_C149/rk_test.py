"""C149, hypothesis (R_kappa): the continuum quantile configurations y^{+-}, defined by
Z_0(y_j) = I_j/L +- kappa0/L with Z_0 the continuum counting function of r_Q at the continuum Fermi point
m(Q) = nu, are a super- and a subsolution: F(y^+) >= 0 >= F(y^-).  By the comparison principle (PROOF Thm 1.1)
this gives |Z_0(x_j) - I_j/L| <= kappa0/L for the true roots, which is also checked.
Continuum: fftsolver.Restricted2 (uniform grid h = eta/hfac, Gregory ends, FFT); Z_0 is the Nystrom counting
function (1/2pi)[p(u) - sum_k w_k r_k theta(u - x_k)], evaluated by FFT on an extended grid and
cubic Hermite interpolation.  Half filling uses r_0 in closed form.  Floating point only.
Usage: python3 rk_test.py eps L1 L2 ...   (writes rk_<eps>.json)"""
import numpy as np, sys, json, time
from scipy.fft import rfft, irfft, next_fast_len
from scipy.interpolate import CubicHermiteSpline
from bethe_ferro import Model, residual, solve
from fftsolver import Kern, solve_Q2

KAPPAS = [0.2, 0.25, 0.3, 0.35, 0.4, 0.45]

def theta(m, u):
    return -2*np.arctan(np.tanh(u)/np.tan(m.eta))

def counting_function(K, m, R, X):
    """Z_0 and s = Z_0' on the grid x extended by X on both sides (same spacing)."""
    U = R.U; h = U.h; N = U.N
    ne = int(np.ceil(X/h))
    u = -R.Q-ne*h+h*np.arange(N+1+2*ne)
    M = len(u)
    # kernel on offsets (i - k) h, i in [0, M), k in [0, N] -> offsets from -(N) - ne ... M - 1 + ne
    offs = np.arange(-(N+ne), M-ne)*h  # u_i - x_k = (i - ne - k) h, i - ne - k in [-(N+ne), M-1-ne]
    g = U.w*R.r
    Lf = next_fast_len(len(offs)+N+1)
    def conv(kern):
        c = irfft(rfft(kern, Lf)*rfft(g, Lf), Lf)
        return c[N:N+M]
    th = theta(m, offs); Tv = K.T(offs)
    Z = (m.p(u)-conv(th))/(2*np.pi)
    s = K.a1(u)+conv(Tv)
    return u, Z, s

def quantiles_from(u, Z, s, targets):
    sp = CubicHermiteSpline(u, Z, s)
    i = np.searchsorted(Z, targets)
    i = np.clip(i, 1, len(u)-1)
    y = u[i-1]+(targets-Z[i-1])/np.maximum(s[i-1], 1e-300)
    for _ in range(30):
        y = np.clip(y, u[0], u[-1])
        f = sp(y)-targets; d = sp(y, 1)
        y = y-f/d
    return y, sp

def case(eps, L, n, hfac=10.0, X=30.0):
    eta = 2*np.arcsin(np.sqrt(eps/2)); m = Model(eta); K = Kern(eta)
    nu = n/L; I = np.arange(1, n+1)-(n+1)/2; t = I/L
    out = dict(eps=eps, eta=eta, L=L, n=n)
    G, x, _ = solve(L, n, eta)
    if 2*n == L:
        g = np.pi-eta; b = np.pi/g
        Zinv = lambda tt: (2/b)*np.arctanh(np.tan(np.pi*tt))
        Z0 = lambda uu: np.arctan(np.tanh(b*uu/2))/np.pi
        out['Q'] = np.inf
        quant = lambda tt: Zinv(tt)
    else:
        R = solve_Q2(K, nu, hfac=hfac, withZ=False)
        u, Z, s = counting_function(K, m, R, X)
        out['Q'] = R.Q; out['lamQ'] = eta*L*R.r[-1]; out['lam0'] = eta*L*R.r[R.U.N//2]
        out['Zmax_room_L'] = float(L*(Z[-1]-nu/2))
        sp = CubicHermiteSpline(u, Z, s)
        Z0 = sp
        quant = lambda tt: quantiles_from(u, Z, s, tt)[0]
    out['maxdisp'] = float(L*np.max(np.abs(Z0(x)-t)))
    res = {}
    for k0 in KAPPAS:
        try:
            yp = quant(t+k0/L); ym = quant(t-k0/L)
            if not (np.all(np.diff(yp) > 0) and np.all(np.diff(ym) > 0) and np.all(np.isfinite(yp)) and np.all(np.isfinite(ym))):
                res[str(k0)] = None; continue
            Fp = residual(m, yp, L)/(2*np.pi*L)*L/eta
            Fm = residual(m, ym, L)/(2*np.pi*L)*L/eta
            res[str(k0)] = [float(Fp.min()), float(Fm.max()), int(np.argmin(Fp)), int(np.argmax(Fm))]
        except Exception as e:
            res[str(k0)] = str(e)
    out['R'] = res
    return out

if __name__ == '__main__':
    eps = float(sys.argv[1]); Ls = [int(a) for a in sys.argv[2:] if not a.startswith('tag=')]
    tag = ''.join(a[4:] for a in sys.argv[2:] if a.startswith('tag='))
    rows = []
    for L in Ls:
        for n in sorted(set([max(2, L//32), L//8, L//4, 3*L//8, L//2-L//16, L//2-8, L//2-2, L//2-1, L//2])):
            t0 = time.time()
            try:
                r = case(eps, L, n)
            except Exception as e:
                r = dict(eps=eps, L=L, n=n, error=str(e))
            r['sec'] = round(time.time()-t0, 1)
            rows.append(r); print(json.dumps(r), flush=True)
    json.dump(rows, open('rk_%g%s.json' % (eps, tag), 'w'), indent=1)
