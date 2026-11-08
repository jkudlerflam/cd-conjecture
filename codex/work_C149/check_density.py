"""C149: consequences of the sandwich on the true roots. Checks, for Bethe roots x,
  (a) L sup_v |N_x(v) - N_0(v)|  <=  kappa0 + 1/2   (cumulative bound, Lemma 3.1),
  (b) eta L sup |sigma_L - s_nu|  vs  (2 kappa0 + 1) eta cot(eta)/pi   (Theorem 3.3 (i)),
  (c) eta^2 L sup |sigma_L' - s_nu'| vs (2 kappa0 + 1) 2 eta^2 T_A(0) T_B(0)  (Theorem 3.3 (ii)),
  (d) the sign of Q_L - Q_nu and L (Z_0(Q_L) - nu/2),
  (e) the factorization T = T_A * T_B at sample points.
s_nu is the Nystrom counting density of r_Q at m(Q) = nu (fftsolver).  Floating point only."""
import numpy as np, sys, json
from scipy.optimize import brentq
from bethe_ferro import Model, solve
from fftsolver import Kern, solve_Q2
from rk_test import counting_function

def TA(K, x):
    g = K.gamma; ep = np.pi*K.eta/g
    return np.sin(ep)/(g*(np.cosh(2*np.pi*x/g)-np.cos(ep)))
def TB(K, x):
    return np.sin(K.eta)/(np.pi*(np.cosh(2*x)-np.cos(K.eta)))

def case(eps, L, n, k0=0.4):
    eta = 2*np.arcsin(np.sqrt(eps/2)); m = Model(eta); K = Kern(eta)
    nu = n/L; G, x, _ = solve(L, n, eta)
    R = solve_Q2(K, nu, hfac=10)
    u, Z, s = counting_function(K, m, R, 6.0)
    # derivative s' by convolution with T' (Nystrom)
    h = u[1]-u[0]
    sp = np.gradient(s, h, edge_order=2)
    sel = np.abs(u) <= R.Q+3*eta
    U = u[sel]
    sig = K.a1(U)+np.concatenate([K.T(U[i:i+4000][:, None]-x[None, :]).sum(1) for i in range(0, len(U), 4000)])/L
    a1p = np.gradient(K.a1(U), h)
    def Tp(d):
        t = np.tanh(d); return -2*K.te*t*(1-t*t)*(1+K.te**2)/(np.pi*(t*t+K.te**2)**2)
    sigp = a1p+np.concatenate([Tp(U[i:i+4000][:, None]-x[None, :]).sum(1) for i in range(0, len(U), 4000)])/L
    out = dict(eps=eps, L=L, n=n, Q=R.Q)
    # (a) cumulative: N_0 on grid vs N_x
    N0 = np.clip(Z+nu/2, 0, nu)
    Nx = np.searchsorted(x, u, side='right')/L
    out['L_cum'] = float(L*np.max(np.abs(Nx-N0))); out['cum_bound'] = k0+0.5
    out['D'] = float(eta*L*np.max(np.abs(sig-s[sel]))); out['D_bound'] = (2*k0+1)*eta/np.tan(eta)/np.pi
    out['D1'] = float(eta**2*L*np.max(np.abs(sigp-sp[sel]))); out['D1_bound'] = (2*k0+1)*2*eta**2*TA(K, 0.0)*TB(K, 0.0)
    if 2*n < L:
        f = lambda X: m.p(X)-np.sum(m.phi(X-x))/L
        QL = brentq(f, x[-1]+1e-12, x[-1]+60, xtol=1e-15)
        Zq = np.interp(QL, u, Z)
        out['QL_minus_Qnu'] = float(QL-R.Q); out['L_Z0QL_minus'] = float(L*(Zq-nu/2))
    # (e) factorization check at a few points via FFT-free quadrature
    yy = np.linspace(-40*eta, 40*eta, 200001); dy = yy[1]-yy[0]
    for X in [0.0, eta, 5*eta]:
        conv = np.sum(TA(K, X-yy)*TB(K, yy))*dy
        out.setdefault('fact_relerr', []).append(float(abs(conv-K.T(X))/K.T(X)))
    return out

if __name__ == '__main__':
    eps = float(sys.argv[1]); rows = []
    for L in [int(a) for a in sys.argv[2:]]:
        for n in [L//8, L//4, 3*L//8, L//2-8, L//2-1]:
            r = case(eps, L, n); rows.append(r); print(json.dumps(r), flush=True)
    json.dump(rows, open('check_density_%g.json' % eps, 'w'), indent=1)
