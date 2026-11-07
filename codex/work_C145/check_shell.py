"""C145: the remaining shell inequality (G). Floating point only. J=1.
lambda(x) = eta L r_Q(x). K0 = {lambda >= M0}, M0 = 2; K1 = {lambda >= M1}, M1 = 4 C_D(eta) (proved comparability there);
shell S = K0 minus K1. Reports min sigma/r on K0, eta*Var_K0(r/sigma), eta*Var_S(r/sigma), and the proved
log bound for eta*Var_K1(r/sigma). Usage: python3 check_shell.py eps [eps ...]"""
import json, sys
import numpy as np
from scipy.optimize import brentq
from bethe_ferro import solve, Model
from continuum import Kernel, Restricted, grading
from check_attractive import bounds

def sigma_on(ker, X, x, L, chunk=4000):
    out = np.empty_like(X)
    for i in range(0, len(X), chunk):
        Xi = X[i:i+chunk]
        out[i:i+chunk] = ker.a1(Xi)+ker.T(Xi[:, None]-x[None, :]).sum(axis=1)/L
    return out

def case(eta, L, n, M0, M1):
    mod = Model(eta); ker = Kernel(eta); g = np.pi-eta; b = np.pi/g
    G, x, d = solve(L, n, eta)
    if 2*n == L:
        QL = np.inf; rfun = lambda X: 1/(2*g*np.cosh(b*X))
    else:
        f = lambda Q: mod.p(Q)-np.sum(mod.phi(Q-x))/L
        QL = brentq(f, x[-1]+1e-12, x[-1]+60, xtol=1e-15)
        R = Restricted(ker, QL, N=800, c=grading(ker, QL, 800))
        rfun = lambda X: R.interp(R.r, X)
    # x_c: lambda(x_c) = M0 on [0, Q]
    lam = lambda X: eta*L*rfun(X)
    top = QL if np.isfinite(QL) else x[-1]+10
    if lam(0.0) < M0:
        return dict(QL=QL, empty=True)
    xc0 = top if lam(top) >= M0 else brentq(lambda X: lam(X)-M0, 0, top, xtol=1e-14)
    xc1 = 0.0 if lam(0.0) < M1 else (top if lam(top) >= M1 else brentq(lambda X: lam(X)-M1, 0, top, xtol=1e-14))
    sig_max = sigma_on(ker, np.array([0.0]), x, L)[0]
    h = min(eta, 1/(L*sig_max))/25
    X = np.linspace(0, xc0, int(np.ceil(xc0/h))+2)
    sig = sigma_on(ker, X, x, L); r = rfun(X); gq = r/sig
    dv = np.abs(np.diff(gq))
    inS = X[1:] > xc1
    return dict(QL=float(QL), xc0=float(xc0), xc1=float(xc1), npts=len(X),
                min_sig_over_r_K0=float(np.min(sig/r)), max_sig_over_r_K0=float(np.max(sig/r)),
                etaVar_K0=float(2*eta*dv.sum()), etaVar_S=float(2*eta*dv[inS].sum()),
                etaVar_K1=float(2*eta*dv[~inS].sum()))

if __name__ == '__main__':
    M0 = 2.0; rows = []
    EPS = [float(a) for a in sys.argv[1:]]
    for eps in EPS:
        eta = 2*np.arcsin(np.sqrt(eps/2)); bd = bounds(eta); M1 = 4*bd['CD']
        for L in [64, 256, 1024, 2048]:
            for n in sorted(set([L//8, L//4, 3*L//8, L//2-2, L//2])):
                try:
                    c = case(eta, L, n, M0, M1)
                except Exception as exc:
                    c = dict(error=str(exc))
                c.update(eps=eps, eta=eta, L=L, n=n, etaL=eta*L, M1=M1)
                rows.append(c); print(json.dumps(c), flush=True)
    json.dump(rows, open('check_shell_%s.json' % '_'.join(sys.argv[1:]), 'w'), indent=1)
