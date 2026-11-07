"""C142 Task 2 checks near Delta=-1 (attractive side). Floating point only. J=1, eta=pi-gamma.

PROOF.md Part II proves, for every eta in (0,pi/2], L, n and x in [-Q_L,Q_L],
  (D_log)  eta L |sigma_L - r_Q| <= C_D(eta) := 2(1+q) eta R(0)  <= (4/pi)(log(1/eta)+pi/2),
  (D1_log) eta^2 L |sigma_L' - r_Q'| <= C_D1(eta) := eta^2 [4 M1T + 4(2+q) T(0) R(0)],
where R(0) = (1/pi) int_0^inf That/(1-That) dw is the peak of the full-line dressing kernel.
This script (a) evaluates R(0) by quadrature and checks R(0) <= (log(1/eta)+pi/2)/(pi eta) and that the
full-line kernel is unimodal, (b) evaluates the left sides for Bethe roots and compares.
"""
import json, sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from bethe_ferro import solve, Model
from continuum import Kernel, Restricted, grading

def Rhat(w, eta):
    w = abs(w)
    if w < 1e-12:
        q = 1-2*eta/np.pi
        return q/(1-q)
    g = np.pi-eta
    # That/(1-That) = sinh((pi-2eta)w/2) / (2 cosh(g w/2) sinh(eta w/2)), written stably
    A = (np.pi-2*eta)*w/2
    num = 0.5*(1-np.exp(-2*A))
    den = (1+np.exp(-g*w))*0.5*(1-np.exp(-eta*w))
    return num*np.exp(A-g*w/2-eta*w/2)/den

def R0(eta):
    f = lambda w: Rhat(w, eta)
    v = 0.
    edges = [0, 1, 1/eta, 10/eta, 60/eta]
    for a, b in zip(edges[:-1], edges[1:]):
        v += quad(f, a, b, limit=400, epsabs=1e-13, epsrel=1e-11)[0]
    return v/np.pi

def kernel_profile(eta, xs):
    # R_inf(x) = (1/pi) int_0^inf Rhat(w) cos(wx) dw, by quadrature on a fine w grid
    W = np.concatenate([np.linspace(0, 1, 4001), np.geomspace(1+1e-9, 80/eta, 40001)])
    W = np.unique(W)
    F = np.array([Rhat(w, eta) for w in W])
    out = []
    for x in xs:
        out.append(np.trapezoid(F*np.cos(W*x), W)/np.pi)
    return np.array(out)

def bounds(eta):
    q = 1-2*eta/np.pi
    r0 = R0(eta)
    T0 = 1/(np.pi*np.tan(eta))
    s = np.sin(eta)
    M1T = (2*np.cos(eta)/np.pi)*(0.32476/s**2+0.25/s)
    return dict(R0=r0, R0_bound=(np.log(1/eta)+np.pi/2)/(np.pi*eta),
                CD=2*(1+q)*eta*r0, CD_closed=(4/np.pi)*(np.log(1/eta)+np.pi/2),
                CD1=eta**2*(4*M1T+4*(2+q)*T0*r0))

def case(eta, L, n, N=700):
    mod = Model(eta); ker = Kernel(eta)
    G, x, d = solve(L, n, eta)
    if 2*n == L:
        QL = np.inf
        Xmax = x[-1]+8
        X = np.concatenate([np.linspace(0, Xmax, 8000), x[x >= 0]])
        g = np.pi-eta; b = np.pi/g
        r = 1/(2*g*np.cosh(b*X)); dr = -b*np.tanh(b*X)*r
        inner = np.ones_like(X, bool)
    else:
        f = lambda Q: mod.p(Q)-np.sum(mod.phi(Q-x))/L
        QL = brentq(f, x[-1]+1e-12, x[-1]+60, xtol=1e-15)
        R = Restricted(ker, QL, N=N, c=grading(ker, QL, N))
        X = np.concatenate([np.linspace(0, QL, 6000), x[x >= 0]])
        r = R.interp(R.r, X)
        gr = np.gradient(R.r, R.y); gr[0] = 0.0
        dr = np.interp(X, R.y, gr)
        inner = X < QL-min(eta/4, QL/4)
    sig = ker.a1(X)+ker.T(X[:, None]-x[None, :]).sum(axis=1)/L
    hh = 1e-6*max(eta, 1e-3)
    Tp = (ker.T(X[:, None]-x[None, :]+hh)-ker.T(X[:, None]-x[None, :]-hh))/(2*hh)
    dsig = (ker.a1(X+hh)-ker.a1(X-hh))/(2*hh)+Tp.sum(axis=1)/L
    return dict(QL=QL, Dmax=float(eta*L*np.max(np.abs(sig-r))),
                D1max=float(eta**2*L*np.max(np.abs(dsig-dr)[inner])))

if __name__ == '__main__':
    rows = []; bnd = {}
    EPS = [float(a) for a in sys.argv[1:]] or [0.5, 0.1, 1e-2, 1e-4, 1e-6]
    tag = '_'.join(sys.argv[1:])
    for eps in EPS:
        eta = 2*np.arcsin(np.sqrt(eps/2))
        bd = bounds(eta); bnd[eps] = bd
        xs = np.linspace(0, 6, 241)
        prof = kernel_profile(eta, xs)
        bd["unimodal_on_grid_tol1e-5"] = bool(np.all(np.diff(prof) <= 1e-5*prof[0]))
        bd['R0_quadrature_vs_profile0'] = float(prof[0]/bd['R0'])
        print(eps, json.dumps(bd), flush=True)
        for L in [16, 64, 256, 1024]:
            ns = [2, L//8, L//4, 3*L//8, L//2-2, L//2] if L < 1024 else [L//4, L//2-2, L//2]
            for n in sorted(set(ns)):
                try:
                    c = case(eta, L, n)
                except Exception as exc:
                    c = dict(error=str(exc))
                c.update(eps=eps, eta=eta, L=L, n=n, etaL=eta*L, CD=bd['CD'], CD1=bd['CD1'])
                c['ok'] = bool(('Dmax' in c) and c['Dmax'] <= bd['CD'] and c['D1max'] <= bd['CD1'])
                rows.append(c); print(json.dumps(c), flush=True)
    summ = dict(bounds={str(k): v for k, v in bnd.items()},
                max_Dmax=max(r['Dmax'] for r in rows if 'Dmax' in r),
                max_D1max=max(r['D1max'] for r in rows if 'D1max' in r),
                min_ratio_CD_over_Dmax=min(r['CD']/r['Dmax'] for r in rows if 'Dmax' in r),
                min_ratio_CD1_over_D1max=min(r['CD1']/r['D1max'] for r in rows if 'D1max' in r),
                all_ok=all(r['ok'] for r in rows), n_cases=len(rows))
    json.dump(dict(rows=rows, summary=summ), open('check_attractive%s.json' % ('_'+tag if tag else ''), 'w'), indent=1)
    print(json.dumps(summ, indent=1))
