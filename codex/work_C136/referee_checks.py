"""Independent referee checks for C136 (floating point, not enclosures).

Run with PYTHONDONTWRITEBYTECODE=1 from this folder.
Uses the local Bethe solver bethe_scan.solve only for the finite roots.
The continuum objects (r_Q, Z_Q, m, E, h, epsilon_Q) are recomputed here by
an independent composite Gauss-Legendre Nystrom discretization.
"""
import json
import sys
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
from scipy.special import roots_legendre

from bethe_scan import solve


def setup(delta):
    g = np.arccos(delta)
    s = np.sin(g)
    a, b = np.pi*s, np.pi/g
    q = 2*g/np.pi-1
    B = 1/(1-q)
    a1 = lambda x: s/(np.pi*(np.cosh(2*x)-delta))
    a1p = lambda x: -2*s*np.sinh(2*x)/(np.pi*(np.cosh(2*x)-delta)**2)
    T = lambda x: -np.sin(2*g)/(2*np.pi*(np.sinh(x)**2+s**2))
    Tp = lambda x: np.sin(2*g)*np.sinh(2*x)/(2*np.pi*(np.sinh(x)**2+s**2)**2)
    e = lambda x: -a*a1(x)
    r0 = lambda x: 1/(2*g*np.cosh(b*x))
    return dict(g=g, s=s, a=a, b=b, q=q, B=B, a1=a1, a1p=a1p, T=T, Tp=Tp, e=e, r0=r0)


def nodes(Q, panels=None, order=40):
    if panels is None:
        panels = max(4, int(np.ceil(2*Q/0.15)))
    x, w = roots_legendre(order)
    edges = np.linspace(-Q, Q, panels+1)
    ys, ws = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        ys.append((hi-lo)/2*x+(hi+lo)/2)
        ws.append((hi-lo)/2*w)
    return np.concatenate(ys), np.concatenate(ws)


class Cont:
    """Restricted continuum at boundary Q."""
    def __init__(self, P, Q):
        self.P, self.Q = P, Q
        y, w = nodes(Q)
        self.y, self.w = y, w
        A = np.eye(len(y))-P['T'](y[:, None]-y[None, :])*w
        self.r, self.Z = np.linalg.solve(A, np.stack([P['a1'](y), np.ones_like(y)], 1)).T
        self.m = float(w@self.r)
        self.E = float(np.cos(P['g'])/4+w@(P['e'](y)*self.r))
        self.rQ = float(self.rf(Q))
        self.ZQ = float(self.Zf(Q))
        self.h = -P['a']*self.rQ/self.ZQ

    def rf(self, x):
        x = np.asarray(x, float)
        return self.P['a1'](x)+self.P['T'](x[..., None]-self.y)@(self.w*self.r)

    def rpf(self, x):
        x = np.asarray(x, float)
        return self.P['a1p'](x)+self.P['Tp'](x[..., None]-self.y)@(self.w*self.r)

    def Zf(self, x):
        x = np.asarray(x, float)
        return 1+self.P['T'](x[..., None]-self.y)@(self.w*self.Z)

    def Zpf(self, x):
        x = np.asarray(x, float)
        return self.P['Tp'](x[..., None]-self.y)@(self.w*self.Z)

    def epsf(self, x):
        return -self.P['a']*self.rf(x)-self.h*self.Zf(x)

    def epspf(self, x):
        return -self.P['a']*self.rpf(x)-self.h*self.Zpf(x)


_bulk_cache = {}


def ebulk(P, delta, nu):
    key = (delta, nu)
    if key in _bulk_cache:
        return _bulk_cache[key]
    if nu == 0:
        return delta/4, None
    if nu >= .5:
        g, s = P['g'], P['s']
        val = delta/4+quad(lambda x: P['e'](x)*P['r0'](x), -np.inf, np.inf, epsabs=1e-15, epsrel=1e-14)[0]
        _bulk_cache[key] = (val, None)
        return val, None
    hi = 1.
    while Cont(P, hi).m < nu:
        hi *= 1.6
    Qn = brentq(lambda Q: Cont(P, Q).m-nu, 1e-9, hi, xtol=1e-14, rtol=1e-14)
    C = Cont(P, Qn)
    _bulk_cache[key] = (C.E, C)
    return C.E, C


def all_constant(gap):
    alpha = np.arccos(1-gap); s = np.sqrt(gap*(2-gap)); B = np.pi/(2*alpha)
    mb = np.sqrt(2)*np.sin(np.pi/4-alpha/2)
    C = 4/(np.pi*s*(1-mb))
    H1 = 2/s**2; H2 = 8/s**4+4/s**2
    c = 2/3*np.log(max(1., 16*C/3))
    t = s/(16*np.pi*np.cosh(c+1)**2)
    kap = min(.5, s/(np.cosh(2*(c+1))+1), t/2*np.exp(-2*c))
    P0 = 8*C/3; P1 = P0*H1; P2 = P0*H2
    W = (P2+4*P0+2*P1+H1*(P1+2*P0))/kap
    Vv = np.pi*(4+4*W/3+32*C*H1/(3*kap))
    D = 1+3*B/(np.pi*s); D1 = 1.5*B*H2
    Vclip = 2*Vv+8*np.pi*(H1+D1/D)+8*np.pi
    K = max((4*np.pi*D/kap+1)**2, Vclip/8+np.pi*(128*D**2/kap**2+72*D*B/kap+B**2/2))
    return dict(K=K, kappa=kap, C0=C, Vv=Vv, c=c)


def finite_case(P, delta, L, n):
    E, diag, x = solve(L, n, delta, return_roots=True)
    E = float(E)
    nu = n/L
    eb, Cnu = ebulk(P, delta, nu)
    out = dict(L=L, n=n, nu=nu, Lambda=L*(E-L*eb))
    if not (0 < nu < .5):
        return out
    T, a1 = P['T'], P['a1']
    p = lambda y: 2*np.arctan(np.tanh(y)/np.tan(P['g']/2))
    th = lambda y: 2*np.arctan(np.tanh(y)/np.tan(P['g']))
    z = lambda y: (p(y)-np.sum(th(np.asarray(y)[..., None]-x), -1)/L)/(2*np.pi)
    sig = lambda y: a1(y)+np.sum(T(np.asarray(y)[..., None]-x), -1)/L
    hi = float(x[-1])+.5
    while z(hi) <= nu/2:
        hi += 1.
    QL = brentq(lambda y: z(y)-nu/2, float(x[-1]), hi, xtol=1e-14)
    C = Cont(P, QL)
    out['QL'] = QL
    out['mass_mismatch_ratio'] = L*abs(nu-C.m)/(P['B']-1) if P['B'] > 1 else 0.
    # exact identity (4.3)
    yq, wq = nodes(QL)
    rhs = np.sum(C.epsf(x))-L*(wq@(C.epsf(yq)*sig(yq)))
    out['identity_err'] = abs(E-L*(C.E+C.h*(nu-C.m))-rhs)
    # cumulative G sup vs 3B/(2L)
    grid = np.linspace(-QL, QL, 4001)
    rg = C.rf(grid)
    cum = np.concatenate([[0], np.cumsum((rg[1:]+rg[:-1])/2*np.diff(grid))])
    before = np.searchsorted(x, grid, 'left')/L
    after = np.searchsorted(x, grid, 'right')/L
    Gsup = max(np.max(abs(before-cum)), np.max(abs(after-cum)), abs(nu-C.m))
    out['G_ratio'] = Gsup/(1.5*P['B']/L)
    # kappa comparison and velocity variation
    out['min_rQ_over_r0'] = float(np.min(rg/P['r0'](grid)))
    v = C.epspf(grid)/rg
    out['Var_v'] = float(np.sum(abs(np.diff(v))))
    out['curv'] = float(C.epspf(QL)/(2*C.rQ*C.ZQ**2))/P['a']
    # midpoint quadrature V_L/8 actual
    out['L_quad'] = L*abs(rhs)
    return out


def main():
    res = dict()
    # (6.2) exponential moment and R_infinity decay
    mom = {}
    for delta in [-0.999, -0.99, -0.9, -0.5, -0.1]:
        P = setup(delta)
        num = quad(lambda t: np.exp(.5*t)*P['T'](t), -80, 80, epsabs=1e-14, limit=500, points=[0.])[0]
        claim = np.sqrt(2)*np.sin(P['g']/2-np.pi/4)
        # full-line resolvent kernel by Fourier: Rhat = 1/(1+Khat)-1
        g = P['g']
        def Rhat(w):
            w = np.where(abs(w) < 1e-12, 1e-12, w)
            return np.sinh(np.pi*w/2)/(2*np.sinh((np.pi-g)*w/2)*np.cosh(g*w/2))-1
        xs = np.array([0., 2., 5., 10.])
        Rx = [quad(lambda w: Rhat(w)*np.cos(w*xx)/np.pi, 0, 200, limit=2000)[0] for xx in xs]
        T0 = abs(delta)/(np.pi*P['s'])
        C0pt = 4*T0/(1-claim)
        mom[delta] = dict(numeric=num, claim=claim,
                          R_inf_over_C0exp=[float(R/(C0pt*np.exp(-.5*xx))) for R, xx in zip(Rx, xs)],
                          R_inf_exp_b=[float(R*np.exp(P['b']*xx)) for R, xx in zip(Rx, xs)])
    res['moment'] = mom
    if '--moment-only' in sys.argv:
        print(json.dumps(mom, indent=1)); return
    # (6.13) majorant
    ds = np.logspace(-6, 0, 400)
    res['Kall_ratio_max'] = float(max(all_constant(d)['K']*d**(25/3)/1e15 for d in ds))
    # finite checks
    finite = {}
    for delta in [-0.5, -0.9]:
        P = setup(delta)
        rows = []
        for L in [16, 17, 32, 33, 64, 65, 128, 129, 256, 257, 400, 401]:
            ns = sorted(set([1, 2, 3, int(np.sqrt(L)), L//8, L//4, 3*L//8,
                             L//2-int(np.sqrt(L)), L//2-3, L//2-2, L//2-1, L//2]))
            for n in ns:
                if not 1 <= n <= L//2:
                    continue
                row = finite_case(P, delta, L, n)
                if 2*L <= 520:
                    E2 = float(solve(2*L, 2*n, delta)[0])
                    E1 = float(solve(L, n, delta)[0])
                    row['defect'] = L*(E2-2*E1)
                rows.append(row)
                print(delta, L, n, {k: (round(v, 6) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
        ca = all_constant(1+delta)
        P0 = setup(delta)
        finite[delta] = dict(rows=rows, kappa_compact=ca['kappa'], Vv_compact=ca['Vv'],
                             K_all=ca['K'])
    res['finite'] = finite
    with open('referee_checks.json', 'w') as f:
        json.dump(res, f, indent=1, default=float)
    # summary
    for delta, d in finite.items():
        rows = d['rows']
        print('Delta', delta)
        print('  max |Lambda|', max(abs(r['Lambda']) for r in rows))
        for key in ['mass_mismatch_ratio', 'G_ratio', 'Var_v', 'curv', 'L_quad', 'identity_err']:
            vals = [r[key] for r in rows if key in r]
            print('  max', key, max(vals))
        print('  min r_Q/r0', min(r['min_rQ_over_r0'] for r in rows if 'min_rQ_over_r0' in r), 'kappa', d['kappa_compact'])
        print('  max defect', max(r.get('defect', -1e9) for r in rows))
    print(json.dumps(res['moment'], indent=1))
    print('Kall ratio', res['Kall_ratio_max'])


if __name__ == '__main__':
    main()
