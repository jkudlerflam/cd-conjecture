"""C147 referee: repulsive side, every new inequality of C142 Part I on Bethe roots, own code.
Restricted continuum by Gauss-Legendre Nystrom. J=1."""
import json, sys
import numpy as np
from scipy.optimize import brentq
from ref_core import Rep, bethe_rep, energy_rep

def consts(g):
    pi = np.pi; q = 1-2*g/pi; qJ = (pi-2*g)/(2*(pi-g)); s = 1-g/pi
    a = pi*np.sin(g); b = pi/g; ab = a*b
    gG = (2+q)*(1+qJ)/(1-qJ)*s*(1+qJ*q/2)
    D = 2*gG/(pi*g); sn = np.sin(g); M1 = (2/pi)*(9/(16*np.sqrt(3))/sn**2+0.25/sn); D1 = 4*M1*gG
    cr = 1+4/pi; Vv = (2+40/pi)*ab; Vf = ab*(4+4*D1/(b*D))+2*Vv
    a2 = ab*(1+q)**2/2
    K = Vf/8 + 2*s*ab*gG*(2*cr/pi+2) + a2*(qJ*s)**2/2
    return dict(q=q, qJ=qJ, s=s, g=gG, D=D, D1=D1, K=K, theta1=1-2*(2-np.log(2))*(pi-g)/pi**2)

class Restr:
    def __init__(self, M, Q, N=300):
        self.M, self.Q = M, Q
        t, w = np.polynomial.legendre.leggauss(N)
        self.y = Q*t; self.w = Q*w
        A = np.eye(N) + M.K(self.y[:, None]-self.y[None, :])*self.w[None, :]
        self.r = np.linalg.solve(A, M.a1(self.y))
        self.m = float(self.w @ self.r)
        self.E = float(self.w @ (M.e(self.y)*self.r))
    def rf(self, x):
        return self.M.a1(x) - (self.M.K(x[:, None]-self.y[None, :])*self.w[None, :]) @ self.r

def ebulk(M, nu):
    if abs(nu-0.5) < 1e-15:
        x = np.linspace(-60/M.b, 60/M.b, 200001)
        return M.D/4 + np.trapezoid(M.e(x)*M.r0(x), x)
    Qhi = 1.0
    while Restr(M, Qhi).m < nu: Qhi *= 2
    Q = brentq(lambda Q: Restr(M, Q).m-nu, 1e-9, Qhi, xtol=1e-14)
    return M.D/4 + Restr(M, Q).E

def analyse(Delta, L, n):
    g = np.arccos(Delta); M = Rep(g); c = consts(g); nu = n/L
    x, mh = bethe_rep(L, n, g)
    F = energy_rep(L, n, g, x)
    sig = lambda y: M.a1(y) - M.K(y[:, None]-x[None, :]).sum(1)/L
    zL = lambda y: M.p(y)/(2*np.pi) - M.th(y[:, None]-x[None, :]).sum(1)/(2*np.pi*L)
    xn = x[-1]; o = {}
    o['Lr0xn_b'] = L*M.r0(xn)/M.b; o['need_I2'] = 1/(2*(np.pi-g))
    if n >= 2:
        ys = np.linspace(x[0], xn, 4001); o['min_sig_r0'] = float((sig(ys)/M.r0(ys)).min())
    o['theta1'] = c['theta1']
    far = np.linspace(xn, xn+80/M.b, 80001); sf = sig(far)
    o['Lneg_side'] = float(L*np.trapezoid(np.clip(-sf, 0, None), far)); o['need_I4'] = (np.pi-2*g)/(2*np.pi)
    if 2*n == L:
        Q = np.inf
    else:
        zf = zL(far); idx = np.where(zf >= nu/2)[0]
        Q = brentq(lambda t: zL(np.array([t]))[0]-nu/2, far[idx[0]-1], far[idx[0]], xtol=1e-15)
    o['Q'] = Q
    if np.isfinite(Q):
        R = Restr(M, Q); m = R.m
        grid = np.linspace(-Q, Q, 60001); rg = R.rf(grid)
    else:
        m = 0.5; grid = np.linspace(-xn-80/M.b, xn+80/M.b, 60001); rg = M.r0(grid)
    Cmu = np.searchsorted(x, grid, side='right')/L - nu/2
    cumr = np.r_[0, np.cumsum((rg[1:]+rg[:-1])/2*np.diff(grid))]
    if not np.isfinite(Q): cumr = np.arctan(np.exp(M.b*grid))/np.pi
    G = Cmu + nu/2 - cumr
    CumR = Cmu - zL(grid)
    o['L_CumR'] = float(L*np.abs(CumR).max()); o['s'] = c['s']
    o['L_mism'] = L*abs(nu-m); o['need_mism'] = c['qJ']*c['s']
    o['L_G'] = float(L*np.abs(G).max()); o['g'] = c['g']
    o['L_sig_r'] = float(L*np.abs(sig(grid)-rg).max()); o['D'] = c['D']
    o['L_F_Le'] = float(L*abs(F-L*ebulk(M, nu))); o['K'] = c['K']
    o['minhess_over_L'] = mh/L
    o['ok'] = bool(o['Lr0xn_b'] >= o['need_I2'] and o.get('min_sig_r0', 9) >= o['theta1']
                   and o['Lneg_side'] <= o['need_I4']+1e-9 and o['L_CumR'] <= max(o['s'], 0.5)+1e-6
                   and o['L_mism'] <= o['need_mism']+1e-6 and o['L_G'] <= o['g'] and o['L_sig_r'] <= o['D']
                   and o['L_F_Le'] <= o['K'])
    return o

if __name__ == '__main__':
    rows = []
    for Delta in [0.5, 0.9, 0.99, 0.999]:
        for L in [6, 13, 32, 64, 129, 256]:
            for n in sorted(set([1, 2, L//5, L//3, L//2-1, L//2])):
                if not 1 <= n <= L//2: continue
                o = analyse(Delta, L, n); o.update(Delta=Delta, L=L, n=n); rows.append(o)
                print(json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in o.items()}), flush=True)
    keys = ['Lr0xn_b', 'min_sig_r0', 'Lneg_side', 'L_CumR', 'L_mism', 'L_G', 'L_sig_r', 'L_F_Le', 'minhess_over_L']
    summ = {k: dict(min=min(r[k] for r in rows if k in r), max=max(r[k] for r in rows if k in r)) for k in keys}
    summ['all_ok'] = all(r['ok'] for r in rows); summ['n_cases'] = len(rows)
    summ['max_LFLe_over_K'] = max(r['L_F_Le']/r['K'] for r in rows)
    json.dump(dict(rows=rows, summary=summ), open('check_rep.json', 'w'), indent=1)
    print(json.dumps(summ, indent=1))
