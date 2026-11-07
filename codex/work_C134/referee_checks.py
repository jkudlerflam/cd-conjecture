"""Referee checks for PROOF.md Sections 2-8 (Delta in (-1,0]). Floating point only.

Run from this folder:  OPENBLAS_NUM_THREADS=1 python3 referee_checks.py
"""
import itertools
import json
import sys
import numpy as np
from scipy.integrate import quad
from bethe_scan import solve, functions
from half_scan import bulk

out = {}

# ---------- A. kernel sign and L1 norm, independent quadrature ----------
kern = []
for Delta in [-0.999, -0.99, -0.9, -0.5, -0.1, 0.0]:
    g = np.arccos(Delta)
    K = lambda x: np.sin(2*g)/(2*np.pi*(np.sinh(x)**2+np.sin(g)**2))
    l1 = 2*quad(lambda x: abs(K(x)), 0, 60, limit=400)[0]
    xs = np.linspace(-30, 30, 20001)
    kern.append(dict(Delta=Delta, L1=l1, q=2*g/np.pi-1, maxK=float(np.max(K(xs)))))
out['kernel'] = kern

# ---------- B. independent ED: min over flux of H (+XY) equals Bethe at K_L periodic ----------
def ed_min_flux(L, n, Delta, nphi=181):
    states = [sum(1 << x for x in c) for c in itertools.combinations(range(L), n)]
    idx = {s: i for i, s in enumerate(states)}
    best = np.inf
    for phi in np.linspace(0, 2*np.pi, nphi):
        H = np.zeros((len(states), len(states)), complex)
        for i, s in enumerate(states):
            for x in range(L):
                y = (x+1) % L
                bx, by = (s >> x) & 1, (s >> y) & 1
                H[i, i] += Delta*(0.25 if bx == by else -0.25)
                if bx != by:
                    t = s ^ (1 << x) ^ (1 << y)
                    # flux phi on bond (L-1,0); hop amplitude +1/2
                    ph = np.exp(1j*phi*(1 if bx else -1)) if x == L-1 else 1.0
                    H[idx[t], i] += 0.5*ph
        best = min(best, np.linalg.eigvalsh(H)[0])
    return best

edrows = []
for Delta in [-0.9, -0.5]:
    for L in [6, 7, 8, 9]:
        n = L//2
        E_ed = ed_min_flux(L, n, Delta)
        E_ba = float(solve(L, n, Delta)[0])
        edrows.append(dict(Delta=Delta, L=L, ed_minflux=E_ed, bethe=E_ba, diff=E_ed-E_ba))
out['ed_flux'] = edrows

# ---------- C. constants ----------
def consts(Delta):
    g = np.arccos(Delta); s = np.sin(g); q = 2*g/np.pi-1; B = 1/(1-q)
    V1 = 2*abs(1/np.tan(g))/np.pi; V2 = q*(8/s**4+4/s**2)
    D = 1+B*V1; D1 = B*V2; b = np.pi/g; a = np.pi*s
    Kl = a*(1.5*b+D1/(2*D)+32*D**2/b+8*D*B)
    return dict(g=g, s=s, q=q, B=B, D=D, D1=D1, b=b, a=a, Klarge=Kl,
                Ksmall=(4*g*D+1)**2, K=max((4*g*D+1)**2, Kl))

def Kdelta(d):
    al = np.arccos(1-d); s = np.sqrt(d*(2-d)); Bb = np.pi/(2*al)
    Db = 1+1/(al*s); D1b = Bb*(8/s**4+4/s**2)
    return (4*np.pi*Db+1)**2+np.pi*(3+D1b/2+32*Db**2+8*Db*Bb)

worst_ratio, worst_dom = 0., 0.
for d in np.geomspace(1e-4, 1, 400):
    Kd = Kdelta(d)
    worst_ratio = max(worst_ratio, Kd*d**3/1200)
    for Delta in np.linspace(-1+d, 0, 60):
        worst_dom = max(worst_dom, consts(Delta)['K']/Kd)
out['Kdelta_over_1200_per_delta3_max'] = worst_ratio
out['K_over_Kdelta_max'] = worst_dom
out['point_K'] = {str(D): consts(D)['K'] for D in [-0.99, -0.9, -0.5, 0.0]}

# ---------- D. Bethe intermediate inequalities ----------
def analyse(L, Delta, e0, E, x):
    c = consts(Delta); g, b, a, D, D1, B = c['g'], c['b'], c['a'], c['D'], c['D1'], c['B']
    p, th, dp, dt, bare, _ = functions(Delta)
    n = len(x); nu = n/L
    rho = lambda y: 1/(2*g*np.cosh(b*y))
    eps = lambda y: -a*rho(y)
    sig = lambda y: (dp(y)-dt(np.subtract.outer(np.atleast_1d(y), x)).sum(axis=1)/L)/(2*np.pi)
    z = lambda y: (p(y)-th(np.subtract.outer(np.atleast_1d(y), x)).sum(axis=1)/L)/(2*np.pi)
    # second derivative of theta, numerically from closed form of K'
    def dsig(y):
        y = np.atleast_1d(y)
        def kp(u):  # derivative of dt(u)
            return -np.sin(2*g)*np.sinh(2*u)/(np.sinh(u)**2+np.sin(g)**2)**2
        def ap(u):
            return -np.sin(g)*np.sinh(2*u)/(np.sinh(u)**2+np.sin(g/2)**2)**2
        return (ap(y)-kp(np.subtract.outer(y, x)).sum(axis=1)/L)/(2*np.pi)
    R = {}
    R['casimir_lower'] = float(L*(L*e0-E))  # must be <= K
    xm = max(10., abs(x[-1])+3)
    xs = np.unique(np.r_[np.linspace(-xm, xm, 6001), x, x-1e-9, x+1e-9])
    cum_rho = .25+np.arctan(np.sinh(b*xs))/(2*np.pi)
    cnt = np.searchsorted(x, xs, side='right')/L
    R['LG_over_B'] = float(L*np.max(np.abs(cnt-cum_rho))/B)
    rh = rho(xs)
    R['Ldsig_over_D'] = float(L*np.max(np.abs(sig(xs)-rh))/D)
    drh = -b*np.tanh(b*xs)*rh
    R['Lddsig_over_D1'] = float(L*np.max(np.abs(dsig(xs)-drh))/D1) if D1 > 0 else 0.
    if L >= 4*g*D:
        A = np.arccosh(L/(4*g*D))/b
        epsA = eps(A)
        gfun = lambda y: np.where(np.abs(y) <= A, eps(y), epsA)
        # int_{-A}^{A} eps sigma
        core = quad(lambda y: float(eps(y)*sig(y)[0]), -A, A, limit=800, epsabs=1e-13, epsrel=1e-12)[0] if A > 0 else 0.
        zA = float(z(A)[0]); cc = (p(50.)-th(50.-x).sum()/L)/(2*np.pi)
        tails_sigma_mass = 2*cc-(2*zA)
        int_g_sigma = core+epsA*tails_sigma_mass
        clipped = float(np.sum(gfun(x))-L*int_g_sigma)
        # midpoint error on central interval only (odd L: add back strips, where f=epsA)
        w = 2*cc-nu
        mid = clipped+L*epsA*w if L % 2 else clipped
        # tail term L int (eps-g) d(mu - sigma)
        tr = quad(lambda y: float((eps(y)-epsA)*sig(y)[0]), A, 60, limit=800, epsabs=1e-15)[0]
        tail = float(np.sum((eps(x)-gfun(x)))-2*L*tr)
        R['A'] = float(A)
        R['L_mid_over_bound'] = float(L*abs(mid)/(a*(1.5*b+D1/(2*D))))
        R['L_tail_over_bound'] = float(L*abs(tail)/(a*(32*D**2/b+8*D*B)))
        R['identity_residual'] = float((E-L*e0)-(clipped+tail))
        # sampled variation of f' in x on [-A,A] including two jumps
        gr = np.linspace(-A, A, 20001)
        fp = a*b*np.tanh(b*gr)*rho(gr)/sig(gr)
        var = abs(fp[0])+np.sum(np.abs(np.diff(fp)))+abs(fp[-1])
        R['var_over_bound'] = float(var/(12*a*b+4*a*D1/D))
        R['min_sigma_over_rho_core'] = float(np.min(sig(gr)/rho(gr)))
    return R

rows = []
Ls = [int(a) for a in sys.argv[1:]] or [8, 9, 16, 17, 64, 65, 128, 129, 256, 257, 512, 513]
for Delta in [-0.5, -0.9, -0.99]:
    e0 = bulk(Delta)
    K = consts(Delta)['K']
    cache = {}
    def get(L):
        if L not in cache:
            cache[L] = solve(L, L//2, Delta, return_roots=True)
        return cache[L]
    for L in Ls:
        E, diag, x = get(L)
        E = float(E)
        R = analyse(L, Delta, e0, E, x)
        E2 = float(get(2*L)[0])
        R.update(Delta=Delta, L=L, K=K, Q_J=float(L*(E2-2*E)), Q_bound_J=2.5*K,
                 upper_casimir=float(L*(E-L*e0)), residual=diag['residual'])
        rows.append(R)
        print(json.dumps(R), flush=True)
out['bethe'] = rows
json.dump(out, open('referee_checks.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != 'bethe'}, indent=1))
