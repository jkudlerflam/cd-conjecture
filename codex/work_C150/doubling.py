"""C150: direct comparison of the (2L, 2n) and (L, n) Bethe root systems.

J = 1, eta = pi - gamma, Delta = -cos(eta). Conventions of C136/C138 (bethe_ferro.py).

For a trial configuration yt of the 2L system (2n points) the residual in
counting units is rho_k = z_2L[yt](yt_k) - t'_k, t'_k = I'_k/(2L), and the
proved homotopy bound (PROOF Theorem 2) is |E_2L(y*) - E(yt)| <= 2a sum|rho_k|.

Two trials are built from the L roots x_j:
  refined: xh_{2j-1}, xh_{2j} = z_L^{-1}(t_j -+ 1/(4L))         (every cell)
  mixed:   refined on a central block of cells, doubled roots (x_j, x_j) outside.
D = F_2L(2n) - 2F_L(n) = A + B with A = E(yt) - 2E_L(x) and B = E_2L(y*) - E(yt).
Floating point only.
"""
import json
import sys
import time
import numpy as np
from bethe_ferro import Model, solve
import constants_C150 as CC

_KC = {}


def KC(eta):
    if eta not in _KC:
        _KC[eta] = CC.consts(eta, N=20000)
    return _KC[eta]


def Tker(m, u):
    return np.sin(2*m.eta)/(2*np.pi*(np.sinh(u)**2+np.sin(m.eta)**2))


def a1(m, x):
    return np.sin(m.eta)/(2*np.pi*(np.sinh(x)**2+np.cos(m.eta/2)**2))


def theta(m, u):
    return -2*np.arctan(np.tanh(u)/m.te)


def zfun(m, x, roots, w):
    """counting function p/(2pi) - w sum_k theta(x - roots_k)/(2pi), w = 1/L."""
    x = np.atleast_1d(x)
    out = np.empty(len(x))
    for i0 in range(0, len(x), 512):
        xx = x[i0:i0+512]
        out[i0:i0+512] = (m.p(xx)-w*theta(m, xx[:, None]-roots[None, :]).sum(1))/(2*np.pi)
    return out


def sigma(m, x, roots, w):
    x = np.atleast_1d(x)
    out = np.empty(len(x))
    for i0 in range(0, len(x), 512):
        xx = x[i0:i0+512]
        out[i0:i0+512] = a1(m, xx)+w*Tker(m, xx[:, None]-roots[None, :]).sum(1)
    return out


def invert_z(m, roots, L, targets, lo, hi, it=80):
    """solve z_L(x) = target on brackets [lo, hi] (vectorized bisection + Newton polish)."""
    lo = lo.copy(); hi = hi.copy()
    for _ in range(it):
        mid = (lo+hi)/2
        f = zfun(m, mid, roots, 1.0/L)-targets
        lo = np.where(f < 0, mid, lo); hi = np.where(f < 0, hi, mid)
        if np.max(hi-lo) < 1e-15*max(1.0, np.max(np.abs(hi))):
            break
    x = (lo+hi)/2
    for _ in range(3):
        f = zfun(m, x, roots, 1.0/L)-targets
        x = np.clip(x-f/sigma(m, x, roots, 1.0/L), lo, hi)
    return x


def refined_grid(m, x, L):
    n = len(x)
    I = np.arange(1, n+1)-(n+1)/2
    t = I/L
    delta = np.arctanh(min(np.tan(m.eta), 1-1e-15))
    eps = 1e-14
    xin = invert_z(m, x, L, t-0.25/L, x-delta-eps, x)
    xout = invert_z(m, x, L, t+0.25/L, x, x+delta+eps)
    # symmetrize (exact symmetry x -> -x)
    xin_s = (xin-xout[::-1])/2
    xout_s = (xout-xin[::-1])/2
    return xin_s, xout_s


def residual_2L(m, yt, L):
    """rho_k = z_2L[yt](yt_k) - t'_k."""
    N = len(yt)
    Ip = np.arange(1, N+1)-(N+1)/2
    return zfun(m, yt, yt, 1.0/(2*L))-Ip/(2*L)


def jac_2L(m, y, L):
    d = y[:, None]-y[None, :]
    Tm = Tker(m, d)/(2*L)
    np.fill_diagonal(Tm, 0.0)
    J = -Tm
    J[np.diag_indices_from(J)] = a1(m, y)+Tm.sum(1)
    return J


def eprime(m, y):
    # e = -a a1, a1 = sin(eta)/(2pi (sinh^2 + cos^2(eta/2)))
    s2, c2 = np.sin(m.eta)**2, np.cos(m.eta/2)**2
    return s2*np.sinh(y)*np.cosh(y)/(np.sinh(y)**2+c2)**2


def r_half(eta, x):
    g = np.pi-eta
    return 1/(2*g*np.cosh(np.pi*x/g))


def case(eta, L, n, cont=True, Nmesh=500):
    m = Model(eta)
    a = np.pi*np.sin(eta)
    t0 = time.time()
    GL, x, _ = solve(L, n, eta)
    G2, y, _ = solve(2*L, 2*n, eta)
    D = G2-2*GL
    out = dict(eta=eta, eps=float(m.eps), L=L, n=n, etaL=eta*L, LD=L*D,
               LD_logas=L*m.eps*2*n*(L-n)/((2*L-1)*(L-1)) if L > 1 else 0.0)
    if 2*n == L:
        g = np.pi-eta
        v = (np.pi/2)*np.sin(g)/g
        out['cft_piv_over4'] = np.pi*v/4
    # discrete local density at roots and cell edges
    xin, xout = refined_grid(m, x, L)
    lam_root = eta*L*sigma(m, x, x, 1.0/L)
    I = np.arange(1, n+1)-(n+1)/2
    tl, tr = I/L-0.5/L, I/L+0.5/L
    delta = np.arctanh(min(np.tan(eta), 1-1e-15))
    # cell edges: z_L^{-1}(t_j -+ 1/(2L)); brackets from neighbouring roots
    xl = np.concatenate([[x[0]-50.0], x[:-1]])
    xr = np.concatenate([x[1:], [x[-1]+50.0]])
    el = invert_z(m, x, L, tl, xl, x)
    er = invert_z(m, x, L, tr, x, xr)
    pts = np.stack([el, xin, x, xout, er], 1)
    lam_pts = eta*L*sigma(m, pts.ravel(), x, 1.0/L).reshape(pts.shape)
    lam_cell_min = lam_pts.min(1)
    out['pair_halfwidth_over_delta'] = float(max(np.max(xout-x), np.max(x-xin))/delta)
    # refined trial
    yr = np.empty(2*n); yr[0::2] = xin; yr[1::2] = xout
    rho = residual_2L(m, yr, L)
    A_ref = float(np.sum(m.e(xin)+m.e(xout)-2*m.e(x)))
    B_ref = D-A_ref
    J = jac_2L(m, yr, L)
    w = np.linalg.solve(J, eprime(m, yr))
    out.update(LA_ref=L*A_ref, LB_ref=L*B_ref, L_bound_ref=L*2*a*np.sum(np.abs(rho)),
               LB_lin_ref=-L*float(w@rho), max_w_over_2a=float(np.max(np.abs(w))/(2*a)),
               L_sumrho_ref=L*float(np.sum(np.abs(rho))))
    # the true 2L roots: check w bound at y* too
    w2 = np.linalg.solve(jac_2L(m, y, L), eprime(m, y))
    out['max_w_over_2a_at_ystar'] = float(np.max(np.abs(w2))/(2*a))
    # discrete hypothesis quantities
    for Lam in (1.0, 2.0, 4.0):
        out['etaNsp_root_%g' % Lam] = eta*int(np.sum(lam_root < Lam))
        out['etaNsp_cell_%g' % Lam] = eta*int(np.sum(lam_cell_min < Lam))
    # mixed trials with discrete central block: cells with min lambda >= Lam_d, central run
    for Lam_d in (1.0, 2.0):
        dense = lam_cell_min >= Lam_d
        # central block: maximal run of dense cells containing the middle
        mid = n//2
        blk = np.zeros(n, bool)
        if n > 0 and dense[mid] and dense[n-1-mid]:
            i = mid
            while i >= 0 and dense[i] and dense[n-1-i]:
                i -= 1
            blk[i+1:n-1-i] = True
        res = mixed(m, x, xin, xout, el, er, blk, L, D, a)
        out.update({k+'_d%g' % Lam_d: v for k, v in res.items()})
    # continuum block [-X2, X2] with lambda_Q(X2) = 2 (theorem's choice)
    if cont:
        try:
            out.update(continuum_part(m, x, xin, xout, el, er, lam_pts, L, n, D, a, Nmesh))
        except Exception as exc:
            out['cont_error'] = str(exc)
    out['time'] = time.time()-t0
    return out


def mixed(m, x, xin, xout, el, er, blk, L, D, a):
    n = len(x)
    yt = np.empty(2*n)
    yt[0::2] = np.where(blk, xin, x)
    yt[1::2] = np.where(blk, xout, x)
    rho = residual_2L(m, yt, L)
    A = float(np.sum(np.where(blk, m.e(xin)+m.e(xout)-2*m.e(x), 0.0)))
    B = D-A
    Nsp = int(np.sum(~blk))
    eta = m.eta
    # I1 = int over block images of dy/lambda_L
    if blk.any():
        ys = np.linspace(el[blk][0], er[blk][-1], 20001)
        lam = eta*L*sigma(m, ys, x, 1.0/L)
        I1 = float(np.trapezoid(1/lam, ys))
    else:
        I1 = 0.0
    K = KC(eta)
    thm = a*(Nsp+(K['H']+1)/4)+(a/eta)*(9/512)*(2*K['C3']*K['H']+K['C3e'])*I1
    # pieces of sum|rho|: sparse members, and |zeta~| everywhere
    sp = np.repeat(~blk, 2)
    shift = np.zeros(2*n); shift[0::2] = 0.25/L; shift[1::2] = -0.25/L
    zt = rho-np.where(sp, shift, 0.0)
    return dict(LA_mix=L*A, LB_mix=L*B, L_bound_mix=L*2*a*np.sum(np.abs(rho)), Nsp_mix=Nsp,
                etaNsp_mix=eta*Nsp, I1_mix=I1, LD_check=L*(A+B), L_thm_bound_actual=thm,
                L_2a_zeta_dense=L*2*a*float(np.sum(np.abs(zt))))


def continuum_part(m, x, xin, xout, el, er, lam_pts, L, n, D, a, Nmesh):
    eta = m.eta
    if 2*n == L:
        QL = np.inf
        rfun = lambda X: r_half(eta, X)
    else:
        from scipy.optimize import brentq
        from continuum import Kernel, Restricted, grading
        f = lambda Q: zfun(m, np.array([Q]), x, 1.0/L)[0]-n/(2*L)
        QL = brentq(f, x[-1]+1e-14, x[-1]+60, xtol=1e-15)
        ker = Kernel(eta)
        R = Restricted(ker, QL, N=Nmesh, c=grading(ker, QL, Nmesh))
        rfun = lambda X: R.interp(R.r, X)
    lamQ_root = eta*L*rfun(x)
    lamQ_pts = eta*L*rfun(np.stack([el, xin, x, xout, er], 1))
    CD = (4/np.pi)*(np.log(1/eta)+np.pi/2)
    shell = (lamQ_pts >= 2) & (lamQ_pts < 2*CD)
    deep = lamQ_pts >= 2*CD
    out = dict(QL=float(QL), CD=CD)
    out['min_lamL_over_lamQ_shell'] = float(np.min(lam_pts[shell]/lamQ_pts[shell])) if shell.any() else None
    out['min_lamL_over_lamQ_deep'] = float(np.min(lam_pts[deep]/lamQ_pts[deep])) if deep.any() else None
    out['n_shell_pts'] = int(shell.sum())
    blk = np.all(lamQ_pts >= 2, axis=1)
    res = mixed(m, x, xin, xout, el, er, blk, L, D, a)
    out.update({k+'_Q': v for k, v in res.items()})
    out['etaN_tail_Q'] = eta*int(np.sum(lamQ_root < 2))
    return out


if __name__ == '__main__':
    tag = sys.argv[1]
    EPS = [float(s) for s in sys.argv[2].split(',')]
    LS = [int(s) for s in sys.argv[3].split(',')]
    rows = []
    for eps in EPS:
        eta = 2*np.arcsin(np.sqrt(eps/2))
        for L in LS:
            fills = sorted(set([max(1, L//64), L//16, L//8, L//4, 3*L//8, L//2-2, L//2-1, L//2]))
            for n in fills:
                if n < 1:
                    continue
                try:
                    r = case(eta, L, n)
                except Exception as exc:
                    r = dict(eta=eta, L=L, n=n, error=str(exc))
                rows.append(r)
                print(json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()}), flush=True)
    json.dump(rows, open('doubling_%s.json' % tag, 'w'), indent=1)
