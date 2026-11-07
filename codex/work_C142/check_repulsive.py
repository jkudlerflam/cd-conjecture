"""C142 Task 1 checks on the repulsive side 0<Delta<1. Floating point, no enclosures. J=1.

For each (Delta, L, n) it tests the new inequalities of PROOF.md Part I:
 (1) L r0(x_n)/b >= 1/(2(pi-gamma))                         [Lemma I.2]
 (2) min_{[x_1,x_n]} sigma_L/r0 >= theta_1 = 1-2(2-ln2)(pi-gamma)/pi^2   [Thm I.3]
 (3) per-side negative mass of sigma_L beyond x_n <= (pi-2gamma)/(2 pi L) [Lemma I.4]
 (4) L sup|Cum R| on [-Q_L,Q_L] <= s = 1-gamma/pi            [Lemma I.5]
 (5) L|nu-m(Q_L)| <= qJ s                                     [Lemma I.6]
 (6) L||G|| <= g(gamma)                                        [Lemma I.6]
 (7) L sup|sigma_L-r_Q| <= 2K(0) L||G|| (true G) and <= D=2g/(pi gamma)  [Lemma I.7]
 (8) L|F_L(n)-L e_bulk(n/L)| <= K_large(gamma)                [Theorem I.9]
"""
import json, sys
import numpy as np
from bethe_scan import solve, functions
from continuum_rep import Model, Restricted, e_bulk
from constants_rep import consts

def analyse(delta, L, n):
    M = Model(delta); g = M.g; b = M.b; qJ = M.qJ
    E, d, x = solve(L, n, delta, return_roots=True)
    x = np.atleast_1d(x)
    p, th, dp, dt, bare, seed = functions(delta)
    nu = n/L
    zL = lambda y: p(y)/(2*np.pi) - th(np.asarray(y)[:, None]-x[None, :]).sum(1)/(2*np.pi*L)
    sig = lambda y: M.a1(y) - M.K(np.asarray(y)[:, None]-x[None, :]).sum(1)/L
    xn = x[-1]
    out = {}
    out['Lr0xn_over_b'] = L*M.r0(xn)/b
    out['need1'] = 1/(2*(np.pi-g))
    # (2)
    if n >= 2:
        ys = np.linspace(0, xn, 3001)
        out['min_sig_over_r0'] = float((sig(ys)/M.r0(ys)).min())
    out['theta1'] = 1-2*(2-np.log(2))*(np.pi-g)/np.pi**2
    # beyond x_n: fine grid in units of gamma out to xn + 40 gamma/pi*... use b-scale
    far = xn + 60/b
    ys = np.linspace(xn, far, 60001)
    s_ = sig(ys)
    neg = np.trapezoid(np.clip(-s_, 0, None), ys)
    out['L_negmass_side'] = L*neg
    out['need3'] = (np.pi-2*g)/(2*np.pi)
    # Q_L: first y>=x_n with z_L(y)=nu/2
    zvals = zL(np.array([xn]))[0] + np.r_[0, np.cumsum((s_[1:]+s_[:-1])/2*np.diff(ys))]
    idx = np.where(zvals >= nu/2)[0]
    if len(idx) > 0:
        QL = float(ys[idx[0]])
    elif 2*n == L:
        QL = np.inf
    else:
        raise RuntimeError('Q_L not found')
    out['QL'] = QL
    # cumulative R on [-QL, QL]; by symmetry use [0, QL]
    if np.isfinite(QL):
        R = Restricted(M, QL)
        m = R.m
        grid = np.linspace(-QL, QL, 40001)
        rg = R.s_full(grid)
        Cr = np.r_[0, np.cumsum((rg[1:]+rg[:-1])/2*np.diff(grid))]
        rprime = np.gradient(rg, grid)
    else:
        m = 0.5
        grid = np.linspace(-far, far, 40001)
        rg = M.r0(grid)
        Cr = np.arctan(np.exp(b*grid))/np.pi
        rprime = -b*np.tanh(b*grid)*rg
    Cmu = np.searchsorted(x, grid, side='right')/L
    G = Cmu - Cr
    zg = zL(grid)
    CumR = (Cmu - nu/2) - zg
    s = 1-g/np.pi
    out['L_supCumR'] = float(L*np.abs(CumR).max()); out['s'] = s
    out['L_mismatch'] = L*abs(nu-m); out['need5'] = qJ*s
    c = consts(g)
    out['L_G'] = float(L*np.abs(G).max()); out['g_bound'] = c['g']
    sg = sig(grid)
    out['L_sig_minus_r'] = float(L*np.abs(sg-rg).max())
    out['twoK0_LG'] = 2*M.K(0)*out['L_G']; out['D_bound'] = c['D']
    sgp = np.gradient(sg, grid)
    out['L_sigp_minus_rp'] = float(L*np.abs(sgp-rprime)[100:-100].max()); out['D1_bound'] = c['D1']
    eb = e_bulk(M, nu)
    out['L_F_minus_Le'] = float(L*abs(float(E)-L*eb)); out['K_large'] = c['K_large']
    ok = (out['Lr0xn_over_b'] >= out['need1'] and out.get('min_sig_over_r0', 1) >= out['theta1']
          and out['L_negmass_side'] <= out['need3']+1e-9 and out['L_supCumR'] <= max(s, .5)+1e-6
          and out['L_mismatch'] <= out['need5']+1e-6 and out['L_G'] <= c['g']
          and out['L_sig_minus_r'] <= out['D_bound'] and out['L_F_minus_Le'] <= c['K_large'])
    out['ok'] = bool(ok)
    return out

if __name__ == '__main__':
    rows = []
    deltas = [0.5, 0.9, 0.99, 0.999]
    Ls = [4, 8, 17, 32, 65, 128, 256]
    for delta in deltas:
        for L in Ls:
            for n in sorted(set([1, 2, L//8, L//4, (3*L)//8, L//2-1, L//2])):
                if not 1 <= n <= L//2: continue
                try:
                    o = analyse(delta, L, n)
                except Exception as exc:
                    o = dict(error=str(exc))
                o.update(delta=delta, L=L, n=n)
                rows.append(o)
                print(json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in o.items()}), flush=True)
    json.dump(rows, open('check_repulsive.json', 'w'), indent=0)
    keys = ['Lr0xn_over_b', 'min_sig_over_r0', 'L_negmass_side', 'L_supCumR', 'L_mismatch', 'L_G', 'L_sig_minus_r', 'L_sigp_minus_rp', 'L_F_minus_Le']
    summ = {}
    for k in keys:
        v = [r[k] for r in rows if k in r]
        summ[k] = dict(min=min(v), max=max(v))
    summ['all_ok'] = all(r.get('ok', False) for r in rows)
    summ['failures'] = [(r['delta'], r['L'], r['n']) for r in rows if not r.get('ok', False)]
    json.dump(summ, open('check_repulsive_summary.json', 'w'), indent=1)
    print(json.dumps(summ, indent=1))
