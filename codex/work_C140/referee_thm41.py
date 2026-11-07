"""Referee checks of C140 Theorem 4.1 and its lemmas. Floating point diagnostics, J=1.
Run from codex/work_C140 with OPENBLAS_NUM_THREADS=1. Output referee_thm41.json."""
import json
import numpy as np
from scipy.integrate import quad
from bethe_scan import solve
from continuum_rep import Model, Restricted

out = dict(kernel=[], cutoff={}, cases=[], continuation=[])
ln2 = np.log(2)

# --- 63/64 cutoff, exact arithmetic in floating point
g_c = np.pi - np.pi**3/(8*(2-ln2))
g64 = np.arccos(63/64)
theta = lambda g: 1-8*(2-ln2)*(np.pi-g)/np.pi**3
out['cutoff'] = dict(gamma_threshold=g_c, Delta_threshold=np.cos(g_c), gamma_63_64=g64,
                     theta_at_63_64=theta(g64))
print(out['cutoff'])

# --- kernel constants J(0), sum_{m>=2} J^{*m}(0)
for D in [0.0, 0.25, 0.5, 0.9, 0.98, 63/64, 0.999]:
    g = np.arccos(D)
    def Kh(w):
        return np.exp(-g*w)*(-np.expm1(-(np.pi-2*g)*w))/(-np.expm1(-np.pi*w)) if w > 0 else (np.pi-2*g)/np.pi
    Jh = lambda w: Kh(w)/(1+Kh(w))
    J0 = quad(Jh, 0, np.inf, limit=500)[0]/np.pi
    S2 = quad(lambda w: Jh(w)**2/(1-Jh(w)), 0, np.inf, limit=500)[0]/np.pi
    row = dict(Delta=D, J0=J0, J0_bound=ln2/(np.pi*g), S2=S2, S2_bound=(1-ln2)/(np.pi*g),
               qJ=(np.pi-2*g)/(2*(np.pi-g)))
    out['kernel'].append(row); print(row, flush=True)


def analyse(L, n, D, x):
    M = Model(D); g, b = M.g, M.b
    sig = lambda y: M.a1(y)-M.K(np.asarray(y)[:, None]-x[None, :]).sum(1)/L
    xn = x[-1]
    # sigma on [x_1,x_n]
    ys = np.linspace(0, xn, 6001)
    s = sig(ys)
    lb_proved = M.r0(ys)-(2-ln2)/(np.pi*g*L)
    # gap masses via z_L
    def z(y):
        th = 2*np.arctan(np.tanh(np.asarray(y)[:, None]-x[None, :])/np.tan(g))
        return 2*np.arctan(np.tanh(y)/np.tan(g/2))/(2*np.pi)-th.sum(1)/(2*np.pi*L)
    I = np.arange(1, n+1)-(n+1)/2
    zres = float(np.max(np.abs(z(x)-I/L)))
    # tail beyond x_n
    yt = np.linspace(xn, xn+60*g+20, 200001); st = sig(yt); dy = yt[1]-yt[0]
    negmass = 2*float(np.sum(np.clip(-st, 0, None))*dy)
    tail = float(np.trapezoid(st, yt))
    c = (np.pi-g)/(2*np.pi)-(n/L)*(np.pi-2*g)/(2*np.pi)
    tail_pred = c-n/(2*L)+1/(2*L)
    # Hessian
    dp = lambda y: 2*np.tan(g/2)*(1-np.tanh(y)**2)/(np.tan(g/2)**2+np.tanh(y)**2)
    dt = lambda y: 2*np.tan(g)*(1-np.tanh(y)**2)/(np.tan(g)**2+np.tanh(y)**2)
    W = dt(x[:, None]-x[None, :])
    H = np.diag(L*dp(x)-W.sum(1))+W
    hmin = float(np.linalg.eigvalsh(H).min())
    # upper comparison sigma <= s_Q beyond Q=x_n (Corollary 3.2) and Lemma 3.3 chain
    R = Restricted(M, xn, 8)
    ytt = np.linspace(xn, xn+12*g+6, 801)
    sQ = R.s_full(ytt)
    up_viol = float(np.max(sig(ytt)-sQ))
    hQ_half = float(np.trapezoid(R.s_full(np.linspace(xn, xn+80*g+30, 40001)), np.linspace(xn, xn+80*g+30, 40001)))
    return dict(L=L, n=n, Delta=D,
                min_sig_over_r0_core=float((s/M.r0(ys)).min()),
                min_sig_minus_lb_proved=float((s-lb_proved).min()),
                theta=theta(g),
                Lr0xn=float(L*M.r0(xn)), Lr0xn_bound=np.pi**2/(8*g*(np.pi-g)),
                gap_residual=zres, negmass_out=negmass, negmass_bound=1/L,
                tail_mass=tail, tail_pred=tail_pred,
                hess_min=hmin, sigma_min_roots=float(sig(x).min()),
                upper_comp_violation=up_viol, half_hQ=hQ_half, half_hQ_needed=1/(2*L))


for D in [0.25, 0.5, 0.9, 0.98, 63/64]:
    for L in list(range(2, 13))+[16, 17, 32, 33, 64, 65, 128, 129, 256, 257]:
        for n in sorted(set([2, 3, L//4, L//3, L//2-1, L//2])):
            if n < 2 or n > L//2:
                continue
            E, d, x = solve(L, n, D, return_roots=True)
            row = analyse(L, n, D, x)
            out['cases'].append(row)
    rows = [r for r in out['cases'] if r['Delta'] == D]
    print(D, 'min sig/r0 core', min(r['min_sig_over_r0_core'] for r in rows),
          'min margin vs proved lb', min(r['min_sig_minus_lb_proved'] for r in rows),
          'min Lr0xn/bound', min(r['Lr0xn']/r['Lr0xn_bound'] for r in rows),
          'max negmass*L', max(r['negmass_out']*r['L'] for r in rows),
          'max |tail-pred|', max(abs(r['tail_mass']-r['tail_pred']) for r in rows),
          'max gap res', max(r['gap_residual'] for r in rows),
          'min hess', min(r['hess_min'] for r in rows),
          'max upper viol', max(r['upper_comp_violation'] for r in rows),
          'min 2L*int_Q s_Q', min(2*r['L']*r['half_hQ'] for r in rows), flush=True)

# --- continuation in gamma from the free point, compared with the direct solver
for L, n in [(7, 3), (8, 4), (12, 6), (13, 6), (33, 16), (64, 32)]:
    gs = np.linspace(np.pi/2, g64, 400)
    I = np.arange(1, n+1)-(n+1)/2
    x = None; worst = 0.; minsig = np.inf
    for g in gs:
        D = np.cos(g)
        if x is None:
            x = np.arctanh(np.tan(np.pi*I/L))  # gamma=pi/2: p=2atan(tanh x)=2pi I/L
        tg, t2 = np.tan(g/2), np.tan(g)
        p = lambda y: 2*np.arctan(np.tanh(y)/tg)
        th = lambda y: 2*np.arctan(np.tanh(y)/t2) if abs(np.cos(g)) > 1e-15 else 0*y
        dp = lambda y: 2*tg*(1-np.tanh(y)**2)/(tg**2+np.tanh(y)**2)
        dt = lambda y: 2*t2*(1-np.tanh(y)**2)/(t2**2+np.tanh(y)**2) if abs(np.cos(g)) > 1e-15 else 0*y
        for it in range(60):
            r = p(x)-th(x[:, None]-x[None, :]).sum(1)/L-2*np.pi*I/L
            if np.max(np.abs(r)) < 1e-14:
                break
            W = dt(x[:, None]-x[None, :])/L
            Hm = np.diag(dp(x)-W.sum(1))+W
            x = x-np.linalg.solve(Hm, r)
        M = Model(max(D, 0.0))
        sr = M.a1(x)-M.K(x[:, None]-x[None, :]).sum(1)/L
        minsig = min(minsig, float((L*sr).min()))
    E, dd, xd = solve(L, n, np.cos(g64), return_roots=True)
    row = dict(L=L, n=n, max_root_diff_to_direct=float(np.max(np.abs(x-xd))), min_L_sigma_roots_along=minsig)
    out['continuation'].append(row); print(row, flush=True)

json.dump(out, open('referee_thm41.json', 'w'), indent=1)
