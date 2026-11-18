"""Theorem 4.2 ratios, Theorem 4.3 TV(D) rho/eta, (I4) Var(v), monotonicity of r_Q, own continuum solver."""
import numpy as np
from c204lib import *
from scipy.integrate import cumulative_trapezoid
print('eta        Q     sup|r\'|rho/r  sup|r\'\'|rho^2/r  max r\' on (0,Q]  sup TV(D)rho/eta   Var(v)  24a+4aB')
for x in [1e-2, 1e-3, 1e-4, 1e-5]:
    eta = float(np.arccos(1-x)); K = Kern(eta); a = PI*np.sin(eta); Bc = PI/(2*eta)
    for Q in [0.02, 0.1, 0.3, 1.0, 3.0, 6.0]:
        c = Continuum(eta, Q, hfac=40)
        v, r, h = c.v, c.r, c.h
        r1 = np.gradient(r, h, edge_order=2); r2 = np.gradient(r1, h, edge_order=2)
        delta = Q - np.abs(v); rho = np.minimum(1, eta + delta)
        inner = delta > 3*h
        R1 = np.max(np.abs(r1[inner])*rho[inner]/r[inner]); R2 = np.max(np.abs(r2[inner])*rho[inner]**2/r[inner])
        mono = np.max(r1[(v > 0) & inner])
        # TV(D) at sample points y in [0,Q)
        C = cumulative_trapezoid(r, v, initial=0)
        worst = 0
        for y in np.concatenate([Q*np.linspace(0, 0.99, 25), Q - eta*np.array([0.5, 1, 2, 4, 8, 16])]):
            if y < 0 or y >= Q: continue
            Cy = np.interp(y, v, C)
            mmax = max(C[-1]-Cy, Cy)
            m = np.linspace(0, mmax, 40001)
            # w_+(m): v such that C(v) - Cy = m, capped
            wp = np.interp(Cy + m, C, v, right=np.inf) - y
            wp = np.where(Cy + m > C[-1], np.inf, wp)
            wm = y - np.interp(Cy - m, C, v, left=-np.inf)
            wm = np.where(Cy - m < 0, np.inf, wm)
            th = lambda w, cap: np.where(np.isinf(w), K.Theta(cap), K.Theta(np.where(np.isinf(w), 0, w)))
            D = th(wp, Q - y) - th(wm, Q + y)
            TV = np.abs(np.diff(D)).sum()
            rh = min(1, eta + Q - y)
            worst = max(worst, TV*rh/eta)
        # (I4)
        Ra1p = c.solve(np.gradient(K.a1(v), h) * 0 + (-2*np.sinh(2*v)*K.a1(v)/(np.cosh(2*v)+np.cos(eta))))
        u = -a*Ra1p; vv = u/r
        Var = np.abs(np.diff(vv)).sum()
        print('%.2e  %5.2f   %8.3f      %8.3f        %+.2e       %8.3f        %7.3f  %7.2f' % (eta, Q, R1, R2, mono, worst, Var, 24*a+4*a*Bc))
