"""C200: sign of S(y) at nodes of the right half near the Fermi edge (kappa = 1/2 trial, nodes t(y_j) = j, y_n = Q).
Prints max over nodes with 0 < y < Q, delta = Q - y < 1 of S/eta, and of [S - B(delta)/2]/eta, with lambda_Q.  Floating point."""
import numpy as np, sys
from scipy.interpolate import CubicHermiteSpline
from bethe_ferro import Model
from fftsolver import Kern, solve_Q2
from gen_trial import counting_function, Fmap
def case(eps, L, n, X=25.0):
    eta = np.arccos(1-eps); m = Model(eta); K = Kern(eta); nu = n/L
    Th = lambda w: np.arctan(np.tanh(w)/np.tan(eta))/np.pi; q = 1-2*eta/np.pi
    I = np.arange(1, n+1)-(n+1)/2
    R = solve_Q2(K, nu, hfac=10.0, withZ=False); Q = R.Q
    u, Z, s = counting_function(K, m, R, X); sp = CubicHermiteSpline(u, Z, s)
    tg = I/L+0.5/L
    i = np.clip(np.searchsorted(Z, tg), 1, len(u)-1)
    y = u[i-1]+(tg-Z[i-1])/np.maximum(s[i-1], 1e-300)
    for _ in range(40):
        y = np.clip(y, u[0], u[-1]); y = y-(sp(y)-tg)/sp(y, 1)
    y[-1] = Q
    F = Fmap(m, y, L, 0)/(2*np.pi)          # = 1/2 - L T*G(y_j), L T*G = q_in/2 + S
    qin = Th(Q-y)+Th(Q+y)
    S = 0.5-F-0.5*qin
    d = Q-y; sel = (y > 0) & (d < 1) & (d > 0)
    Bd = q/2-Th(d)
    if not sel.any(): return eta*L*R.r[-1], float("nan"), float("nan"), float("nan"), 0
    return eta*L*R.r[-1], (S[sel]/eta).max(), ((S-Bd/2)[sel]/eta).max(), ((S-Bd/2)[sel]/Bd[sel]).max(), int(sel.sum())
if __name__ == '__main__':
    for eps in [1e-2, 1e-3, 1e-4]:
        for L in [64, 256, 1024]:
            for n in [L//8, L//4, 3*L//8, L//2-L//16, L//2-2]:
                lamQ, a, b, c, k = case(eps, L, n)
                print(f'eps={eps:g} L={L} n={n} lamQ={lamQ:.3f} nodes={k}: max S/eta={a:.4f}  max (S-B/2)/eta={b:.4f}  max (S-B/2)/B={c:.4f}', flush=True)
