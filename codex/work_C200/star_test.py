"""C200: the interior inequality (*) at general filling.
X_j = [ int_{-Q}^{Q} T(y_j - v) s(t(v) - t_j) dv - (1/2) int_{|v|>Q} T(y_j - v) dv ] / eta   for trial nodes y_j in [-Q, Q],
obtained from the exactly evaluated Bethe map of the trial (kappa0 in 1/2 + N), the nodes beyond Q and the phantom levels.
(*) asks X_j <= (2/pi) kappa0.  Floating point only."""
import numpy as np, sys
from scipy.interpolate import CubicHermiteSpline
from bethe_ferro import Model
from fftsolver import Kern, solve_Q2
from gen_trial import counting_function, Fmap

def case(eps, L, n, k0, X=25.0):
    eta = np.arccos(1-eps); m = Model(eta); K = Kern(eta); nu = n/L
    Th = lambda w: np.arctan(np.tanh(w)/np.tan(eta))/np.pi
    q = 1-2*eta/np.pi
    I = np.arange(1, n+1)-(n+1)/2
    R = solve_Q2(K, nu, hfac=10.0, withZ=False); Q = R.Q
    u, Z, s = counting_function(K, m, R, X); sp = CubicHermiteSpline(u, Z, s)
    def quant(tg):
        i = np.clip(np.searchsorted(Z, tg), 1, len(u)-1)
        y = u[i-1]+(tg-Z[i-1])/np.maximum(s[i-1], 1e-300)
        for _ in range(40):
            y = np.clip(y, u[0], u[-1]); y = y-(sp(y)-tg)/sp(y, 1)
        return y
    lev = I/L+k0/L; fin = lev < Z[-1]
    y = quant(lev[fin]); ninf = n-fin.sum()
    F = Fmap(m, y, L, ninf)/(2*np.pi)          # = kappa0 - L T*G_+(y_j)
    K0 = int(round(k0-0.5))
    tph = (k0-0.5-np.arange(K0))               # phantom levels in t
    vph = quant(tph/L-nu/2)
    inside = y <= Q+1e-12
    yb = y[~inside]                             # finite nodes beyond Q
    Xs = []
    for j in np.where(inside)[0]:
        yj = y[j]
        qin = Th(Q-yj)+Th(Q+yj); qout = q-qin
        ph = np.sum(Th(yj+Q)-Th(yj-vph))
        bey = np.sum(Th(yb-yj)-Th(Q-yj))+ninf*(q/2-Th(Q-yj))
        ints = k0-F[j]-k0*qin+ph-bey
        Xs.append((ints-0.5*qout)/eta)
    Xs = np.array(Xs); jm = int(np.argmax(Xs))
    return Xs.max(), jm+1, int(inside.sum()), Q-y[jm], eta*L*R.r[-1], Xs[-3:]

if __name__ == '__main__':
    eps = float(sys.argv[1]); Ls = [int(a) for a in sys.argv[2:]]
    for L in Ls:
        for n in sorted(set([max(2, L//32), L//8, L//4, 3*L//8, L//2-L//16, L//2-8, L//2-2, L//2-1])):
            if n < 8: continue
            row = []
            for k0 in [1.5, 2.5, 5.5]:
                xm, jm, nin, d, lamQ, last = case(eps, L, n, k0)
                row.append(f'k0={k0}: maxX={xm:.3f}@{jm}/{nin} d={d:.3g} last3={np.round(last,3)}')
            print(f'eps={eps:g} L={L} n={n} lamQ={lamQ:.2f} | '+' | '.join(row), flush=True)
