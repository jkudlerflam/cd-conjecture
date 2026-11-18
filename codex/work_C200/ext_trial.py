"""C200: the trial of PROOF Theorem 4.11 (Astra call D): interior nodes t(y_j) = j + k (j <= n - k), M exterior nodes at
Q + i A d, d = max(eta, ell), ell = width of the last continuum cell, the other k - M nodes at +infinity.
Evaluates the Bethe map exactly and prints min Ftilde/(2 pi eta) over interior nodes and over exterior nodes.  Floating point."""
import numpy as np, sys
from scipy.interpolate import CubicHermiteSpline
from bethe_ferro import Model
from fftsolver import Kern, solve_Q2
from gen_trial import counting_function, Fmap
def case(eps, L, n, k, M, A, X=25.0):
    eta = np.arccos(1-eps); m = Model(eta); K = Kern(eta); nu = n/L
    R = solve_Q2(K, nu, hfac=10.0, withZ=False); Q = R.Q
    u, Z, s = counting_function(K, m, R, X); sp = CubicHermiteSpline(u, Z, s)
    def quant(tlev):
        tg = tlev/L-nu/2
        i = np.clip(np.searchsorted(Z, tg), 1, len(u)-1)
        y = u[i-1]+(tg-Z[i-1])/np.maximum(s[i-1], 1e-300)
        for _ in range(40):
            y = np.clip(y, u[0], u[-1]); y = y-(sp(y)-tg)/sp(y, 1)
        return y
    yin = quant(np.arange(1, n-k+1)+k); yin[-1] = Q
    ell = Q-float(quant(np.array([n-1.0]))[0]); d = max(eta, ell)
    ok5 = M*A*d <= min(1, Q)/4
    yext = Q+A*d*np.arange(1, M+1) if (ok5 or FORCE) else np.zeros(0)
    y = np.concatenate([yin, yext]); ninf = n-len(y)
    F = Fmap(m, y, L, ninf)/(2*np.pi*eta)
    Fi = F[:len(yin)]; Fe = F[len(yin):]
    return eta*L*R.r[-1], ell/eta, ok5, Fi.min(), int(np.argmin(Fi))+1, len(yin), (Fe.min() if len(Fe) else float('nan'))
FORCE = False
if __name__ == '__main__':
    k, M, A = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]); FORCE = len(sys.argv) > 4
    for eps in [1e-3, 1e-4, 1e-5]:
        for L in [64, 256, 1024]:
            for n in [L//8, L//4, 3*L//8, L//2-L//16, L//2-2]:
                if n <= k+2: continue
                lamQ, le, ok5, fi, ji, nin, fe = case(eps, L, n, k, M, A)
                print(f'eps={eps:g} L={L} n={n} lamQ={lamQ:.3f} ell/eta={le:.1f} case(5)={ok5}: interior min={fi:.3f}@{ji}/{nin}  exterior min={fe:.3f}', flush=True)
