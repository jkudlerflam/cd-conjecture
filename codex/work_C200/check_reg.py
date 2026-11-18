"""C200: numerical look at the regularity estimate (F3): sup |r_Q'| rho/r_Q and sup |r_Q''| rho^2/r_Q, rho = min(1, eta + Q - |y|).
Also the identity P_{pi/2-eta} * g = a1 with g = 1/(2 eta cosh(pi x/eta)).  Floating point (FFT solver, spline derivatives)."""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from fftsolver import Kern, Restricted2
for eps in [1e-2, 1e-3, 1e-4]:
    eta = np.arccos(1-eps); K = Kern(eta)
    # identity check at a few points
    err = 0
    for x in [0.0, 0.3*eta, 2*eta, 0.5, 2.0]:
        Pk = lambda v: np.sin(2*eta)/(np.pi*(np.cosh(2*(x-v))+np.cos(2*eta)))
        g = lambda v: 1/(2*eta*np.cosh(np.pi*v/eta))
        val = quad(lambda v: Pk(v)*g(v), -40*eta, 40*eta, points=[0, x] if abs(x) < 40*eta else [0], limit=400, epsabs=1e-15, epsrel=1e-12)[0]
        err = max(err, abs(val-K.a1(x))/K.a1(x))
    for Q in [0.3, 1.0, 3.0]:
        R = Restricted2(K, Q, hfac=14.0, withZ=False)
        x = R.x; r = R.r; cs = CubicSpline(x, r)
        sel = slice(8, len(x)-8)
        rho = np.minimum(1, eta+Q-np.abs(x))
        d1 = np.abs(cs(x, 1))*rho/r; d2 = np.abs(cs(x, 2))*rho**2/r
        print(f'eps={eps:g} Q={Q}: rel err of P*g=a1 {err:.1e};  sup |r\'|rho/r = {d1[sel].max():.3f}  sup |r\'\'|rho^2/r = {d2[sel].max():.3f}  r(Q)/r(0) = {r[-1]/r[len(r)//2]:.4f}', flush=True)
