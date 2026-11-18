"""C200: the continuum quantity TV(Dt_y) for r_Q (general filling), L-independent.
Dt_y(m) = Theta(min(w_+(m), d)) - Theta(min(w_-(m), d')), matched masses of r_Q on [-Q, Q], d = Q - y, d' = Q + y.
Reports sup over y (with d >= A eta) of TV(Dt_y) * min(1, d)/eta, and whether Dt_y restricted to the paired range is unimodal.
Floating point (FFT solver, cubic interpolation of the cumulative function)."""
import numpy as np, sys
from scipy.interpolate import CubicSpline, PchipInterpolator
from fftsolver import Kern, Restricted2
def run(eps, Q, hfac=12.0, ny=160):
    eta = np.arccos(1-eps); K = Kern(eta)
    R = Restricted2(K, Q, hfac=hfac, withZ=False)
    x = R.x; r = R.r
    rs = CubicSpline(x, r); Ns = rs.antiderivative()
    Nx = Ns(x); Ninv = PchipInterpolator(Nx, x)
    Th = lambda w: np.arctan(np.tanh(w)/np.tan(eta))/np.pi
    out = []
    ys = np.concatenate([Q-np.geomspace(0.3*eta, Q, ny//2), -Q+np.geomspace(0.3*eta, Q, ny//2)])
    for y in ys:
        d = Q-y; dp = Q+y; Ny = float(Ns(y)); mR = Nx[-1]-Ny; mL = Ny-Nx[0]
        w = np.geomspace(1e-3*eta, 2*Q, 3000)
        mm = np.unique(np.concatenate([Ny-Ns(np.clip(y-w, -Q, Q)), Ns(np.clip(y+w, -Q, Q))-Ny, [mR, mL]]))
        mm = mm[mm > 0]
        def inv(target):
            v = Ninv(target)
            for _ in range(6):
                v = np.clip(v-(Ns(v)-target)/np.maximum(rs(v), 1e-300), -Q, Q)
            return v
        wp = np.where(mm < mR, inv(np.minimum(Ny+mm, Nx[-1]))-y, d)
        wm = np.where(mm < mL, y-inv(np.maximum(Ny-mm, Nx[0])), dp)
        D = np.concatenate([[0.0], Th(wp)-Th(wm)])
        tv = np.abs(np.diff(D)).sum()
        mp_ = min(mR, mL); sel = np.concatenate([[True], mm <= mp_])
        Dp = D[sel]; tvp = np.abs(np.diff(Dp)).sum(); supp = np.abs(Dp).max()
        out.append((y, d if y > 0 else dp, tv, tvp, supp))
    return eta, np.array(out)
if __name__ == '__main__':
    for eps in [1e-2, 1e-3, 1e-4]:
        for Q in [0.3, 1.0, 2.5, 5.0]:
            eta, o = run(eps, Q)
            dd = o[:, 1]; Rr = np.minimum(1, dd)
            for A in [1, 5, 20]:
                s = dd >= A*eta
                if not s.any(): continue
                # paired-range TV scaled by min(1,d)/eta ; ratio TVpaired/(2 sup) (1 = unimodal)
                print(f'eps={eps:g} Q={Q} A={A}: max TVpaired*R/eta = {np.max(o[s,3]*Rr[s]/eta):.3f}  max TVpaired/(2 sup) = {np.max(o[s,3]/(2*o[s,4]+1e-300)):.3f}  max TVfull*R/eta = {np.max(o[s,2]*Rr[s]/eta):.3f}', flush=True)
