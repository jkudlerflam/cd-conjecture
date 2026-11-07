"""C147 referee: attractive side. Own Nystrom and Bethe code. J=1, eta = pi - gamma, eps = 1+Delta.
(1) calR(0) by quadrature vs the Lemma II.3 bound, calR(0) >= T(0).
(2) Lemma II.4/II.5: Var_y calR_Q(x,.) <= 2(1+q)calR(0), Var_y d_x calR_Q(x,.) <= 4M1T + 4(2+q)T(0)calR(0).
(3) C145 Theorem B.1 (reflection) and Lemma B.2; spot values for (W1'), (W3').
(4) C142 Theorem II.6 on Bethe roots, up to the Fermi point (derivative of the Nystrom interpolant)."""
import json, sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from ref_core import Att, bethe_att

def Rhat2(w, eta):
    """That/(1-That) = sinh(A)/(2 cosh(A+c) sinh(c)), A=(pi-2eta)w/2, c=eta w/2, written without overflow."""
    if w < 1e-12:
        return (1-2*eta/np.pi)/(2*eta/np.pi)
    A = (np.pi-2*eta)*w/2; c = eta*w/2
    return np.exp(-c)*(1-np.exp(-2*A))/((1+np.exp(-2*A-2*c))*2*np.sinh(c)) if c < 300 else 0.0

def calR0(eta):
    f = lambda w: Rhat2(w, eta)
    edges = [0, 1e-3, 1, 1/eta, 10/eta, 80/eta]
    return sum(quad(f, a, b, limit=500, epsabs=1e-13, epsrel=1e-12)[0] for a, b in zip(edges[:-1], edges[1:]))/np.pi

class Nys:
    def __init__(self, M, Q, N):
        t, w = np.polynomial.legendre.leggauss(N)
        self.M, self.Q, self.y, self.w = M, Q, Q*t, Q*w
        Tm = M.T(self.y[:, None]-self.y[None, :])
        self.Tm = Tm
        self.A = np.eye(N) - Tm*self.w[None, :]
        self.r = np.linalg.solve(self.A, M.a1(self.y))
        self.Z = np.linalg.solve(self.A, np.ones(N))
    def interp(self, f, x, src):  # u(x) = src(x) + sum_k w_k T(x-y_k) u_k
        return src(x) + (self.M.T(x[:, None]-self.y[None, :])*self.w[None, :]) @ f
    def dinterp(self, f, x, dsrc):
        return dsrc(x) + (self.M.Tp(x[:, None]-self.y[None, :])*self.w[None, :]) @ f

def var(f):
    return float(np.abs(np.diff(f)).sum())

def kernel_checks(eps, Qs, N):
    eta = np.arccos(1-eps); M = Att(eta); q = M.q
    R0 = calR0(eta); T0 = 1/(np.pi*np.tan(eta)); s = np.sin(eta)
    M1T = (2*np.cos(eta)/np.pi)*(9/(16*np.sqrt(3))/s**2+0.25/s)
    o = dict(eps=eps, eta=eta, calR0=R0, calR0_bound=(np.log(1/eta)+np.pi/2)/(np.pi*eta), T0=T0,
             bound_II4=2*(1+q)*R0, bound_II5=4*M1T+4*(2+q)*T0*R0, cases=[])
    for Q in Qs:
        S = Nys(M, Q, N)
        U = np.linalg.solve(S.A, S.Tm)  # calR_Q at nodes
        yf = np.linspace(-Q, Q, int(40*2*Q/eta)+2001)
        TF = M.T(yf[:, None]-S.y[None, :])*S.w[None, :]  # fine x nodes
        worst4 = worst5 = 0.
        for x in [0.0, 0.5*Q, Q-min(eta, Q/4), Q]:
            # u(y) = calR_Q(x,y) = calR_Q(y,x): source T(y-x)
            ux = np.linalg.solve(S.A, M.T(S.y-x))
            uf = M.T(yf-x) + TF @ ux
            worst4 = max(worst4, var(uf))
            c = S.w*M.Tp(x-S.y)
            vf = M.Tp(x-yf) + M.T(yf[:, None]-S.y[None, :]) @ c + TF @ (U @ c)
            worst5 = max(worst5, var(vf))
        # Theorem B.1 and Lemma B.2 at nodes
        r0 = M.r0(S.y); h = r0 - S.r
        refl = M.r0(2*Q-S.y)+M.r0(2*Q+S.y)
        d = Q-np.abs(S.y); kap = 1-4*np.exp(-2*M.b*d)
        epsp = np.linalg.solve(S.A, -M.a*(-2*np.sinh(2*S.y)*M.a1(S.y)/(np.cosh(2*S.y)+np.cos(eta))))
        # edge values by interpolation
        Zq = S.interp(S.Z, np.array([Q]), lambda z: np.ones_like(z))[0]
        rq = S.interp(S.r, np.array([Q]), M.a1)[0]
        da1 = lambda z: -2*np.sinh(2*z)*M.a1(z)/(np.cosh(2*z)+np.cos(eta))
        epq = S.interp(epsp, np.array([Q]), lambda z: -M.a*da1(z))[0]
        xs = np.linspace(0, Q, 400); vq = S.interp(epsp, xs, lambda z: -M.a*da1(z))/S.interp(S.r, xs, M.a1)
        o['cases'].append(dict(Q=Q, var_calRQ=worst4, var_dcalRQ=worst5,
                               B1_min_slack=float((refl-h).min()/r0.max()), h_min=float(h.min()/r0.max()),
                               B1_ratio_ok=bool(np.all(S.r/r0 >= kap-1e-9)),
                               B2_max=float(np.max(np.abs(epsp)/(2*M.a*np.tanh(2*Q)*S.r))),
                               W1p_monotone=bool(np.all(np.diff(vq) >= -1e-9*np.abs(vq).max())),
                               W3p=float(eta*Zq**2/np.tanh(2*Q)), ebulk2_over_eta2=float(epq/(2*rq*Zq**2)/eta**2)))
    return o

def bethe_case(eps, L, n, N):
    eta = np.arccos(1-eps); M = Att(eta)
    x = bethe_att(L, n, eta)
    sig = lambda z: M.a1(z) + M.T(z[:, None]-x[None, :]).sum(1)/L
    dsig = lambda z: -2*np.sinh(2*z)*M.a1(z)/(np.cosh(2*z)+np.cos(eta)) + M.Tp(z[:, None]-x[None, :]).sum(1)/L
    da1 = lambda z: -2*np.sinh(2*z)*M.a1(z)/(np.cosh(2*z)+np.cos(eta))
    if 2*n == L:
        X = np.linspace(0, x[-1]+10, 20001); X = np.unique(np.r_[X, x[x >= 0]])
        r = M.r0(X); dr = -M.b*np.tanh(M.b*X)*r; Q = np.inf
    else:
        f = lambda Q: M.p(Q)-M.phi(Q-x).sum()/L
        Q = brentq(f, x[-1]+1e-13, x[-1]+80, xtol=1e-15)
        S = Nys(M, Q, N)
        X = np.unique(np.r_[np.linspace(0, Q, 20001), x[x >= 0]])
        r = S.interp(S.r, X, M.a1); dr = S.dinterp(S.r, X, da1)
    return dict(eps=eps, L=L, n=n, etaL=eta*L, Q=Q,
                Dmax=float(eta*L*np.abs(sig(X)-r).max()), D1max=float(eta**2*L*np.abs(dsig(X)-dr).max()))

if __name__ == '__main__':
    out = dict(kernel=[], bethe=[])
    for eps, Qs, N in [(0.5, [0.3, 1, 3], 400), (0.1, [0.3, 1, 3], 600), (1e-2, [0.3, 1, 3], 1000), (1e-3, [0.3, 1, 2], 1600)]:
        o = kernel_checks(eps, Qs, N); out['kernel'].append(o); print(json.dumps(o), flush=True)
    for eps, N in [(0.5, 400), (0.1, 600), (1e-2, 1000), (1e-3, 1600)]:
        eta = np.arccos(1-eps)
        q = 1-2*eta/np.pi; R0 = calR0(eta); T0 = 1/(np.pi*np.tan(eta)); s = np.sin(eta)
        M1T = (2*np.cos(eta)/np.pi)*(9/(16*np.sqrt(3))/s**2+0.25/s)
        CD = 2*(1+q)*eta*R0; CD1 = eta**2*(4*M1T+4*(2+q)*T0*R0)
        for L in [16, 64, 256]:
            for n in sorted(set([2, L//8, L//4, 3*L//8, L//2-2, L//2])):
                c = bethe_case(eps, L, n, N); c.update(CD=CD, CD1=CD1, ok=bool(c['Dmax'] <= CD and c['D1max'] <= CD1))
                out['bethe'].append(c); print(json.dumps(c), flush=True)
    out['summary'] = dict(all_bethe_ok=all(c['ok'] for c in out['bethe']),
                          max_Dmax=max(c['Dmax'] for c in out['bethe']), max_D1max=max(c['D1max'] for c in out['bethe']))
    json.dump(out, open('check_att.json', 'w'), indent=1)
    print(json.dumps(out['summary']))
