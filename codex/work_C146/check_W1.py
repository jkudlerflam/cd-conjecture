"""C146: check of the velocity identity and the variation bound for v_Q = eps_Q'/r_Q. Floating point only. J=1.
Identity (PROOF Lemma 3): v' = F1 - (r(Q)/r) [ v(Q) (G+H) - v (G-H) ],
  F1 = a [ (R a1')^2 - (R a1)(R a1'') ]/(R a1)^2, G+H = R J_+, G-H = R J_0,
  J_+(x) = T(x-Q)+T(x+Q), J_0(x) = T(x-Q)-T(x+Q).
Bound (PROOF Theorem 4): Var_[-Q,Q](v) <= 2[ v(Q) + 2 int_0^Q (F1)_- + 3 a tanh(2Q)(Z(Q)-1) * 2 ].
Usage: python3 check_W1.py eps [eps ...]"""
import json, sys
import numpy as np
from continuum import Kernel, mesh, grading, _far_or_near

def build_pm(K, Q, N, c):
    y = mesh(K, Q, N, c); y[0] = 0.0; N1 = len(y)
    Ae, Ao = np.eye(N1), np.eye(N1)
    yl, yr = y[:-1], y[1:]
    for i in range(N1):
        xi = np.full(N1-1, y[i])
        m0, m1 = _far_or_near(K, xi, yl, yr)
        r0, r1 = _far_or_near(K, -xi, yl, yr)
        for A, s in ((Ae, 1.0), (Ao, -1.0)):
            c0, c1 = m0+s*r0, m1+s*r1
            A[i, :-1] -= c0-c1
            A[i, 1:] -= c1
    return y, Ae, Ao

def a1d(K, x, k):
    h = 1e-4
    if k == 1:
        return (K.a1(x+h)-K.a1(x-h))/(2*h)
    return (K.a1(x+h)-2*K.a1(x)+K.a1(x-h))/h**2

def case(eps, Q, N=700):
    eta = 2*np.arcsin(np.sqrt(eps/2)); K = Kernel(eta); a = K.a
    y, Ae, Ao = build_pm(K, Q, N, grading(K, Q, N))
    r = np.linalg.solve(Ae, K.a1(y)); Z = np.linalg.solve(Ae, np.ones_like(y))
    R1 = np.linalg.solve(Ao, a1d(K, y, 1)); R2 = np.linalg.solve(Ae, a1d(K, y, 2))
    Jp = np.linalg.solve(Ae, K.T(y-Q)+K.T(y+Q)); J0 = np.linalg.solve(Ao, K.T(y-Q)-K.T(y+Q))
    v = -a*R1/r
    F1 = a*(R1**2-r*R2)/r**2
    rhs = F1-(r[-1]/r)*(v[-1]*Jp-v*J0)
    dv = np.gradient(v, y)
    # compare away from the two mesh ends
    sl = slice(2, -3)
    err = float(np.max(np.abs(dv[sl]-rhs[sl]))/max(np.max(np.abs(dv[sl])), 1e-300))
    trap = lambda f: float(np.sum((f[1:]+f[:-1])*np.diff(y))/2)
    F1m = trap(np.maximum(-F1, 0))
    Var = 2*float(np.sum(np.abs(np.diff(v))))
    vQ = float(v[-1])
    edge = trap((r[-1]/r)*np.abs(v[-1]*Jp-v*J0))
    bound = 2*(vQ+2*F1m+2*3*a*np.tanh(2*Q)*(Z[-1]-1))
    return dict(eps=eps, eta=eta, Q=Q, a=a, identity_relerr=err, Var_v=Var, Var_v_over_a=Var/a,
                vQ_over_a=vQ/a, intF1minus_over_a=F1m/a, edge_term_over_a=edge/a,
                a_ZQm1=float(a*(Z[-1]-1)), bound=bound, bound_over_Var=bound/max(Var, 1e-300),
                min_F1_over_a=float(F1.min()/a), int_GplusH_0Q=trap(Jp)/1.0, ZQm1=float(Z[-1]-1))

if __name__ == '__main__':
    rows = []
    for eps in [float(s) for s in sys.argv[1:]]:
        for Q in [0.01, 0.1, 0.5, 1.0, 2.0, 4.0]:
            c = case(eps, Q); rows.append(c); print(json.dumps(c), flush=True)
    json.dump(rows, open('check_W1_%s.json' % '_'.join(sys.argv[1:]), 'w'), indent=1)
