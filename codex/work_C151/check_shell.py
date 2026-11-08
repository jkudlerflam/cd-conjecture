"""C151 shell and Harnack checks on Bethe roots.
H1  Harnack: Psi(x) = 2 pi Im theta(x + i eta/2) >= 2 pi h(1/2) eta L sigma(x), h(1/2) = 2 log 1.5 - 1/2.
H2  margins of the shell hypotheses on S = {M0 <= lambda < 4 C_D(eta)}, lambda = eta L r_Q:
    min lambda_sigma/lambda (needed c, (G_c) asks 3/4), min lambda_sigma/log(lambda), min Psi.
Floating point only."""
import json, sys, time
import numpy as np
from bethe_ferro import solve
from continuum import Kernel, Restricted, grading
from ddv import DdV
H = 2*np.log(1.5)-0.5
def run(eps, L, n, M0=2.0, N=600):
    eta = 2*np.arcsin(np.sqrt(eps/2)); t0 = time.time()
    G, x, _ = solve(L, n, eta)
    D = DdV(eta, L, x); K = Kernel(eta)
    g = np.pi-eta; b = np.pi/g
    if 2*n == L:
        Q = np.inf; Xm = x[-1]+2.0
        rfun = lambda X: 1/(2*g*np.cosh(b*X))
    else:
        Q = D.counting_boundary(); Xm = Q
        R = Restricted(K, Q, N=N, c=grading(K, Q, N))
        rfun = lambda X: R.interp(R.r, X)
    smax = rfun(0.0)
    hgrid = min(eta/12, 1/(12*L*max(smax, 1e-12)))
    X = np.arange(0, Xm, hgrid)
    lam = eta*L*rfun(X)
    lsig = eta*L*D.sigma(X).real
    Psi = D.Psi(X, eta/2)
    CD = (4/np.pi)*(np.log(1/eta)+np.pi/2)
    M1 = 4*CD
    S = (lam >= M0) & (lam < M1)
    K0 = lam >= M0
    o = dict(eps=eps, eta=eta, L=L, n=n, Q=float(Q), lam0=float(lam[0]), M1=M1,
             harnack_min=float(np.min(Psi/(2*np.pi*H*lsig))), shell_pts=int(S.sum()))
    if S.any():
        o.update(c_min_shell=float(np.min(lsig[S]/lam[S])), logratio_min_shell=float(np.min(lsig[S]/np.log(lam[S]))),
                 Psi_min_shell=float(np.min(Psi[S])), Psi_over_pilam_min=float(np.min(Psi[S]/(np.pi*lam[S]))))
    if K0.any():
        o.update(c_min_K0=float(np.min(lsig[K0]/lam[K0])), Psi_min_K0=float(np.min(Psi[K0])))
    o['dev_max_all'] = float(np.max(np.abs(lsig-lam)))
    o['sec'] = round(time.time()-t0, 1)
    return o
if __name__ == '__main__':
    tag = sys.argv[1]; eps = float(sys.argv[2]); Ls = [int(v) for v in sys.argv[3].split(',')]
    res = []
    for L in Ls:
        for n in sorted(set([L//2, L//2-2, 3*L//8, L//4, L//8])):
            o = run(eps, L, n); print(o, flush=True); res.append(o)
    json.dump(res, open(f'check_shell_{tag}.json', 'w'), indent=1)
