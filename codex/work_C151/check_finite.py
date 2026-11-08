"""C151 check at a finite Fermi boundary: the DdV contour identity for the energy, the Fermi-edge
(CFT) pieces, and the C136 identity. f = eps_ext, an analytic extension of the restricted dressed
energy built from a node representation of eps_Q (exactly analytic, real on R)."""
import json, sys, time
import numpy as np
from numpy.polynomial.legendre import leggauss
from bethe_ferro import solve
from continuum import Kernel, Restricted, grading, bulk_at
from ddv import DdV, path_integral, real_integral, gauss_nodes

def run(eta, L, n, N=500):
    t0 = time.time()
    G, x, _ = solve(L, n, eta)
    D = DdV(eta, L, x); K = Kernel(eta)
    Q = D.counting_boundary()
    R = Restricted(K, Q, N=N, c=grading(K, Q, N))
    # nodes for int_{-Q}^{Q} T(z-y) eps_Q(y) dy, panels <= eta/2
    yn, wn = gauss_nodes(-Q, Q, eta/2)
    yn = yn.real; wn = wn.real
    W = wn*R.interp(R.epsQ, yn)
    h = R.h
    s2, c2 = np.sin(eta)**2, np.cos(eta/2)**2
    def f(z):
        z = np.atleast_1d(z)
        out = -s2/(2*(np.sinh(z)**2+c2))-h
        for i in range(0, len(z), 2000):
            zz = z[i:i+2000]
            out[i:i+2000] = out[i:i+2000] + (D.T(zz[:, None]-yn[None, :])*W[None, :]).sum(axis=1)
        return out
    # direct side
    xr, wr = gauss_nodes(-Q, Q, eta/8)
    xr = xr.real; wr = wr.real
    fx = f(xr+0j).real
    direct = np.sum(f(x+0j).real) - L*np.sum(wr*fx*D.sigma(xr).real)
    # DdV: rectangle top at height ep, f ell' form (boundary terms real)
    out = dict(eta=eta, L=L, n=n, Q=Q, direct_L=L*direct)
    hstep = min(eta/6, 0.3/(L*R.r[0]+1))
    F_ = lambda z: f(z)*D.dell(z)
    for ep in [eta/4, eta/2]:
        val = np.imag(path_integral(F_, [Q+0j, Q+1j*ep, -Q+1j*ep, -Q+0j], hstep))/np.pi
        out[f'ddv_L_{ep/eta:.2f}'] = L*val
    # Fermi edge pieces with f' ell form: -(1/pi) Im int f' ell on the two vertical sides
    ep = eta/2
    hd = 1e-6
    fp = lambda z: (f(z+hd)-f(z-hd))/(2*hd)
    Fv = lambda z: fp(z)*D.ell(z)
    edge = -np.imag(path_integral(Fv, [Q+0j, Q+1j*ep], ep/40))/np.pi - np.imag(path_integral(Fv, [-Q+1j*ep, -Q+0j], ep/40))/np.pi
    sigQ = D.sigma(np.array([Q+0j]))[0].real
    vF = fp(np.array([Q+0j]))[0].real/(2*np.pi*sigQ)
    out.update(edge_L=L*edge, cft_L=-np.pi*vF/6, lamQ=eta*L*R.r[-1], lam0=eta*L*R.r[0])
    # C136 identity: F - L[E + h(nu - m)] = L int eps_Q dR_L
    nu = n/L
    F = K.Delta*L/4+G
    out['C136_L'] = L*(F - L*(R.E + R.h*(nu-R.m)))
    Rb = bulk_at(K, nu, Q0=Q, N=N)
    out['Fdev_L'] = L*(F - L*Rb.E)
    out['sec'] = time.time()-t0
    return out

if __name__ == '__main__':
    res = []
    for eta in [0.3, 0.1]:
        for L, n in [(64, 16), (64, 24), (128, 48), (256, 64)]:
            o = run(eta, L, n); print({k: (round(v, 9) if isinstance(v, float) else v) for k, v in o.items()}, flush=True); res.append(o)
    json.dump(res, open('check_finite.json', 'w'), indent=1)
