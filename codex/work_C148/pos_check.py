"""Sufficient positivity inequality of PROOF Theorem 4.3 on the window.
Pi(lam,Q) = [K(0) + (K-J)(d_Q)] / Loc(Q),  d_Q = min(Q, pi-2Q),
Loc(Q) = max(4 r0(Q), (1-omega_Q) r0(Q) / (2 int_Q^{pi/2} r0)),  omega_Q = int_{|u|<pi/2-Q} J.
Also the direct bound margin [tanh(lam/2) - coth(lam)/2]/pi. Exact Fourier/periodized evaluations."""
import numpy as np, json
from circle import Circ
def fourier(lam, coef, x, kmax=None):
    kmax = kmax or int(80/lam)+20
    k = np.arange(1, kmax+1)
    return (coef(0)+2*(coef(k)[None,:]*np.cos(2*np.outer(x,k))).sum(1))/np.pi
def sup_Pi(lam, nQ=3000):
    M = Circ(np.cosh(lam)); t = lambda k: np.exp(-2*np.abs(k)*lam)
    K0 = 1/(np.tanh(lam)*np.pi)
    Q = np.linspace(1e-4, np.pi/2-1e-4, nQ)
    d = np.minimum(Q, np.pi-2*Q)
    KJ = fourier(lam, lambda k: t(k)**2/(1+t(k)), d)
    # omega_Q = int_{-(pi/2-Q)}^{pi/2-Q} J = (1/pi)[(1/2)(pi-2Q) + 2 sum Jhat_k sin(2k(pi/2-Q))/k]
    kmax = int(80/lam)+20; k = np.arange(1, kmax+1); Jk = t(k)/(1+t(k))
    h = np.pi/2-Q
    om = (h+ (2*Jk[None,:]*np.sin(2*np.outer(h,k))/k[None,:]).sum(1))/np.pi
    r0Q = M.r0(Q)
    tail = M.Z0(np.full_like(Q, np.pi/2))-M.Z0(Q)
    Loc = np.maximum(4*r0Q, (1-om)*r0Q/(2*tail))
    P = (K0+KJ)/Loc
    i = int(np.argmax(P))
    return float(P[i]), float(Q[i]), float(om.max())
out = {}
for lam in np.r_[np.linspace(0.02, 0.5, 13), np.linspace(0.5, 1.5, 51), [2, 3, 5]]:
    P, Qs, om = sup_Pi(lam)
    direct = (np.tanh(lam/2)-0.5/np.tanh(lam))/np.pi
    out[f'{lam:.4f}'] = dict(sup_Pi=P, at_Q=Qs, max_omega=om, direct_margin=float(direct))
mx = max(v['sup_Pi'] for k, v in out.items() if float(k) <= 1.5)
print('max sup_Pi on lam<=1.5:', mx)
for k in list(out)[::6]: print(k, out[k])
json.dump(out, open('pos_check.json','w'), indent=1)
