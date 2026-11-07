"""C138 diagnostics of the C136 finite-boundary identity near Delta=-1.

For Bethe roots at (L,n), nu=n/L<1/2, compute the counting boundary Q_L
(z_L(Q_L)=nu/2), the restricted continuum at Q=Q_L (Richardson in the mesh),
and the quantities entering C136 (4.2)-(4.4), (5.4):

  dN    = L(nu-m(Q_L))                      [C136 bound: B-1]
  curv  = e_bulk''(m) / eta^2               [C136 bound: a/eta^2 ~ pi/eta]
  TI    = L^2 [e_bulk(nu)-E-h(nu-m)] ~ L^2 e'' (nu-m)^2/2   (mismatch cost x L)
  TII   = L [F_L - L(E+h(nu-m))]            (quadrature term x L)
  dev   = TII - TI  ~ L(F_L - L e_bulk)
  Gcum  = L sup|mu(-inf,x] - int_{-Q}^x r_Q|  [C136 bound 3B/2]
  reldens(M) = max |sigma_L/r_Q - 1| over roots with L eta r_Q(x_j) >= M
Usage: python3 diag_identity.py K   -> diag_<K>.json
"""
import json
import sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from bethe_ferro import solve, Model
from continuum import Kernel, Restricted, grading

K = int(sys.argv[1])
NBASE = int(sys.argv[2]) if len(sys.argv) > 2 else 700
eps = 10.0**(-K)
eta = 2*np.arcsin(np.sqrt(eps/2))
mod, ker = Model(eta), Kernel(eta)
B = np.pi/(2*eta)
cache = json.loads(Path(f'energies_{K}.json').read_text())


def counting_boundary(L, x):
    f = lambda Q: mod.p(Q)-np.sum(mod.phi(Q-x))/L
    lo = x[-1]+1e-12
    hi = lo+1.0
    while f(hi) <= 0:
        hi *= 2
    return brentq(f, lo, hi, xtol=1e-15, rtol=1e-15)


def sigma_L(L, x, pts):
    return ker.a1(pts)+ker.T(pts[:, None]-x[None, :]).sum(axis=1)/L


def rich(Q):
    c = grading(ker, Q, NBASE)
    R1 = Restricted(ker, Q, N=NBASE, c=c)
    R2 = Restricted(ker, Q, N=2*NBASE, c=c)
    ext = lambda a, b: (4*b-a)/3
    out = dict(m=ext(R1.m, R2.m), E=ext(R1.E, R2.E), h=ext(R1.h, R2.h),
               curv=R2.curvature(), Z0=R2.Z[0], ZQ=R2.Z[-1], rQ=R2.r[-1],
               mesh_change_m=abs(R2.m-R1.m), mesh_change_E=abs(R2.E-R1.E))
    return out, R2


rows = []
Ls = [32, 64, 128, 256, 512, 1024, 2048]
for L in Ls:
    ns = sorted(set([int(round(L*f)) for f in [1/16, 1/8, 1/4, 3/8, 7/16]]
                    + [L//2-1, int(np.sqrt(L))]))
    for n in ns:
        if not 2 <= n < L/2:
            continue
        G, x, d = solve(L, n, eta)
        nu = n/L
        QL = counting_boundary(L, x)
        c, R = rich(QL)
        m, E, h = c['m'], c['E'], c['h']
        F = mod.Delta*L/4+G
        TII = L*(F-L*(E+h*(nu-m)))
        curv = c['curv']
        TI = L*L*curv*(nu-m)**2/2
        rQ = R.interp(R.r, x)
        sig = sigma_L(L, x, x)
        rel = np.abs(sig/rQ-1)
        dens = L*eta*rQ
        reld = {str(M): (float(np.max(rel[dens >= M])) if np.any(dens >= M) else None)
                for M in [1, 4, 16, 64]}
        # cumulative error at roots, one-sided limits
        yy = R.y
        Cr = np.concatenate([[0], np.cumsum((R.r[1:]+R.r[:-1])/2*np.diff(yy))])
        Ix = np.sign(x)*np.interp(np.abs(x), yy, Cr)
        Cnu = Ix   # centered: int_0^x r_Q (continuum mass m, not nu)
        j = np.arange(1, n+1)
        Cmu_minus = (j-1)/L-nu/2
        Cmu_plus = j/L-nu/2
        Gcum = L*max(np.max(np.abs(Cmu_minus-Cnu)), np.max(np.abs(Cmu_plus-Cnu)))
        Ln = L*(nu-m)
        rows.append(dict(L=L, n=n, nu=nu, etaL=eta*L, QL=QL, xmax=float(x[-1]),
                         dN=Ln, dN_over_B=Ln/(B-1), curv_over_eta2=curv/eta**2,
                         TI=TI, TII=TII, dev=TII-TI, dev_over_eta=(TII-TI)/eta,
                         TII_over_eta=TII/eta, Gcum=Gcum, Gcum_over_B=Gcum/B,
                         reldens=reld, ZQ=c['ZQ'], Z0=c['Z0'],
                         ZQ2_eta=c['ZQ']**2*eta,
                         mesh_change_E_times_L2=c['mesh_change_E']*L*L,
                         mesh_change_m_times_L=c['mesh_change_m']*L))
        print(json.dumps(rows[-1]), flush=True)
    Path(f'diag_{K}.json').write_text(json.dumps(dict(K=K, eta=eta, B=B, rows=rows), indent=1))
