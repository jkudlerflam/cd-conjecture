"""Central-switch cell with a Peierls phase phiT on carrier-T hops: Fermi points
and Fermi-sea centres of the order-lambda^2 logical chain at mu = 0."""
import sys, json, os
import numpy as np, scipy.linalg as la
from onebody2 import logical_h_precise
from onebody import blocks, max_c

def fermi_points(T, M, nk=200000):
    s = max(abs(T).max(), abs(M).max()); T, M = T / s, M / s
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    z = np.exp(1j * ks)[:, None, None]
    H = M[None] + T[None] * z + T.conj().T[None] / z
    E = np.linalg.eigvalsh(H)
    pts = []
    for b in range(E.shape[1]):
        e = E[:, b]; sgn = np.sign(e)
        idx = np.where(sgn != np.roll(sgn, 1))[0]
        pts.append([float(ks[i]) for i in idx])
    return pts

if __name__ == '__main__':
    m = int(sys.argv[1]); tau = float(sys.argv[2]); out = []
    base = dict(wD=-0.3, tauA=tau, tauB=tau)
    import os
    if os.environ.get('CHARGED_PHASE'):
        base.update(eT=1.2)
    for phi in [float(x) for x in sys.argv[3:]]:
        p = dict(base)
        if os.environ.get('CHARGED_PHASE'):
            p.update(phiH=phi, phiD=phi)
        else:
            p.update(phiT=phi)
        d0, h = logical_h_precise('switch', m, p)
        T, M = blocks(h, m)
        pts = fermi_points(T, M)
        centres = []
        for p in pts:
            if len(p) == 2:
                a, b = p
                c = 0.5 * (a + b)
                if b - a > np.pi:        # the occupied arc wraps around +-pi
                    c = c + np.pi if c < 0 else c - np.pi
                centres.append(c)
        r = dict(m=m, tau=tau, phiT=phi, c_max=max_c(T, M, nk=4000)[0],
                 T_sym_defect=float(abs(T - T.T).max() / abs(T).max()),
                 T_herm_defect=float(abs(T - T.conj().T).max() / abs(T).max()),
                 fermi_points_mu0=pts, centres=centres)
        print(json.dumps(r), flush=True); out.append(r)
    tag = 'charged' if os.environ.get('CHARGED_PHASE') else 'neutral'
    json.dump(out, open(f'switch_phase_{tag}_m{m}.json', 'w'), indent=1)
