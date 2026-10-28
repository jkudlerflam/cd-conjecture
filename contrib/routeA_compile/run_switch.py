"""ED checks and logical band analysis for the central-switch cell
P L^m M R^m P (model2.LaneChain, layout 'switch').

Usage: python3 run_switch.py ed      -> switch_ed_results.json
       python3 run_switch.py bands   -> switch_band_results.json
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import sys
import json
import time
import numpy as np
import scipy.linalg as la
from model import assemble, DEFAULT
from model2 import LaneChain
from checks import second_third, gap, effective, one_body
from onebody import blocks, max_c, bands, fermi_count
import onebody2

PARAMS = dict(eT=1.0, eH=1.0, eD=1.0, wT=0.3, wH=0.3, wD=-0.3, phiT=0.0,
              tauA=0.3, tauB=0.3, kappa=1.0, muc=0.0)


def ed(m, M, K, lams, p):
    t0 = time.time()
    ch = LaneChain('switch', m, M, K=K, params=p)
    H0, V = assemble(ch.lane_parts(), ch.p)
    out = dict(m=m, cells=M, K=K, dim=ch.dim, code_dim=len(ch.code), params=ch.p)
    r, c = H0.nonzero()
    out['number_changing'] = int(np.max(np.abs(ch.number[r] - ch.number[c])))
    r, c = V.nonzero()
    out['number_changing'] = max(out['number_changing'], int(np.max(np.abs(ch.number[r] - ch.number[c]))))
    out['code_residual'] = float(abs(H0[:, ch.code]).max()) if H0[:, ch.code].nnz else 0.0
    VPP = V[ch.code][:, ch.code]
    out['first_order'] = float(abs(VPP).max()) if VPP.nnz else 0.0
    out['gap'] = gap(ch, H0)
    H2, H3 = second_third(ch, H0, V)
    E0, h, G = one_body(ch, H2)
    tr = np.trace(H2).real / len(ch.code)
    out['H2_norm_traceless'] = float(la.norm(H2 - tr * np.eye(len(ch.code)), 2))
    out['nonGaussian_H2'] = float(la.norm(H2 - G, 2))
    rows = []
    for lam in lams:
        Heff, nxt = effective(ch, H0, V, lam)
        E = la.norm(Heff - lam ** 2 * H2, 2)
        d = Heff - lam ** 2 * H2
        Etr = la.norm(d - np.trace(d).real / len(d) * np.eye(len(d)), 2)
        E3 = la.norm(Heff - lam ** 2 * H2 - lam ** 3 * H3, 2)
        rows.append(dict(lam=lam, E=float(E), E_over_lam3=float(E / lam ** 3),
                         E_traceless=float(Etr), E_after_H3_over_lam4=float(E3 / lam ** 4),
                         next_level=float(nxt)))
    out['ed'] = rows
    out['seconds'] = time.time() - t0
    print(json.dumps({k: v for k, v in out.items() if k != 'params'}), flush=True)
    return out


def band_report(m, p, mu_frac=0.3):
    """Exact order-lambda^2 one-body bond, symmetry of T, c_max, and the Fermi
    points at a chemical potential mu = mu_frac * (smallest |eigenvalue| of T)."""
    E0, h, g = onebody2.logical_h('switch', m, p)
    T, M = blocks(h, m)
    sv = la.svdvals(T)
    ev = la.eigvals(T)
    cmax, mustar = max_c(T, M, nk=4000)
    mu = mu_frac * float(np.min(np.abs(ev)))
    E = bands(T, M + 0 * M, 8000)
    ks = np.linspace(-np.pi, np.pi, 8000, endpoint=False)
    # Fermi points for H = M + T e^{ik} + h.c. - mu
    pts = []
    occ = (E < mu).sum(axis=1)
    for i in range(len(ks)):
        if occ[i] != occ[i - 1]:
            pts.append(float(ks[i]))
    return dict(m=m, gap=g, T_abs_over_max=(np.abs(T) / np.abs(T).max()).tolist(),
                T_asym=float(np.abs(T - T.T).max() / np.abs(T).max()),
                M_over_maxT=float(np.abs(M).max() / np.abs(T).max()),
                T_singular_values=sv.tolist(), c_max=cmax, mu_at_cmax=mustar,
                mu_test=mu, fermi_points_at_mu_test=pts, n_fermi_points=len(pts))


if __name__ == '__main__':
    what = sys.argv[1]
    if what == 'ed':
        res = []
        for (m, M, K) in [(1, 2, 2), (1, 3, 2), (2, 2, 2)]:
            res.append(ed(m, M, K, [0.08, 0.04, 0.02], PARAMS))
            json.dump(res, open('switch_ed_results.json', 'w'), indent=1, default=str)
    elif what == 'ed3':
        res = [ed(2, 3, 2, [0.08, 0.04, 0.02], PARAMS)]
        json.dump(res, open('switch_ed3_results.json', 'w'), indent=1, default=str)
    elif what == 'bands':
        res = []
        for (m, tau, phi) in [(1, 0.1, 0.0), (2, 0.1, 0.0), (2, 0.1, 0.05), (3, 0.01, 0.0), (3, 0.01, 0.002)]:
            p = dict(PARAMS, tauA=tau, tauB=tau, phiT=phi)
            r = band_report(m, p)
            r.update(tau=tau, phiT=phi)
            print(json.dumps(r), flush=True)
            res.append(r)
            json.dump(res, open('switch_band_results.json', 'w'), indent=1)
