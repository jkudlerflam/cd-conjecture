#!/usr/bin/env python3
"""Decompose the cluster-ED logical bond h2(lam) (data/heff_T.npz) into fermion couplings.

With global fermions c_{j,a} (JW across cells) on the two cells of a bond, the parity-even
traceless part of h2 is written as
    sum_{k,l} T_kl c_k^dag c_l + sum_{k<l} (P_kl c_k c_l + h.c.) + (Majorana degree >= 4),
k, l in {(j,1), (j,2), (j+1,1), (j+1,2)}.  Reported in units of J = J0 lam^2:
  t_a, theta_a      hopping c^dag_{j,a} c_{j+1,a}: -t e^{-i theta}
  mu_a              on-site (both cells of the bond summed)
  t_12              inter-flavour hopping (largest |.|)
  D_aa              same-flavour pairing c_{j,a} c_{j+1,a}
  D_12_onsite, D_12_nn   inter-flavour pairing on one cell and across the bond
  inter             operator norm of the degree >= 4 part
plus the dressed Fermi points of the number-conserving part, delta_2 of the dressed points, and
the smallest |p + p'| mod 2pi over pairs of Fermi points (pairing resonance; 0 means resonant).
usage: analyze_heff.py data/heff_INC.npz
Written by Claude (Opus 5.5) for contrib/routeA_incommensurate, 2026-10-28.
"""
import sys, json
import numpy as np
from scipy import linalg
import ed_check as EC


def coeffs(h, n):
    c, _ = EC.majoranas_two_cells(n)
    dim = h.shape[0]
    h = h - np.trace(h) / dim * np.eye(dim)
    m = len(c)
    T = np.zeros((m, m), complex)
    for k in range(m):
        for l in range(m):
            O = c[k].conj().T @ c[l]
            if k == l:
                O = O - 0.5 * np.eye(dim)
            T[k, l] = np.trace(O.conj().T @ h) / np.trace(O.conj().T @ O)
    Pm = np.zeros((m, m), complex)
    for k in range(m):
        for l in range(k + 1, m):
            O = c[k] @ c[l]
            Pm[k, l] = np.trace(O.conj().T @ h) / np.trace(O.conj().T @ O)
    parts, _ = EC.majorana_decompose(h, n)
    return T, Pm, parts


def main(path):
    z = np.load(path)
    J0 = float(z['J0'])
    lams = [float(k.split('_')[1]) for k in z.files if k.startswith('h2_')]
    n = 2
    out = []
    for lam in sorted(lams):
        h = z[f'h2_{lam}'] / (J0 * lam ** 2)
        T, Pm, parts = coeffs(h, n)
        rec = dict(lam=lam)
        for a in range(n):
            hop = T[a, n + a]                      # coefficient of c^dag_{j,a} c_{j+1,a}
            rec[f't{a+1}'] = float(abs(hop)); rec[f'theta{a+1}'] = float(np.angle(-hop))
            rec[f'mu{a+1}'] = float(-(T[a, a] + T[n + a, n + a]).real)
            rec[f'D{a+1}{a+1}'] = float(abs(Pm[a, n + a]))
        rec['t12'] = float(max(abs(T[0, n + 1]), abs(T[1, n + 0]), abs(T[0, 1]) * 0 + abs(T[0, 1])))
        rec['D12_onsite'] = float(max(abs(Pm[0, 1]), abs(Pm[n, n + 1])))
        rec['D12_nn'] = float(max(abs(Pm[0, n + 1]), abs(Pm[1, n])))
        rec['inter'] = float(np.linalg.norm(parts[4], 2))
        cons = sum(T[k, l] * (EC.majoranas_two_cells(n)[0][k].conj().T @ EC.majoranas_two_cells(n)[0][l]) for k in range(2 * n) for l in range(2 * n))
        pts, _, _ = EC.band_from_bond(cons, n)
        allp = np.array(sorted(sum(pts, [])))
        rec['fermi_points'] = allp.round(4).tolist()
        pm = np.concatenate([allp, -allp]) % (2 * np.pi)
        ps = np.sort(pm); gaps = np.diff(np.concatenate([ps, [ps[0] + 2 * np.pi]]))
        rec['delta2'] = float(gaps.min())
        s = [(allp[i] + allp[j]) % (2 * np.pi) for i in range(len(allp)) for j in range(i, len(allp))]
        rec['pair_resonance'] = float(min(min(x, 2 * np.pi - x) for x in s))
        out.append(rec)
        print(json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in rec.items()}))
    json.dump(out, open(path.replace('.npz', '_coeffs.json'), 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1])
