#!/usr/bin/env python3
"""Exact checks for the 13-letter two-flavour bond (task 1 at n = 2).

1. One cell (P_L B P_R, 324 states): spectrum, exact 4-fold code, gap, both flavour charges of the code,
   compression of every port operator (all 3 x 3 matrix units on x_L and on x_R, and on y).
2. Two cells: the second-order bond (exact resolvent), decomposed into flavour hopping (with phase),
   intra-flavour density couplings, fields and inter-flavour couplings; exact des Cloizeaux bond at
   several lam with the traceless residual / lam^3 and the size of inter-flavour terms.
3. Control kap = 0: the second-order bond is a scalar (C18).
4. Symmetries of h: both U(1) charges.
Writes data/ed_check2.json.  Claude (Opus 5.5), 2026-10-28.
"""
import os, json
import numpy as np
from scipy import linalg
import cell2 as C
import effective2 as E

HERE = os.path.dirname(os.path.abspath(__file__))
PAR = dict(E.DEFAULT)
LAMS = [0.2, 0.15, 0.1, 0.05, 0.025]


def one_cell():
    H0, ports, code, qc = C.cell_ops(PAR['DN'], PAR['DC'], PAR['kap'], PAR['sg'], PAR['DCd'])
    w = linalg.eigvalsh(H0)
    comp = 0.0
    for pos in ('L', 'R'):
        for a in range(3):
            for b in range(3):
                for fac in ('x', 'y'):
                    M = np.zeros((3, 3)); M[a, b] = 1
                    if a == b:
                        M = M - np.eye(3) / 3          # traceless part
                    op9 = np.kron(M, np.eye(3)) if fac == 'x' else np.kron(np.eye(3), M)
                    O = np.kron(op9, np.eye(36)) if pos == 'L' else np.kron(np.eye(36), op9)
                    comp = max(comp, float(np.abs(code.conj().T @ O @ code).max()))
    qcode = [np.real(np.diag(code.conj().T @ np.diag(qc[:, k]) @ code)).round(12).tolist() for k in range(2)]
    return dict(low=np.round(w[:12], 10).tolist(), degeneracy=int((np.abs(w) < 1e-10).sum()), gap=float(w[4]),
                code_residual=float(np.abs(H0 @ code).max()), max_traceless_port_compression=comp, code_charges=qcode)


def summarize(H):
    c, r = E.fit(H)
    out = dict(fit_residual=r)
    for a in (1, 2):
        t = c[f'hop{a}'] + 1j * c[f'hopi{a}']
        out[f't{a}'] = abs(t); out[f'phase{a}_over_pi'] = float(np.angle(t) / np.pi)
        out[f'nn{a}'] = c[f'nn{a}']; out[f'field{a}'] = c[f'n{a}_1'] + c[f'n{a}_2']
    out['interflavour'] = max(abs(c[k]) for k in ('n1n2_12', 'n2n1_12', 'n1n2_11', 'n1n2_22'))
    out['coef'] = c
    return out


if __name__ == '__main__':
    res = dict(params={k: (list(v) if isinstance(v, tuple) else v) for k, v in PAR.items()})
    res['one_cell'] = one_cell()
    print('one cell:', {k: v for k, v in res['one_cell'].items() if k != 'low'})
    first, H2, leak, gap = E.second_order(PAR)
    res['first_order_max'] = float(np.abs(first).max())
    res['H2'] = summarize(H2)
    print('first order', res['first_order_max'], ' H2:', {k: v for k, v in res['H2'].items() if k != 'coef'})
    p0 = dict(PAR); p0['kap'] = (0.0, 0.0)
    _, H20, _, _ = E.second_order(p0)
    res['kappa0_traceless_norm'] = float(np.linalg.norm(H20 - np.trace(H20) / 16 * np.eye(16), 2))
    print('kappa = 0: traceless norm of H2 =', res['kappa0_traceless_norm'])
    res['exact'] = []
    for lam in LAMS:
        He, w, sv = E.exact_eff(PAR, lam)
        s = summarize(He)
        r = He - lam ** 2 * H2
        r = r - np.trace(r) / 16 * np.eye(16)
        row = dict(lam=lam, gap_16_17=float(w[16] - w[15]), min_overlap=sv, res_over_lam3=float(np.linalg.norm(r, 2) / lam ** 3),
                   **{k: v for k, v in s.items() if k != 'coef'})
        res['exact'].append(row)
        print(f"lam={lam}: t1/lam2={row['t1'] / lam ** 2:.5f} (phase {row['phase1_over_pi']:.4f} pi) t2/lam2={row['t2'] / lam ** 2:.5f} "
              f"(phase {row['phase2_over_pi']:.4f} pi) inter/lam2={row['interflavour'] / lam ** 2:.2e} res/lam3={row['res_over_lam3']:.4f} "
              f"fitres={row['fit_residual']:.1e}")
    os.makedirs(os.path.join(HERE, 'data'), exist_ok=True)
    json.dump(res, open(os.path.join(HERE, 'data', 'ed_check2.json'), 'w'), indent=1, default=float)
