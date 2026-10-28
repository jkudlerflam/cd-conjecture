#!/usr/bin/env python3
"""Exact checks of the cell and of the second-order logical bond (task 1).

1. One cell (P B P, dimension 32): spectrum, exact 2-fold code at zero energy, gap, charge of the code,
   and the compression G^dag O G of every port operator (all 2 x 2 Hermitian O on pL and on pR).
2. Two cells: second-order bond H2 (closed form Jxy = Jp Jz kap^2 / (4 D (D^2 - kap^2)) for D0 = D1 = D),
   exact des Cloizeaux effective bond for several lam, residual / lam^3.
3. Control kap = 0 (pure Bell checks, frustration free on port and partner only): H2 is a scalar,
   as the C18 obstruction requires.
4. Three cells: exact connected three-cell term h3 = Heff(3) - h2 (x) 1 - 1 (x) h2 at each lam.
Writes data/ed_check.json.  Claude (Opus 5.5), 2026-10-27.
"""
import os, json
import numpy as np
from scipy import linalg
import scipy.sparse as sps
from scipy.sparse.linalg import eigsh
import cell as C
import effective as E

HERE = os.path.dirname(os.path.abspath(__file__))
P = dict(D0=1.0, D1=1.0, kap=0.6, Jp=2.0, Jz=-1.0)
LAMS = [0.2, 0.15, 0.1, 0.05, 0.025, 0.0125]


def one_cell():
    H0, ports, code, qc = C.cell_ops(P['D0'], P['D1'], P['kap'])
    w = linalg.eigvalsh(H0)
    paulis = [C.X, np.array([[0, -1j], [1j, 0]]), C.Z]
    comp = []
    for side in ('L', 'R'):
        for s in paulis:
            if side == 'L':
                O = np.kron(s, np.eye(16))
            else:
                O = np.kron(np.eye(16), s)
            comp.append(float(np.abs(code.T @ O @ code).max()))
    qcode = code.T @ np.diag(qc) @ code
    return dict(spectrum=np.round(w, 12).tolist(), code_residual=float(np.abs(H0 @ code).max()),
                degeneracy=int((np.abs(w) < 1e-10).sum()), gap=float(w[2]),
                max_port_compression=max(comp), code_charge=np.round(np.diag(qcode), 12).tolist())


def two_cells():
    first, H2, gap = E.second_order(P['D0'], P['D1'], P['kap'], P['Jp'], P['Jz'])
    s2 = E.summary(H2)
    D, k = P['D0'], P['kap']
    analytic = P['Jp'] * P['Jz'] * k ** 2 / (4 * D * (D ** 2 - k ** 2))
    rows = []
    for lam in LAMS:
        He, w, smin = E.exact_eff(lam, P['D0'], P['D1'], P['kap'], P['Jp'], P['Jz'])
        s = E.summary(He)
        r = He - lam ** 2 * H2
        rows.append(dict(lam=lam, Jxy=s['Jxy'], Jzz=s['Jzz'], a0=s['a0'], other=s['other'], h=s['h'], anis=s['anis'],
                         res=float(np.linalg.norm(r, 2)), res_over_lam3=float(np.linalg.norm(r, 2) / lam ** 3),
                         res_traceless_over_lam3=float(np.linalg.norm(r - np.trace(r) / 4 * np.eye(4), 2) / lam ** 3),
                         min_overlap=float(smin), low=w[:6].tolist()))
    f0, H2k0, _ = E.second_order(P['D0'], P['D1'], 0.0, P['Jp'], P['Jz'])
    s0 = E.summary(H2k0)
    return dict(first_order_max=float(np.abs(first).max()), H2=s2, Jxy_analytic=analytic, pair_gap=gap,
                exact=rows, kappa0=dict(H2=s0, traceless_norm=float(np.linalg.norm(H2k0 - np.trace(H2k0) / 4 * np.eye(4), 2))))


def three_cells(lam):
    H0, ports, code, qc = C.cell_ops(P['D0'], P['D1'], P['kap'])
    I = sps.identity(32, format='csr')
    S = lambda op: sps.csr_matrix(op)
    H0s = S(H0)
    H = sps.kron(sps.kron(H0s, I), I) + sps.kron(sps.kron(I, H0s), I) + sps.kron(sps.kron(I, I), H0s)

    def V(a, b):
        R, L = ports['R'], ports['L']
        ops = [(P['Jp'], R['sp'], L['sm']), (P['Jp'], R['sm'], L['sp']), (P['Jz'], R['z'], L['z'])]
        out = 0
        for c, A, B in ops:
            mats = [I, I, I]
            mats[a] = S(A); mats[b] = S(B)
            out = out + c * sps.kron(sps.kron(mats[0], mats[1]), mats[2])
        return out
    H = H + lam * (V(0, 1) + V(1, 2))
    w, U = eigsh(H.tocsc(), k=10, which='SA', tol=1e-13)
    o = np.argsort(w); w, U = w[o], U[:, o]
    G = np.kron(np.kron(code, code), code)
    Phi = G.T @ U[:, :8]
    O = Phi @ np.real(linalg.inv(linalg.sqrtm(Phi.T @ Phi)))
    H3 = O @ np.diag(w[:8]) @ O.T
    He2, _, _ = E.exact_eff(lam, P['D0'], P['D1'], P['kap'], P['Jp'], P['Jz'])
    h3 = H3 - np.kron(He2, np.eye(2)) - np.kron(np.eye(2), He2)
    # Pauli decomposition of h3
    PA = dict(I=np.eye(2), X=C.X, Y=np.array([[0, -1j], [1j, 0]]), Z=C.Z)
    coef = {}
    for a in 'IXYZ':
        for b in 'IXYZ':
            for c in 'IXYZ':
                v = float(np.real(np.trace(np.kron(np.kron(PA[a], PA[b]), PA[c]).conj().T @ h3)) / 8)
                if abs(v) > 1e-12:
                    coef[a + b + c] = v
    return dict(lam=lam, gap_8_9=float(w[8] - w[7]), min_overlap=float(linalg.svdvals(Phi).min()),
                h3_norm=float(np.linalg.norm(h3 - np.trace(h3) / 8 * np.eye(8), 2)), h3_pauli=coef)


if __name__ == '__main__':
    out = dict(params=P, one_cell=one_cell(), two_cells=two_cells())
    print(json.dumps(out['one_cell'], indent=0)[:600])
    t = out['two_cells']
    print('first order', t['first_order_max'], 'H2', t['H2'], 'analytic Jxy', t['Jxy_analytic'])
    print('kappa = 0 control:', t['kappa0'])
    for r in t['exact']:
        print(f"lam={r['lam']:.4f} Jxy/lam2={r['Jxy'] / r['lam'] ** 2:.5f} Jzz/lam2={r['Jzz'] / r['lam'] ** 2:.5f} "
              f"other={r['other']:.1e} res/lam3={r['res_over_lam3']:.4f} traceless res/lam3={r['res_traceless_over_lam3']:.4f}")
    out['three_cells'] = []
    for lam in (0.15, 0.1, 0.05):
        r = three_cells(lam)
        out['three_cells'].append(r)
        big = sorted(r['h3_pauli'].items(), key=lambda kv: -abs(kv[1]))[:8]
        print(f"three cells lam={lam}: gap above 8 = {r['gap_8_9']:.4f}, |h3 traceless| = {r['h3_norm']:.3e} "
              f"(/lam^3 = {r['h3_norm'] / lam ** 3:.3f}), largest:", [(k, round(v / lam ** 3, 4)) for k, v in big])
    os.makedirs(os.path.join(HERE, 'data'), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, 'data', 'ed_check.json'), 'w'), indent=1)
