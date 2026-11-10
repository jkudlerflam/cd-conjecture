#!/usr/bin/env python3
"""C193b validation of meas193.measure_state on small chains (numpy).  Run: python3 -E -s val193.py [L]

For p153 (complex bond), sym40 with the staggered gauge (tests the gauge phases), the XX chain and the two-leg XX
control, in both parity sectors (two-leg: even only, the odd sector is degenerate), the ground state is found by
z2dmrg.dmrg at a bond dimension that is exact at this L, and every correlator of measure_state is compared with
(a) the same correlator evaluated on the dense vector of that MPS (tests the contraction alone) and (b) sparse
exact diagonalization of the unrotated Hamiltonian (tests the operators, the gauge and the state).
"""
import os
import sys
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import meas193 as MS            # noqa: E402
Z, Mo = MS.Z, MS.Mo


def site_op(op, j, L):
    return sps.kron(sps.kron(sps.identity(2 ** j), sps.csr_matrix(op)), sps.identity(2 ** (L - j - 1)), format='csr')


def ed_ham(L, h2s, h1s, nnn=None):
    H = sps.csr_matrix((2 ** L, 2 ** L), dtype=complex)
    for j in range(L - 1):
        H = H + sps.kron(sps.kron(sps.identity(2 ** j), sps.csr_matrix(h2s[j])), sps.identity(2 ** (L - j - 2)))
    for j in range(L):
        H = H + site_op(h1s[j], j, L)
    if nnn is not None:
        for j in range(L - 2):
            a, b = nnn[j]
            H = H + a * site_op(Z.SP, j, L) @ site_op(Z.SM, j + 2, L) + b * site_op(Z.SM, j, L) @ site_op(Z.SP, j + 2, L)
    return H.tocsr()


def ed_sector(H, L, P):
    idx = np.array([i for i in range(2 ** L) if bin(i).count('1') % 2 == P])
    e, v = spla.eigsh(H[idx][:, idx], k=1, which='SA', tol=1e-15)
    psi = np.zeros(2 ** L, dtype=complex)
    psi[idx] = v[:, 0]
    return e[0], psi


def mps_vector(psi):
    v = psi.M[0][0]                                       # s b
    for A in psi.M[1:]:
        v = np.tensordot(v, A, axes=([-1], [0]))
    return v.reshape(-1)


def mps_from_vector(vec, L, P, gph):
    """Exact right-canonical Z2 MPS (centre at site 0, labels sorted even first) of a dense vector given in the
    unrotated frame, rotated into the gauge frame of the DMRG (U = prod u_j, u = diag(e^{-i th/2}, e^{i th/2}))."""
    v = vec.reshape([2] * L).astype(complex)
    for j in range(L):
        th = -np.angle(gph[j])
        u = np.array([np.exp(-1j * th / 2), np.exp(1j * th / 2)])
        sh = [1] * L
        sh[j] = 2
        v = v * u.reshape(sh)
    rest = v.reshape(2 ** L, 1)                  # rows: configs of sites 0..j, cols: right bond label
    qb = np.array([P])                           # parity labels of the right bond
    Ms = [None] * L
    n0 = [None] * (L + 1)
    n0[L] = 1 if P == 0 else 0
    for j in range(L - 1, 0, -1):
        T = rest.reshape(2 ** j, 2, -1)           # left configs, s_j, b
        lpar = np.array([bin(c).count('1') % 2 for c in range(2 ** j)])
        blocks, labels = [], []
        for q in (0, 1):
            rows = np.where(lpar == q)[0]
            cols = [(s, b) for s in (0, 1) for b in range(len(qb)) if (q ^ s) == qb[b]]
            Mb = np.array([[T[r, s, b] for (s, b) in cols] for r in rows])
            U, S, Vh = np.linalg.svd(Mb, full_matrices=False)
            k = int(np.sum(S > 1e-14 * max(S.max(), 1e-300))) if S.size else 0
            blocks.append((q, rows, cols, U[:, :k] * S[:k], Vh[:k]))
        Da = sum(b[3].shape[1] for b in blocks)
        B = np.zeros((Da, 2, len(qb)), dtype=complex)
        newrest = np.zeros((2 ** j, Da), dtype=complex)
        off = 0
        for q, rows, cols, US, Vh in blocks:
            k = Vh.shape[0]
            for c, (s, b) in enumerate(cols):
                B[off:off + k, s, b] = Vh[:, c]
            newrest[rows, off:off + k] = US
            off += k
        n0[j] = blocks[0][4].shape[0]
        qb = np.array([0] * n0[j] + [1] * (Da - n0[j]))
        Ms[j] = B
        rest = newrest
    M0 = np.zeros((1, 2, len(qb)), dtype=complex)
    for s in (0, 1):
        M0[0, s, :] = rest[s, :]
    Ms[0] = M0
    n0[0] = 1
    return Z.MPS(Ms, n0, P)


def dense_obs(vec, L, sp, sm, out):
    """The quantities of measure_state on a dense vector with site operators sp[j], sm[j]."""
    n = np.diag([0.0, 1.0])
    ex = lambda O: complex(np.vdot(vec, O @ vec))
    so = lambda op, j: site_op(op, j, L)
    res = {}
    dens = np.array([ex(so(n, j)).real for j in range(L)])
    res['dens'] = dens
    S = []
    for l in range(1, L):
        s = np.linalg.svd(vec.reshape(2 ** l, -1), compute_uv=False) ** 2
        S.append(Z.entropy(s))
    res['S'] = np.array(S)
    gb = {}
    for i_s, row in out['Gb_rows'].items():
        i = int(i_s)
        for t in range(len(row[0])):
            j = i + 1 + t
            gb[(i, j)] = ex(so(sm[i], i) @ so(sp[j], j))
    res['Gb'] = gb
    cnn = {}
    for i0 in out['REFS']:
        cnn[i0] = [ex(so(n, i0) @ so(n, j)).real - dens[i0] * dens[j] for j in range(i0 + 1, L)]
    res['Cnn'] = cnn
    dp, mm, pp = {}, {}, {}
    for d in (2, 3):
        for j in range(L - d):
            pp[(d, j)] = ex(so(sp[j], j) @ so(sp[j + d], j + d))
        for i0 in out['REFS']:
            if i0 + d >= L:
                continue
            A = so(sm[i0], i0) @ so(sm[i0 + d], i0 + d)
            mm[(d, i0)] = ex(A)
            dp[(d, i0)] = [ex(A @ so(sp[j], j) @ so(sp[j + d], j + d)) for j in range(i0 + d + 1, L - d)]
    res['Dpair'], res['mm'], res['pp'] = dp, mm, pp
    return res


def compare(out, ref):
    e = {}
    e['dens'] = np.max(np.abs(np.array(out['dens']) - ref['dens']))
    e['S'] = np.max(np.abs(np.array(out['S']) - ref['S']))
    e['Gb'] = 0.0
    for i_s, row in out['Gb_rows'].items():
        i = int(i_s)
        for t, (a, b) in enumerate(zip(*row)):
            e['Gb'] = max(e['Gb'], abs(a + 1j * b - ref['Gb'][(i, i + 1 + t)]))
    e['Cnn'] = max(np.max(np.abs(np.array(out['Cnn'][str(i0)]) - np.array(ref['Cnn'][i0]))) for i0 in out['REFS'])
    e['Dpair'] = 0.0
    for key, row in out['Dpair'].items():
        d, i0 = map(int, key.split('_'))
        if row:
            e['Dpair'] = max(e['Dpair'], np.max(np.abs(np.array([a + 1j * b for a, b in row]) - np.array(ref['Dpair'][(d, i0)]))))
    e['mm'] = max(abs(complex(*v) - ref['mm'][tuple(map(int, k.split('_')))]) for k, v in out['mm'].items())
    e['pp'] = max(abs(complex(*v) - ref['pp'][(int(d), int(j))]) for d, dd in out['pp'].items() for j, v in dd.items())
    ndp = sum(len(r) for r in out['Dpair'].values())
    return e, ndp


def run_case(label, mpo, gph, cfg, h2o, h1o, L, nnn=None, sectors=(0, 1)):
    H = ed_ham(L, h2o, h1o, nnn)
    worst = {}
    for P in sectors:
        e_ed, v_ed = ed_sector(H, L, P)
        dt = complex if mpo.dtype is complex else float
        psi = Z.random_mps(L, P, D0=4, seed=11 + P, dtype=dt)
        chi = 2 ** (L // 2)
        sw = [dict(chi=16, alpha=1e-3, res_tol=1e-8), dict(chi=chi, alpha=1e-5, res_tol=1e-10)] + \
             [dict(chi=chi, alpha=0.0, res_tol=1e-14, maxit=2000)] * 10
        psi, hist = Z.dmrg(psi, mpo, sw)
        c = dict(cfg, L=L, P=P, _n0=list(psi.n0))
        out = MS.measure_state(psi.M, c, np, log=lambda s: None)
        E = out['E_mpo'] + mpo.const
        # (a) the MPS's own dense vector, rotated-frame operators
        sp = [gph[j] * Z.SP for j in range(L)]
        sm = [np.conj(gph[j]) * Z.SM for j in range(L)]
        ea, ndp = compare(out, dense_obs(mps_vector(psi), L, sp, sm, out))
        # (b) ED of the unrotated Hamiltonian, plain operators
        eb, _ = compare(out, dense_obs(v_ed, L, [Z.SP] * L, [Z.SM] * L, out))
        eb['E'] = abs(E - e_ed)
        # (c) the ED vector itself as an exact MPS, through measure_state (operators, gauge, contraction)
        pe = mps_from_vector(v_ed, L, P, gph)
        assert pe.check(1e-12)
        oe = MS.measure_state(pe.M, dict(cfg, L=L, P=P, _n0=list(pe.n0)), np, log=lambda s: None)
        ec, _ = compare(oe, dense_obs(v_ed, L, [Z.SP] * L, [Z.SM] * L, oe))
        ec['E'] = abs(oe['E_mpo'] + mpo.const - e_ed)
        print('   ED-MPS vs ED:  ' + ' '.join(f'{k}={v:.1e}' for k, v in ec.items()))
        for k, v in ec.items():
            worst['EDMPS_' + k] = max(worst.get('EDMPS_' + k, 0.0), v)
        Smps = Z.entropies(psi)
        ea['S_vs_entropies'] = np.max(np.abs(Smps - np.array(out['S'])))
        ea['rc_err'] = out['rc_err']
        ea['norm-1'] = abs(out['norm'] - 1)
        print(f'{label} L={L} P={P} (pair values {ndp})')
        print('   vs MPS vector: ' + ' '.join(f'{k}={v:.1e}' for k, v in ea.items()))
        print('   vs ED:         ' + ' '.join(f'{k}={v:.1e}' for k, v in eb.items()))
        for k, v in list(ea.items()) + [('ED_' + k, v) for k, v in eb.items()]:
            worst[k] = max(worst.get(k, 0.0), v)
    return worst


if __name__ == '__main__':
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    allw = {}
    for name, gauge in (('p153', None), ('sym40', 'stag')):
        b = Mo.bond_family(name)
        hz = b['hz0']
        mpo, gph, shift, info = Mo.candidate(name, L, hz, gauge=gauge)
        cfg = dict(model=name, hz=hz, gauge=gauge)
        allw[f'{name}/{gauge}'] = run_case(f'{name}/{gauge}', mpo, gph, cfg, [info['h2']] * (L - 1),
                                           [hz * Z.ZZ_ + shift * Z.I2] * L, L)
    n_op = np.diag([0.0, 1.0])
    for legs in (1, 2):
        mpo, gph = Mo.xx_control(L, mu=-1.2, legs=legs)
        h2o = [-(np.kron(Z.SP, Z.SM) + np.kron(Z.SM, Z.SP))] * (L - 1) if legs == 1 else [np.zeros((4, 4))] * (L - 1)
        nnn = None if legs == 1 else [(-1.0, -1.0)] * (L - 2)
        cfg = dict(model=f'xx{legs}', mu=-1.2)
        allw[f'xx{legs}'] = run_case(f'xx{legs}', mpo, gph, cfg, h2o, [1.2 * n_op] * L, L, nnn=nnn,
                                     sectors=(0, 1) if legs == 1 else (0,))
    print('worst per case:')
    for k, w in allw.items():
        print(' ', k, ' '.join(f'{a}={v:.1e}' for a, v in w.items()))
