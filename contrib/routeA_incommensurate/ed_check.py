#!/usr/bin/env python3
"""Two- and three-cell exact diagonalization of the full route-A cell chain, and the logical
effective Hamiltonian at finite lam.

Cluster (open, m cells):  H_m = sum_j H0_j + lam sum_{j<m} sum_a FR_a(j) FL_a(j+1),  Delta = 1.
Low band: the N^m lowest eigenvectors Psi (energies E).  With B = (G^{(x)m})^dag Psi and its polar
factor W, the des Cloizeaux (canonical Schrieffer-Wolff) Hamiltonian is H_dc = W E W^dag; it is
linked-cluster additive, so h2(lam) = H_dc(2 cells) and h3(lam) = H_dc(3) - h2 (x) 1 - 1 (x) h2.

Reported: ||h2(lam) - lam^2 h2_pred||, ||H_dc(3) - lam^2 (h2_pred (x) 1 + 1 (x) h2_pred)|| and ||h3||
(after removing the trace), with their power laws in lam; and a Majorana-degree decomposition of
h2(lam) (Gaussian part, its pairing (U(1)-breaking) piece, and the interaction part of degree >= 4),
plus the Fermi points of the Gaussian part.

usage: python3 ed_check.py [--target INC|SYM] [--lams ...] [--lams3 ...] [--out data/heff_INC.npz]
Written by Claude (Opus 5.5) for contrib/routeA_incommensurate, 2026-10-27.
"""
import os, sys, json, time, argparse, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from scipy import linalg
from scipy.sparse.linalg import LinearOperator, eigsh
import cell as CE


def eigbasis(cell):
    """Cell data in the eigenbasis of H0 (code = first N states, in the logical basis of G)."""
    H0, G, N = cell['H0'], cell['G'], cell['N']
    udiag = np.real(np.diag(cell['U']))           # U is diagonal in the physical basis
    cols, ws = [], []
    for sgn in (1, -1):                           # diagonalize H0 inside each parity sector
        idx = np.where(udiag == sgn)[0]
        w, V = linalg.eigh(H0[np.ix_(idx, idx)])
        Vf = np.zeros((H0.shape[0], len(idx)), complex); Vf[idx] = V
        cols.append(Vf[:, w > 1e-9]); ws.append(w[w > 1e-9])
    exc = np.hstack(cols)
    o = np.argsort(np.concatenate(ws), kind='stable')
    exc = exc[:, o]
    Vb = np.hstack([G, exc])                      # code columns are exactly G
    e = np.concatenate([np.zeros(N), np.concatenate(ws)[o]])
    assert np.allclose(Vb.conj().T @ Vb, np.eye(len(e)), atol=1e-10)
    FL = [Vb.conj().T @ F @ Vb for F in cell['FL']]
    FR = [Vb.conj().T @ F @ Vb for F in cell['FR']]
    U = np.diag(np.rint(np.real(np.diag(Vb.conj().T @ cell['U'] @ Vb))))
    return dict(e=e, FL=FL, FR=FR, U=U, N=N, D=len(e))


def cluster_op(cb, m, lam):
    e, FR, FL, D = cb['e'], cb['FR'], cb['FL'], cb['D']
    dim = D ** m

    def mv(x):
        x = x.reshape([D] * m + [-1])
        out = np.zeros_like(x)
        for j in range(m):
            sh = [1] * (m + 1); sh[j] = D
            out += e.reshape(sh) * x
        for j in range(m - 1):
            for a in range(3):
                t = np.moveaxis(np.tensordot(FR[a], x, ([1], [j])), 0, j)
                t = np.moveaxis(np.tensordot(FL[a], t, ([1], [j + 1])), 0, j + 1)
                out += lam * t
        return out.reshape(dim, -1)
    return LinearOperator((dim, dim), matvec=lambda v: mv(v)[:, 0], matmat=mv, dtype=complex), dim


def low_band(cb, m, lam, extra=4):
    N, D = cb['N'], cb['D']
    k = N ** m
    op, dim = cluster_op(cb, m, lam)
    if dim <= 5000:
        Hd = op.matmat(np.eye(dim, dtype=complex))
        w, U = linalg.eigh((Hd + Hd.conj().T) / 2)
        return w[:k], U[:, :k], w[k:k + extra]
    # block Davidson from the code space: Rayleigh-Ritz on [Q, T], T = preconditioned residual
    # with (diag(H0 sum) - w)^{-1} on excited configurations; the band is ~Delta below the rest.
    idx = np.array(list(itertools.product(range(N), repeat=m)))
    Q = np.zeros((dim, k), complex)
    flat = np.ravel_multi_index(tuple(idx.T), [D] * m)
    Q[flat, np.arange(k)] = 1.0
    diag = np.zeros([D] * m)
    for j in range(m):
        sh = [1] * m; sh[j] = D
        diag = diag + cb['e'].reshape(sh)
    diag = diag.ravel()
    exc = diag > 0.5
    nxt = None
    for it in range(300):
        HQ = op.matmat(Q)
        Hs = Q.conj().T @ HQ
        w, s = linalg.eigh(0.5 * (Hs + Hs.conj().T))
        Q = Q @ s; HQ = HQ @ s
        R = HQ - Q * w[None, :]
        res = float(np.linalg.norm(R, axis=0).max())
        if res < 1e-12:
            break
        den = np.where(exc[:, None], diag[:, None] - w[None, :], 1.0)
        T = R / den
        for _ in range(2):
            T = T - Q @ (Q.conj().T @ T)
        Qb, _ = np.linalg.qr(np.hstack([Q, T]))
        HQb = op.matmat(Qb)
        Hs = Qb.conj().T @ HQb
        w2, s2 = linalg.eigh(0.5 * (Hs + Hs.conj().T))
        nxt = w2[k]
        Q = Qb @ s2[:, :k]
    HQ = op.matmat(Q)
    Hs = Q.conj().T @ HQ
    w, s = linalg.eigh(0.5 * (Hs + Hs.conj().T))
    return w, Q @ s, (None if nxt is None else np.array([nxt])), res


def heff(cb, m, lam):
    N, D = cb['N'], cb['D']
    out = low_band(cb, m, lam)
    E, Psi, above = out[0], out[1], out[2]
    res = out[3] if len(out) > 3 else 0.0
    B = Psi.reshape([D] * m + [-1])[tuple([slice(0, N)] * m)].reshape(N ** m, -1)
    u, s, vh = linalg.svd(B)
    W = u @ vh
    H = (W * E) @ W.conj().T
    H = 0.5 * (H + H.conj().T)
    info = dict(band_top=float(E.max()), band_bottom=float(E.min()),
                next=None if above is None else float(np.min(above)), min_overlap_sv=float(s.min()),
                eig_res=float(res))
    return H, info


# --------------------------------------------------------------------------- Majorana analysis
def majoranas_two_cells(n):
    """2n complex modes on two cells (JW order: cell j flavours, then cell j+1), as 4n Majoranas."""
    a, P = CE.fermion_ops(n)
    N = 2 ** n
    I = np.eye(N)
    c = [np.kron(x, I) for x in a] + [np.kron(P, x) for x in a]
    g = []
    for x in c:
        g.append(x + x.conj().T)
        g.append(1j * (x.conj().T - x))
    return c, g


def majorana_decompose(h, n):
    """Split a two-cell parity-even operator into Majorana degree 0, 2, >=4 parts."""
    _, g = majoranas_two_cells(n)
    M = len(g)
    dim = h.shape[0]
    parts = {0: np.trace(h) / dim * np.eye(dim)}
    h2 = np.zeros_like(h)
    A = np.zeros((M, M), complex)       # h2 = (i/4) sum A_kl g_k g_l, A real antisymmetric
    for k in range(M):
        for l in range(k + 1, M):
            O = 1j * g[k] @ g[l]
            coef = np.trace(O @ h) / dim
            h2 += coef * O
            A[k, l] = coef
            A[l, k] = -coef
    parts[2] = h2
    parts[4] = h - parts[0] - h2
    return parts, A


def gaussian_split(h, n):
    """Number-conserving and pairing parts of the Gaussian (degree-2) piece."""
    c, _ = majoranas_two_cells(n)
    dim = h.shape[0]
    cons = np.zeros_like(h)
    pair = np.zeros_like(h)
    m = len(c)
    # orthonormal-ish basis: c_k^dag c_l (k != l and k = l shifted), c_k c_l, c_k^dag c_l^dag
    ops_c, ops_p = [], []
    for k in range(m):
        for l in range(m):
            O = c[k].conj().T @ c[l]
            if k == l:
                O = O - 0.5 * np.eye(dim)
            ops_c.append(O)
    for k in range(m):
        for l in range(k + 1, m):
            ops_p.append(c[k] @ c[l]); ops_p.append(c[l].conj().T @ c[k].conj().T)
    for O in ops_c:
        cons += np.trace(O.conj().T @ h) / np.trace(O.conj().T @ O) * O
    for O in ops_p:
        pair += np.trace(O.conj().T @ h) / np.trace(O.conj().T @ O) * O
    return cons, pair


def band_from_bond(h, n):
    """Single-particle bands of the translation-invariant chain whose Gaussian bond is h (number-
    conserving part).  Returns Fermi points (zeros of the bands) per band."""
    c, _ = majoranas_two_cells(n)
    dim = h.shape[0]
    # t_{kl} with H = sum t_kl c_k^dag c_l: on-site (both cells) and hopping blocks
    T = np.zeros((2 * n, 2 * n), complex)
    for k in range(2 * n):
        for l in range(2 * n):
            O = c[k].conj().T @ c[l]
            if k == l:
                O = O - 0.5 * np.eye(dim)
            T[k, l] = np.trace(O.conj().T @ h) / np.trace(O.conj().T @ O)
    eps0 = T[:n, :n] + T[n:, n:]           # on-site, both cells of the bond land on one site per bond
    t1 = T[:n, n:]                          # c^dag_j c_{j+1}
    ks = np.linspace(-np.pi, np.pi, 20001)
    bands = []
    for kk in ks:
        Hk = eps0 + t1 * np.exp(1j * kk) + t1.conj().T * np.exp(-1j * kk)
        bands.append(linalg.eigvalsh(0.5 * (Hk + Hk.conj().T)))
    bands = np.array(bands)
    pts = []
    for b in range(n):
        s = np.sign(bands[:, b])
        idx = np.where(np.diff(s) != 0)[0]
        pts.append(sorted((ks[idx] % (2 * np.pi)).tolist()))
    return pts, eps0, t1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--target', default='INC')
    ap.add_argument('--lams', type=float, nargs='+', default=[0.0025, 0.005, 0.01, 0.02, 0.04, 0.05, 0.08, 0.1, 0.15, 0.2])
    ap.add_argument('--lams3', type=float, nargs='*', default=[0.005, 0.01, 0.02, 0.04, 0.05, 0.1, 0.2])
    ap.add_argument('--out', default=None)
    a = ap.parse_args()
    par = getattr(CE, a.target)
    cell = CE.build(par, extra=6 if a.target == 'INC3' else 0)
    n, N = cell['n'], cell['N']
    cb = eigbasis(cell)
    h2p = CE.expected_h2(cell)
    J0 = cell['pred']['J0']
    out = a.out or os.path.join(HERE, 'data', f'heff_{a.target}.npz')
    res, log = {}, []
    print(json.dumps(dict(target=a.target, J0=J0, fermi=CE.fermi_data(par))), flush=True)
    for lam in a.lams:
        t0 = time.time()
        H2, i2 = heff(cb, 2, lam)
        r = H2 - lam ** 2 * h2p
        r = r - np.trace(r) / N ** 2 * np.eye(N ** 2)          # a scalar shift is irrelevant
        parts, _ = majorana_decompose(H2, n)
        cons, pair = gaussian_split(H2, n)
        pts, eps0, t1 = band_from_bond(cons, n)
        rec = dict(lam=lam, res2=float(np.linalg.norm(r, 2)), res2_over_J=float(np.linalg.norm(r, 2) / (J0 * lam ** 2)),
                   gauss=float(np.linalg.norm(parts[2], 2) / (J0 * lam ** 2)),
                   pair=float(np.linalg.norm(pair, 2) / (J0 * lam ** 2)),
                   inter=float(np.linalg.norm(parts[4], 2) / (J0 * lam ** 2)),
                   fermi_points=pts, info=i2, t=time.time() - t0)
        res[f'h2_{lam}'] = H2
        log.append(rec)
        print(json.dumps(rec), flush=True)
    for lam in a.lams3:
        t0 = time.time()
        H2 = res.get(f'h2_{lam}')
        if H2 is None:
            H2, _ = heff(cb, 2, lam); res[f'h2_{lam}'] = H2
        H3, i3 = heff(cb, 3, lam)
        I = np.eye(N)
        h3 = H3 - np.kron(H2, I) - np.kron(I, H2)
        h3 = h3 - np.trace(h3) / N ** 3 * np.eye(N ** 3)
        r3 = H3 - lam ** 2 * (np.kron(h2p, I) + np.kron(I, h2p))
        r3 = r3 - np.trace(r3) / N ** 3 * np.eye(N ** 3)
        rec = dict(lam3=lam, res3=float(np.linalg.norm(r3, 2)), h3=float(np.linalg.norm(h3, 2)),
                   h3_over_J=float(np.linalg.norm(h3, 2) / (J0 * lam ** 2)), info=i3, t=time.time() - t0)
        res[f'h3_{lam}'] = h3
        log.append(rec)
        print(json.dumps(rec), flush=True)
    # power laws
    L2 = [x for x in log if 'res2' in x]
    lams = np.array([x['lam'] for x in L2]); r2 = np.array([x['res2'] for x in L2])
    sm = lams <= 0.02
    fit = dict(res2_power=float(np.polyfit(np.log(lams[sm]), np.log(r2[sm]), 1)[0]))
    L3 = [x for x in log if 'res3' in x]
    if len(L3) >= 3:
        l3 = np.array([x['lam3'] for x in L3]); r3 = np.array([x['res3'] for x in L3]); hh = np.array([x['h3'] for x in L3])
        sm3 = l3 <= 0.02
        fit.update(res3_power=float(np.polyfit(np.log(l3[sm3]), np.log(r3[sm3]), 1)[0]),
                   h3_power=float(np.polyfit(np.log(l3[sm3]), np.log(hh[sm3]), 1)[0]))
    print(json.dumps(fit), flush=True)
    np.savez(out, lams=np.array(a.lams), J0=J0, h2pred=h2p, **res)
    json.dump(dict(log=log, fit=fit, J0=J0, target=a.target), open(out.replace('.npz', '.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
