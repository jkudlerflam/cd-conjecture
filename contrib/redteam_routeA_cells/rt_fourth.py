"""Fourth-order Schrieffer-Wolff term of one inter-cell bond (two cells), with
P V P = 0:

  H4 = P V R V R V R V P - 1/2 {P V R^2 V P, P V R V P},   R = Q (E0 - H0)^{-1} Q,

applied in the Kronecker eigenbasis of H0 = H_A + H_B (full two-cell legal
space, all excited and vacuum-walker states kept).  We compare its
data-dependent quadratic part (one-cell matrix and hopping) with the
second-order hopping T, and check it against two-cell ED.
Usage: python3 rt_fourth.py n td g
"""
import sys, json
import numpy as np
import scipy.sparse as sp
from rt_bond import P
from rt_cell import cell, code_basis, port_factors, port_V, second_order, one_particle_T, hs_quadratic


def fourth_order(n, p):
    Hd, info = cell(n, p)
    cb = code_basis(n, Hd, info)
    G, E0 = cb['G'], cb['E0']
    w, U = np.linalg.eigh(Hd)
    D = len(w)
    FM, FX = port_factors(n, info)
    terms = port_V(p, FM, FX)
    Vh = None
    for (A, B, c) in terms:
        t = c * sp.kron(sp.csr_matrix(A), sp.csr_matrix(B))
        Vh = t if Vh is None else Vh + t
    herm = abs(Vh - Vh.conj().T).max()
    Vs = Vh.tocsr()
    assert herm < 1e-14, herm
    den = w[:, None] + w[None, :] - 2 * E0               # Kronecker eigenvalues of H0 - 2E0
    N = G.shape[1]
    code_mask = np.abs(den) < 1e-9
    assert code_mask.sum() == N * N, (code_mask.sum(), N * N)
    inv = np.where(code_mask, 0.0, -1.0 / np.where(code_mask, 1.0, den))   # R = Q (E0 - H0)^{-1}

    def toeig(X):   # X: (D*D, k) -> coefficients in U x U basis
        k = X.shape[1]
        Xr = X.reshape(D, D, k)
        return np.einsum('ai,bj,abk->ijk', U.conj(), U.conj(), Xr, optimize=True).reshape(D * D, k)

    def fromeig(Y):
        k = Y.shape[1]
        Yr = Y.reshape(D, D, k)
        return np.einsum('ai,bj,ijk->abk', U, U, Yr, optimize=True).reshape(D * D, k)

    def Rapply(X, power=1):
        Y = toeig(X)
        Y = Y * (inv.reshape(-1)[:, None] ** power)
        return fromeig(Y)
    GG = np.einsum('ia,jb->ijab', G, G).reshape(D * D, N * N)
    X1 = Vs @ GG                     # V P
    RX1 = Rapply(X1)
    H2 = GG.conj().T @ (Vs @ RX1)    # P V R V P
    R2X1 = Rapply(X1, 2)
    H2b = GG.conj().T @ (Vs @ R2X1)  # P V R^2 V P
    Y = Rapply(Vs @ RX1)             # R V R V P
    Z = Rapply(Vs @ Y)               # R V R V R V P
    H4a = GG.conj().T @ (Vs @ Z)
    H4 = H4a - 0.5 * (H2b @ H2 + H2 @ H2b)
    return dict(E0=E0, H2=H2, H4=H4, G=G, Hd=Hd, info=info, terms=terms, Vs=Vs, D=D)


def ed_effective(n, p, lam, F):
    """Two-cell ED in each total-charge sector: canonical (polar) effective
    Hamiltonian on the code."""
    import scipy.linalg as sla
    Hd, G, E0, Vs, D, info = F['Hd'], F['G'], F['E0'], F['Vs'], F['D'], F['info']
    H = sp.kron(sp.csr_matrix(Hd), sp.identity(D)) + sp.kron(sp.identity(D), sp.csr_matrix(Hd)) + lam * Vs
    H = H.tocsr()
    Qc = info['Q']
    Qt = (Qc[:, None] + Qc[None, :]).reshape(-1)
    N = G.shape[1]
    GG = np.einsum('ia,jb->ijab', G, G).reshape(D * D, N * N)
    qcode = np.array([bin(a).count('1') + bin(b).count('1') for a in range(N) for b in range(N)])
    Heff = np.zeros((N * N, N * N), complex)
    gapn = np.inf
    for q in np.unique(qcode):
        si = np.where(Qt == q)[0]
        ci = np.where(qcode == q)[0]
        Hs = H[si][:, si].toarray()
        w, Vv = sla.eigh(Hs, subset_by_index=[0, len(ci)], driver='evr')
        W = GG[si][:, ci].conj().T @ Vv[:, :len(ci)]
        u, s_, vh = np.linalg.svd(W)
        Uc = u @ vh
        Heff[np.ix_(ci, ci)] = Uc @ np.diag(w[:len(ci)]) @ Uc.conj().T
        gapn = min(gapn, w[len(ci)] - w[len(ci) - 1])
    return Heff, gapn


if __name__ == '__main__':
    n, td, g = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    p = P(td=td, g=g)
    F = fourth_order(n, p)
    H2, H4 = F['H2'], F['H4']
    ref = second_order(n, p)['H2']
    out = dict(n=n, td=td, g=g, H2_vs_sparse_solve=float(np.abs(H2 - ref).max()))
    T2, mA2, mB2 = one_particle_T(H2, n)
    T4, mA4, mB4 = one_particle_T(H4, n)
    e2, h2, r2, _ = hs_quadratic(H2, n)
    e4, h4, r4, nc4 = hs_quadratic(H4, n)
    sT = np.abs(T2).max()
    out.update(T2_scale=float(sT), H4_norm=float(np.linalg.norm(H4 - np.trace(H4) / H4.shape[0] * np.eye(H4.shape[0]), 2)),
               H4_scalar=float(np.real(np.trace(H4)) / H4.shape[0]),
               T4_over_T2=float(np.abs(h4[:n, n:]).max() / sT),
               mu4_hs_over_T2=float(np.abs(h4[:n, :n] + h4[n:, n:]).max() / sT),
               Re_mu4_hs_over_T2=float(np.abs((h4[:n, :n] + h4[n:, n:]).real).max() / sT),
               Re_mu4_diag_over_T2=float(np.abs(np.diag(h4[:n, :n] + h4[n:, n:]).real).max() / sT),
               NG4_over_T2=float(r4 / sT), N_commutator=nc4)
    # lambda at which fourth-order one-body terms reach 10 percent of the narrowest window
    from rt_cell import herm_part_analysis
    ev, cr, asym = herm_part_analysis(h2[:n, n:])
    out['Rs_min'] = float(ev[0])
    worst4 = max(np.abs(h4[:n, :n] + h4[n:, n:]).max(), np.abs(h4[:n, n:]).max())
    out['lambda_max_for_window'] = float(np.sqrt(0.1 * 2 * ev[0] * sT / worst4))
    if n <= 2:
        errs = []
        for lam in (0.1, 0.05, 0.025):
            He, gapn = ed_effective(n, p, lam, F)
            pred = 2 * F['E0'] * np.eye(H2.shape[0]) + lam ** 2 * H2 + lam ** 4 * H4
            D_ = He - pred
            D_ -= np.trace(D_) / D_.shape[0] * np.eye(D_.shape[0])
            D2 = He - 2 * F['E0'] * np.eye(H2.shape[0]) - lam ** 2 * H2
            D2 -= np.trace(D2) / D2.shape[0] * np.eye(D2.shape[0])
            errs.append(dict(lam=lam, traceless_err_after_H4=float(np.linalg.norm(D_, 2)),
                             traceless_err_after_H2=float(np.linalg.norm(D2, 2)), next_gap=float(gapn)))
        out['ed'] = errs
    print(json.dumps(out, indent=1, default=float))
    json.dump(out, open('fourth_n%d_td%g_g%g.json' % (n, td, g), 'w'), indent=1, default=float)
