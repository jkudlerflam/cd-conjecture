"""Exact diagonalization of M = 2, 3 open chains of folded-bus cells, compared
with the exact second-order bond.  Usage: python3 ed_check.py n M lam [lam ...]
Writes ed_n{n}_M{M}.json.  numpy/scipy only."""
import sys, json, itertools
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from ycell import params, build_cell, code, second_order, logical, fit_quadratic, FMAP

ED_PARAMS = dict(delta=0.0, Delta=0.5, Dx=2.0, DM=2.0, gp=0.3, gm=0.3, td=0.05,
                 phd={'0': 0.9, '+': 0.0, '-': 0.0})


def chain(cell, M, lam):
    D = len(cell['S'])
    H0 = sp.csr_matrix(cell['H'])
    I = sp.identity(D, format='csr')

    def emb(ops):  # ops: dict position -> matrix
        out = None
        for m in range(M):
            A = sp.csr_matrix(ops[m]) if m in ops else I
            out = A if out is None else sp.kron(out, A, format='csr')
        return out
    H = sum(emb({m: H0}) for m in range(M))
    for m in range(M - 1):
        for (fm, fx, c) in cell['terms']:
            H = H + lam * c * emb({m: cell['FM'][fm], m + 1: cell['FX'][fx]})
    Q = cell['Q']
    Qt = np.zeros(D ** M, int)
    for m in range(M):
        Qt += np.kron(np.kron(np.ones(D ** m, int), Q), np.ones(D ** (M - m - 1), int))
    return H.tocsr(), Qt


def main():
    global ED_PARAMS
    tag = ''
    if sys.argv[1] == 'ir':      # IR working point of logical_chain.py, td = 1e-3
        from logical_chain import IR_PARAMS
        ED_PARAMS = dict(IR_PARAMS, td=1e-3)
        tag = '_ir'
        sys.argv.pop(1)
    n, M = int(sys.argv[1]), int(sys.argv[2])
    lams = [float(x) for x in sys.argv[3:]]
    P = params(**ED_PARAMS)
    cell = build_cell(n, P)
    E0, G, gap = code(cell)
    H2, kpoles = second_order(cell, G, E0)
    N = 2 ** n
    Gt = G
    for m in range(M - 1):
        Gt = np.kron(Gt, G)
    Qcode = np.array([sum(int(b) for b in format(i, '0%db' % (n * M))) for i in range(N ** M)])
    pred2 = np.zeros((N ** M, N ** M), complex)
    for m in range(M - 1):
        pred2 += np.kron(np.kron(np.eye(N ** m), H2), np.eye(N ** (M - m - 2)))
    rec = dict(n=n, M=M, params={k: v for k, v in ED_PARAMS.items()}, E0=E0, cell_gap=gap,
               cell_dim=len(cell['S']), chain_dim=len(cell['S']) ** M, poles=kpoles, runs=[])
    Lg = logical(H2, n)
    for lam in lams:
        H, Qt = chain(cell, M, lam)
        Heff = np.zeros((N ** M, N ** M), complex)
        maxres, nextgap = 0.0, np.inf
        for q in range(n * M + 1):
            ci = np.where(Qcode == q)[0]
            si = np.where(Qt == q)[0]
            mq = len(ci)
            Hs = H[si][:, si]
            kk = mq + 2
            if Hs.shape[0] <= 2500:
                w, V = np.linalg.eigh(Hs.toarray())
                w, V = w[:kk], V[:, :kk]
            else:   # shift-invert just below the unperturbed code energy
                w, V = spla.eigsh(Hs.tocsc(), k=kk, sigma=M * E0 - 0.05, which='LM', tol=1e-14)
                o = np.argsort(w); w, V = w[o], V[:, o]
            res = np.linalg.norm(Hs @ V - V * w, axis=0).max()
            maxres = max(maxres, res)
            nextgap = min(nextgap, w[mq] - w[mq - 1])
            Wm = Gt[si][:, ci].conj().T @ V[:, :mq]
            u, s, vh = np.linalg.svd(Wm)
            U = u @ vh
            Heff[np.ix_(ci, ci)] = U @ np.diag(w[:mq]) @ U.conj().T
        pred = M * E0 * np.eye(N ** M) + lam ** 2 * pred2
        err = np.linalg.norm(Heff - pred, 2)
        D_ = Heff - pred
        D_ -= np.trace(D_) / D_.shape[0] * np.eye(D_.shape[0])
        err_tr = np.linalg.norm(D_, 2)
        run = dict(lam=lam, err=err, err_tr=err_tr, err_over_lam4=err / lam ** 4,
                   max_residual=maxres, min_gap_to_next=nextgap)
        if M == 2:
            Le = logical((Heff - M * E0 * np.eye(N ** M)) / lam ** 2, n)
            run['T_pred'] = [[str(x) for x in r] for r in Lg['T']]
            run['T_ed'] = [[str(x) for x in r] for r in Le['T']]
            run['T_rel_err'] = float(np.abs(Le['T'] - Lg['T']).max() / np.abs(Lg['T']).max())
        print(json.dumps(run), flush=True)
        rec['runs'].append(run)
    with open('ed_n%d_M%d%s.json' % (n, M, tag), 'w') as f:
        json.dump(rec, f, indent=1, default=str)


if __name__ == '__main__':
    main()
