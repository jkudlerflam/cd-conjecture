"""(Copied from work_C144, C155 adds the edge_field option.)
Compact two-site DMRG (NumPy/SciPy) for an open chain H = sum_{j=1}^{L-1} h_{j,j+1}, h a 4x4 bond on qubits.

Conventions: MPS tensors M[j] of shape (Dl, d, Dr); MPO W[j] of shape (wl, wr, d_out, d_in);
environments Lenv[j] of shape (D_bra, w, D_ket). Lanczos via scipy eigsh on the two-site effective H.
Entanglement S(l) for the cut after site l (l = 1..L-1) is recorded on the last left-to-right sweep.
"""
import numpy as np
import scipy.sparse.linalg as spla


def bond_mpo(h, L):
    d = int(round(np.sqrt(h.shape[0])))
    h4 = h.reshape(d, d, d, d)  # (s1, s2, t1, t2)
    hm = h4.transpose(0, 2, 1, 3).reshape(d * d, d * d)  # (s1 t1), (s2 t2)
    U, S, Vh = np.linalg.svd(hm)
    keep = S > 1e-12 * S[0]
    U = U[:, keep] * np.sqrt(S[keep]); Vh = (Vh[keep].T * np.sqrt(S[keep])).T
    nk = U.shape[1]
    w = nk + 2
    dt = np.result_type(h.dtype, float)
    W = np.zeros((w, w, d, d), dtype=dt)
    I = np.eye(d)
    W[0, 0] = I; W[w - 1, w - 1] = I
    for k in range(nk):
        W[0, 1 + k] = U[:, k].reshape(d, d)
        W[1 + k, w - 1] = Vh[k].reshape(d, d)
    Ws = []
    for j in range(L):
        if j == 0:
            Ws.append(W[0:1].copy())
        elif j == L - 1:
            Ws.append(W[:, w - 1:w].copy())
        else:
            Ws.append(W.copy())
    return Ws


def _lupdate(Le, M, W):
    # Lnew[c,v,d] = sum Le[a,w,b] conj(M[a,s,c]) W[w,v,s,t] M[b,t,d]
    T = np.tensordot(Le, M, axes=([2], [0]))            # a w t d
    T = np.tensordot(T, W, axes=([1, 2], [0, 3]))        # a d v s
    T = np.tensordot(np.conj(M), T, axes=([0, 1], [0, 3]))  # c d v
    return T.transpose(0, 2, 1)


def _rupdate(Re, M, W):
    # Rnew[a,w,b] = sum conj(M[a,s,c]) W[w,v,s,t] M[b,t,d] Re[c,v,d]
    T = np.tensordot(M, Re, axes=([2], [2]))             # b t c v
    T = np.tensordot(T, W, axes=([3, 1], [1, 3]))        # b c w s
    T = np.tensordot(np.conj(M), T, axes=([1, 2], [3, 1]))  # a b w
    return T.transpose(0, 2, 1)


def _heff(Le, W1, W2, Re, shape, dtype):
    def mv(x):
        th = x.reshape(shape)                              # b t1 t2 d
        T = np.tensordot(Le, th, axes=([2], [0]))          # a w t1 t2 d
        T = np.tensordot(T, W1, axes=([1, 2], [0, 3]))     # a t2 d v s1
        T = np.tensordot(T, W2, axes=([3, 1], [0, 3]))     # a d s1 u s2
        T = np.tensordot(T, Re, axes=([1, 3], [2, 1]))     # a s1 s2 c
        return T.reshape(-1)
    n = int(np.prod(shape))
    return spla.LinearOperator((n, n), matvec=mv, dtype=dtype)


def entropy(s):
    p = s ** 2; p = p[p > 1e-16]; p = p / p.sum()
    return float(-np.sum(p * np.log(p)))


def dmrg(h, L, chis=(16, 32, 64, 64), cutoff=1e-10, seed=0, verbose=False, lanczos_tol=1e-9, edge_field=None):
    d = int(round(np.sqrt(h.shape[0])))
    dt = np.result_type(h.dtype, float)
    Ws = bond_mpo(h, L)
    if edge_field is not None:  # C155: extra on-site term at the two end sites
        Ws[0][0, -1] = Ws[0][0, -1] + edge_field
        Ws[L - 1][0, 0] = Ws[L - 1][0, 0] + edge_field
    rng = np.random.default_rng(seed)
    # random right-canonical start, bond dim 4
    M = []
    D = [1]
    for j in range(L):
        D.append(min(4, d ** (j + 1), d ** (L - j - 1)))
    for j in range(L):
        M.append(rng.standard_normal((D[j], d, D[j + 1])).astype(dt))
    # right-canonicalize
    for j in range(L - 1, 0, -1):
        A = M[j].reshape(M[j].shape[0], -1)
        Q, R = np.linalg.qr(A.T.conj())
        M[j] = Q.T.conj().reshape(-1, d, M[j].shape[2])
        M[j - 1] = np.tensordot(M[j - 1], R.T.conj(), axes=([2], [0]))
    M[0] /= np.linalg.norm(M[0])
    Lenv = [None] * (L + 1); Renv = [None] * (L + 1)
    Lenv[0] = np.ones((1, 1, 1), dtype=dt)
    Renv[L] = np.ones((1, 1, 1), dtype=dt)
    for j in range(L - 1, 0, -1):
        Renv[j] = _rupdate(Renv[j + 1], M[j], Ws[j])
    E = None; hist = []
    S_cut = np.zeros(L - 1); trunc_max = 0.0
    for sweep, chi in enumerate(chis):
        trunc_max = 0.0
        # left to right
        for j in range(L - 1):
            th = np.tensordot(M[j], M[j + 1], axes=([2], [0]))
            shape = th.shape
            op = _heff(Lenv[j], Ws[j], Ws[j + 1], Renv[j + 2], shape, dt)
            if th.size <= 16:
                Hd = np.array([op.matvec(e) for e in np.eye(th.size, dtype=dt)]).T
                w, v = np.linalg.eigh(Hd); E = w[0]; th = v[:, 0]
            else:
                w, v = spla.eigsh(op, k=1, which="SA", v0=th.reshape(-1), tol=lanczos_tol, ncv=min(20, th.size - 1))
                E = w[0]; th = v[:, 0]
            th = th.reshape(shape[0] * d, d * shape[3])
            U, s, Vh = np.linalg.svd(th, full_matrices=False)
            keep = min(chi, max(1, int(np.sum(s > cutoff * s[0]))))
            tr = float(np.sum(s[keep:] ** 2) / np.sum(s ** 2)); trunc_max = max(trunc_max, tr)
            s = s[:keep] / np.linalg.norm(s[:keep])
            S_cut[j] = entropy(s)
            M[j] = U[:, :keep].reshape(shape[0], d, keep)
            M[j + 1] = (s[:, None] * Vh[:keep]).reshape(keep, d, shape[3])
            Lenv[j + 1] = _lupdate(Lenv[j], M[j], Ws[j])
        # right to left
        for j in range(L - 2, -1, -1):
            th = np.tensordot(M[j], M[j + 1], axes=([2], [0]))
            shape = th.shape
            op = _heff(Lenv[j], Ws[j], Ws[j + 1], Renv[j + 2], shape, dt)
            if th.size <= 16:
                Hd = np.array([op.matvec(e) for e in np.eye(th.size, dtype=dt)]).T
                w, v = np.linalg.eigh(Hd); E = w[0]; th = v[:, 0]
            else:
                w, v = spla.eigsh(op, k=1, which="SA", v0=th.reshape(-1), tol=lanczos_tol, ncv=min(20, th.size - 1))
                E = w[0]; th = v[:, 0]
            th = th.reshape(shape[0] * d, d * shape[3])
            U, s, Vh = np.linalg.svd(th, full_matrices=False)
            keep = min(chi, max(1, int(np.sum(s > cutoff * s[0]))))
            tr = float(np.sum(s[keep:] ** 2) / np.sum(s ** 2)); trunc_max = max(trunc_max, tr)
            s = s[:keep] / np.linalg.norm(s[:keep])
            M[j + 1] = Vh[:keep].reshape(keep, d, shape[3])
            M[j] = (U[:, :keep] * s[None, :]).reshape(shape[0], d, keep)
            Renv[j + 1] = _rupdate(Renv[j + 2], M[j + 1], Ws[j + 1])
        hist.append((chi, float(np.real(E)), trunc_max))
        if verbose:
            print(f"  sweep {sweep} chi={chi} E={np.real(E):.12f} trunc={trunc_max:.2e}", flush=True)
    # final left-to-right pass to record entanglement and leave M left-canonical up to the last site
    for j in range(L - 1):
        th = np.tensordot(M[j], M[j + 1], axes=([2], [0]))
        shape = th.shape
        th = th.reshape(shape[0] * d, d * shape[3])
        U, s, Vh = np.linalg.svd(th, full_matrices=False)
        keep = min(chis[-1], max(1, int(np.sum(s > cutoff * s[0]))))
        s = s[:keep] / np.linalg.norm(s[:keep])
        S_cut[j] = entropy(s)
        M[j] = U[:, :keep].reshape(shape[0], d, keep)
        M[j + 1] = (s[:, None] * Vh[:keep]).reshape(keep, d, shape[3])
    return dict(E=float(np.real(E)), S=S_cut, hist=hist, M=M)


def expect_bonds(M, ops2):
    """<O1_j O2_{j+1}> for every bond j, M left-canonical with the orthogonality center at the last site.
    ops2: list of (O1, O2). Returns array (n_ops, L-1)."""
    L = len(M)
    # right environments of identity
    R = [None] * (L + 1)
    R[L] = np.ones((1, 1))
    for j in range(L - 1, -1, -1):
        A = M[j]
        R[j] = np.einsum("asc,cd,bsd->ab", np.conj(A), R[j + 1], A)
    out = np.zeros((len(ops2), L - 1), dtype=complex)
    Lf = np.ones((1, 1))
    for j in range(L - 1):
        A = M[j]; B = M[j + 1]
        for k, (O1, O2) in enumerate(ops2):
            T = np.einsum("ab,asc,st,btd->cd", Lf, np.conj(A), O1, A)
            T = np.einsum("ab,asc,st,btd->cd", T, np.conj(B), O2, B)
            out[k, j] = np.einsum("cd,cd->", T, R[j + 2])
        Lf = np.einsum("ab,asc,bsd->cd", Lf, np.conj(A), A)
    return out
