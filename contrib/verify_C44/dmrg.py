"""Dense two-site finite DMRG (open chain) and one-site VUMPS for a
translation-invariant nearest-neighbour chain, numpy/scipy only.

Conventions: MPO tensors W[wl, wr, s_out, s_in]; left boundary picks row D-1,
right boundary picks column 0; H = sum_i h1_i + sum_i h2_(i,i+1).
"""
import numpy as np
import scipy.linalg as sla
import scipy.sparse.linalg as spla


def mpo_from_bond(h2, h1, d, tol=1e-14):
    """Operator-Schmidt decomposition of the bond into a lower-triangular MPO."""
    M = h2.reshape(d, d, d, d).transpose(0, 2, 1, 3).reshape(d * d, d * d)   # (s1' s1),(s2' s2)
    U, S, Vh = np.linalg.svd(M)
    r = int((S > tol * S[0]).sum())
    A = [(U[:, a] * np.sqrt(S[a])).reshape(d, d) for a in range(r)]
    B = [(Vh[a] * np.sqrt(S[a])).reshape(d, d) for a in range(r)]
    D = r + 2
    W = np.zeros((D, D, d, d), complex)
    W[0, 0] = np.eye(d)
    W[D - 1, D - 1] = np.eye(d)
    W[D - 1, 0] = h1
    for a in range(r):
        W[D - 1, a + 1] = A[a]
        W[a + 1, 0] = B[a]
    # check
    rec = sum(np.kron(A[a], B[a]) for a in range(r))
    err = np.abs(rec - h2).max()
    return W, err


def lanczos(matvec, v0, k=40, tol=1e-12, abs_tol=None):
    """Lowest eigenpair by Lanczos with full reorthogonalization."""
    v = v0 / np.linalg.norm(v0)
    V = [v]
    alph, beta = [], []
    w = matvec(v)
    e_old = None
    for j in range(k):
        a = np.vdot(V[j], w).real
        alph.append(a)
        w = w - a * V[j] - (beta[-1] * V[j - 1] if j > 0 else 0)
        for u in V:
            w = w - np.vdot(u, w) * u
        b = np.linalg.norm(w)
        T = np.diag(alph) + np.diag(beta, 1) + np.diag(beta, -1)
        ev, evec = np.linalg.eigh(T)
        thr = abs_tol if abs_tol is not None else tol * max(1, abs(ev[0]))
        if e_old is not None and abs(ev[0] - e_old) < thr:
            break
        e_old = ev[0]
        if b < 1e-14:
            break
        beta.append(b)
        V.append(w / b)
        w = matvec(V[-1])
    T = np.diag(alph) + np.diag(beta[:len(alph) - 1], 1) + np.diag(beta[:len(alph) - 1], -1)
    ev, evec = np.linalg.eigh(T)
    x = sum(evec[i, 0] * V[i] for i in range(len(alph)))
    return ev[0], x / np.linalg.norm(x)


def apply_two(L, W1, W2, R, th):
    # L (a,w,a'), th (a',s1,s2,b'), W (w,w',s',s), R (b,w,b')
    t = np.tensordot(L, th, axes=(2, 0))                  # a w s1 s2 b'
    t = np.tensordot(t, W1, axes=([1, 2], [0, 3]))        # a s2 b' w2 s1'
    t = np.tensordot(t, W2, axes=([3, 1], [0, 3]))        # a b' s1' w3 s2'
    t = np.tensordot(t, R, axes=([3, 1], [1, 2]))         # a s1' s2' b
    return t


def grow_left(L, A, W):
    t = np.tensordot(L, A, axes=(2, 0))                   # a w s b'
    t = np.tensordot(t, W, axes=([1, 2], [0, 3]))         # a b' w' s'
    t = np.tensordot(A.conj(), t, axes=([0, 1], [0, 3]))  # b b' w'
    return t.transpose(0, 2, 1)


def grow_right(R, B, W):
    t = np.tensordot(B, R, axes=(2, 2))                   # a s b w   (R index b' contracted with B's right)
    t = np.tensordot(t, W, axes=([1, 3], [3, 1]))         # a b w s'
    t = np.tensordot(t, B.conj(), axes=([3, 1], [1, 2]))  # a_in w a_out
    return t.transpose(2, 1, 0)                          # (out, w, in)


def svd_trunc(th, chi, cut=1e-14, gap_window=0.0):
    """Truncated SVD. With gap_window>0 the cut is placed at the largest
    relative Schmidt gap among keep in [(1-gap_window) chi, chi], so that a
    near-degenerate multiplet is never split."""
    a, d1, d2, b = th.shape
    U, S, Vh = np.linalg.svd(th.reshape(a * d1, d2 * b), full_matrices=False)
    keep = min(chi, int((S > cut * S[0]).sum()))
    if gap_window > 0 and keep == chi and len(S) > chi:
        lo = max(1, int(np.floor((1 - gap_window) * chi)))
        ratios = S[lo - 1:chi] / np.maximum(S[lo:chi + 1], 1e-300)
        keep = lo + int(np.argmax(ratios))
    trunc = float((S[keep:] ** 2).sum() / (S ** 2).sum())
    S = S[:keep]
    S = S / np.linalg.norm(S)
    return U[:, :keep].reshape(a, d1, keep), S, Vh[:keep].reshape(keep, d2, b), trunc


def entropy(S):
    p = S ** 2
    p = p[p > 1e-30]
    return float(-(p * np.log(p)).sum())


# ------------------------------------------------------------------ finite DMRG
def finite_dmrg(W, L, chi_schedule, d, ops=None, seed=1, init=None, lanczos_k=40, verbose=False):
    D = W.shape[0]
    rng = np.random.default_rng(seed)
    if init is None:
        Ms = [rng.standard_normal((1, d, 1)) + 0j for _ in range(L)]
    else:
        Ms = [init(i).reshape(1, d, 1).astype(complex) for i in range(L)]
    # right-canonicalize
    for i in range(L - 1, 0, -1):
        a, dd, b = Ms[i].shape
        Q, Rm = np.linalg.qr(Ms[i].reshape(a, dd * b).T)
        Ms[i] = Q.T.reshape(-1, dd, b)
        Ms[i - 1] = np.tensordot(Ms[i - 1], Rm.T, axes=(2, 0))
    Ms[0] /= np.linalg.norm(Ms[0])
    Lb = np.zeros((1, D, 1), complex); Lb[0, D - 1, 0] = 1
    Rb = np.zeros((1, D, 1), complex); Rb[0, 0, 0] = 1
    Lenv = [None] * (L + 1)
    Renv = [None] * (L + 1)
    Lenv[0] = Lb
    Renv[L] = Rb
    for i in range(L - 1, -1, -1):
        Renv[i] = grow_right(Renv[i + 1], Ms[i], W)
    hist = []
    S_prof = np.zeros(L - 1)
    for sweep, chi in enumerate(chi_schedule):
        maxtr = 0.0
        # left-to-right
        for i in range(L - 1):
            th = np.tensordot(Ms[i], Ms[i + 1], axes=(2, 0))
            shp = th.shape
            mv = lambda v: apply_two(Lenv[i], W, W, Renv[i + 2], v.reshape(shp)).ravel()
            E, x = lanczos(mv, th.ravel(), k=lanczos_k)
            A, S, B, tr = svd_trunc(x.reshape(shp), chi)
            maxtr = max(maxtr, tr)
            Ms[i] = A
            Ms[i + 1] = np.tensordot(np.diag(S), B, axes=(1, 0))
            Lenv[i + 1] = grow_left(Lenv[i], A, W)
            S_prof[i] = entropy(S)
        # right-to-left
        for i in range(L - 2, -1, -1):
            th = np.tensordot(Ms[i], Ms[i + 1], axes=(2, 0))
            shp = th.shape
            mv = lambda v: apply_two(Lenv[i], W, W, Renv[i + 2], v.reshape(shp)).ravel()
            E, x = lanczos(mv, th.ravel(), k=lanczos_k)
            A, S, B, tr = svd_trunc(x.reshape(shp), chi)
            maxtr = max(maxtr, tr)
            Ms[i + 1] = B
            Ms[i] = np.tensordot(A, np.diag(S), axes=(2, 0))
            Renv[i + 1] = grow_right(Renv[i + 2], B, W)
            S_prof[i] = entropy(S)
        hist.append(dict(sweep=sweep, chi=chi, E=float(E), max_trunc=maxtr))
        if verbose:
            print(hist[-1], flush=True)
    # measurements (MPS is left-normalized up to site 0 center after the last right-to-left sweep)
    meas = {}
    if ops:
        for name, o in ops.items():
            vals = np.zeros(L)
            # Ms[0] carries the norm; bring the center along with QR to the right
            Mc = [m.copy() for m in Ms]
            for i in range(L):
                vals[i] = np.einsum('asb,ts,atb->', Mc[i].conj(), o, Mc[i]).real
                if i < L - 1:
                    a, dd, b = Mc[i].shape
                    Q, Rm = np.linalg.qr(Mc[i].reshape(a * dd, b))
                    Mc[i] = Q.reshape(a, dd, -1)
                    Mc[i + 1] = np.tensordot(Rm, Mc[i + 1], axes=(1, 0))
            meas[name] = vals
    return dict(E=float(E), S=S_prof.copy(), hist=hist, meas=meas)


# ------------------------------------------------------------------ VUMPS (one-site unit cell)
def _polar_left(M):
    U, S, Vh = np.linalg.svd(M, full_matrices=False)
    return U @ Vh


def _tl(GL, A, W):
    """left transfer of a (bra, w, ket) environment through one site."""
    return grow_left(GL, A, W)


def vumps(W, chi, d, iters=400, tol=1e-9, seed=7, ops=None, verbose=False, A0=None):
    D = W.shape[0]
    rng = np.random.default_rng(seed)
    if A0 is None:
        A = rng.standard_normal((chi, d, chi)) + 1j * rng.standard_normal((chi, d, chi))
        AL = _polar_left(A.reshape(chi * d, chi)).reshape(chi, d, chi)
        AR = _polar_left(A.reshape(chi, d * chi).T).T.reshape(chi, d, chi)
        C = np.diag(rng.random(chi) + 0.1).astype(complex)
    else:
        AL, AR, C = A0                      # mixed-gauge initial state
    C = C / np.linalg.norm(C)
    AC = np.tensordot(AL, C, axes=(2, 0))
    I = np.eye(chi)
    hist = []

    def T_L(X, Wop):          # X (bra, ket) -> sum A*^(s') X A^s Wop[s',s]
        Y = np.tensordot(X, AL, axes=(1, 0))                       # bra s b
        Y = np.tensordot(Y, Wop, axes=(1, 1))                      # bra b s'
        return np.tensordot(AL.conj(), Y, axes=([0, 1], [0, 2]))   # b* b

    def T_R(X, Wop):          # X (bra, ket) on the right
        Y = np.tensordot(AR, X, axes=(2, 1))                       # a s bra
        Y = np.tensordot(Y, Wop, axes=(1, 1))                      # a bra s'
        return np.tensordot(AR.conj(), Y, axes=([1, 2], [2, 1]))    # (bra a*, ket a)

    def solve_left(Y, Rfix):
        e = np.trace(Y @ Rfix.T) if False else np.sum(Y * Rfix.T)
        rhs = (Y - e * I).ravel()

        def mv(x):
            X = x.reshape(chi, chi)
            return (X - T_L(X, np.eye(d)) + np.sum(X * Rfix.T) * I).ravel()
        op = spla.LinearOperator((chi * chi, chi * chi), matvec=mv, dtype=complex)
        x, info = spla.gmres(op, rhs, rtol=1e-12, atol=0, restart=60, maxiter=200)
        return x.reshape(chi, chi), e

    def solve_right(Y, Lfix):
        e = np.sum(Y * Lfix)
        rhs = (Y - e * I).ravel()

        def mv(x):
            X = x.reshape(chi, chi)
            return (X - T_R(X, np.eye(d)) + np.sum(X * Lfix) * I).ravel()
        op = spla.LinearOperator((chi * chi, chi * chi), matvec=mv, dtype=complex)
        x, info = spla.gmres(op, rhs, rtol=1e-12, atol=0, restart=60, maxiter=200)
        return x.reshape(chi, chi), e

    for it in range(iters):
        # environments; left fixed point of T_L is I, right fixed point is R = C C^dag
        Rf = C @ C.conj().T
        Lf = C.conj().T @ C
        GL = [None] * D
        GL[D - 1] = I.astype(complex)
        for w in range(D - 2, 0, -1):
            GL[w] = sum(T_L(GL[v], W[v, w]) for v in range(w + 1, D) if np.any(W[v, w]))
        Y = sum(T_L(GL[v], W[v, 0]) for v in range(1, D) if np.any(W[v, 0]))
        GL[0], eL = solve_left(Y, Rf)
        GR = [None] * D
        GR[0] = I.astype(complex)
        for w in range(1, D - 1):
            GR[w] = sum(T_R(GR[v], W[w, v]) for v in range(0, w) if np.any(W[w, v]))
        Y = sum(T_R(GR[v], W[D - 1, v]) for v in range(0, D - 1) if np.any(W[D - 1, v]))
        GR[D - 1], eR = solve_right(Y, Lf)
        Lenv = np.stack(GL, axis=1)        # (bra, w, ket)
        Renv = np.stack(GR, axis=1)        # (bra, w, ket)

        def hac(v):
            X = v.reshape(chi, d, chi)
            t = np.tensordot(Lenv, X, axes=(2, 0))                  # a w s b'
            t = np.tensordot(t, W, axes=([1, 2], [0, 3]))            # a b' w' s'
            t = np.tensordot(t, Renv, axes=([2, 1], [1, 2]))         # a s' b
            return t.ravel()

        def hc(v):
            X = v.reshape(chi, chi)
            t = np.tensordot(Lenv, X, axes=(2, 0))                  # a w b'
            t = np.tensordot(t, Renv, axes=([1, 2], [1, 2]))         # a b
            return t.ravel()
        eac, xac = lanczos(hac, AC.ravel(), k=40, abs_tol=1e-14)
        ec, xc = lanczos(hc, C.ravel(), k=40, abs_tol=1e-14)
        AC = xac.reshape(chi, d, chi)
        C = xc.reshape(chi, chi)
        UL = _polar_left(AC.reshape(chi * d, chi))
        UC = _polar_left(C)
        AL = (UL @ UC.conj().T).reshape(chi, d, chi)
        URt = _polar_left(AC.reshape(chi, d * chi).T)          # (d chi, chi)
        UCt = _polar_left(C.T)
        AR = (URt @ UCt.conj().T).T.reshape(chi, d, chi)
        err = max(np.linalg.norm(AC - np.tensordot(AL, C, axes=(2, 0))),
                  np.linalg.norm(AC - np.tensordot(C, AR, axes=(1, 0))))
        e_site = float(np.real(eL))
        hist.append((it, e_site, err))
        if verbose and it % 10 == 0:
            print(it, e_site, err, flush=True)
        if err < tol:
            break
    S = np.linalg.svd(C, compute_uv=False)
    S = S / np.linalg.norm(S)

    def tmv(xv):
        X = xv.reshape(chi, chi)
        return T_L(X, np.eye(d)).ravel()
    T = spla.LinearOperator((chi * chi, chi * chi), matvec=tmv, dtype=complex)
    vals = spla.eigs(T, k=6, which='LM', return_eigenvectors=False, tol=1e-12, maxiter=50000)
    vals = np.sort(np.abs(vals))[::-1]
    out = dict(iters=it + 1, e_site=e_site, err=float(err), S=entropy(S), xi=float(-1 / np.log(vals[1] / vals[0])),
               tm_eigs=[float(v) for v in vals], schmidt_min=float(S[-1]))
    if ops:
        for name, o in ops.items():
            out[name] = float(np.einsum('asb,ts,atb->', AC.conj(), o, AC).real / np.vdot(AC, AC).real)
    out['AL'], out['AR'], out['C'] = AL, AR, C
    out['hist'] = hist
    return out
