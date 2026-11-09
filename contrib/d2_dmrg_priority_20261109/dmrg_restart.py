# C183 solver copy: only optional initial MPS and sweep diagnostics added, 2026-11-09.
"""C177 (Claude Opus 5.5 agent, 2026-11-09): dmrg.py of C166 with two changes, everything else identical.

1. Bond-dependent MPO: `h` may be a list of 4x4 bonds (h[j] on bond (j, j+1)), so that a staggered real gauge
   can be used (the parity-class bond at J_x + J_y = 0 becomes real after the site rotation
   exp(-i theta_j Z_j/2), theta_j = (-1)^j pi/4; see PROOF.md 1.2).  dtype is real when all bonds are real.
2. Local eigensolver: Davidson with the diagonal of the two-site effective Hamiltonian as preconditioner
   (instead of ARPACK Lanczos).  At the Fibonacci cusp the two-site problem has a spectral spread of order V and
   low-energy gaps of order t2 = Delta^2/(2V), and ARPACK needed about 500 matvecs per solve at V/Delta = 8.
   The diagonal carries the large V n_j n_{j+1} terms exactly, which removes most of that stiffness.
   Convergence criterion: residual norm ||H x - theta x|| < res_tol in the last two sweeps, and a loose
   loose_tol (at most 60 matvecs) in the earlier sweeps, where chi still grows.

Validated against sparse ED in ed_check2.out.
"""
import numpy as np
import scipy.sparse.linalg as spla
from dmrg import _lupdate, _rupdate, entropy


def bond_mpo_list(hs):
    d = int(round(np.sqrt(hs[0].shape[0])))  # C183: local dimension read from the bond (d = 4 control)
    L = len(hs) + 1
    mats = []
    nk_max = 0
    for h in hs:
        h4 = h.reshape(d, d, d, d)
        hm = h4.transpose(0, 2, 1, 3).reshape(d * d, d * d)
        U, S, Vh = np.linalg.svd(hm)
        keep = S > 1e-12 * S[0]
        U = U[:, keep] * np.sqrt(S[keep]); Vh = (Vh[keep].T * np.sqrt(S[keep])).T
        mats.append((U, Vh)); nk_max = max(nk_max, U.shape[1])
    dt = np.result_type(*[h.dtype for h in hs], float)
    w = nk_max + 2
    I = np.eye(d)
    Ws = []
    for j in range(L):
        W = np.zeros((w, w, d, d), dtype=dt)
        W[0, 0] = I; W[w - 1, w - 1] = I
        if j < L - 1:  # opens bond (j, j+1)
            U, _ = mats[j]
            for k in range(U.shape[1]):
                W[0, 1 + k] = U[:, k].reshape(d, d)
        if j > 0:      # closes bond (j-1, j)
            _, Vh = mats[j - 1]
            for k in range(Vh.shape[0]):
                W[1 + k, w - 1] = Vh[k].reshape(d, d)
        if j == 0:
            W = W[0:1].copy()
        elif j == L - 1:
            W = W[:, w - 1:w].copy()
        Ws.append(W)
    return Ws, dt


def _heff_apply(Le, W1, W2, Re, shape):
    def mv(x):
        th = x.reshape(shape)
        T = np.tensordot(Le, th, axes=([2], [0]))
        T = np.tensordot(T, W1, axes=([1, 2], [0, 3]))
        T = np.tensordot(T, W2, axes=([3, 1], [0, 3]))
        T = np.tensordot(T, Re, axes=([1, 3], [2, 1]))
        return T.reshape(-1)
    return mv


def _heff_diag(Le, W1, W2, Re):
    l = np.einsum('awa->aw', Le)                  # (a, w)
    r = np.einsum('cvc->cv', Re)                  # (c, v)   Re[bra, w, ket]
    w1 = np.einsum('wvss->wvs', W1)
    w2 = np.einsum('vuss->vus', W2)
    t = np.einsum('aw,wvs->avs', l, w1)
    t = np.einsum('avs,vut->astu', t, w2)
    return np.einsum('astu,cu->astc', t, r).reshape(-1)


def davidson(mv, x0, diag, res_tol=1e-8, maxit=400, mmax=24, keep=3):
    n = x0.size
    dt = np.result_type(x0.dtype, diag.dtype)
    x = x0 / np.linalg.norm(x0)
    V = np.zeros((mmax, n), dtype=dt); AV = np.zeros((mmax, n), dtype=dt)
    V[0] = x; AV[0] = mv(x); m = 1
    diag = diag.real
    nmv = 1
    for it in range(maxit):
        Hs = V[:m].conj() @ AV[:m].T
        Hs = (Hs + Hs.conj().T) / 2
        th, Y = np.linalg.eigh(Hs)
        y = Y[:, 0]; theta = th[0]
        u = y @ V[:m]; Au = y @ AV[:m]
        r = Au - theta * u
        rn = np.linalg.norm(r)
        if rn < res_tol:
            break
        den = diag - theta
        small = np.abs(den) < 1e-3
        den[small] = np.where(den[small] >= 0, 1e-3, -1e-3)
        t = r / den
        if m == mmax:  # restart with the lowest `keep` Ritz vectors
            Vn = Y[:, :keep].T @ V[:m]; AVn = Y[:, :keep].T @ AV[:m]
            V[:keep] = Vn; AV[:keep] = AVn; m = keep
        for _ in range(2):
            t = t - V[:m].T @ (V[:m].conj() @ t)
        tn = np.linalg.norm(t)
        if tn < 1e-14:
            t = r - V[:m].T @ (V[:m].conj() @ r); tn = np.linalg.norm(t)
            if tn < 1e-14:
                break
        V[m] = t / tn; AV[m] = mv(V[m]); m += 1; nmv += 1
    return theta, u / np.linalg.norm(u), nmv, rn


def add_parity_string(Ws, lam):
    """Add lam * prod_j Z_j to the MPO (a string channel).  P = prod Z = (-1)^N commutes with the parity-class H,
    so lam > 0 selects the odd sector and lam < 0 the even one, with E = E' + |lam| in the selected sector."""
    L = len(Ws)
    Z = np.diag([1.0, -1.0])
    out = []
    for j, W in enumerate(Ws):
        a, b = W.shape[:2]
        # channel map: old 0..w-2 keep, old w-1 ('done') -> w, new string channel at w-1
        w = b if j == 0 else a
        def m(i, n):
            return i if i < n - 1 else n
        if j == 0:
            Wn = np.zeros((1, b + 1, 2, 2), dtype=W.dtype)
            for q in range(b):
                Wn[0, m(q, b)] = W[0, q]
            Wn[0, b - 1] = Z if L > 1 else lam * Z
        elif j == L - 1:
            Wn = np.zeros((a + 1, 1, 2, 2), dtype=W.dtype)
            for p in range(a):
                Wn[m(p, a), 0] = W[p, 0]
            Wn[a - 1, 0] = lam * Z
        else:
            Wn = np.zeros((a + 1, b + 1, 2, 2), dtype=W.dtype)
            for p in range(a):
                for q in range(b):
                    Wn[m(p, a), m(q, b)] = W[p, q]
            Wn[a - 1, b - 1] = Z
        out.append(Wn)
    return out


def dmrg(h, L, chis=(16, 32, 64, 64), cutoff=1e-10, seed=0, verbose=False, res_tol=1e-8, edge_field=None,
         edge_fields=None, loose_tol=1e-6, parity_lam=None, initial_M=None):
    """h: one 4x4 bond (uniform) or a list of L-1 bonds.  edge_field adds the same operator at both ends;
    edge_fields=(left, right) adds different ones."""
    hs = list(h) if isinstance(h, (list, tuple)) else [h] * (L - 1)
    d = int(round(np.sqrt(hs[0].shape[0])))  # C183: was d = 2
    Ws, dt = bond_mpo_list(hs)
    if edge_field is not None:
        edge_fields = (edge_field, edge_field)
    if edge_fields is not None:
        Ws[0][0, -1] = Ws[0][0, -1] + edge_fields[0]
        Ws[L - 1][0, 0] = Ws[L - 1][0, 0] + edge_fields[1]
    if parity_lam is not None:
        Ws = add_parity_string(Ws, parity_lam)
    Ws = [W.astype(np.result_type(W.dtype, dt)) for W in Ws]
    dt = np.result_type(*[W.dtype for W in Ws])
    rng = np.random.default_rng(seed)
    M = []
    D = [1]
    for j in range(L):
        D.append(min(4, d ** (j + 1), d ** (L - j - 1)))
    for j in range(L):
        M.append(rng.standard_normal((D[j], d, D[j + 1])).astype(dt))
    if initial_M is not None:
        if len(initial_M) != L:
            raise ValueError("initial MPS length mismatch")
        M = [np.asarray(A,dtype=dt).copy() for A in initial_M]
        if any(A.shape[1] != d for A in M):
            raise ValueError("initial MPS physical dimension mismatch")
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
    S_cut = np.zeros(L - 1)
    stats = []

    def local(j, rt=None, mi=400):
        th = np.tensordot(M[j], M[j + 1], axes=([2], [0]))
        shape = th.shape
        mv = _heff_apply(Lenv[j], Ws[j], Ws[j + 1], Renv[j + 2], shape)
        if th.size <= 64:
            Hd = np.array([mv(e) for e in np.eye(th.size, dtype=dt)]).T
            w, v = np.linalg.eigh((Hd + Hd.conj().T) / 2)
            return w[0], v[:, 0], shape, 0, 0.0
        dg = _heff_diag(Lenv[j], Ws[j], Ws[j + 1], Renv[j + 2])
        x0 = th.reshape(-1)
        if not np.any(x0):
            x0 = rng.standard_normal(x0.size).astype(dt)
        e, v, nmv, rn = davidson(mv, x0, dg, res_tol=(res_tol if rt is None else rt), maxit=mi)
        return e, v, shape, nmv, rn

    for sweep, chi in enumerate(chis):
        trunc_max = 0.0; nmvs = []; rmax = 0.0
        # loose local solves while chi still grows (all but the last two sweeps), tight ones at the end
        loose = sweep < len(chis) - 2
        rt, mi = (max(res_tol, loose_tol), 60) if loose else (res_tol, 400)
        for j in range(L - 1):
            E, th, shape, nmv, rn = local(j, rt, mi); nmvs.append(nmv); rmax = max(rmax, rn)
            th = th.reshape(shape[0] * d, d * shape[3])
            U, s, Vh = np.linalg.svd(th, full_matrices=False)
            keep = min(chi, max(1, int(np.sum(s > cutoff * s[0]))))
            tr = float(np.sum(s[keep:] ** 2) / np.sum(s ** 2)); trunc_max = max(trunc_max, tr)
            s = s[:keep] / np.linalg.norm(s[:keep])
            S_cut[j] = entropy(s)
            M[j] = U[:, :keep].reshape(shape[0], d, keep)
            M[j + 1] = (s[:, None] * Vh[:keep]).reshape(keep, d, shape[3])
            Lenv[j + 1] = _lupdate(Lenv[j], M[j], Ws[j])
        for j in range(L - 2, -1, -1):
            E, th, shape, nmv, rn = local(j, rt, mi); nmvs.append(nmv); rmax = max(rmax, rn)
            th = th.reshape(shape[0] * d, d * shape[3])
            U, s, Vh = np.linalg.svd(th, full_matrices=False)
            keep = min(chi, max(1, int(np.sum(s > cutoff * s[0]))))
            tr = float(np.sum(s[keep:] ** 2) / np.sum(s ** 2)); trunc_max = max(trunc_max, tr)
            s = s[:keep] / np.linalg.norm(s[:keep])
            M[j + 1] = Vh[:keep].reshape(keep, d, shape[3])
            M[j] = (U[:, :keep] * s[None, :]).reshape(shape[0], d, keep)
            Renv[j + 1] = _rupdate(Renv[j + 2], M[j + 1], Ws[j + 1])
        hist.append((chi, float(np.real(E)), trunc_max))
        stats.append(dict(chi=chi, mv_mean=float(np.mean(nmvs)), mv_max=int(np.max(nmvs)), res_max=float(rmax), local_solves=len(nmvs), local_maxit_hits=sum(n >= mi+1 for n in nmvs)))
        if verbose:
            print(f"  sweep {sweep} chi={chi} E={np.real(E):.12f} trunc={trunc_max:.2e} mv/solve={np.mean(nmvs):.1f} "
                  f"(max {np.max(nmvs)}) res_max={rmax:.1e}", flush=True)
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
    return dict(E=float(np.real(E)), S=S_cut, hist=hist, M=M, stats=stats)
