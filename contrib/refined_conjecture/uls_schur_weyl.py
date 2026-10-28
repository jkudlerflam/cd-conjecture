"""Spectrum of the SU(d) Uimin-Lai-Sutherland ring H = sum_i P_{i,i+1} for all d at once.

Schur-Weyl duality: (C^d)^{(x)L} = sum_lambda V_lambda(GL_d) (x) S^lambda, and H acts as the
group-algebra element X = sum_{i=1}^{L} (i, i+1 mod L) on the S_L irrep S^lambda. So the
spectrum of H is the union over partitions lambda of L (at most d rows) of spec(X|S^lambda),
each eigenvalue with multiplicity dim V_lambda(d) = prod_cells (d + content)/hook.

X|S^lambda is built in Young's orthogonal form. Small irreps are diagonalized exactly; large
ones are treated by the finite-temperature Lanczos method (FTLM), which returns a list of
(Ritz energy, weight) pairs whose weights sum to dim S^lambda.

Usage:  python3 uls_schur_weyl.py L [maxrows] [dense_max] [R] [M] [ntasks] [task]
Writes data/uls_L{L}.npz (or data/uls_L{L}_r{maxrows}_part{task}.npz when split over
ntasks array tasks) with, per partition, energies and weights.
"""
import sys, time
import numpy as np
import scipy.sparse as sp
from math import factorial


def partitions(n, maxpart=None, maxrows=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    if maxrows == 0:
        return
    for k in range(min(n, maxpart), 0, -1):
        for p in partitions(n - k, k, None if maxrows is None else maxrows - 1):
            yield (k,) + p


def hook_dim(lam):
    n = sum(lam)
    conj = [sum(1 for r in lam if r > j) for j in range(lam[0])]
    h = 1
    for i, r in enumerate(lam):
        for j in range(r):
            h *= (r - j - 1) + (conj[j] - i - 1) + 1
    return factorial(n) // h


def gl_mult(lam, d):
    """dim of the GL(d) irrep with highest weight lam (0 if more than d rows), as float."""
    if len(lam) > d:
        return 0.0
    conj = [sum(1 for r in lam if r > j) for j in range(lam[0])]
    m = 1.0
    for i, r in enumerate(lam):
        for j in range(r):
            hook = (r - j - 1) + (conj[j] - i - 1) + 1
            m *= (d + j - i) / hook
    return m


def syt(lam):
    """Standard Young tableaux as Yamanouchi words (row of entry k) and contents."""
    L, nr = sum(lam), len(lam)
    words = np.zeros((1, 0), dtype=np.int8)
    rowlen = np.zeros((1, nr), dtype=np.int16)
    conts = np.zeros((1, 0), dtype=np.int16)
    lam_a = np.array(lam)
    for k in range(L):
        W, R, C = [], [], []
        for r in range(nr):
            ok = rowlen[:, r] < lam_a[r]
            if r > 0:
                ok &= rowlen[:, r] < rowlen[:, r - 1]
            idx = np.nonzero(ok)[0]
            if len(idx) == 0:
                continue
            w = np.concatenate([words[idx], np.full((len(idx), 1), r, np.int8)], axis=1)
            c = np.concatenate([conts[idx], (rowlen[idx, r] - r)[:, None].astype(np.int16)], axis=1)
            rl = rowlen[idx].copy()
            rl[:, r] += 1
            W.append(w); C.append(c); R.append(rl)
        words, conts, rowlen = np.concatenate(W), np.concatenate(C), np.concatenate(R)
    return words, conts


def encode(words, base):
    code = np.zeros(len(words), dtype=np.int64)
    for k in range(words.shape[1]):
        code = code * base + words[:, k].astype(np.int64)
    return code


def yor_generators(lam):
    """Sparse matrices of s_k = (k, k+1), k = 0..L-2, in Young's orthogonal form."""
    words, conts = syt(lam)
    D, L = words.shape
    base = max(len(lam), 2)
    codes = encode(words, base)
    order = np.argsort(codes)
    codes_sorted = codes[order]
    gens = []
    ar = np.arange(D)
    for k in range(L - 1):
        a = (conts[:, k + 1] - conts[:, k]).astype(float)  # axial distance
        rows, cols, vals = [ar], [ar], [1.0 / a]
        sw = np.abs(a) >= 2
        if sw.any():
            w2 = words[sw].copy()
            w2[:, [k, k + 1]] = w2[:, [k + 1, k]]
            pos = np.searchsorted(codes_sorted, encode(w2, base))
            j = order[pos]
            assert np.all(codes[j] == encode(w2, base))
            rows.append(ar[sw]); cols.append(j)
            vals.append(np.sqrt(1.0 - 1.0 / a[sw] ** 2))
        S = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(D, D))
        gens.append(S)
    return gens


class XOp:
    """X = sum_{k} s_k + (1 L), with (1 L) = s_1 ... s_{L-2} s_{L-1} s_{L-2} ... s_1 (1-based)."""

    def __init__(self, lam):
        self.g = yor_generators(lam)
        self.D = self.g[0].shape[0] if self.g else 1
        self.Ssum = sum(self.g)
        n = len(self.g)
        self.wrap_seq = list(range(n)) + list(range(n - 2, -1, -1))  # apply right to left

    def matvec(self, v):
        out = self.Ssum @ v
        w = v
        for k in self.wrap_seq[::-1]:
            w = self.g[k] @ w
        return out + w

    def dense(self):
        return self.matvec(np.eye(self.D))


def ftlm(op, D, R=20, M=120, seed=0):
    rng = np.random.default_rng(seed)
    E_all, w_all = [], []
    for r in range(R):
        v0 = rng.choice([-1.0, 1.0], size=D)
        v0 /= np.linalg.norm(v0)
        V = np.zeros((M + 1, D))
        V[0] = v0
        al, be = [], []
        v, vprev, b = v0, np.zeros(D), 0.0
        m_eff = M
        for j in range(M):
            w = op.matvec(v) - b * vprev
            a = w @ v
            w -= a * v
            # full reorthogonalization keeps the Ritz values clean
            w -= V[: j + 1].T @ (V[: j + 1] @ w)
            al.append(a)
            b = np.linalg.norm(w)
            if b < 1e-10 or j == M - 1:
                m_eff = j + 1
                break
            be.append(b)
            vprev, v = v, w / b
            V[j + 1] = v
        T = np.diag(al) + np.diag(be[: m_eff - 1], 1) + np.diag(be[: m_eff - 1], -1)
        ev, U = np.linalg.eigh(T)
        E_all.append(ev)
        w_all.append(U[0] ** 2 * D / R)
    return np.concatenate(E_all), np.concatenate(w_all)


def main():
    L = int(sys.argv[1])
    maxrows = int(sys.argv[2]) if len(sys.argv) > 2 else L
    dense_max = int(sys.argv[3]) if len(sys.argv) > 3 else 6000
    R = int(sys.argv[4]) if len(sys.argv) > 4 else 20
    M = int(sys.argv[5]) if len(sys.argv) > 5 else 120
    ntasks = int(sys.argv[6]) if len(sys.argv) > 6 else 1
    task = int(sys.argv[7]) if len(sys.argv) > 7 else 0
    out = {}
    t0 = time.time()
    lams = sorted(partitions(L, maxrows=maxrows), key=hook_dim, reverse=True)
    for i, lam in enumerate(lams):
        if i % ntasks != task:
            continue
        D = hook_dim(lam)
        op = XOp(lam)
        key = "_".join(map(str, lam))
        if D <= dense_max:
            ev = np.linalg.eigvalsh(op.dense())
            out["E_" + key] = ev
            out["w_" + key] = np.ones_like(ev)
            out["exact_" + key] = np.array(1)
        else:
            E, w = ftlm(op, D, R=R, M=M, seed=hash(lam) % 2**32)
            out["E_" + key] = E
            out["w_" + key] = w
            out["exact_" + key] = np.array(0)
        print(f"L={L} lam={lam} D={D} {'dense' if D <= dense_max else 'FTLM'} t={time.time()-t0:.1f}s", flush=True)
    fn = f"data/uls_L{L}.npz" if ntasks == 1 else f"data/uls_L{L}_r{maxrows}_part{task}.npz"
    np.savez_compressed(fn, **out)


if __name__ == "__main__":
    main()
