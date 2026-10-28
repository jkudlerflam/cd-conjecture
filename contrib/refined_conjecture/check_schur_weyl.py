"""Check uls_schur_weyl.py against brute-force exact diagonalization on small rings."""
import numpy as np
from uls_schur_weyl import partitions, XOp, gl_mult, hook_dim, yor_generators

def brute(L, d):
    D = d**L
    idx = np.arange(D)
    digits = np.array([(idx // d**(L-1-i)) % d for i in range(L)]).T
    H = np.zeros((D, D))
    for i in range(L):
        j = (i+1) % L
        sw = digits.copy(); sw[:, [i, j]] = sw[:, [j, i]]
        tgt = (sw * d**np.arange(L-1, -1, -1)).sum(1)
        H[tgt, idx] += 1.0
    return np.sort(np.linalg.eigvalsh(H))

for L, d in [(6, 3), (8, 2), (5, 4), (6, 2)]:
    ev = []
    for lam in partitions(L):
        m = gl_mult(lam, d)
        if m == 0: continue
        e = np.linalg.eigvalsh(XOp(lam).dense())
        ev += list(np.repeat(e, int(round(m))))
    ev = np.sort(ev); b = brute(L, d)
    print(L, d, len(ev), len(b), np.max(np.abs(ev - b)))
# Coxeter relations on one irrep
g = yor_generators((3, 2, 2, 1))
I = np.eye(g[0].shape[0])
print("s^2=1:", max(np.abs((s @ s).toarray() - I).max() for s in g))
print("braid:", max(np.abs((g[k] @ g[k+1] @ g[k] - g[k+1] @ g[k] @ g[k+1]).toarray()).max() for k in range(len(g)-1)))
print("far commute:", np.abs((g[0] @ g[3] - g[3] @ g[0]).toarray()).max())
