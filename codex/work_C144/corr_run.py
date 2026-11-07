"""Connected two-point functions of bond operators a_j = X_j Z_{j+1}, b_j = Z_j X_{j+1} and of Y_j in the DMRG
ground state of an open chain, from the centre-left site j0 = L/4 outwards. Distinguishes a gapped A sector
(exponential a-a decay) from a critical one (power law)."""
import sys, json
import numpy as np
from xzlib import term, SX, SY, SZ, I2
from dmrg import dmrg

def right_envs(M):
    L = len(M); R = [None] * (L + 1); R[L] = np.ones((1, 1))
    for j in range(L - 1, -1, -1):
        R[j] = np.einsum("asc,cd,bsd->ab", np.conj(M[j]), R[j + 1], M[j])
    return R

def apply(T, A, O):
    return np.einsum("ab,asc,st,btd->cd", T, np.conj(A), O, A)

def corr(M, R, ops, j0, rmax):
    """ops: list of single-site operators (len 1 or 2) forming the local operator starting at site j."""
    L = len(M); n = len(ops)
    T = np.ones((1, 1))
    for j in range(j0):
        T = apply(T, M[j], I2)
    # one-point functions
    def onept(j):
        U = np.ones((1, 1))
        for k in range(j):
            U = apply(U, M[k], I2)
        for k in range(n):
            U = apply(U, M[j + k], ops[k])
        return np.einsum("cd,cd->", U, R[j + n]).real
    T1 = T
    for k in range(n):
        T1 = apply(T1, M[j0 + k], ops[k])
    out = []
    Tc = T1
    pos = j0 + n
    for r in range(n, rmax + 1):
        # operator at j = j0 + r
        while pos < j0 + r:
            Tc = apply(Tc, M[pos], I2); pos += 1
        U = Tc
        for k in range(n):
            U = apply(U, M[j0 + r + k], ops[k])
        val = np.einsum("cd,cd->", U, R[j0 + r + n]).real
        out.append((r, val))
    o0 = onept(j0)
    return [(r, v - o0 * onept(j0 + r)) for r, v in out]

L = int(sys.argv[1]); al, be = float(sys.argv[2]), float(sys.argv[3])
h = term("xz") + al * term("zx") + be * term("yy")
r = dmrg(h, L, chis=(16, 32, 48, 64, 64))
R = right_envs(r["M"])
j0 = L // 4; rmax = L // 2
res = {}
for name, ops in [("a", [SX, SZ]), ("b", [SZ, SX]), ("y", [SY])]:
    res[name] = corr(r["M"], R, ops, j0, rmax)
    rs = [x[0] for x in res[name]]; vs = [x[1] for x in res[name]]
    sel = [i for i, rr in enumerate(rs) if rr in (2, 3, 4, 6, 8, 12, 16, 24, 32)]
    print(f"alpha={al} beta={be} L={L} {name}: " + " ".join(f"r={rs[i]}:{vs[i]:+.2e}" for i in sel), flush=True)
json.dump(dict(L=L, alpha=al, beta=be, corr=res), open(f"corr_L{L}_a{al}_b{be}.json", "w"))
