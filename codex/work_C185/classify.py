"""C185 R1, R2: extensive ground spaces among classical and frustration-free qubit bonds.

R1. Classical bond J ZZ + b Z (per bond (b/2)(Z1+Z2)): exact count of open-chain ground states by transfer matrix.
R2. Frustration-free: ground space of sum_j P_{j,j+1} (P = projector onto range of the bond, kernel S) on the open
chain, built iteratively G_{L+1} = (G_L (x) C^2) cap (C^{2^{L-1}} (x) S). Only S matters.
"""
import numpy as np, itertools
rng = np.random.default_rng(185)
np.set_printoptions(linewidth=160)

def classical_count(J, b, L):
    s = np.array([1, -1])
    e = np.array([[J*s1*s2 + b/2*(s1+s2) for s2 in s] for s1 in s])
    # count minimum-energy strings: enumerate via DP of (min energy, multiplicity)
    best = {s1: (0.0, 1) for s1 in (0, 1)}
    for _ in range(L-1):
        nb = {}
        for s2 in (0, 1):
            cands = [(best[s1][0] + e[s1, s2], best[s1][1]) for s1 in (0, 1)]
            m = min(c[0] for c in cands)
            nb[s2] = (m, sum(c[1] for c in cands if abs(c[0]-m) < 1e-9))
        best = nb
    m = min(v[0] for v in best.values())
    return sum(v[1] for v in best.values() if abs(v[0]-m) < 1e-9)

def ff_dims(Pbond, Lmax=14, cap=4000):
    """Pbond: 4x4 projector (range = forbidden). Returns [dim G_L for L=2..Lmax]."""
    w, v = np.linalg.eigh(Pbond)
    K = v[:, w < 0.5]  # kernel basis 4 x k
    B = K.copy()  # basis of G_2, 4 x d
    dims = [B.shape[1]]
    Pi = Pbond
    for L in range(3, Lmax+1):
        d = B.shape[1]
        if d == 0:
            dims.append(0); continue
        if d > cap:
            dims.append(-1); break
        # candidates: B (x) I2, shape 2^L x 2d
        C = np.kron(B, np.eye(2))
        # constraint on last two sites: (I_{2^{L-2}} (x) Pi) C = 0
        Cr = C.reshape(2**(L-2), 4, 2*d)
        M = np.einsum('ab,xbk->xak', Pi, Cr).reshape(2**L, 2*d)
        # null space of M in coefficient space
        u, s, vh = np.linalg.svd(M, full_matrices=False)
        tol = 1e-9 * max(1.0, s.max() if s.size else 1.0)
        null = vh[np.sum(s > tol):].conj().T
        B = C @ null
        # re-orthonormalize
        q, r = np.linalg.qr(B)
        B = q[:, np.abs(np.diag(r)) > 1e-10] if B.shape[1] else B
        dims.append(B.shape[1])
    return dims

def proj(vecs):
    q, _ = np.linalg.qr(np.array(vecs, complex).T)
    return q @ q.conj().T

def ket(a): return np.array(a, complex)/np.linalg.norm(a)

if __name__ == "__main__":
    print("R1: classical bond J ZZ + b Z, J = 1. Open-chain ground-state count for L = 4, 8, 12, 16.")
    for J in (1.0, -1.0):
        for b in (0.0, 0.5, 1.0, 1.5, 1.9, 2.0, 2.1, 3.0, -2.0):
            print(f"  J = {J:+.0f}, b = {b:+.2f}: ", [classical_count(J, b, L) for L in (4, 8, 12, 16)])
    print("\nR2: frustration-free bonds, dim G_L for L = 2..14 (open chain)")
    up, dn = ket([1, 0]), ket([0, 1])
    def rank1(psi): return proj([psi])
    th = 0.7
    a = ket([np.cos(th), np.sin(th)])
    b2 = ket([np.cos(1.9), np.exp(0.4j)*np.sin(1.9)])
    cases = {
        "rank1 |11> (Fibonacci)": rank1(np.kron(dn, dn)),
        "rank1 |aa>, a generic": rank1(np.kron(a, a)),
        "rank1 |10> (a perp b)": rank1(np.kron(dn, up)),
        "rank1 |a b>, generic a,b": rank1(np.kron(a, b2)),
        "rank1 |1 a>, <1|a> != 0": rank1(np.kron(dn, a)),
        "rank1 singlet (FM Heisenberg)": rank1(ket([0, 1, -1, 0])),
        "rank1 Bell |00>+|11>": rank1(ket([1, 0, 0, 1])),
        "rank1 |01>+i|10>": rank1(ket([0, 1, 1j, 0])),
        "rank2 S^perp={|11>,|10>}": proj([np.kron(dn, dn), np.kron(dn, up)]),
        "rank2 S^perp={|11>,|00>}": proj([np.kron(dn, dn), np.kron(up, up)]),
        "rank2 S^perp={|11>,|01>-|10>}": proj([np.kron(dn, dn), ket([0, 1, -1, 0])]),
        "rank2 S^perp={|11>,|01>+|10>}": proj([np.kron(dn, dn), ket([0, 1, 1, 0])]),
        "rank2 S^perp={|01>,|10>}": proj([np.kron(up, dn), np.kron(dn, up)]),
        "rank3 S={|00>}": proj([np.kron(dn, dn), np.kron(up, dn), np.kron(dn, up)]),
        "rank3 S={singlet}": proj([ket([1, 0, 0, 0]), ket([0, 0, 0, 1]), ket([0, 1, 1, 0])]),
    }
    for name, P in cases.items():
        print(f"  {name:34s}", ff_dims(P, 14))
    print("\n  random kernels (20 each): max over samples of dim G_L at L = 6, 10, 14")
    for r in (1, 2, 3):
        mx = np.zeros(3, int)
        for _ in range(20):
            M = rng.normal(size=(4, r)) + 1j*rng.normal(size=(4, r))
            P = proj(list(M.T))
            d = ff_dims(P, 14)
            mx = np.maximum(mx, [d[4], d[8], d[12]])
        print(f"  rank {r}: {mx.tolist()}")
    print("\n  random product rank-1 |a b> (20 samples): dims at L = 6, 10, 14; and |a a> (5 samples)")
    out = []
    for _ in range(20):
        a_ = ket(rng.normal(size=2) + 1j*rng.normal(size=2)); b_ = ket(rng.normal(size=2) + 1j*rng.normal(size=2))
        d = ff_dims(rank1(np.kron(a_, b_)), 14); out.append((d[4], d[8], d[12]))
    print("   |ab>:", sorted(set(out)))
    out = []
    for _ in range(5):
        a_ = ket(rng.normal(size=2) + 1j*rng.normal(size=2))
        d = ff_dims(rank1(np.kron(a_, a_)), 14); out.append((d[4], d[8], d[12]))
    print("   |aa>:", sorted(set(out)))
