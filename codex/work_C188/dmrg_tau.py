"""C188 (Claude Opus 5.5; einsum optimize=True added 18:40 on 2026-11-09, results unchanged, much faster): open-chain DMRG of the cyclic stratum in the tau frame (C154 Prop. 2), where
h = a XZ + b YX + g ZY becomes the real period-3 XYZ chain with bond j couplings
(Jx, Jy, Jz) = (a, b, g), (g, a, b), (b, g, a) for j = 0, 1, 2 mod 3. The map is a product of on-site unitaries,
so open-chain energies and entanglement equal those of the sigma-frame chain.
Usage: python3 dmrg_tau.py L chis out.jsonl a,b,g [a,b,g ...]
Records E, E/L, S(l), bond energies, single-site <X>, <Z>, parity-resolved and averaged CC fits, truncation."""
import sys, json, time
import numpy as np
from dmrg3 import dmrg

X = np.array([[0., 1.], [1., 0.]]); iY = np.array([[0., 1.], [-1., 0.]]); Z = np.diag([1., -1.]); I2 = np.eye(2)


def couplings(a, b, g, j):
    return [(a, b, g), (g, a, b), (b, g, a)][j % 3]


def mpo_tau(a, b, g, L):
    # Y (x) Y = -(iY) (x) (iY)
    ops = [X, iY, Z]; sgn = [1., -1., 1.]
    Ws = []
    for j in range(L):
        W = np.zeros((5, 5, 2, 2))
        W[0, 0] = I2; W[4, 4] = I2
        if j < L - 1:
            J = couplings(a, b, g, j)
            for k in range(3):
                W[0, 1 + k] = sgn[k] * J[k] * ops[k]
        for k in range(3):
            W[1 + k, 4] = ops[k]
        if j == 0:
            W = W[0:1]
        elif j == L - 1:
            W = W[:, 4:5]
        Ws.append(W)
    return Ws


def site_expect(M, O):
    L = len(M)
    R = [None] * (L + 1); R[L] = np.ones((1, 1))
    for j in range(L - 1, -1, -1):
        A = M[j]; R[j] = np.einsum("asc,cd,bsd->ab", np.conj(A), R[j + 1], A, optimize=True)
    out = np.zeros(L); Lf = np.ones((1, 1))
    for j in range(L):
        A = M[j]
        T = np.einsum("ab,asc,st,btd->cd", Lf, np.conj(A), O, A, optimize=True)
        out[j] = float(np.real(np.einsum("cd,cd->", T, R[j + 1])))
        Lf = np.einsum("ab,asc,bsd->cd", Lf, np.conj(A), A, optimize=True)
    return out


def bond_expect(M, O1, O2):
    L = len(M)
    R = [None] * (L + 1); R[L] = np.ones((1, 1))
    for j in range(L - 1, -1, -1):
        A = M[j]; R[j] = np.einsum("asc,cd,bsd->ab", np.conj(A), R[j + 1], A, optimize=True)
    out = np.zeros(L - 1); Lf = np.ones((1, 1))
    for j in range(L - 1):
        A = M[j]; B = M[j + 1]
        T = np.einsum("ab,asc,st,btd->cd", Lf, np.conj(A), O1, A, optimize=True)
        T = np.einsum("ab,asc,st,btd->cd", T, np.conj(B), O2, B, optimize=True)
        out[j] = float(np.real(np.einsum("cd,cd->", T, R[j + 2])))
        Lf = np.einsum("ab,asc,bsd->cd", Lf, np.conj(A), A, optimize=True)
    return out


def cfit(S, L, lo, hi):
    S = np.asarray(S); l = np.arange(1, L); x = np.log((2 * L / np.pi) * np.sin(np.pi * l / L))
    m = (l >= lo * L) & (l <= hi * L)
    par = [float(6 * np.polyfit(x[m & (l % 2 == p)], S[m & (l % 2 == p)], 1)[0]) for p in (0, 1)]
    Sa = 0.5 * (S[:-1] + S[1:]); la = l[:-1] + 0.5; xa = np.log((2 * L / np.pi) * np.sin(np.pi * la / L)); ma = (la >= lo * L) & (la <= hi * L)
    return par + [float(6 * np.polyfit(xa[ma], Sa[ma], 1)[0])]


if __name__ == "__main__":
    L = int(sys.argv[1]); chis = tuple(int(c) for c in sys.argv[2].split(",")); fout = sys.argv[3]
    for spec in sys.argv[4:]:
        a, b, g = [float(x) for x in spec.split(",")]
        t0 = time.time()
        r = dmrg(None, L, chis=chis, Ws=mpo_tau(a, b, g, L))
        M = r["M"]
        xx = bond_expect(M, X, X); yy = -bond_expect(M, iY, iY); zz = bond_expect(M, Z, Z)
        e = np.array([np.dot(couplings(a, b, g, j), (xx[j], yy[j], zz[j])) for j in range(L - 1)])
        mx = site_expect(M, X); mz = site_expect(M, Z)
        S = np.asarray(r["S"])
        rec = dict(a=a, b=b, g=g, L=L, chis=list(chis), E=r["E"], E_per_site=r["E"] / L,
                   trunc=r["hist"][-1][2], dE_last=(r["hist"][-1][1] - r["hist"][-2][1]) if len(r["hist"]) > 1 else None,
                   c_fit=cfit(S, L, 0.25, 0.75), c_fit_wide=cfit(S, L, 0.15, 0.85), S=[round(float(s), 7) for s in S],
                   S_mid4=float(np.mean(S[L // 2 - 3:L // 2 + 1])),
                   e=[round(float(x), 8) for x in e], mx=[round(float(x), 5) for x in mx], mz=[round(float(x), 5) for x in mz],
                   sec=round(time.time() - t0, 1))
        with open(fout, "a") as f:
            f.write(json.dumps(rec) + "\n")
        print(f"({a},{b},{g}) L={L}: E/L={rec['E_per_site']:.10f} trunc={rec['trunc']:.1e} dE={rec['dE_last']} c_fit={np.round(rec['c_fit'],3).tolist()} "
              f"S_mid4={rec['S_mid4']:.4f} |m|max={max(np.max(np.abs(mx)), np.max(np.abs(mz))):.3f} ({rec['sec']}s)", flush=True)
