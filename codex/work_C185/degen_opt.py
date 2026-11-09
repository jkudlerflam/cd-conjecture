"""C185 R3: unbiased search for extensively degenerate ground spaces in B_2.

Bond h = sum_ab J_ab s^a (x) s^b + sum_a B_a (s^a (x) 1 + 1 (x) s^a)/2, ring of L sites (so the field per site is B).
Parameters p = (J (9), B (3)) normalized to |p| = 1. Objective, for a target count (K8, K10):
   W(p) = max_L (E_{K_L - 1} - E_0)(L) / L       over L = 8, 10   (width of the lowest K_L states per site)
normalized by the bond norm (fixed to 1), so W = 0 exactly at a point whose ground space has >= K_L states.
Also reported: the gap above the manifold (E_{K_L} - E_{K_L - 1}), which distinguishes a manifold from a dense band.
Nelder-Mead from random starts, then a polish. Usage: python3 degen_opt.py K8 K10 seed nstart out.jsonl
"""
import numpy as np, sys, json, time
from scipy.optimize import minimize
from scipy.sparse import kron, identity, csr_matrix
s = [csr_matrix(np.array(m, complex)) for m in ([[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]])]

def site_op(L, j, o):
    M = identity(1, format='csr', dtype=complex)
    for k in range(L):
        M = kron(M, o if k == j else identity(2, format='csr', dtype=complex), format='csr')
    return M

def build_ops(L):
    S = [[site_op(L, j, s[a]) for a in range(3)] for j in range(L)]
    JJ = [[sum(S[j][a] @ S[(j+1) % L][b] for j in range(L)).toarray() for b in range(3)] for a in range(3)]
    BB = [sum(S[j][a] for j in range(L)).toarray() for a in range(3)]
    return [JJ[a][b] for a in range(3) for b in range(3)] + BB

OPS = {}
def spectrum(p, L):
    if L not in OPS:
        OPS[L] = build_ops(L)
    H = sum(c*O for c, O in zip(p, OPS[L]))
    return np.linalg.eigvalsh(H)

def W(p, Ks):
    p = np.asarray(p)/np.linalg.norm(p)
    out = 0.0
    for L, K in Ks.items():
        e = spectrum(p, L)
        out = max(out, (e[K-1]-e[0])/L)
    return out

def run(K8, K10, seed, nstart):
    rng = np.random.default_rng(seed)
    Ks = {8: K8, 10: K10}
    res = []
    for i in range(nstart):
        p0 = rng.normal(size=12); p0 /= np.linalg.norm(p0)
        t0 = time.time()
        r = minimize(W, p0, args=(Ks,), method='Nelder-Mead', options={'maxfev': 4000, 'xatol': 1e-7, 'fatol': 1e-10, 'adaptive': True})
        r = minimize(W, r.x, args=(Ks,), method='Nelder-Mead', options={'maxfev': 3000, 'xatol': 1e-9, 'fatol': 1e-12, 'adaptive': True})
        p = r.x/np.linalg.norm(r.x)
        info = {}
        for L, K in Ks.items():
            e = spectrum(p, L)
            info[L] = dict(width=float((e[K-1]-e[0])/L), gap_above=float(e[K]-e[K-1]), K=K)
        J = p[:9].reshape(3, 3); B = p[9:]
        sv = np.linalg.svd(J, compute_uv=False)
        res.append(dict(K8=K8, K10=K10, seed=seed, start=i, W=float(r.fun), p=p.tolist(), Jsv=sv.tolist(),
                        Jsym_eig=np.linalg.eigvalsh((J+J.T)/2).tolist(), Janti=[(J-J.T)[1, 2]/2, (J-J.T)[2, 0]/2, (J-J.T)[0, 1]/2],
                        Bnorm=float(np.linalg.norm(B)), info={str(k): v for k, v in info.items()}, sec=time.time()-t0))
    return res

if __name__ == "__main__":
    K8, K10, seed, nstart, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    for r in run(K8, K10, seed, nstart):
        with open(out, 'a') as f:
            f.write(json.dumps(r) + "\n")
