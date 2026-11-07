"""Search for qubit bonds with a two-site local integral of motion q_{j,j+1} (each one separately conserved).

Such operators are the mechanism of exact extensive degeneracy at the classical point (q = Z_j) and at the
XZ point (q = Z_j X_{j+1}). For a bond x = (J, B), let C(x) be the real linear map q -> [H_L, q_{0,1}] on the
15-dimensional space of traceless Hermitian two-site operators, on a ring of L = 6 sites. The cost
sigma_min(C)/||H|| vanishes iff a nonzero q exists (on this ring). Nelder-Mead from random starts.
Usage: python3 local_ioms.py NSTARTS SEED OUT.json
"""
import sys, json
import numpy as np
from scipy.optimize import minimize
from ed import bond, ham_sparse, PAULI, I2, params_to_JB

L = 6
P4 = [I2] + PAULI


def embed(op2, L):
    return np.kron(op2, np.eye(2 ** (L - 2)))


QBASIS = []
for a in range(4):
    for b in range(4):
        if a == 0 and b == 0:
            continue
        QBASIS.append(((a, b), embed(np.kron(P4[a], P4[b]), L)))


def cmat(x):
    J, B = params_to_JB(x)
    H = ham_sparse(bond(J, B), L).toarray()
    cols = []
    for _, Q in QBASIS:
        C = H @ Q - Q @ H
        cols.append(np.concatenate([C.real.ravel(), C.imag.ravel()]))
    return np.array(cols).T, np.linalg.norm(H)


def comm3(x):
    """Gauge-invariant three-body part of [h_12, h_23], relative to |J|^2: zero iff J = lambda n n^T (classical)."""
    J, B = params_to_JB(x)
    tot = 0.0
    for a in range(3):
        for d in range(3):
            tot += np.sum(np.cross(J[a, :], J[:, d]) ** 2)
    return np.sqrt(tot) / max(np.sum(J ** 2), 1e-12)


def cost(x):
    x = x / np.linalg.norm(x)
    C, nH = cmat(x)
    s = np.linalg.svd(C, compute_uv=False)
    # exclude the classical set (J = lambda n n^T, which includes J = 0) by a penalty
    J, B = params_to_JB(x)
    return s[-1] / nH + 10.0 * max(0.0, 0.25 - comm3(x)) + 10.0 * max(0.0, 0.5 - np.linalg.norm(J))


def iom(x):
    x = x / np.linalg.norm(x)
    C, nH = cmat(x)
    u, s, vt = np.linalg.svd(C)
    v = vt[-1]
    labs = "1xyz"
    terms = {labs[a] + labs[b]: float(v[i]) for i, ((a, b), _) in enumerate(QBASIS) if abs(v[i]) > 1e-3}
    return s[-1] / nH, s[-2] / nH, terms


if __name__ == "__main__":
    n = int(sys.argv[1]); seed = int(sys.argv[2]); out = sys.argv[3]
    rng = np.random.default_rng(seed)
    res = []
    for i in range(n):
        x0 = rng.normal(size=12)
        r = minimize(cost, x0, method="Nelder-Mead", options=dict(maxfev=6000, xatol=1e-8, fatol=1e-12, adaptive=True))
        x = r.x / np.linalg.norm(r.x)
        s1, s2, terms = iom(x)
        J, B = params_to_JB(x)
        res.append({"x": x.tolist(), "sigma_min": s1, "sigma_next": s2, "q": terms,
                    "J_sv": np.linalg.svd(J, compute_uv=False).tolist(), "B": B.tolist(),
                    "Jsym_minus_Jsym_T": float(np.linalg.norm(J - J.T)), "comm3": float(comm3(x))})
        print(i, f"sigma_min {s1:.2e} next {s2:.2e}", "J_sv", np.round(np.linalg.svd(J, compute_uv=False), 4),
              "|B|", round(float(np.linalg.norm(B)), 4), "comm3", round(float(comm3(x)), 3), "q", {k: round(v, 3) for k, v in terms.items()}, flush=True)
    json.dump(res, open(out, "w"), indent=1)
