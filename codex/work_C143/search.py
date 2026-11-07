"""Search B_2 for singular qubit bonds.

Criterion (PROOF section 1): the soft count of low ring states at the CFT scale,
    S_L(h) = sum_n exp(-L (E_n - E_0)/(kappa |x|)),
is bounded in L at regular bonds (finite CFT or gapped with finite degeneracy) and grows at singular
bonds (z > 1 bands: S_L grows like a power of L; extensive degeneracy: exponentially).
Objective: g(h) = min_L [log S_{L+2} - log S_L] over L1 <= L < L2 (even), maximized by Nelder-Mead from random starts.
|x| is the Euclidean norm of the 12 coefficients (J, B), an equivalent norm to M on the quotient.

Usage: python3 search.py SEED NSTARTS L1 L2 OUTFILE [maxfev]
"""
import sys, json, time
import numpy as np
from scipy.optimize import minimize
from ed import bond, spectrum, all_levels, params_to_JB, ham_sparse, PAULI, I2, gauge_M, gauge_lmin

KAPPA = 3.0


def levels(x, L):
    J, B = params_to_JB(x)
    w, k = all_levels(spectrum(bond(J, B), L))
    return w, k


def softcount(x, L):
    nx = np.linalg.norm(x)
    w, _ = levels(x, L)
    e = (w - w.min()) * L / (KAPPA * nx)
    return np.sum(np.exp(-e[e < 60]))


def objective(x, L1, L2):
    """min over consecutive even lengths L1, L1+2, ..., L2 of log S_{L+2} - log S_L.
    (Version 1 used the single increment log S_L2 - log S_L1 and was fooled by ground-state level crossings
    at one length in incommensurate phases; see PROGRESS.md.)"""
    nx = np.linalg.norm(x)
    if nx < 1e-9:
        return 0.0
    logs = [np.log(softcount(x, L)) for L in range(L1, L2 + 1, 2)]
    return float(np.min(np.diff(logs)))


def u1_defect(x, L=6, stag=False):
    """min over unit n of ||[H_L, sum_j s_j n.sigma_j]||/||H_L|| with s_j = 1 or (-1)^j."""
    J, B = params_to_JB(x)
    H = ham_sparse(bond(J, B), L).toarray()
    ops = []
    for a in range(3):
        Q = np.zeros_like(H)
        for j in range(L):
            M = np.array([[1.0 + 0j]])
            for i in range(L):
                M = np.kron(M, PAULI[a] if i == j else I2)
            Q += ((-1) ** j if stag else 1) * M
        ops.append(Q)
    C = [H @ Q - Q @ H for Q in ops]
    G = np.array([[np.vdot(C[a], C[b]).real for b in range(3)] for a in range(3)])
    lam, vec = np.linalg.eigh(G)
    return np.sqrt(max(lam[0], 0)) / np.linalg.norm(H), vec[:, 0]


def diagnostics(x, Ls=(8, 10, 12)):
    x = np.asarray(x, float) / np.linalg.norm(x)
    J, B = params_to_JB(x)
    h = bond(J, B)
    out = {"x": x.tolist()}
    out["J_sv"] = np.linalg.svd(J, compute_uv=False).tolist()
    out["J_sym_sv"] = np.linalg.svd((J + J.T) / 2, compute_uv=False).tolist()
    Dvec = [J[1, 2] - J[2, 1], J[2, 0] - J[0, 2], J[0, 1] - J[1, 0]]
    out["DM"] = list(map(float, Dvec))
    out["B"] = B.tolist()
    M, _ = gauge_M(h)
    out["M"] = float(M)
    lmin, _ = gauge_lmin(h)
    out["u1_uniform"], n1 = u1_defect(x, stag=False)
    out["u1_stag"], n2 = u1_defect(x, stag=True)
    rows = []
    for L in Ls:
        w, k = levels(x, L)
        o = np.argsort(w); w = w[o]
        e = w - w[0]
        nlow = [int(np.sum(e < t * M / L)) for t in (0.5, 1.0, 2.0)]
        deg = int(np.sum(e < 1e-8))
        gaps = e[e > 1e-8][:3].tolist()
        rows.append({"L": L, "E0_per_L": float(w[0] / L), "frustration": float(w[0] / L - lmin),
                     "deg": deg, "n_low(0.5,1,2 M/L)": nlow, "gaps": gaps,
                     "S": float(softcount(x, L))})
    out["rings"] = rows
    return out


def run(seed, nstarts, L1, L2, outfile, maxfev=800):
    rng = np.random.default_rng(seed)
    res = []
    for s in range(nstarts):
        x0 = rng.normal(size=12); x0 /= np.linalg.norm(x0)
        t0 = time.time()
        f = lambda x: -objective(x / np.linalg.norm(x), L1, L2)
        r = minimize(f, x0, method="Nelder-Mead",
                     options=dict(maxfev=maxfev, xatol=1e-4, fatol=1e-6, initial_simplex=None, adaptive=True))
        x = r.x / np.linalg.norm(r.x)
        d = {"seed": seed, "start": s, "g_start": -f(x0), "g_final": -r.fun, "nfev": r.nfev,
             "secs": time.time() - t0, "x0": x0.tolist()}
        d.update(diagnostics(x))
        res.append(d)
        with open(outfile, "a") as fh:
            fh.write(json.dumps(d) + "\n")
        print(f"seed {seed} start {s}: g {d['g_start']:.3f} -> {d['g_final']:.3f} ({d['secs']:.0f}s)", flush=True)
    return res


if __name__ == "__main__":
    seed = int(sys.argv[1]); n = int(sys.argv[2]); L1 = int(sys.argv[3]); L2 = int(sys.argv[4])
    outfile = sys.argv[5]
    maxfev = int(sys.argv[6]) if len(sys.argv) > 6 else 800
    run(seed, n, L1, L2, outfile, maxfev)
