"""Packed Typhon search: one random start per worker at (L1, L2) = (10, 14), objective min over the increments 10->12->14, diagnostics at 10, 12, 14,
and a truncated soft count at L = 16 (lowest NEV levels per momentum sector, a lower bound on S_16)."""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
import search
from ed import bond, spectrum, all_levels, params_to_JB

OUT = sys.argv[1]
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 48
MAXFEV = int(sys.argv[3]) if len(sys.argv) > 3 else 900

def s16(x, nev=40):
    J, B = params_to_JB(x)
    w, _ = all_levels(spectrum(bond(J, B), 16, nev=nev))
    e = (w - w.min()) * 16 / (search.KAPPA * np.linalg.norm(x))
    return float(np.sum(np.exp(-e[e < 60])))

def task(seed):
    rng = np.random.default_rng(1000 + seed)
    x0 = rng.normal(size=12); x0 /= np.linalg.norm(x0)
    t0 = time.time()
    f = lambda x: -search.objective(x / np.linalg.norm(x), 10, 14)
    r = minimize(f, x0, method="Nelder-Mead", options=dict(maxfev=MAXFEV, xatol=1e-4, fatol=1e-6, adaptive=True))
    x = r.x / np.linalg.norm(r.x)
    d = {"seed": seed, "g_start": -f(x0), "g_final": -r.fun, "nfev": r.nfev, "secs_opt": time.time() - t0, "x0": x0.tolist()}
    d.update(search.diagnostics(x, Ls=(10, 12, 14)))
    try:
        d["S16_trunc"] = s16(x)
    except Exception as ex:
        d["S16_trunc"] = str(ex)
    d["secs"] = time.time() - t0
    with open(OUT, "a") as fh:
        fh.write(json.dumps(d) + "\n")
    return seed

if __name__ == "__main__":
    with Pool(NW) as p:
        for s in p.imap_unordered(task, range(NW)):
            print("done", s, flush=True)
