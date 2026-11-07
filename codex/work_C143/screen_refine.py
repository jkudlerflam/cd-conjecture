"""Two-stage search for singular qubit bonds (version 3).

Stage 1 (screen): sample random bonds uniformly on the unit sphere of the 12 coefficients (J, B) and compute
the v2 indicator g = min over increments of log S_L at L = 8 -> 10 -> 12 (kappa = 3).
Stage 2 (refine): Nelder-Mead from the NTOP best screened bonds, maximizing g at the lengths LREF
(8,10,12 locally, 10,12,14 on Typhon), then full diagnostics.
Random Nelder-Mead starts (version 2) stall on the flat g = 0 plateau of gapped phases, hence the screen.

Usage: python3 screen_refine.py OUTPREFIX NWORKERS NSCREEN_PER_WORKER NTOP MAXFEV LREF(comma) [seed]
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize
import search

OUT = sys.argv[1]; NW = int(sys.argv[2]); NS = int(sys.argv[3]); NTOP = int(sys.argv[4])
MAXFEV = int(sys.argv[5]); LREF = tuple(int(a) for a in sys.argv[6].split(",")); SEED = int(sys.argv[7]) if len(sys.argv) > 7 else 0


def gL(x, Ls):
    logs = [np.log(search.softcount(x, L)) for L in Ls]
    return float(np.min(np.diff(logs))), [float(np.exp(v)) for v in logs]


def screen(w):
    rng = np.random.default_rng(SEED * 1000 + w)
    out = []
    for i in range(NS):
        x = rng.normal(size=12); x /= np.linalg.norm(x)
        g, S = gL(x, (8, 10, 12))
        out.append({"x": x.tolist(), "g": g, "S": S})
    return out


def refine(item):
    x0 = np.array(item["x"])
    t0 = time.time()
    f = lambda x: -gL(x / np.linalg.norm(x), LREF)[0]
    r = minimize(f, x0, method="Nelder-Mead", options=dict(maxfev=MAXFEV, xatol=1e-4, fatol=1e-6, adaptive=True))
    x = r.x / np.linalg.norm(r.x)
    d = {"g_screen": item["g"], "g_start_ref": -f(x0), "g_final": -r.fun, "nfev": r.nfev, "LREF": LREF, "x0": x0.tolist()}
    d.update(search.diagnostics(x, Ls=(8, 10, 12) if max(LREF) <= 12 else (10, 12, 14)))
    d["secs"] = time.time() - t0
    with open(OUT + "_refined.jsonl", "a") as fh:
        fh.write(json.dumps(d) + "\n")
    return d["g_final"]


if __name__ == "__main__":
    with Pool(NW) as p:
        scr = [r for chunk in p.map(screen, range(NW)) for r in chunk]
    scr.sort(key=lambda r: -r["g"])
    with open(OUT + "_screen.jsonl", "w") as fh:
        for r in scr:
            fh.write(json.dumps(r) + "\n")
    gs = np.array([r["g"] for r in scr])
    print("screened", len(scr), "quantiles of g (50, 90, 99, max):", np.quantile(gs, [0.5, 0.9, 0.99]).round(3), gs.max().round(3), flush=True)
    top = scr[:NTOP]
    with Pool(NW) as p:
        for g in p.imap_unordered(refine, top):
            print("refined g =", round(g, 3), flush=True)
