"""Spectra of selected landing points at L = 7..14: soft count S_L, number of states below M/L, low gaps times L,
and the momenta of the lowest states (to diagnose incommensurate level crossings)."""
import json, sys
import numpy as np
import search
from ed import bond, spectrum, all_levels, params_to_JB, gauge_M

def probe(x, Ls):
    x = np.asarray(x) / np.linalg.norm(x)
    J, B = params_to_JB(x)
    M = gauge_M(bond(J, B))[0]
    for L in Ls:
        w, k = all_levels(spectrum(bond(J, B), L)); o = np.argsort(w); w = w[o]; k = k[o]
        e = w - w[0]
        print(f"   L={L:2d} S={search.softcount(x, L):7.3f}  n(<M/L)={int(np.sum(e < M / L)):3d}  "
              f"L*gaps={np.round(L * e[1:5], 3).tolist()}  k0..k4={k[:5].tolist()}  e0/L={w[0]/L:.6f}", flush=True)

if __name__ == "__main__":
    fn, idx = sys.argv[1], int(sys.argv[2])
    Ls = [int(a) for a in sys.argv[3].split(",")]
    d = [json.loads(l) for l in open(fn)][idx]
    print(fn, idx)
    probe(d["x"], Ls)
