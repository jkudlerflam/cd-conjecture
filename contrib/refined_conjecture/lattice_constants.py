"""Microscopic constants of the benchmark chains: ||h|| (half spread), ||j|| (energy current
i[h12,h23]), the Casimir coefficient A = sum_a c_a v_a from E_0(L)/L = e_0 - pi A/(6 L^2) + b/L^4,
and, for ULS, the f-sum ratio k v dim g / (2 pi <-sum_a [t^a,[t^a,h]]>) = (d^2-1)/(d (1 - d e_0)).
Output: printed table and data/lattice_constants.json.
"""
import glob
import json
import sys
import numpy as np
from scipy.special import digamma
from models_ed import models, half_spread, current_norm, ground_energy
from thermo import load_uls

Ls = {"Ising": range(8, 19, 2), "XX": range(8, 19, 4), "XXZ_D0.5": range(8, 19, 2),
      "Potts3": range(6, 13), "TB_spin1": range(6, 13, 2)}


def fit_casimir(L, E):
    L = np.asarray(L, float); y = np.asarray(E) / L
    X = np.vstack([np.ones_like(L), -np.pi / (6 * L ** 2), 1 / L ** 4]).T
    coef = np.linalg.lstsq(X, y, rcond=None)[0]
    return coef[0], coef[1]


def main():
    only = sys.argv[1] if len(sys.argv) > 1 else ""
    try:
        out = json.load(open("data/lattice_constants.json"))
    except FileNotFoundError:
        out = {}
    for name, (h, d, c, v) in models().items():
        if only and not name.startswith(only):
            continue
        row = dict(d=d, c=c, v_exact=v, h_norm=half_spread(h), j_norm=current_norm(h, d))
        if name.startswith("ULS"):
            dd = d
            avail = sorted(int(f.split("uls_L")[1].split(".")[0].split("_")[0]) for f in glob.glob("data/uls_L*.npz"))
            LL = sorted(set(L for L in avail if L % dd == 0 and L >= max(6, 2 * dd)))
            E = []
            for L in LL:
                spec = load_uls(L)
                E.append(min(e.min() for lam, (e, w, ex) in spec.items() if len(lam) <= dd))
            e0_exact = 1 + 2 * (digamma(1 / dd) + np.euler_gamma) / dd
            row.update(e0_exact=e0_exact, fsum_ratio=(dd * dd - 1) / (dd * (1 - dd * e0_exact)))
        else:
            LL = list(Ls[name]) if name in Ls else []
            E = [ground_energy(h, d, L)[0] for L in LL]
        if len(LL) >= 3:
            e0, A = fit_casimir(LL, E)
            row.update(L_fit=list(map(int, LL)), e0_fit=e0, casimir_fit=A, v_fit=A / c)
        vv = v if v is not None else row.get("v_fit")
        row["casimir_ratio_2pi"] = c * vv / (2 * np.pi * row["h_norm"])
        out[name] = row
        print(name, {k: (round(x, 5) if isinstance(x, float) else x) for k, x in row.items() if k != "L_fit"}, flush=True)
    json.dump(out, open("data/lattice_constants.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
