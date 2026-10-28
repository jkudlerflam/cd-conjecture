"""Conformal window of the SU(d) ULS ring from the Schur-Weyl spectra.

For each d and ring length L = 0 mod d available in data/, compares the lattice entropy per site
with (a) the SU(d)_1 torus entropy at the same L (finite-size CFT, v = 2 pi/d), and (b) the
thermodynamic CFT line pi c T/(3v). Reports, for eps in EPS:
  T2  two-sided window: |s_lat/s_ref - 1| <= eps on [T_lo, T]
  T1  one-sided window: s_lat >= (1 - eps) s_ref on [T_lo, T]   (what the entropy bound uses)
with T_lo = 2 v/L (below that the ring is in its finite-size regime), and the derived
sigma* = s(T*)/log d, l* = v/T*, R_T = ||h||/T* (||h|| = 1).
Usage: python3 analyze_uls.py [dmax]
"""
import sys, glob, re
import numpy as np
from thermo import load_uls, uls_levels, entropy_from_levels, sun1_torus_entropy

EPS = (0.1, 0.2)
Ts = np.round(np.arange(0.01, 3.0, 0.005), 4)


def lengths():
    Ls = set()
    for f in glob.glob("data/uls_L*.npz"):
        Ls.add(int(re.search(r"uls_L(\d+)", f).group(1)))
    return sorted(Ls)


def first_fail(ok):
    bad = np.nonzero(~ok)[0]
    return None if len(bad) == 0 else bad[0]


def win(Ts, ok):
    if not ok[0]:
        return np.nan
    b = first_fail(ok)
    return Ts[-1] if b is None else Ts[b - 1]


def main():
    dmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    rows = []
    for L in lengths():
        spec = load_uls(L)
        exact = all(e for (_, _, e) in spec.values())
        rmax = max(len(l) for l in spec)
        for d in range(2, dmax + 1):
            if L % d or d > rmax:
                continue
            v, c = 2 * np.pi / d, d - 1
            E, g = uls_levels(spec, d)
            assert abs(g.sum() / d ** L - 1) < 1e-6, (L, d, g.sum() / d ** L)
            s_lat = entropy_from_levels(E, g, Ts) / L
            s_tor, _ = sun1_torus_entropy(d, L, v, Ts)  # finite-size CFT reference
            s_inf = np.pi * c * Ts / (3 * v)
            m = Ts >= 2 * v / L
            T = Ts[m]
            out = dict(d=d, L=L, exact=exact, Tlo=T[0])
            for eps in EPS:
                for name, ref in (("tor", s_tor[m]), ("inf", s_inf[m])):
                    r = s_lat[m] / ref
                    out[f"T2_{name}_{eps}"] = win(T, np.abs(r - 1) <= eps)
                    out[f"T1_{name}_{eps}"] = win(T, r >= 1 - eps)
            rows.append(out)
            Tst = out["T1_tor_0.1"]
            sst = np.interp(Tst, Ts, s_lat) if np.isfinite(Tst) else np.nan
            print(f"d={d} L={L:2d} exact={int(exact)} Tlo={T[0]:.3f} | "
                  + " ".join(f"{k}={out[k]:.3f}" for k in out if k.startswith("T") and k != "Tlo")
                  + f" | T1tor0.1: sigma*={sst/np.log(d):.3f} l*={v/Tst:.2f} R_T={1/Tst:.2f}", flush=True)
    np.save("data/uls_windows.npy", rows, allow_pickle=True)


if __name__ == "__main__":
    main()
