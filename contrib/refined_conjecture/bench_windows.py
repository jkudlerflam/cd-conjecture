"""Conformal windows of the benchmark chains from data/bench_entropy.npz.
Uses the larger ring and only temperatures where the two ring lengths agree to 1%."""
import json
import numpy as np
z = np.load("data/bench_entropy.npz"); T = z["T"]
C = json.load(open("data/lattice_constants.json"))
specs = {"Ising": ("Ising_inf", None), "XX": ("XX_inf", None), "XXZ_D0.5": ("XXZ_L14", "XXZ_L12"),
         "Potts3": ("Potts_L9", "Potts_L8"), "TB_spin1": ("TB_L9", "TB_L8")}
out = {}
for name, (k1, k0) in specs.items():
    row = C[name]; c, d, hn = row["c"], row["d"], row["h_norm"]
    v = row["v_exact"] if row["v_exact"] else row["v_fit"]
    s = z[k1]
    conv = np.ones_like(T, bool) if k0 is None else np.abs(z[k1] / z[k0] - 1) < 0.01
    i0 = np.argmax(conv)  # first converged temperature
    r = s / (np.pi * c * T / (3 * v))
    for eps in (0.1,):
        ok = np.abs(r[i0:] - 1) <= eps
        if not ok[0]:
            Ts = np.nan
        else:
            bad = np.nonzero(~ok)[0]
            Ts = T[i0:][bad[0] - 1] if len(bad) else T[-1]
        ok1 = r[i0:] >= 1 - eps
        bad1 = np.nonzero(~ok1)[0]
        T1 = T[i0:][bad1[0] - 1] if len(bad1) and ok1[0] else np.nan
    sT = np.interp(Ts, T, s)
    out[name] = dict(T_conv=float(T[i0]), Tstar=float(Ts), Tstar_onesided=float(T1), sigma=float(sT / np.log(d)), ell=float(v / Ts), R_T=float(hn / Ts), max_ratio=float(r[i0:].max()))
    print(name, {k: round(x, 3) for k, x in out[name].items()})
json.dump(out, open("data/bench_windows.json", "w"), indent=1)
