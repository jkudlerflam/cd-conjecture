"""Single-quasiparticle dispersion on the gapped side of a landing point, along the cut x(s) = (J, (1+s) B):
lowest excitation in every momentum sector, E(k) - E0, at L = 12 and 14. A Pokrovsky-Talapov onset shows a
minimum at incommensurate k0 that closes quadratically as s decreases to s_c."""
import json, sys
import numpy as np
from ed import bond, spectrum, params_to_JB
fn, idx = sys.argv[1], int(sys.argv[2])
svals = [float(a) for a in sys.argv[3].split(",")]
d = [json.loads(l) for l in open(fn)][idx]
J, B = params_to_JB(np.array(d["x"]))
for s in svals:
    for L in [12, 14]:
        sp_ = spectrum(bond(J, (1 + s) * B), L, nev=3)
        E0 = min(w[0] for w in sp_.values())
        m0 = [m for m, w in sp_.items() if abs(w[0] - E0) < 1e-12][0]
        disp = []
        for m in range(L):
            w = sp_[m]
            e = w[1] if m == m0 else w[0]
            disp.append(e - E0)
        disp = np.array(disp)
        mm = int(np.argmin(disp))
        print(f"s={s:+.3f} L={L}: ground k-index {m0}; lowest excitation {disp.min():.4f} at k = 2pi*{mm}/{L} (k/pi={2*mm/L:.3f}); E(k)-E0 = {np.round(disp,3).tolist()}", flush=True)
