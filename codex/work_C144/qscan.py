"""ED of the doubling ratio Q_L(h) = L (E_{2L} - 2 E_L) / M(h) near XZ, for L = 2..Lmax (rings L and 2L).
Directions are bond matrices added to h0 = X(x)Z with amplitude s, normalized so that M(direction) = 1.
Usage: python3 qscan.py Lmax label:spec ...   spec = comma list of name=coef (names as in xzlib.term)."""
import sys, time, json
import numpy as np
from xzlib import term, gauge_M, ground_energy

def build(spec):
    h = np.zeros((4, 4), complex)
    for item in spec.split(","):
        n, c = item.split("=")
        if n == "rand":  # random direction in the 12-dim (J, B) tangent space, seed c
            from xzlib import bond
            g = np.random.default_rng(int(c)).standard_normal(12)
            h += bond(g[:9].reshape(3, 3), g[9:])
        else:
            h += float(c) * term(n)
    return h

Lmax = int(sys.argv[1])
h0 = term("xz")
for arg in sys.argv[2:]:
    label, spec, amps = arg.split(":")
    u = build(spec) if spec != "none" else np.zeros((4, 4), complex)
    Mu = gauge_M(u) if spec != "none" else 1.0
    if spec != "none":
        u = u / Mu
    for s in [float(a) for a in amps.split("/")]:
        h = h0 + s * u
        M = gauge_M(h)
        t0 = time.time()
        E = {}
        Qs = []
        for L in range(2, Lmax + 1):
            for LL in (L, 2 * L):
                if LL not in E:
                    E[LL] = float(ground_energy(h, LL))
            Qs.append(L * (E[2 * L] - 2 * E[L]) / M)
        print(json.dumps(dict(label=label, s=s, M=M, Q={L: round(q, 5) for L, q in zip(range(2, Lmax + 1), Qs)},
                              E={k: E[k] for k in sorted(E)}, sec=round(time.time() - t0, 1))), flush=True)
