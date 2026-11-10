"""C191: KT flow of (K_s, y1) from the bare values in predict.json / predict.py.
Flow (Giamarchi normalization, x(O1) = 2 K_s, y1 = g1/(pi u)):  dK/dl = -K^2 y^2 / 2,  dy/dl = (2 - 2K) y.
Bare scale: l = 0 at a0 = 1/k_F (lattice units).  Reports K_s at L = 128, 256 and, if the flow runs away,
the scale l* where y = 1 (gap Delta ~ E_F e^{-l*}, xi ~ a0 e^{l*})."""
import numpy as np
from predict import run

def flow(K0, y0, lmax=40, dl=1e-3):
    K, y, l = K0, y0, 0.0
    traj = []
    while l < lmax:
        traj.append((l, K, y))
        if abs(y) > 1.0:
            return np.array(traj), l
        dK = -0.5 * K * K * y * y; dy = (2 - 2 * K) * y
        K += dK * dl; y += dy * dl; l += dl
    return np.array(traj), None

for n in (0.05, 1/16, 0.0703, 3/32, 1/8, 5/32, 3/16):
    o = run(n)
    kF = np.pi * n / 2; a0 = 1 / kF; EF = None
    for r in o['rows']:
        if r['method'] != 'hybrid-ABA': continue
        traj, lstar = flow(r['Ks'], o['y1'])
        out = f"n={n:.4f} g4={r['g4']:6s} K_s0={r['Ks']:.4f} y1={o['y1']:.4f} "
        if lstar is None:
            def K_at(L):
                l = np.log(L / a0); i = min(int(l / 1e-3), len(traj) - 1); return traj[i, 1]
            out += f"-> fixed line K_s*={traj[-1,1]:.4f}; K_s(128)={K_at(128):.5f} K_s(256)={K_at(256):.5f} drift={K_at(256)-K_at(128):+.1e}"
        else:
            out += f"-> runs away at l*={lstar:.2f}: xi ~ {a0*np.exp(lstar):.3g} sites"
        print(out)
