"""C186 route A: parity-class cusp bond beyond second order. des Cloizeaux coefficients of the full open chain
(L = 11) in units of the measured NNN hop: v = V2/|t2| + 2, g = -V3/|t2|, V4/|t2|, y3 = (dN=2, separation 3)/|t2|."""
import numpy as np
from sw_general import des_cloizeaux, measure
L = 11
for VD in (8.0, 10.0, 14.0, 20.0, 40.0):
    V = VD; p = 1.0
    t2 = 1.0 / (2 * V)
    w = 1j * 0.5 * t2           # u = 0.5, phi = pi/2
    mu = -2.1 * t2 + 3.0 / V
    He, pos, gap, smin = des_cloizeaux(L, V, w, p, 0.0, 0.0, mu, 0.0, False)
    m = measure(He, pos, L, False)
    T2 = abs(m["hop2"])
    print(f"V/Delta={VD:5.1f}  |t2|meas/t2={T2/t2:.5f}  v={m['V2'].real/T2+2:+.5f}  g={-m['V3'].real/T2:.5f}"
          f"  V4/|t2|={m['V4'].real/T2:+.5f}  y3={abs(m['dN2 sep3'])/T2:.5f} (4/3 Delta/V = {4/3/VD:.5f})"
          f"  hop1/|t2|={complex(m['hop1 (b+_{j+1} b_j)'])/T2:.4f}  gap={gap:.2f}", flush=True)
