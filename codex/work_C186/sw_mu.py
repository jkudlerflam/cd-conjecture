"""C186: chemical-potential dependence of the all-orders (des Cloizeaux) cusp coefficients in the parity class,
open L = 11, u = 1, phi = pi/2. Compares mu_eff = -12 (C177's two-body fit) with mu_eff = -2.25 (the dilute onset).
Units: the measured NNN hop. D = 2 V2 V3 + 4 V2 + V3 u^2 + 4 V3 + 8 (C175's resonance function, u = |hop1|/|hop2|)."""
import numpy as np
from sw_general import des_cloizeaux, measure
L = 11
for VD in (8.0, 16.0, 40.0):
    for mu_eff in (-12.0, -2.25):
        V = VD; t2 = 1.0 / (2 * V)
        w = 1j * t2
        mu = mu_eff * t2 + 3.0 / V
        He, pos, gap, smin = des_cloizeaux(L, V, w, 1.0, 0.0, 0.0, mu, 0.0, False)
        m = measure(He, pos, L, False)
        T2 = abs(m["hop2"]); u = abs(m["hop1 (b+_{j+1} b_j)"]) / T2
        V2 = m["V2"].real / T2; V3 = m["V3"].real / T2; V4 = m["V4"].real / T2
        D = 2 * V2 * V3 + 4 * V2 + V3 * u * u + 4 * V3 + 8
        print(f"V/Delta={VD:4.0f} mu_eff={mu_eff:6.2f} |t2|/t2={T2/t2:.4f} u={u:.4f} V2={V2:+.4f} (v={V2+2:+.4f}, v(V/D)^2={(V2+2)*VD**2:6.2f}) "
              f"V3={V3:+.4f} V4={V4:+.4f} D={D:+.3f} y3={abs(m['dN2 sep3'])/T2:.4f} gap={gap:.2f}", flush=True)
