"""C191: Fermi-surface two-body data at n = 1/8, 5/32, 3/16 (rho = n/2 per valley), |w| = 1, phi = pi/2."""
import numpy as np
from scat191 import Model, unlike_amplitudes, like_amplitude
m = Model(w=1j)
ka, kb = m.valleys()[:2]
print(f"# valleys ka={ka:.6f} kb={kb:.6f} ka+kb={ka+kb:.6f} eps_min={m.eps(ka):.6f}")
print("# n  pL  pR  vR  vL  eF-emin | like: delta_s (fermion ref) | unlike: |r| arg(t) arg(tp) |t| arg(-t)  unit")
for n in (0.0625, 0.09375, 0.125, 0.140625, 5/32, 0.171875, 3/16, 0.203125, 0.21875, 0.25):
    rho = n/2
    pL, pR, eF = m.fermi_points(ka, rho)
    vR, vL = m.deps(ka+pR), m.deps(ka+pL)
    Lk = like_amplitude(m, ka, pL, pR)
    U = unlike_amplitudes(m, ka, kb, pL, pR)
    print(f"{n:.5f} {pL:+.4f} {pR:+.4f} {vR:.4f} {vL:+.4f} {eF-m.eps(ka):.4f} | {Lk['delta']:+.4f} nopen={Lk['nopen']} | "
          f"{abs(U['r']):.4f} {abs(U['rp']):.4f} {np.angle(U['t']):+.4f} {np.angle(U['tp']):+.4f} {abs(U['t']):.4f} {np.angle(-U['t']):+.4f} nopen={U['nopen']} {U['unit']:.1e}")
