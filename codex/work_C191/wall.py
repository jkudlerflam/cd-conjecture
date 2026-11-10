"""C191: single-particle reflection at an open end of model (3.2) (free hops w = i, NNN -1), energy near eps_F.
Half-line x >= 1; free modes z^x; boundary: psi(0) = psi(-1) = 0 for the free extension.
Incoming at a left wall = left movers (v<0); outgoing = right movers. Prints |R| in the valley basis."""
import numpy as np
from scat191 import Model
m = Model(w=1j)
ka, kb = m.valleys()[:2]
for n in (1/8, 3/16):
    pL, pR, eF = m.fermi_points(ka, n/2)
    E = eF
    # single-particle: sum_d a_d psi(x-d) = E psi(x)  -> modes z^x with sum_d a_d z^{-d} = E
    coeffs = {d: ad for d, ad in m.a.items()}   # term a_d z^{-d}
    # polynomial in z: multiply by z^2: sum_d a_d z^{2-d} - E z^2
    poly = np.zeros(5, complex)  # powers 4..0
    for d, ad in coeffs.items():
        poly[4 - (2 - d)] += ad
    poly[4 - 2] -= E
    z = np.roots(poly)
    ks = np.angle(z)
    v = np.array([m.deps(k) for k in ks])
    lab = []
    for k in ks:
        dists = {'aR': ka + pR, 'aL': ka + pL, 'bR': kb - pL, 'bL': kb - pR}
        lab.append(min(dists, key=lambda s: abs(np.angle(np.exp(1j*(k - dists[s]))))))
    inc = [i for i in range(4) if v[i] < 0]; out = [i for i in range(4) if v[i] > 0]
    B = lambda i: np.array([z[i]**0, z[i]**(-1)])   # psi(0), psi(-1)
    Mo = np.array([B(i) for i in out]).T
    R = np.zeros((2, 2), complex)
    for j, ii in enumerate(inc):
        sol = np.linalg.solve(Mo, -B(ii))
        R[:, j] = sol * np.sqrt(np.abs(v[out]) / abs(v[ii]))
    print(f"n={n}: incoming {[lab[i] for i in inc]}, outgoing {[lab[i] for i in out]}")
    print("  |R| =\n", np.round(np.abs(R), 4), " unitarity", np.linalg.norm(R.conj().T @ R - np.eye(2)))
