"""Numerical check of the f-sum (Drude) chain of inequalities used in the proof:

  NC:  sum_F |v_F|  <=  2 pi <-h_hop>_GS  <=  2 pi ||B||_1 ,  so  A <= pi ||B||_1 <= pi M.
  MJ:  sum_cross |v| <= 4 pi <-h_hop>_GS  <=  2 pi ||B||_1 ,  so  A <= (pi/2) ||B||_1 <= pi M.

<h_hop>_GS is the ground-state expectation of the hopping part of one bond, computed
from the Bloch negative-energy projector:  NC: int dk/2pi tr[P(k)(zB+z^*B^dag)];
MJ: (1/2) int dk/2pi tr[P(k) i(zB - z^* B^T)].
"""
import numpy as np
from ff_bound import velocities, trace_norm, random_nc, random_mj, bloch_nc, bloch_mj

def hop_energy(kind, A0, B, order=64):
    """GS hopping energy per bond, integrated piecewise between Fermi momenta (P(k) is smooth there)."""
    from ff_bound import fermi_points
    h = bloch_nc(A0, B) if kind == "nc" else bloch_mj(A0, B)
    kf = sorted(k for k, _ in fermi_points(kind, A0, B))
    if not kf:
        kf = [-np.pi]
    edges = kf + [kf[0] + 2 * np.pi]
    xg, wg = np.polynomial.legendre.leggauss(order)
    tot = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        for x, w in zip(xg, wg):
            k = (a + b) / 2 + (b - a) / 2 * x
            e, U = np.linalg.eigh(h(k))
            P = U[:, e < 0] @ U[:, e < 0].conj().T
            z = np.exp(1j * k)
            hop = (B * z + B.conj().T / z) if kind == "nc" else 1j * (B * z - B.T / z)
            tot += w * (b - a) / 2 * np.trace(P @ hop).real
    tot /= 2 * np.pi
    return tot if kind == "nc" else tot / 2

rng = np.random.default_rng(3)
worst = {"nc": (0, 0), "mj": (0, 0)}
for kind, sizes, pref in (("nc", (1, 2, 3), 2 * np.pi), ("mj", (1, 2, 3, 4), 4 * np.pi)):
    for n in sizes:
        for _ in range(300):
            A0, B = random_nc(n, rng) if kind == "nc" else random_mj(n, rng)
            sv = velocities(kind, A0, B).sum()
            K = -hop_energy(kind, A0, B)
            r1 = sv / (pref * K) if K > 1e-12 else 0.0
            r2 = sv / (2 * np.pi * trace_norm(B))
            worst[kind] = (max(worst[kind][0], r1), max(worst[kind][1], r2))
    print(f"{kind}: max of  sum|v| / (prefactor * <-h_hop>) = {worst[kind][0]:.5f}  (must be <= 1);"
          f"  max of sum|v| / (2 pi ||B||_1) = {worst[kind][1]:.5f}  (must be <= 1; XX gives 2/pi = {2/np.pi:.4f})")

# exactness for decoupled bands (no interband term): XX at generic filling
A0 = np.array([[0.7]], dtype=complex); B = -np.ones((1, 1), dtype=complex)
sv = velocities("nc", A0, B).sum(); K = -hop_energy("nc", A0, B)
print("XX at mu=0.7: sum|v| =", round(sv, 6), "  2 pi <-h_hop> =", round(2 * np.pi * K, 6), " (equal: no interband term)")
