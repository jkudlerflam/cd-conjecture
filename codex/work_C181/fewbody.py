"""C181: one-, two- and three-body data of model (3.2) off the symmetric line phi = pi/2.

For each |w| and phi: valley minima (k_a, eps_a) < (k_b, eps_b), offset delta = eps_b - eps_a, valley inverse
masses eps''; the lowest two-body level over total momentum K (the ab pair) with its binding relative to
eps_a + eps_b and its inverse mass; the lowest three-body level near K_pair + k_a (the aab trimer) on the infinite
line in gap coordinates with hard walls at gap R (nbody_box.py, C175), and its binding B3 relative to pair + a atom.
Usage: python3 fewbody.py <wa> <R3> <phi1,phi2,...>
"""
import numpy as np, scipy.linalg as sl, sys, time
from scat import valleys, eps, rel_coeffs, Vpot
import nbody_box


def E2(w, K, Lr=1500, V2=-2.0, V3=-2/3):
    c = rel_coeffs(K, w); n = Lr
    ab = np.zeros((3, n), complex)
    ab[0, :] = c[0] + np.array([{2: V2, 3: V3}.get(r, 0.0) for r in range(2, n + 2)])
    ab[1, :-1] = c[-1]; ab[2, :-2] = c[-2]
    return sl.eig_banded(ab, lower=True, select='i', select_range=(0, 0), eigvals_only=True)[0]


def epp(k, w, h=1e-4):
    return (eps(k + h, w) - 2 * eps(k, w) + eps(k - h, w)) / h ** 2


def pair_min(w, K0, Lr=1500):
    # golden-section minimisation of E2(K) around K0, then curvature
    a, b = K0 - 0.05, K0 + 0.05
    g = (np.sqrt(5) - 1) / 2
    for _ in range(40):
        c1 = b - g * (b - a); c2 = a + g * (b - a)
        if E2(w, c1, Lr) < E2(w, c2, Lr):
            b = c2
        else:
            a = c1
    K = (a + b) / 2
    h = 2e-3
    e0 = E2(w, K, Lr)
    curv = (E2(w, K + h, Lr) - 2 * e0 + E2(w, K - h, Lr)) / h ** 2
    return K, e0, curv


def trimer(w, Kc, R, dK=0.01):
    es = [nbody_box.lowest(3, Kc + s * dK, w, R)[0] for s in (-1, 0, 1)]
    # parabolic minimum
    a = (es[0] - 2 * es[1] + es[2]) / (2 * dK ** 2)
    bq = (es[2] - es[0]) / (2 * dK)
    Kmin = Kc - bq / (2 * a) if a > 0 else Kc
    emin = es[1] - bq ** 2 / (4 * a) if a > 0 else min(es)
    return Kmin, emin, 2 * a


if __name__ == "__main__":
    wa = float(sys.argv[1]); R3 = int(sys.argv[2])
    phis = [float(x) for x in sys.argv[3].split(',')]
    print("# |w| phi k_a k_b delta 1/m_a 1/m_b | K_p E_pair Eb=(ea+eb-Ep) 1/M_p | K_3 E_3 B3=(Ep+ea-E3) 1/M_3 | secs", flush=True)
    for phi in phis:
        t0 = time.time()
        w = wa * np.exp(1j * phi)
        (ka, ea), (kb, eb) = valleys(w)[:2]
        Kp, Ep, Mp = pair_min(w, ka + kb)
        if R3 > 0:
            K3, E3, M3 = trimer(w, Kp + ka, R3)
        else:
            K3, E3, M3 = np.nan, np.nan, np.nan
        print(f"{wa:.2f} {phi:.4f} {ka:+.5f} {kb:+.5f} {eb-ea:.5f} {epp(ka,w):.4f} {epp(kb,w):.4f} | "
              f"{Kp:+.5f} {Ep:.6f} {ea+eb-Ep:.6f} {Mp:.4f} | {K3:+.5f} {E3:.6f} {Ep+ea-E3:.6f} {M3:.4f} | {time.time()-t0:.0f}",
              flush=True)
