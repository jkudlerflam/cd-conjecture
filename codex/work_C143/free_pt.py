"""Free-fermion locus F: f(k) = -a sin k + sqrt((c cos k + B)^2 + g^2 sin^2 k), quasiparticle energy E = 2f.
a = J_xy - J_yx, c = J_xx + J_yy, g = |(J_xx - J_yy) - i(J_xy + J_yx)|, B = uniform z field (per site).
main(): Pokrovsky-Talapov (PT) points a_c(c, B, g) where min_k f = 0 at incommensurate k0, with f''(k0) > 0 (z = 2).
main2(): ED gaps across a_c in the ordered region (first output, kept in free_pt.out; the low gap there is the
Kitaev doublet splitting, so the clean dispersion check is free_pt2.py).
"""
import numpy as np
from scipy.optimize import brentq
from ed import bond, spectrum, all_levels


def f(k, a, c, B, g):
    return -a * np.sin(k) + np.sqrt((c * np.cos(k) + B) ** 2 + g ** 2 * np.sin(k) ** 2)


ks = np.linspace(-np.pi, np.pi, 200001)


def fmin(a, c, B, g):
    v = f(ks, a, c, B, g)
    i = np.argmin(v)
    return v[i], ks[i]


def spin_bond(a, c, B, g):
    """J_xx = (c+g)/2, J_yy = (c-g)/2, J_xy = a/2, J_yx = -a/2, field B along z."""
    J = np.zeros((3, 3)); J[0, 0] = (c + g) / 2; J[1, 1] = (c - g) / 2; J[0, 1] = a / 2; J[1, 0] = -a / 2
    return bond(J, (0, 0, B))


def main():
    print("PT points: slice c = 1, g = 0.5, several B")
    for B in [0.0, 0.3, 0.6, 1.5]:
        ac = brentq(lambda a: fmin(a, 1.0, B, 0.5)[0], 1e-6, 5.0)
        m, k0 = fmin(ac, 1.0, B, 0.5)
        h = 1e-3
        f2 = (f(k0 + h, ac, 1, B, .5) - 2 * f(k0, ac, 1, B, .5) + f(k0 - h, ac, 1, B, .5)) / h ** 2
        print(f"  B={B}: a_c={ac:.6f}, k0={k0:.5f} (k0/pi={k0/np.pi:.5f}), f(k0)={m:.1e}, f''(k0)={f2:.4f}")


def main2():
    B, g = 0.3, 0.5
    ac = brentq(lambda a: fmin(a, 1.0, B, g)[0], 1e-6, 5.0)
    for L in [10, 12, 14]:
        row = []
        for a in [0.6 * ac, 0.9 * ac, ac, 1.1 * ac, 1.4 * ac]:
            w, k = all_levels(spectrum(spin_bond(a, 1.0, B, g), L))
            w = np.sort(w)
            row.append(w[1] - w[0])
        print(f"  L={L}: a/a_c = 0.6,0.9,1,1.1,1.4 -> gaps " + ", ".join(f"{x:.4f}" for x in row))


if __name__ == "__main__":
    main(); main2()
