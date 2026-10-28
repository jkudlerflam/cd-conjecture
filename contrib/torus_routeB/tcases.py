#!/usr/bin/env python3
"""The effective-theory test cases of contrib/torus_routeB (task 1).

N = 4 flavours (r = 2), Fermi points symmetric about zero (+-k_a), velocity 1 for every flavour.
  K_DIST  distinct generic momenta in [0.3 pi, 0.7 pi], chosen to maximize the momentum margins
  K_EQ    equal momenta
Cases (each with a sign s = +-1 multiplying the "dangerous" couplings):
  'i'    U(1)-only: classes U V J R P X, s multiplies P (pair hopping), momenta K_DIST
  'ip'   parity-only: case 'i' plus quadratic pairing f_{0a} f_{1b} (breaks total number)
  'ii'   torus: classes U V J R, s multiplies J and R, momenta K_DIST
  'iii'  torus, equal momenta: classes U V J R, s multiplies J and R, momenta K_EQ
  'free' / 'freeeq': decoupled free chains (calibration, exact c = 4)
The random couplings are drawn once per seed, so the cases share the same U, V, J, R values.
no=True subtracts the Hartree shift of the U and V terms at the free densities k_a/pi.

Written by Claude (Opus 5.5) for contrib/torus_routeB, 2026-10-27.
"""
import numpy as np
import fermions as F

K_DIST = np.pi * np.array([0.30, 0.375, 0.41, 0.555])   # margins: all 4-field sums 0.035 pi, torus-allowed 0.07 pi
K_EQ = np.pi * np.array([0.40, 0.40, 0.40, 0.40])


def normal_order_density(cl, nbar):
    """Hartree subtraction: n_a n_b -> (n_a - nbar_a)(n_b - nbar_b) + const for the U and V classes,
    so that the free densities nbar = k/pi stay self-consistent at mean-field level."""
    out = []
    for cf, m in cl['U'] + cl['V']:
        a, b = m[0][2], m[2][2]
        out.append((-cf * nbar[b], ((1, 0, a), (0, 0, a))))
        out.append((-cf * nbar[a], ((1, 0, b), (0, 0, b))))
    return out


def case_monos(case, s=1.0, g=0.3, seed=0, gpair=None, no=False):
    rng = np.random.default_rng(seed)
    cl = F.interaction_monos(rng, g, 'UVJRPX')
    pair = F.pairing_monos(np.random.default_rng(seed + 1000), gpair if gpair is not None else g / 2)
    if case in ('free', 'freeeq'):
        k = K_DIST if case == 'free' else K_EQ
        return k, F.free_monos(k), []
    if case in ('i', 'ip'):
        k = K_DIST
        inter = (cl['U'] + cl['V'] + cl['J'] + cl['R'] + cl['X']
                 + [(s * c, m) for c, m in cl['P']])
        if case == 'ip':
            inter = inter + pair
    elif case in ('ii', 'iii'):
        k = K_DIST if case == 'ii' else K_EQ
        inter = cl['U'] + cl['V'] + [(s * c, m) for c, m in cl['J'] + cl['R']]
    else:
        raise ValueError(case)
    quad = normal_order_density(cl, k / np.pi) if no else []
    return k, F.free_monos(k) + quad + inter, inter


def two_site_term(case, s=1.0, g=0.3, seed=0, no=False):
    k, monos, inter = case_monos(case, s, g, seed, no=no)
    h = F.terms_matrix(monos, 2)
    return k, (h + h.conj().T) / 2, inter
