"""Floating-point band test of C32's mobile core at the leading order g^2,
for n beyond C32's n <= 2, using the validated structured symbol of
mobile_leading.py.

For each (n, D, psi_0): (i) tree positivity, the spectrum of the Hermitian part
of the phase-rotated hopping at t_d = 0; (ii) a loop t_d chosen so that the
loop-induced onsite matrix has norm kappa times the smallest tree eigenvalue
(normalized units), linearity checked by doubling t_d; (iii) common field equal
to half the smallest Hermitian eigenvalue (as in C32 band_checks.py); (iv) the
2n Fermi points, sea centres modulo pi, delta2 (smallest separation of Fermi
points), delta4 (smallest |p1 + p2 - p3 - p4| mod 2 pi over distinct unordered
pairs, repeated points allowed) and the smallest velocity.

All quantities are in units of max|T|. Status: numerical, double precision.
The smallest tree eigenvalue falls by about two decades per flavour, so the
narrow bands lose resolution near n = 7; the run records a precision flag.
numpy and scipy only; about a minute. Output: mobile_bands.json.
"""
import itertools
import json
import os

import numpy as np
from scipy.optimize import brentq

from mobile_leading import leading

KAPPA = 0.25


def analyse(n, D, psi0, m=1.0, grid=20001):
    T0, M0 = leading(n, m=m, D=D, td=0.0, psi0=psi0)
    sc = np.max(abs(T0))
    phase = np.angle(T0[0, 0])
    R0 = T0 * np.exp(-1j * phase) / sc
    imag_resid = float(np.max(abs(R0.imag)))
    gam = np.linalg.eigvalsh((R0 + R0.conj().T) / 2)
    rec = dict(n=n, D=D, psi0=psi0, m=m, T_scale=float(sc), tree_imag_residual=imag_resid,
               tree_onsite_over_T=float(np.max(abs(M0)) / sc),
               tree_eigs=gam.tolist(), tree_positive=bool(gam[0] > 0))
    if gam[0] <= 0:
        return rec
    rec['control_tree'] = bands(T0 / sc, M0 / sc, phase, n, grid)
    if n == 1:
        # one flavour: the onsite matrix is a 1x1 imaginary number, hence zero
        return rec
    tref = 1e-6 * gam[0] * D ** (-(n + 1))
    _, Mr = leading(n, m=m, D=D, td=tref, psi0=psi0)
    _, Mr2 = leading(n, m=m, D=D, td=2 * tref, psi0=psi0)
    lin = float(np.linalg.norm(Mr2 - 2 * Mr, 2) / np.linalg.norm(Mr, 2))
    td = tref * KAPPA * gam[0] / (np.linalg.norm(Mr, 2) / sc)
    T, M = leading(n, m=m, D=D, td=td, psi0=psi0)
    T, M = T / sc, M / sc
    rec.update(loop_td=float(td), loop_linearity=lin,
               loop_T_change=float(np.linalg.norm(T - T0 / sc, 2)))
    rec.update(bands(T, M, phase, n, grid))
    return rec


def bands(T, M, phase, n, grid):
    """Common field = half the smallest Hermitian eigenvalue; Fermi data."""
    R = T * np.exp(-1j * phase)
    g2 = np.linalg.eigvalsh((R + R.conj().T) / 2)
    field = g2[0] / 2
    mu = M + field * np.eye(n)

    def vals(k):
        return np.linalg.eigvalsh(mu + T * np.exp(1j * k) + T.conj().T * np.exp(-1j * k))

    ks = np.linspace(-np.pi, np.pi, grid)
    E = np.array([vals(k) for k in ks])
    roots, centres, vel, counts = [], [], [], []
    for a in range(n):
        rr = [brentq(lambda k: vals(k)[a], ks[j], ks[j + 1], xtol=1e-15)
              for j in range(grid - 1) if E[j, a] * E[j + 1, a] < 0]
        counts.append(len(rr))
        roots.extend(rr)
        if len(rr) == 2:
            centres.append(float(((rr[0] + rr[1]) / 2 + np.pi / 2) % np.pi - np.pi / 2))
        for k in rr:
            e, U = np.linalg.eigh(mu + T * np.exp(1j * k) + T.conj().T * np.exp(-1j * k))
            der = 1j * (T * np.exp(1j * k) - T.conj().T * np.exp(-1j * k))
            vel.append(float(abs(np.vdot(U[:, a], der @ U[:, a]))))

    def dist(x):
        return abs((x + np.pi) % (2 * np.pi) - np.pi)

    rec = dict(onsite_norm=float(np.linalg.norm(M, 2)), field=float(field),
               hermitian_eigs=g2.tolist(), roots_per_band=counts, centres_mod_pi=centres,
               min_velocity=min(vel) if vel else None)
    # a root moves by about (entry error)/(velocity); entries carry ~1e-15
    # relative error of max|T| = 1, so this bounds the resolvable momentum
    rec['root_error_estimate'] = float(20 * np.finfo(float).eps / min(vel)) if vel else None
    if all(c == 2 for c in counts):
        d2 = min(dist(a - b) for a, b in itertools.combinations(roots, 2))
        pairs = list(itertools.combinations_with_replacement(range(len(roots)), 2))
        sums = [roots[a] + roots[b] for a, b in pairs]
        d4 = min(dist(x - y) for x, y in itertools.combinations(sums, 2))
        cs = sorted(centres)
        dc = min((dist(2 * (a - b)) / 2 for a, b in itertools.combinations(cs, 2)), default=None)
        rec.update(delta2=float(d2), delta4=float(d4), min_centre_separation_mod_pi=dc,
                   margins_resolved=bool(d4 > 40 * rec['root_error_estimate']))
    return rec


if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    out = []
    for psi0 in (0.9, np.pi):
        for n in range(1, 8):
            for D in (100.0, 1000.0):
                r = analyse(n, D, psi0)
                keep = {k: r.get(k) for k in ('n', 'D', 'psi0', 'tree_positive', 'roots_per_band',
                                              'delta2', 'delta4', 'min_centre_separation_mod_pi',
                                              'min_velocity', 'root_error_estimate',
                                              'margins_resolved', 'loop_T_change')}
                keep['tree_min'] = r['tree_eigs'][0]
                ct = r.get('control_tree', {})
                keep['tree_control_delta4'] = ct.get('delta4')
                print(json.dumps(keep), flush=True)
                out.append(r)
    with open(os.path.join(here, 'mobile_bands.json'), 'w') as f:
        json.dump({'status': 'numerical, double precision, leading order g^2', 'kappa': KAPPA,
                   'runs': out}, f, indent=1)
        f.write('\n')
