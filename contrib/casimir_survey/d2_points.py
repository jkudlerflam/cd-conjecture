#!/usr/bin/env python3
"""A/(2 pi ||h||) for the d = 2 candidates of ../d2_classification and ../d2_followup_220_1282,
from their stored ring energies E_0(L), L = 8..24(26), and VUMPS e_0 where available.
c = 1 is taken from those folders (floating spirals). For a sliding mode with incommensurate
wavevector q the ring acts as a twist 2 pi q L, so
    E_0(L) - e_0 L = -pi v/(6L) + (2 pi^2 rho/L) min_m (qL - m)^2,
which we fit for (v, rho) [and e_0 when no VUMPS value exists] with q scanned on a grid.
The bound v <= v_LR = 2e||h|| gives A/(2 pi ||h||) <= e/pi = 0.865 for any c = 1 chain.
"""
import json, os
import numpy as np
from core import hnorm
PAULI = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.diag([1.0 + 0j, -1])]
HERE = os.path.dirname(os.path.abspath(__file__))


def bond(J, f):
    h = np.zeros((4, 4), complex)
    for a in range(3):
        h += f[a] * 0.5 * (np.kron(PAULI[a], np.eye(2)) + np.kron(np.eye(2), PAULI[a]))
        for b in range(3):
            h += J[a][b] * np.kron(PAULI[a], PAULI[b])
    return h


def fit(Ls, E, e0=None, qgrid=np.linspace(0, 0.5, 2001)):
    Ls = np.asarray(Ls, float); E = np.asarray(E)
    best = None
    for q in qgrid:
        dist = np.abs(q * Ls - np.round(q * Ls))
        cols = [-np.pi / (6 * Ls), 2 * np.pi ** 2 * dist ** 2 / Ls]
        y = E.copy()
        if e0 is None:
            cols.insert(0, Ls)
        else:
            y = y - e0 * Ls
        X = np.vstack(cols).T
        coef, res, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = np.sum((X @ coef - y) ** 2)
        if best is None or r < best[0]:
            best = (r, q, coef)
    r, q, coef = best
    v = coef[-2]; rho = coef[-1]
    return dict(q=q, v=v, rho=rho, rms=np.sqrt(r / len(Ls)), e0=(coef[0] if e0 is None else e0))


vumps_e0 = {}
for name, fn in (('220', 'fes_220_n11_damped.json'), ('1282', 'fes_1282_n10_damped.json')):
    r = json.load(open(os.path.join(HERE, '..', 'd2_followup_220_1282', 'data', fn)))
    vumps_e0[name] = min(x['e'] for x in r['results'])
params = json.load(open(os.path.join(HERE, '..', 'd2_classification', 'candidates_params.json')))
rows = []
for name, p in params.items():
    path = os.path.join(HERE, '..', 'd2_classification', 'results_cluster', name + '.json')
    if not os.path.exists(path):
        path = os.path.join(HERE, '..', 'd2_classification', 'results', name + '.json')
    if not os.path.exists(path):
        continue
    sizes = json.load(open(path))['sizes']
    if name in ('220', '1282'):
        s2 = json.load(open(os.path.join(HERE, '..', 'd2_followup_220_1282', 'data', 'ed', name + '_k6.json')))['sizes']
        sizes.update(s2)
    Ls = sorted(int(L) for L in sizes if int(L) >= 10)
    E = [sizes[str(L)]['E0'] for L in Ls]
    h = bond(p['J'], p['field'])
    hn = hnorm(h, 2)
    f = fit(Ls, E, vumps_e0.get(name))
    Aeff = [-6 * L * (e - f['e0'] * L) / np.pi for L, e in zip(Ls, E)]
    rows.append((name, hn, f, max(Aeff)))
    print('%-10s ||h||=%.4f (unopt %.4f)  q=%.4f v=%.3f rho=%.3f rms=%.1e e0=%.7f%s  R=v/(2pi||h||)=%.3f  max_L A_eff/(2pi||h||)=%.3f'
          % (name, hn, np.ptp(np.linalg.eigvalsh(h)) / 2, f['q'], f['v'], f['rho'], f['rms'], f['e0'],
             ' (VUMPS)' if name in vumps_e0 else '', f['v'] / (2 * np.pi * hn), max(Aeff) / (2 * np.pi * hn)))
