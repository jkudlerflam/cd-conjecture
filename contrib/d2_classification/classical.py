#!/usr/bin/env python3
"""Classical (unit-vector) ground state on a periodic ring, versus the quantum energy density.

e_cl: minimum of sum_j s_j^T J s_{j+1} + h.s_j over unit vectors on rings of N sites
(N = 60..64, several random starts, L-BFGS on angles), per site.
The ratio e_q/e_cl is 1 for a ferromagnet (product ground state) and grows with
quantum fluctuations (Heisenberg AF: 1.77).  The spiral wave vector q_cl is read
from the Fourier peak of s_j.
"""
import json, sys
import numpy as np
from scipy.optimize import minimize

def energy(x, J, h, N):
    th, ph = x[:N], x[N:]
    s = np.stack([np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)], 1)
    sn = np.roll(s, -1, 0)
    E = np.einsum('ia,ab,ib->', s, J, sn) + (s@h).sum()
    # gradient
    gs = sn@J.T + np.roll(s, 1, 0)@J + h  # dE/ds_i
    dth = np.stack([np.cos(th)*np.cos(ph), np.cos(th)*np.sin(ph), -np.sin(th)], 1)
    dph = np.stack([-np.sin(th)*np.sin(ph), np.sin(th)*np.cos(ph), 0*th], 1)
    g = np.concatenate([(gs*dth).sum(1), (gs*dph).sum(1)])
    return E, g

def classical(J, h, Ns=(60, 61, 62, 63, 64), starts=6, seed=0):
    J = np.asarray(J, float); h = np.asarray(h, float)
    rng = np.random.default_rng(seed)
    best = (np.inf, None, None)
    for N in Ns:
        for _ in range(starts):
            x0 = np.concatenate([np.arccos(rng.uniform(-1, 1, N)), rng.uniform(0, 2*np.pi, N)])
            r = minimize(energy, x0, args=(J, h, N), jac=True, method='L-BFGS-B', options=dict(maxiter=20000, gtol=1e-10))
            if r.fun/N < best[0]:
                th, ph = r.x[:N], r.x[N:]
                s = np.stack([np.sin(th)*np.cos(ph), np.sin(th)*np.sin(ph), np.cos(th)], 1)
                best = (r.fun/N, N, s)
    e, N, s = best
    F = np.abs(np.fft.fft(s, axis=0))**2
    F = F.sum(1); F[0] = F[0]  # include uniform part
    q = np.argmax(F)/N
    m = np.linalg.norm(s.mean(0))
    return e, min(q, 1-q), m

if __name__ == '__main__':
    pars = json.load(open('candidates_params.json'))
    out = {}
    for k, p in pars.items():
        e, q, m = classical(p['J'], p['field'])
        out[k] = dict(e_cl=e, q_cl=q, m_cl=m)
        print(f"{k:10s} e_cl={e:.6f} q_cl={q:.4f} |m|={m:.3f}", flush=True)
    json.dump(out, open('classical.json', 'w'), indent=1)
