#!/usr/bin/env python3
"""Calibration of the ring Casimir estimator on exactly solvable (quadratic) chains.

For a quadratic chain with Bloch Hamiltonian h(k), the grand-canonical ground-state energy of a ring of L
cells with flux theta is E(L, theta) = sum over k = (2 pi m + theta)/L of the negative band energies.
We minimize over theta (64 values, which removes the current zero mode of one pair of Fermi points) and form
    c_Cas(L) = 6 L [L e_inf - min_theta E(L, theta)] / (pi v_sum),   v_sum = sum over Fermi seas of (v_R + v_L)/2,
which tends to the number of Fermi seas only if the density zero modes of all seas vanish together.
usage: ring_casimir_free.py SPEC [SPEC ...]  (SPEC as in run_chain.py; 'hf' or quadratic levels)
Written by Claude (Opus 5.5) for contrib/routeA_n2_chains, 2026-10-28.
"""
import os, sys, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import bonds as B, bands as BD, run_chain as R


def ring_E(T, M, L, th):
    ks = (2 * np.pi * np.arange(L) + th) / L
    E = np.array([np.linalg.eigvalsh(BD.bloch(T, M, k)) for k in ks]).ravel()
    return E[E < 0].sum()


def analyse(spec, Ls=range(6, 41)):
    K, meta = R.build(spec)
    c, h, U, Hq, H4, r = B.normal_order(K)
    T, M = BD.blocks(h)
    pts, _ = BD.fermi_points(T, M, 0.0)
    vsum = sum(abs(p['v']) for p in pts) / 2
    nk = 200000
    ks = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    Eb = np.array([np.linalg.eigvalsh(BD.bloch(T, M, k)) for k in ks])
    einf = Eb[Eb < 0].sum() / nk
    cs = []
    for L in Ls:
        e = min(ring_E(T, M, L, th) for th in np.linspace(0, 2 * np.pi, 64, endpoint=False))
        cs.append(float(6 * L * (L * einf - e) / (np.pi * vsum)))
    cs = np.array(cs)
    return dict(spec=spec, seas=len(pts) // 2, vsum=vsum, L=list(Ls), c=np.round(cs, 3).tolist(),
                median_L20_40=float(np.median(cs[np.array(list(Ls)) >= 20])),
                range_L20_40=[float(cs[np.array(list(Ls)) >= 20].min()), float(cs[np.array(list(Ls)) >= 20].max())])


if __name__ == '__main__':
    for s in sys.argv[1:]:
        o = analyse(s)
        print(json.dumps(o), flush=True)
