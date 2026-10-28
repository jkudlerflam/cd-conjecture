"""(i) |T| against the neutral express phase: the leading-order prefactor should be
proportional to sin(n phi_x) (two channels with conjugate transfer phases).
(ii) larger n at moderate D_x: does a fixed D_x keep R_s > 0?"""
import json, numpy as np
from rt_bond import P
from rt_token import T_token, rs_eigs
out = dict(phase=[], large_n=[])
for n in (3, 5, 8):
    for phx in (0.1, 0.2, 0.3, 0.4, np.pi / n, 0.5, 0.6):
        p = P(td=0.0, Dx=100.0, DM=100.0, ph_x={'0': phx, 'p': 0.0, 'm': 0.0})
        T = T_token(n, p)
        p0 = P(td=0.0, Dx=100.0, DM=100.0, ph_x={'0': np.pi / (2 * n), 'p': 0.0, 'm': 0.0})
        ref = np.abs(T_token(n, p0)).max()
        r = dict(n=n, phx=float(phx), T_rel=float(np.abs(T).max() / ref), sin_n_phi=float(abs(np.sin(n * phx))))
        out['phase'].append(r); print(r, flush=True)
for n in (8, 9, 10):
    for Dx in (10.0, 20.0, 30.0, 50.0, 100.0, 300.0):
        p = P(td=0.0, Dx=Dx, DM=Dx, ph_x={'0': np.pi / (2 * n), 'p': 0.0, 'm': 0.0})
        T = T_token(n, p)
        ev = rs_eigs(T)
        R = T * np.exp(-1j * np.angle(T[0, 0]))
        r = dict(n=n, Dx=Dx, T_scale=float(np.abs(T).max()), lam_min=float(ev[0]), imag=float(np.abs(R.imag).max() / np.abs(R).max()))
        out['large_n'].append(r); print(r, flush=True)
json.dump(out, open('scan_phase.json', 'w'), indent=1)
