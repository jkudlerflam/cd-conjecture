"""Central charge of the projected n=1 chain (referee item i): model setup and open-chain DMRG.

python3 -B run_c.py prepare            (writes kpq.npz from the rebuilt bond)
python3 -B run_c.py MODEL finite CHI L  (open-chain DMRG entanglement profile)
MODEL: M0 (decoupled free reference, c=2 exactly), C0 (payload only, c=1),
       M1 (scheduled point rho=0.05, s/delta->0), M2 (M1 with mixed traceless
       terms x300), C1 (rho=0.10: attractive length interaction beyond -2 j_l).
Reads kpq.npz (written by `python3 -B run_c.py prepare`).
"""
import json
import os
import sys
import time
import numpy as np
from effmodel import build
from dmrg import mpo_from_bond, finite_dmrg

MODELS = {
    'M0': dict(rho=0.05, kw=dict(decouple=True, drop_zz=True, nu=0.30)),
    'C0': dict(rho=0.05, kw=dict(decouple=True, drop_zz=True, nu=-0.5)),   # kappa=+... empties lengths
    'M1': dict(rho=0.05, kw=dict(nu=0.30)),
    'M2': dict(rho=0.05, kw=dict(nu=0.30, stress=300.0)),
    'C1': dict(rho=0.10, kw=dict(nu=0.30)),
}


def model(name):
    z = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kpq.npz'))
    spec = MODELS[name]
    tag = 'rho%s_' % spec['rho']
    Kpq = {k: z[tag + k] for k in ('11', '12', '21', '22')}
    kw = dict(spec['kw'])
    nu = kw.pop('nu')
    if nu < 0:
        h2, h1, info, ops = build(Kpq, nu=0.0, **kw)
        h1 = h1 + 2.0 * ops['nl']                 # length bits fully empty
        info['kappa'] = float(info['kappa'] - 2.0)
    else:
        h2, h1, info, ops = build(Kpq, nu=nu, **kw)
    return h2, h1, info, ops


if __name__ == '__main__':
    if sys.argv[1] == 'prepare':
        from precise import kbar_bonds
        out = {}
        for rho in (0.05, 0.1):
            kb = kbar_bonds(1, rho)
            for k, v in kb.items():
                out['rho%s_%s' % (rho, k)] = v
        np.savez('kpq.npz', **out)
        sys.exit()
    name, method, chi = sys.argv[1], sys.argv[2], int(sys.argv[3])
    h2, h1, info, ops = model(name)
    W, err = mpo_from_bond(h2, h1, 4)
    t0 = time.time()
    os.makedirs('results', exist_ok=True)
    if method == 'finite':
        L = int(sys.argv[4])
        sched = [min(32, chi), min(64, chi)] + [chi] * 8
        o = finite_dmrg(W, L, sched, 4, ops={'nl': ops['nl'], 'nz': ops['nz']})
        o['S'] = list(o['S'])
        o['meas'] = {k: list(v) for k, v in o['meas'].items()}
        fn = 'results/%s_finite_L%d_chi%d.json' % (name, L, chi)
    o['info'] = {k: (v if not isinstance(v, complex) else [v.real, v.imag]) for k, v in info.items()}
    o['mpo_D'] = int(W.shape[0])
    o['mpo_err'] = float(err)
    o['seconds'] = time.time() - t0
    json.dump(o, open(fn, 'w'), indent=1, default=float)
    print(fn, o['E'], o['seconds'])
