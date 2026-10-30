"""Finite-entanglement scaling with one-site VUMPS (verify_C44/dmrg.py) for the n = 2 chain.

python3 -B run_fes.py MODEL CHI [CHI_L]  -> results/MODEL_vumps_chi{CHI}[_l{CHI_L}].json

From random initial states VUMPS repeatedly lands in a metastable state with
payload density 0.75 per cell (the exact density is 1). Runs therefore start
from the product of a payload-only VUMPS state (model C0, bond dimension
CHI/CHI_L, which converges to density 1 from random starts) and a free
length-chain VUMPS state (bond dimension CHI_L), and VUMPS then optimizes the
full model from there. C0 and C1 start from the payload state alone.
"""
import sys
import os
import json
import time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'verify_C44'))
import numpy as np
from model2 import build
from dmrg import mpo_from_bond, vumps

CACHE = 'results/init_cache'


def payload_state(chi_p):
    fn = '%s/payload_chi%d.npz' % (CACHE, chi_p)
    if os.path.exists(fn):
        z = np.load(fn)
        return z['AL'], z['AR'], z['C']
    h2, h1, info, ops, T = build('C0')
    W, err = mpo_from_bond(h2, h1, 8)
    best, log = None, []
    for seed in (7, 11, 13, 17, 19, 23):
        o = vumps(W, chi_p, 8, iters=1500, tol=1e-9, seed=seed, ops={'nz': ops['nz']})
        log.append((seed, o['e_site'], o['nz']))
        print('payload init chi %d seed %d e %.10f nz %.4f' % (chi_p, seed, o['e_site'], o['nz']), flush=True)
        if abs(o['nz'] - 1) < 0.01 and (best is None or o['e_site'] < best['e_site']):
            best = o
    if best is None:
        raise RuntimeError('no half-filled payload state at chi %d' % chi_p)
    AL, AR, C = best['AL'][:, :4, :], best['AR'][:, :4, :], best['C']
    os.makedirs(CACHE, exist_ok=True)
    np.savez(fn, AL=AL, AR=AR, C=C)
    return AL, AR, C


def length_state(chi_l, jl, kappa):
    hop = np.zeros((4, 4))
    hop[1, 2] = hop[2, 1] = -jl
    nl = np.diag([0., 1.])
    W, err = mpo_from_bond(hop.astype(complex), (-kappa * nl).astype(complex), 2)
    o = vumps(W, chi_l, 2, iters=2000, tol=1e-10, seed=3)
    return o['AL'], o['AR'], o['C']


def product(Al, Ap):
    # A[(al ap), 4 l + z, (bl bp)]
    t = np.einsum('alb,czd->aclzbd', Al, Ap)
    a, c, l, z, b, d = t.shape
    return t.reshape(a * c, l * z, b * d)


if __name__ == '__main__':
    name, chi = sys.argv[1], int(sys.argv[2])
    chi_l = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    h2, h1, info, ops, T = build(name)
    W, err = mpo_from_bond(h2, h1, 8)
    if name in ('C0', 'C1'):
        ALp, ARp, Cp = payload_state(chi)
        z = np.zeros((chi, 4, chi), complex)
        A0 = (np.concatenate([ALp, z], 1), np.concatenate([ARp, z], 1), Cp)
    else:
        ALp, ARp, Cp = payload_state(chi // chi_l)
        ALl, ARl, Cl = length_state(chi_l, info['jl'], info['kappa'])
        A0 = (product(ALl, ALp), product(ARl, ARp), np.kron(Cl, Cp))
    out = dict(info=info, mpo_D=int(W.shape[0]), mpo_err=float(err), chi_l=chi_l, runs=[])
    os.makedirs('results', exist_ok=True)
    fn = 'results/%s_vumps_chi%d.json' % (name, chi) if chi_l == 3 else 'results/%s_vumps_chi%d_l%d.json' % (name, chi, chi_l)
    t0 = time.time()
    o = vumps(W, chi, 8, iters=3000, tol=1e-9, A0=A0, ops={'nl': ops['nl'], 'nz': ops['nz']})
    rec = dict(chi=chi, chi_l=chi_l, iters=o['iters'], err=o['err'], e_site=o['e_site'], S=o['S'], xi=o['xi'],
               tm_eigs=o['tm_eigs'], nl=o['nl'], nz=o['nz'], schmidt_min=o['schmidt_min'],
               seconds=time.time() - t0, last=o['hist'][-1])
    print(json.dumps(rec, default=float), flush=True)
    out['runs'].append(rec)
    json.dump(out, open(fn, 'w'), indent=1, default=float)
