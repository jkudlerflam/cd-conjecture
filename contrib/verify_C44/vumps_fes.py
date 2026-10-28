"""Finite-entanglement scaling with one-site VUMPS for the projected n=1 chain.
python3 -B vumps_fes.py MODEL chi1 chi2 ...   -> results/MODEL_vumps.json"""
import sys, json, time
import numpy as np
from dmrg import mpo_from_bond, vumps
from run_c import model

if __name__ == '__main__':
    name = sys.argv[1]
    chis = [int(x) for x in sys.argv[2:]]
    fn = 'results/%s_vumps.json' % name if len(chis) > 1 else 'results/%s_vumps_chi%d.json' % (name, chis[0])
    h2, h1, info, ops = model(name)
    W, err = mpo_from_bond(h2, h1, 4)
    out = dict(info={k: (v if not isinstance(v, complex) else [v.real, v.imag]) for k, v in info.items()},
               mpo_D=int(W.shape[0]), points=[])
    for chi in chis:
        best = None
        for seed in ((7, 11) if chi < 128 else (7,)):
            t0 = time.time()
            o = vumps(W, chi, 4, iters=800 if chi < 128 else 400, tol=1e-9 if chi < 128 else 3e-8, seed=seed,
                      ops={'nl': ops['nl'], 'nz': ops['nz']})
            rec = dict(chi=chi, seed=seed, iters=o['iters'], err=o['err'], e_site=o['e_site'], S=o['S'],
                       xi=o['xi'], tm_eigs=o['tm_eigs'], nl=o['nl'], nz=o['nz'], schmidt_min=o['schmidt_min'],
                       seconds=time.time() - t0)
            if o['err'] < 1e-6 and (best is None or rec['e_site'] < best['e_site'] - 1e-12):
                best = rec
            print(json.dumps(rec), flush=True)
        if best is not None:
            out['points'].append(best)
        json.dump(out, open(fn, 'w'), indent=1)
