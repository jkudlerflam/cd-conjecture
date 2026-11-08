"""Kitaev chain in a [110] field = XZ + B (x + z)/sqrt 2 (C156). Ring ED at 12/16 and 16/20, open DMRG L = 64 (chi 64).
Prints each stage as it finishes. Usage: python3 probe_kitaev.py B1 B2 ..."""
import sys, time, json
import numpy as np
from multiprocessing import Pool


def one(B):
    import est, dmrg_c
    from search import W0, F_BIS, w_to_v
    v = w_to_v(W0 + B * F_BIS)
    out = dict(B=B, v=v.tolist())
    t = time.time()
    for Ls in ((12, 16), (16, 20)):
        e = est.estimate(est.bond_v(v), Ls=Ls, nev=3)
        out[f'ed{Ls}'] = {k: e[k] for k in ('c_S', 'c_fit', 'c_cas', 'crit', 'z1', 'z2', 'ndeg', 'gap1', 'gap2', 'm0')}
        print(f"B={B} ED{Ls}: c_S {e['c_S']:.3f} c_fit {e['c_fit']:.3f} c_cas {e['c_cas']:.3f} crit {e['crit']:.3f} gap1 {np.round(e['gap1'], 5).tolist()} ndeg {e['ndeg']} m0 {e['m0']} ({time.time()-t:.0f}s)", flush=True)
    o = dmrg_c.run(v, Ls=(64,), chis=(16, 32, 64, 64, 64))
    x = o['L'][64]
    out['dmrg64'] = {k: x[k] for k in ('c_prof', 'c_sm', 'trunc', 'dE_last', 'Smid')}
    print(f"B={B} DMRG64: c_prof {x['c_prof']:.3f} c_sm {x['c_sm']:.3f} trunc {x['trunc']:.1e} dE {x['dE_last']:.1e} ({time.time()-t:.0f}s)", flush=True)
    json.dump(out, open(f'kitaev_B{B}.json', 'w'))
    return B


if __name__ == '__main__':
    Bs = [float(a) for a in sys.argv[1:]]
    with Pool(len(Bs)) as p:
        list(p.imap_unordered(one, Bs))
