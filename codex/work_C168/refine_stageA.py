"""C168: refinement of the 28 C155 stage-A uniform points flagged gapless (crit < 0.15) with c_obj > 1.1.
These were never refined by C155 stage D (refine.py collects only search and gap-closing outputs).
For each point: ring ED at L = 16, 20 (C155 est.estimate, nev = 3) and open-chain DMRG at L = 32, 64
(C155 dmrg_c.run). Claude Opus 5.5, 2026-11-08. Usage: python3 refine_stageA.py [nproc]"""
import json, os, sys, time
import numpy as np
from multiprocessing import Pool

KEYS = ('c_S', 'c_fit', 'c_cas', 'c_obj', 'z1', 'z2', 'crit', 'ndeg', 'v', 'm0', 'gap1', 'gap2', 'Shalf')


def work(item):
    import est_c155 as est, dmrg_c_c155 as dc
    fn = f"refineA/pt_{item['idx']}.json"
    if os.path.exists(fn):
        return fn
    t0 = time.time(); out = dict(item)
    try:
        e = est.estimate(est.bond_v(item['v']), Ls=(16, 20), nev=3)
        out['ed_16_20'] = {k: e[k] for k in KEYS}
    except Exception as ex:
        out['ed_16_20'] = dict(error=str(ex))
    out['t_ed'] = time.time() - t0
    try:
        o = dc.run(item['v'], Ls=(32, 64), chis=(16, 32, 64, 96, 96))
        for L in o['L']:
            o['L'][L].pop('S', None)
        out['dmrg'] = o
    except Exception as ex:
        out['dmrg'] = dict(error=str(ex))
    out['t_total'] = time.time() - t0
    json.dump(out, open(fn, 'w'))
    return fn


if __name__ == '__main__':
    os.makedirs('refineA', exist_ok=True)
    items = json.load(open('stageA_flags.json'))
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    with Pool(n) as p:
        for fn in p.imap_unordered(work, items):
            print('done', fn, flush=True)
