"""C159 refinement of scan candidates (Claude Opus agent, 2026-11-07).
Selects from scan/*.jsonl:
  hi_<mode>   : c_obj = min(c_S, c_cas) > THR and crit < 0.15, best TOP per mode
  ising_<T|K> : 0.3 < c_S < 0.7, crit < 0.15, unique ground state at both sizes (ndeg == [1, 1]), best 8 by crit
  dist_<mode> : 4 random gapless points (crit < 0.15) per mode as controls
For each: ring estimate at (16, 20) and open-chain DMRG (dmrg_c.run, L = 32, 64, 128).
usage: refine_tk.py [--thr 1.1] [--top 8] [--procs 48] [--chi 128]"""
import sys, os, json, glob, time, argparse
import numpy as np
from multiprocessing import Pool



def work(arg):
    k, (tag, r), CHI, LS = arg
    os.environ['OMP_NUM_THREADS'] = '1'
    import est as E, dmrg_c as DC
    fn = f'refine/{tag}_{k}.json'
    if os.path.exists(fn):
        return fn
    v = np.array(r['vec'])
    out = dict(tag=tag, scan=r)
    t0 = time.time()
    try:
        J, f = E.J_f(v)
        out['ring16_20'] = E.estimate(E.bond_from_JB(J, f), Ls=(16, 20), nev=2)
    except Exception as ex:
        out['ring16_20'] = str(ex)
    json.dump(out, open(fn, 'w'), default=float)
    try:
        chis = tuple(c for c in (16, 32, 64, 96, CHI, CHI) if c <= CHI)
        out['dmrg'] = DC.run(v, tuple(LS), chis=chis)
    except Exception as ex:
        out['dmrg'] = str(ex)
    out['time_s'] = time.time() - t0
    json.dump(out, open(fn, 'w'), default=float)
    return fn



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--thr', type=float, default=1.1); ap.add_argument('--top', type=int, default=8)
    ap.add_argument('--procs', type=int, default=48); ap.add_argument('--chi', type=int, default=128)
    ap.add_argument('--Ls', type=int, nargs='+', default=[32, 64, 128])
    a = ap.parse_args()
    
    recs = []
    for fn in glob.glob('scan/*.jsonl'):
        for l in open(fn):
            try:
                r = json.loads(l)
            except Exception:
                continue
            if 'c_S' in r:
                recs.append(r)
    modes = sorted({r['mode'] for r in recs})
    rng = np.random.default_rng(7)
    cands = []
    for mo in modes:
        R = [r for r in recs if r['mode'] == mo]
        hi = sorted([r for r in R if r['c_obj'] > a.thr and r['crit'] < 0.15], key=lambda r: -r['c_obj'])[:a.top]
        cands += [('hi_' + mo, r) for r in hi]
        if mo in ('T', 'K', 'TFM', 'KFM'):
            isg = sorted([r for r in R if 0.3 < r['c_S'] < 0.7 and r['crit'] < 0.15 and list(r['ndeg']) == [1, 1]],
                         key=lambda r: r['crit'])[:8]
            cands += [('ising_' + mo, r) for r in isg]
        gl = [r for r in R if r['crit'] < 0.15]
        if gl:
            for j in rng.choice(len(gl), size=min(4, len(gl)), replace=False):
                cands.append(('dist_' + mo, gl[j]))
    print(len(recs), 'records;', len(cands), 'candidates', flush=True)
    os.makedirs('refine', exist_ok=True)
    with Pool(a.procs) as p:
        for fn in p.imap_unordered(work, [(k, c, a.chi, a.Ls) for k, c in enumerate(cands)]):
            print('done', fn, flush=True)


if __name__ == '__main__':
    main()
