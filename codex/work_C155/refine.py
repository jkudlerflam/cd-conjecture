"""C155 stage D: refinement of search candidates (Claude Opus agent, 2026-11-07).

Collects candidate points from results/*.jsonl (c maximization) and results_crit/gap_*.jsonl (gap closings),
clusters them (greedy, distance 0.03 on S^8), and for every cluster representative with screening value
> --thr (plus the --top best representatives regardless of value) runs
  (a) ring ED at L = 16, 20 (est.estimate, nev = 3),
  (b) open-chain DMRG at L = 32, 64, 128 (dmrg_c.run, chi up to --chi),
writing refine/cand_<id>.json (skipped if present, so the stage can be resumed).
"""
import argparse, glob, json, os, time
import numpy as np
from multiprocessing import Pool


def collect(thr, top):
    pts = []
    for fn in glob.glob('results/task_*.jsonl'):
        for line in open(fn):
            try:
                r = json.loads(line)
            except Exception:
                continue
            if r.get('best_v') is not None:
                pts.append((r['best_val'], r['best_v'], 'search:' + r['seed'] + ':' + r['mode']))
            for hgh in r.get('highs', []):
                pts.append((hgh['val'], hgh['point'], 'search-high'))
    for fn in glob.glob('results_crit/gap_*.jsonl'):
        for line in open(fn):
            try:
                r = json.loads(line)
            except Exception:
                continue
            for key, vk in (('est_min', 'vmin'), ('est_slide', 'vslide')):
                e = r.get(key)
                if e:
                    val = e['c_obj'] - 2 * max(0.0, min(e['crit'], 3.0) - 0.2)
                    pts.append((val, r[vk], 'gap:' + key + ':' + ','.join(r.get('labels_' + key[4:], []))))
    if os.path.exists('extra_candidates.json'):   # local-run flags, always refined
        for val, v, src in json.load(open('extra_candidates.json')):
            pts.append((max(val, thr + 1e-3), v, src))
    pts.sort(key=lambda p: -p[0])
    reps = []
    for val, v, src in pts:
        v = np.asarray(v, float); v /= np.linalg.norm(v)
        if any(np.linalg.norm(v - np.asarray(q[1])) < 0.03 for q in reps):
            continue
        reps.append((val, v.tolist(), src))
    sel = [r for r in reps if r[0] > thr]
    sel += [r for r in reps[:top] if r[0] <= thr]
    return sel, len(pts), len(reps)


def work(item):
    cid, (val, v, src), chi = item
    fn = f'refine/cand_{cid}.json'
    if os.path.exists(fn):
        return fn
    import est, dmrg_c
    from crit import labels
    t0 = time.time()
    out = dict(id=cid, screen_val=val, v=v, source=src, labels=labels(v))
    try:
        e = est.estimate(est.bond_v(v), Ls=(16, 20), nev=3)
        out['ed_16_20'] = {k: e[k] for k in ('c_S', 'c_fit', 'c_cas', 'c_obj', 'z1', 'z2', 'crit', 'ndeg', 'v', 'm0', 'gap1', 'gap2')}
    except Exception as ex:
        out['ed_16_20'] = dict(error=str(ex))
    out['t_ed'] = time.time() - t0
    try:
        o = dmrg_c.run(v, Ls=(32, 64, 128), chis=(16, 32, 64, chi, chi, chi))
        for L in o['L']:
            o['L'][L].pop('S', None) if False else None
        out['dmrg'] = o
    except Exception as ex:
        out['dmrg'] = dict(error=str(ex))
    out['t_total'] = time.time() - t0
    json.dump(out, open(fn, 'w'))
    return fn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--thr', type=float, default=1.1)
    ap.add_argument('--top', type=int, default=24)
    ap.add_argument('--chi', type=int, default=128)
    ap.add_argument('--procs', type=int, default=48)
    args = ap.parse_args()
    os.makedirs('refine', exist_ok=True)
    sel, npts, nrep = collect(args.thr, args.top)
    json.dump(dict(n_points=npts, n_clusters=nrep, n_selected=len(sel), selected=sel), open('refine/selection.json', 'w'))
    print(f'{npts} candidate points, {nrep} clusters, {len(sel)} selected', flush=True)
    items = [(i, s, args.chi) for i, s in enumerate(sel)]
    with Pool(args.procs) as p:
        for fn in p.imap_unordered(work, items):
            print('done', fn, flush=True)


if __name__ == '__main__':
    main()
