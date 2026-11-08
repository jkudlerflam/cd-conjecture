"""C155 summary of all stages (Claude Opus agent, 2026-11-07). Run in the job directory."""
import glob, json, sys
import numpy as np

GAPLESS = 0.15


def jl(pattern):
    out = []
    for fn in sorted(glob.glob(pattern)):
        for line in open(fn):
            try:
                out.append(json.loads(line))
            except Exception:
                pass
    return out


def main():
    print('== Stage A: uniform sample of S^8 (rings 12, 16) ==')
    A = [r for r in jl('results_crit/uniform_*.jsonl') if r.get('est')]
    if A:
        crit = np.array([r['est']['crit'] for r in A]); cobj = np.array([r['est']['c_obj'] for r in A])
        cS = np.array([r['est']['c_S'] for r in A]); deg = np.array([r['est']['ndeg'][-1] >= 3 for r in A])
        g = (crit < GAPLESS) & ~deg
        print(f'points {len(A)}; flat-band (ndeg>=3) {deg.sum()}; gapless by crit<{GAPLESS}: {g.sum()} ({g.mean():.4f})')
        for t in (0.05, 0.1, 0.15, 0.25):
            print(f'  crit<{t}: fraction {np.mean((crit < t) & ~deg):.4f}')
        print(f'  c_obj quantiles (all) 50/90/99/max: {np.percentile(cobj, [50, 90, 99]).round(3)} {cobj.max():.3f}')
        if g.any():
            print(f'  c_obj among gapless: median {np.median(cobj[g]):.3f}, max {cobj[g].max():.3f}; c_S among gapless max {cS[g].max():.3f}')
        print(f'  points with c_obj > 1.1: {(cobj > 1.1).sum()}; with c_obj > 1.1 and gapless: {((cobj > 1.1) & g).sum()}')
    print('\n== Stage B: gap-closing minimization ==')
    B = jl('results_crit/gap_*.jsonl')
    if B:
        cm = np.array([r['crit_min'] for r in B])
        print(f'restarts {len(B)}; minima with crit < {GAPLESS}: {(cm < GAPLESS).sum()}')
        cmin = np.array([(r['est_min'] or {}).get('c_obj', np.nan) for r in B])
        csl = np.array([(r['est_slide'] or {}).get('c_obj', np.nan) for r in B])
        print(f'  c_obj at minima: median {np.nanmedian(cmin):.3f}, max {np.nanmax(cmin):.3f}')
        print(f'  c_obj after slide: median {np.nanmedian(csl):.3f}, max {np.nanmax(csl):.3f}')
        labs = {}
        for r in B:
            if r['crit_min'] < GAPLESS:
                key = 'known(U1/free)' if any(l in r['labels_min'] for l in ('U1', 'U1stag', 'free')) else 'other'
                labs[key] = labs.get(key, 0) + 1
        print('  critical minima by stratum:', labs)
        fr3 = [r['ball']['0.03']['frac_gapless'] for r in B if r['crit_min'] < GAPLESS and r['ball']['0.03']['n']]
        fr8 = [r['ball']['0.08']['frac_gapless'] for r in B if r['crit_min'] < GAPLESS and r['ball']['0.08']['n']]
        if fr3:
            print(f'  ball test around critical minima: mean gapless fraction r=0.03 {np.mean(fr3):.3f}, r=0.08 {np.mean(fr8):.3f}')
            print(f'  minima with full gapless ball r=0.03: {sum(f == 1 for f in fr3)}, r=0.08: {sum(f == 1 for f in fr8)}')
        mb = [r['ball']['0.08']['max_c_obj'] for r in B if r['ball']['0.08']['n']]
        print(f'  max c_obj inside any ball: {np.max(mb):.3f}')
        top = sorted(B, key=lambda r: -np.nan_to_num((r['est_slide'] or {}).get('c_obj', -9)))[:8]
        for r in top:
            e = r['est_slide']
            print(f"   slide top: c_obj {e['c_obj']:.3f} c_S {e['c_S']:.3f} c_cas {e['c_cas']:.3f} crit {e['crit']:.3f} labels {r['labels_slide']}")
    print('\n== Stage C: c maximization ==')
    C = jl('results/task_*.jsonl')
    if C:
        bv = np.array([r['best_val'] for r in C])
        nf = sum(r['nfev'] for r in C)
        print(f'restarts {len(C)}, evaluations {nf}; best values 50/90/max: {np.percentile(bv, [50, 90]).round(3)} {bv.max():.3f}')
        for mode in ('min', 'Sfit'):
            sub = [r for r in C if r['mode'] == mode]
            if sub:
                print(f'  mode {mode}: {len(sub)} restarts, max {max(r["best_val"] for r in sub):.3f}, '
                      f'> 1.1: {sum(r["best_val"] > 1.1 for r in sub)}')
        seeds = {}
        for r in C:
            k = r['seed'].split(':')[0]
            seeds.setdefault(k, []).append(r['best_val'])
        for k, vals in seeds.items():
            print(f'  seed {k}: {len(vals)} restarts, max {max(vals):.3f}')
        for r in sorted(C, key=lambda r: -r['best_val'])[:10]:
            e = r['best_est']
            print(f"   top: {r['best_val']:.3f} {r['seed']} {r['mode']} c_S {e['c_S']:.3f} c_fit {e['c_fit']:.3f} c_cas {e['c_cas']:.3f} "
                  f"crit {e.get('crit', np.nan):.3f} m0 {e['m0']} v {np.round(r['best_v'], 3).tolist()}")
    print('\n== Stage E: XZ neighbourhood and Kitaev [110] field ==')
    for fn in sorted(glob.glob('xz/pt_*.json'), key=lambda s: int(s.split('_')[-1].split('.')[0])):
        r = json.load(open(fn))
        e1 = r.get('ed_12_16', {}); e2 = r.get('ed_16_20', {}); d = r.get('dmrg', {})
        dm = ' '.join(f"L{L}:{x['c_prof']:.2f}/{x['c_sm']:.2f}" for L, x in d.get('L', {}).items())
        print(f"{r['dir']:16s} d={r['delta']:.2f} ED12/16 c_S {e1.get('c_S', np.nan):.2f} c_cas {e1.get('c_cas', np.nan):.2f} crit {e1.get('crit', np.nan):.2f} | "
              f"ED16/20 c_S {e2.get('c_S', np.nan):.2f} c_cas {e2.get('c_cas', np.nan):.2f} | DMRG {dm} two {d.get('c_two', np.nan):.2f} "
              f"two128/256 {d.get('c_two_128_256', np.nan):.2f} | {r['labels']}")
    print('\n== Stage D: refinement ==')
    for fn in sorted(glob.glob('refine/cand_*.json'), key=lambda s: int(s.split('_')[-1].split('.')[0])):
        r = json.load(open(fn))
        e = r.get('ed_16_20', {}); d = r.get('dmrg', {})
        Ld = d.get('L', {})
        dm = ' '.join(f"L{L}:{x['c_prof']:.2f}/{x['c_sm']:.2f}(tr{x['trunc']:.0e})" for L, x in Ld.items())
        print(f"cand {r['id']:3d} screen {r['screen_val']:.3f} ED16/20 c_S {e.get('c_S', np.nan):.3f} c_cas {e.get('c_cas', np.nan):.3f} "
              f"crit {e.get('crit', np.nan):.3f} | DMRG {dm} c_two {d.get('c_two', np.nan):.3f} | {r['labels']} {r['source']}")


if __name__ == '__main__':
    main()
