#!/usr/bin/env python3
"""C183 (Claude Opus 5.5 agent, 2026-11-09): the packed Typhon pipeline (stage S screen and stage C confirmation),
run inside one SLURM allocation (never on a login node).  Families, rules and thresholds are those registered in
RESULT.md before any output existed.

Phase 0  tune every (r, |w|, phi) (tune183.tune, ring ED L = 12) -> tunes.json
Phase 1  density ladders, L = 32, chi <= 64, dressed chemical potential x in XS -> lad/
Phase 2  screen at the target dressed densities, L = 32 and 64 -> scr/ ; free controls -> ctrl/
Phase 3  screen summary and selection -> screen_summary.json, screen_summary.txt
Phase 4  confirmation of the selected points: L = 128, 256, parity sectors at L = 64, 128, VUMPS n = 2 -> conf/, vc/
usage: pipeline183.py [--workers 24] [--threads 2] [--phase4-workers 16] [--phase4-threads 3] [--dry]
"""
import os, sys, json, math, time, argparse, subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
PI2 = math.pi / 2
XS = [-0.4, -0.2, -0.1, 0.0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.45, 0.6, 0.8, 1.0, 1.6, 2.4]
CH_LAD = [16, 32, 64, 64]
CH_32 = [32, 64, 96, 128, 128]
CH_64 = [32, 64, 96, 128, 160, 192, 192]
CH_128 = [32, 64, 96, 128, 192, 256, 256, 256]
CH_256 = [32, 64, 96, 128, 192, 256, 320, 320, 320]
CH_P64 = [32, 64, 96, 128, 160, 160]
CH_P128 = [32, 64, 96, 128, 192, 192, 192]
CH_V = [16, 24, 32, 48, 64, 96, 128, 160, 192, 256]
LS = dict(lad=32, s1=32, s2=64, c1=128, c2=256, p1=64, p2=128)


def families():
    F = []
    for r in (4, 6, 8, 10, 12, 16):
        for wa in (0.5, 1.0, 2.0, 3.0):
            tg = [0.08, 0.17, 0.28]
            if wa == 1.0 and r in (8, 12):
                tg = [0.08, 0.28]          # n_d = 0.17 there is C177's
            F.append(dict(id=f"A_r{r}_w{wa:g}", group='A', r=float(r), wa=wa, phi=PI2, extra={}, targets=tg))
    for phi in (1.56, 1.53, 1.48, 1.40):
        for r in (6, 10):
            for wa in (1.0, 2.0):
                F.append(dict(id=f"B_r{r}_w{wa:g}_p{phi:.2f}", group='B', r=float(r), wa=wa, phi=phi, extra={}, targets=[0.08, 0.17]))
    # addendum D (registered in RESULT.md after the r <= 6 rows of stage T, before any screen output)
    for wa in (0.5, 1.0, 1.5):
        F.append(dict(id=f"D_r5_w{wa:g}", group='D', r=5.0, wa=wa, phi=PI2, extra={}, targets=[0.06, 0.12, 0.17, 0.24]))
    for wa in (0.5, 1.0):
        F.append(dict(id=f"D_r4_w{wa:g}", group='D', r=4.0, wa=wa, phi=PI2, extra={}, targets=[0.06, 0.12]))
    for phi in (1.565, 1.55, 1.53):
        F.append(dict(id=f"D_r5_w1_p{phi:.3f}", group='D', r=5.0, wa=1.0, phi=phi, extra={}, targets=[0.06, 0.12]))
    # family G (coordinator correction, registered in RESULT.md before any screen output; replaces addendum E):
    # C186's transverse-DM bond, Kxz = Gamma (XZ - ZX), Gamma = sqrt(v)/4, |w| = 1.035, valleys made degenerate
    # (offset 0) or split by 0.05 |t2| in the full chain (tune183.tune_dm).
    for r in (12, 16):
        for v in (0.3, 0.6):
            for off in ((0.0, 0.05) if v == 0.3 else (0.0,)):
                F.append(dict(id=f"G_r{r}_v{v:g}_o{off:g}", group='G', r=float(r), wa=1.035, phi=PI2, extra={},
                              dm=float(np.sqrt(v) / 4), off=off, targets=[0.06, 0.12, 0.17]))
    # family H (coordinator note, C181's partially polarized point, registered before any screen output):
    # parity class, |w| = 1, phi = 1.555 and 1.53 (and 1.545), r = 12, 16, 24, n = 3/16 from the onset, plus the
    # one-valley control phi = 1.30 at r = 16 on the same chain.
    for r in (12, 16, 24):
        for phi in (1.555, 1.545, 1.53):
            F.append(dict(id=f"H_r{r}_p{phi:.3f}", group='H', r=float(r), wa=1.0, phi=phi, extra={}, targets=[0.1875]))
    F.append(dict(id="H_r16_p1.300", group='H', r=16.0, wa=1.0, phi=1.30, extra={}, targets=[0.1875]))
    for r in (6, 10):
        for wa in (1.0, 2.0):
            F.append(dict(id=f"K_r{r}_w{wa:g}_p1.30", group='K', r=float(r), wa=wa, phi=1.30, extra={}, targets=[0.17]))
    for r in (6, 10):
        for wa in (1.0, 2.0):
            for term in ('Jxz', 'Kxz', 'hx'):
                for lam in (0.1, 0.3):
                    F.append(dict(id=f"C_r{r}_w{wa:g}_{term}{lam:g}", group='C', r=float(r), wa=wa, phi=PI2,
                                  extra={term: lam}, targets=[0.17], parent=f"A_r{r}_w{wa:g}"))
    return F


def _tune(f):
    import tune183 as TU
    if f.get('dm') is not None:
        return tkey(f), TU.tune_dm(f['r'], f['wa'], f['dm'], offset=f['off'])
    return tkey(f), TU.tune(f['r'], f['wa'], f['phi'])


def tkey(f):
    if f.get('dm') is not None:
        return f"dm_{f['r']:g}_{f['wa']:g}_{f['dm']:.6f}_{f['off']:g}"
    return f"{f['r']:g}_{f['wa']:g}_{f['phi']:.6f}"


def run(cmd, threads, log):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS=str(threads), MKL_NUM_THREADS=str(threads))
    t0 = time.time()
    with open(log, 'w') as fh:
        p = subprocess.run(cmd, cwd=HERE, env=env, stdout=fh, stderr=subprocess.STDOUT)
    return p.returncode, time.time() - t0


def fc_cmd(tag, L, cp, chis, outdir, extra_args=()):
    return ['python3', 'fc.py', tag, str(L), '--cp', json.dumps(cp), '--chis', *map(str, chis), '--outdir', outdir, *extra_args]


def load(fn):
    try:
        return json.load(open(fn))
    except Exception:
        return None


def n_central(d):
    Z = np.array(d['Z']); L = len(Z)
    return float(np.mean((1 - Z[L // 4: 3 * L // 4]) / 2))


def choose_x(lad, n_vac, target):
    """lad: list of (x, n_raw, Sbar).  Amendment (RESULT.md, before the job started): the dressed density is counted from
    the gas onset, n_d = n_raw - n_on, where the onset is the first ladder point whose Sbar(32) exceeds the Sbar of the
    lowest x by 0.3, and n_on is n_raw at the ladder point before it (the dressed vacuum just below the onset, which
    carries more virtual pairs than the far-below vacuum n_vac).  Points below the onset count as n_d = 0.
    Linear interpolation of n_d in x; flag a jump > 0.1 across the target."""
    lad = sorted(lad)
    S0 = lad[0][2] if lad else 0.0
    ion = 0 if S0 > 0.6 else next((i for i, p in enumerate(lad) if p[2] > S0 + 0.3), None)   # gas already at the lowest x
    if ion is None:
        return None, 'nogas'
    n_on = lad[ion - 1][1] if ion > 0 else n_vac
    pts = [(p[0], (p[1] - n_on) if i >= ion else 0.0, p[2]) for i, p in enumerate(lad)]
    for (x0, n0, s0), (x1, n1, s1) in zip(pts, pts[1:]):
        if n0 <= target <= n1 and n1 > n0:
            if s1 < s0 - 0.3:          # the upper point is in the dense (blockade-packed) phase: stay on the gas side
                return x0, 'jump'
            if n1 - n0 > 0.15:
                return x1, 'jump'
            return x0 + (target - n0) * (x1 - x0) / (n1 - n0), 'ok'
    if pts and pts[-1][1] < target:
        return None, 'unreached'
    return None, 'none'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workers', type=int, default=24)
    ap.add_argument('--threads', type=int, default=2)
    ap.add_argument('--phase4-workers', type=int, default=16)
    ap.add_argument('--phase4-threads', type=int, default=3)
    ap.add_argument('--dry', action='store_true')
    ap.add_argument('--only', default=None, help='comma list of family ids (testing)')
    ap.add_argument('--test', action='store_true', help='tiny sizes for a local plumbing test')
    a = ap.parse_args()
    global XS, CH_LAD, CH_32, CH_64, CH_128, CH_256, CH_P64, CH_P128, CH_V, LS
    if a.test:
        XS = [0.0, 0.3, 1.0]; CH_LAD = [8, 16]; CH_32 = CH_64 = CH_128 = CH_256 = CH_P64 = CH_P128 = [8, 16]; CH_V = [8, 12]
        LS = dict(lad=12, s1=12, s2=16, c1=16, c2=20, p1=12, p2=16)
    for d in ('lad', 'scr', 'ctrl', 'conf', 'vc', 'logs'):
        os.makedirs(os.path.join(HERE, d), exist_ok=True)
    F = families()
    if a.only:
        keep = set(a.only.split(','))
        F = [f for f in F if f['id'] in keep or f.get('parent') in keep]
    print(f"{len(F)} families", flush=True)
    if a.dry:
        for f in F:
            print(f['id'], f['targets'])
        return

    # Phase 0: tunes
    import tune183 as TU
    tfn = os.path.join(HERE, 'tunes.json')
    tunes = load(tfn) or {}
    todo = sorted({tkey(f): f for f in F}.values(), key=tkey)
    todo = [f for f in todo if tkey(f) not in tunes]

    os.environ.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=min(a.workers, 24)) as ex:
        for k, T in ex.map(_tune, todo):
            tunes[k] = T
    json.dump(tunes, open(tfn, 'w'), indent=1)
    print('phase 0 done', flush=True)

    def cp_of(f, x):
        T = tunes[tkey(f)]
        extra = {k: v * T['t2'] for k, v in f['extra'].items()}
        if f.get('dm') is not None:
            extra = dict(T['extra'])
        return TU.couplings(T, x, extra)

    # Phase 1: ladders
    jobs = []
    for f in F:
        if f['group'] == 'C':
            continue
        for x in XS:
            tag = f"{f['id']}_x{x:g}"
            if not os.path.exists(os.path.join(HERE, 'lad', f"{tag}_{LS['lad']}.json")):
                jobs.append((fc_cmd(tag, LS['lad'], cp_of(f, x), CH_LAD, 'lad', ('--nocorr',)), a.threads, os.path.join(HERE, 'logs', f"lad_{tag}.log")))
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        list(ex.map(lambda j: run(*j), jobs))
    print('phase 1 done', flush=True)

    # Phase 2: choose x, screen
    plan = []
    for f in F:
        if f['group'] == 'C':
            continue
        T = tunes[tkey(f)]
        lad = []
        for x in XS:
            d = load(os.path.join(HERE, 'lad', f"{f['id']}_x{x:g}_{LS['lad']}.json"))
            if d:
                lad.append((x, n_central(d), d['Sbar']))
        f['ladder'] = lad
        for tg in f['targets']:
            x, flag = choose_x(lad, T['n_vac'] if T.get('n_vac') is not None else 0.0, tg)
            plan.append(dict(fam=f['id'], group=f['group'], target=tg, x=x, flag=flag))
    byid = {f['id']: f for f in F}
    for f in F:
        if f['group'] != 'C':
            continue
        par = [p for p in plan if p['fam'] == f['parent'] and p['target'] == 0.17]
        x = par[0]['x'] if par else None
        plan.append(dict(fam=f['id'], group='C', target=0.17, x=x, flag=(par[0]['flag'] if par else 'noparent')))
    json.dump(dict(plan=plan, ladders={f['id']: f.get('ladder') for f in F}), open(os.path.join(HERE, 'plan.json'), 'w'), indent=1)
    jobs = []
    for p in plan:
        if p['x'] is None:
            continue
        f = byid[p['fam']]
        p['tag'] = f"{p['fam']}_n{p['target']:g}"
        p['cp'] = cp_of(f, p['x'])
        for L, ch in ((LS['s2'], CH_64), (LS['s1'], CH_32)):
            if not os.path.exists(os.path.join(HERE, 'scr', f"{p['tag']}_{L}.json")):
                jobs.append((fc_cmd(p['tag'], L, p['cp'], ch, 'scr'), a.threads, os.path.join(HERE, 'logs', f"scr_{p['tag']}_{L}.log")))
    for kind in ('xx', 'leg2'):
        for L in (LS['s1'], LS['s2']):
            jobs.append((['python3', 'ctrl183.py', kind, str(L), '--chis', *map(str, CH_64), '--outdir', 'ctrl'], a.threads,
                         os.path.join(HERE, 'logs', f"ctrl_{kind}_{L}.log")))
    jobs.sort(key=lambda j: -int(j[0][3]) if j[0][1] == 'fc.py' else 0)   # long (L = 64) first
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        list(ex.map(lambda j: run(*j), jobs))
    json.dump(dict(plan=plan, ladders={f['id']: f.get('ladder') for f in F}), open(os.path.join(HERE, 'plan.json'), 'w'), indent=1)
    print('phase 2 done', flush=True)

    # Phase 3: summary and selection (registered rule)
    rows = []
    for p in plan:
        if p.get('tag') is None:
            continue
        d32 = load(os.path.join(HERE, 'scr', f"{p['tag']}_{LS['s1']}.json")); d64 = load(os.path.join(HERE, 'scr', f"{p['tag']}_{LS['s2']}.json"))
        if not (d32 and d64):
            continue
        T = tunes[tkey(byid[p['fam']])]
        c2 = 6 * (d64['Sbar'] - d32['Sbar']) / math.log(LS['s2'] / LS['s1'])
        row = dict(tag=p['tag'], fam=p['fam'], group=p['group'], target=p['target'], x=p['x'], flag=p['flag'],
                   n_d=n_central(d64) - (T['n_vac'] or 0.0), n_on_note='n_d here is relative to the far vacuum; the target is relative to the onset', Sbar32=d32['Sbar'], Sbar64=d64['Sbar'], c2=c2,
                   cprof32=d32['c_prof'], cprof64=d64['c_prof'], trunc64=d64['trunc'], cp=p['cp'])
        row['passed'] = bool(row['cprof64'] >= 1.4 and c2 >= 1.3 and row['trunc64'] < 1e-6)
        row['score'] = min(c2, row['cprof64']) if row['trunc64'] < 1e-6 else -9
        rows.append(row)
    ctrl = {}
    for kind in ('xx', 'leg2'):
        a32 = load(os.path.join(HERE, 'ctrl', f"{kind}_{LS['s1']}.json")); a64 = load(os.path.join(HERE, 'ctrl', f"{kind}_{LS['s2']}.json"))
        if a32 and a64:
            ctrl[kind] = dict(c2=6 * (a64['Sbar'] - a32['Sbar']) / math.log(LS['s2'] / LS['s1']), cprof64=a64['c_prof'])
    par = sorted([r for r in rows if r['group'] in ('A', 'B')], key=lambda r: -r['score'])[:2]   # amended from 4 to 2
    gen = sorted([r for r in rows if r['group'] == 'C'], key=lambda r: -r['score'])[:1]
    dd = sorted([r for r in rows if r['group'] == 'D'], key=lambda r: -r['score'])[:1]   # addendum D
    forced = []   # family G: the best G point at r = 16 and at r = 12 (any v, offset)
    for r_ in (16.0, 12.0):
        cand = [r for r in rows if r['group'] == 'G' and byid[r['fam']]['r'] == r_]
        forced += sorted(cand, key=lambda r: -r['score'])[:1]
    for t_ in ('H_r16_p1.555_n0.1875', 'H_r16_p1.530_n0.1875', 'H_r24_p1.555_n0.1875', 'H_r24_p1.530_n0.1875'):
        forced += [r for r in rows if r['tag'] == t_]      # family H forced confirmations
    sel = par + gen + dd + [f for f in forced if f['tag'] not in {s_['tag'] for s_ in par + gen + dd}]
    json.dump(dict(rows=rows, controls=ctrl, selected=[s['tag'] for s in sel]), open(os.path.join(HERE, 'screen_summary.json'), 'w'), indent=1)
    with open(os.path.join(HERE, 'screen_summary.txt'), 'w') as fh:
        fh.write(f"controls: {json.dumps(ctrl)}\n")
        for r in sorted(rows, key=lambda r: -r['score']):
            fh.write(f"{r['tag']:28s} x={r['x']:+.3f} {r['flag']:8s} n_d={r['n_d']:.3f} S32={r['Sbar32']:.3f} S64={r['Sbar64']:.3f} "
                     f"c2={r['c2']:+.2f} cprof32={r['cprof32']:+.2f} cprof64={r['cprof64']:+.2f} trunc={r['trunc64']:.1e} "
                     f"{'PASS' if r['passed'] else ''}\n")
        fh.write('selected: ' + ' '.join(s['tag'] for s in sel) + '\n')
    print('phase 3 done, selected', [s['tag'] for s in sel], flush=True)

    # Phase 4: confirmation
    jobs = []
    t4 = a.phase4_threads
    for s in sel:
        tag, cp = s['tag'], s['cp']
        jobs.append((fc_cmd(tag, LS['c2'], cp, CH_256, 'conf'), t4, os.path.join(HERE, 'logs', f"conf_{tag}_{LS['c2']}.log")))
        jobs.append((['python3', 'vumps_fc.py', tag, '--cp', json.dumps(cp), '--n', '2', '--mode', 'damped', '--maxit', '800',
                      '--chis', *map(str, CH_V), '--out', f'vc/{tag}_n2.json', '--npz', f'vc/{tag}_n2.npz'], t4,
                     os.path.join(HERE, 'logs', f"vc_{tag}.log")))
        jobs.append((fc_cmd(tag, LS['c1'], cp, CH_128, 'conf'), t4, os.path.join(HERE, 'logs', f"conf_{tag}_{LS['c1']}.log")))
        if s['group'] not in ('C', 'G'):  # parity sectors for parity-class points (A, B, D)
            for L, ch in ((LS['p2'], CH_P128), (LS['p1'], CH_P64)):
                for par_ in ('odd', 'even'):
                    jobs.append((fc_cmd(f"{tag}_{par_}", L, cp, ch, 'conf', ('--parity', par_, '--nocorr')), t4,
                                 os.path.join(HERE, 'logs', f"conf_{tag}_{par_}_{L}.log")))
    for kind in ('xx', 'leg2'):
        for L in (LS['c2'], LS['c1']):
            jobs.append((['python3', 'ctrl183.py', kind, str(L), '--chis', *map(str, CH_128), '--outdir', 'ctrl'], t4,
                         os.path.join(HERE, 'logs', f"ctrl_{kind}_{L}.log")))
    def rank(j):
        c = j[0]
        if c[1] == 'vumps_fc.py':
            return 1
        if c[1] == 'ctrl183.py':
            return 2 if c[3] == str(LS['c2']) else 4
        if '--parity' in c:
            return 5
        return 0 if c[3] == str(LS['c2']) else 3
    jobs.sort(key=rank)
    with ThreadPoolExecutor(max_workers=a.phase4_workers) as ex:
        list(ex.map(lambda j: run(*j), jobs))
    print('phase 4 done', flush=True)


if __name__ == '__main__':
    main()
