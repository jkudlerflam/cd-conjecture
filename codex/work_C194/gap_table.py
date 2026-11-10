"""C194: pair-breaking gap table and the registered extrapolation (fitgap.fit_gap).
Usage (from work_C194): python3 -E -s gap_table.py [outdir ...]   default out_main and typhon_out/out_fin if present.
For each density and L: Delta_pb = E(N+1) + E(N-1) - 2E(N) from the final energies (same schedule for the three
states), with the convergence residuals (last-sweep change, and change from the last chi = 128 sweep to the final
sweep when the schedule goes beyond 128). Densities 1/8, 5/32, 3/16 are read from C176's converged runs
(../work_C176/typhon_out/out/), at the largest chi common to N and N +- 1.
Writes gap_table.json."""
import json, glob, os, sys, re
import numpy as np
from fitgap import fit_gap


def load_dir(d):
    out = {}
    for f in glob.glob(os.path.join(d, '*.jsonl')):
        recs = [json.loads(l) for l in open(f) if l.strip()]
        fin = [r for r in recs if r.get('kind') == 'final']
        if not fin:
            continue
        sw = [r for r in recs if r.get('kind') == 'sweep']
        c = fin[0]['cfg']
        eprev = [r['E'] for r in sw if r['chi'] < sw[-1]['chi']]
        out[c['tag']] = dict(nom=c['nom'], L=c['L'], N=c['N'], E=fin[0]['E'], dlast=sw[-1]['E'] - sw[-2]['E'],
                             dchi=sw[-1]['E'] - eprev[-1], chi=sw[-1]['chi'],
                             trunc=sw[-1]['trunc'])
    return out


def c176():
    base = '../work_C176/typhon_out/out/'
    res = {}
    for key, nom in (('sym18', '1/8'), ('sym532', '5/32'), ('sym316', '3/16')):
        for f in glob.glob(base + key + '_L*.jsonl'):
            L = int(re.search(r'_L(\d+)\.jsonl', f).group(1))
            fm, fp = base + f'{key}m1_L{L}.jsonl', base + f'{key}p1_L{L}.jsonl'
            if not (os.path.exists(fm) and os.path.exists(fp)):
                continue
            R = [[json.loads(l) for l in open(x) if l.strip()] for x in (f, fm, fp)]
            chi = min(r[-1]['chi'] for r in R)
            Es = []
            for recs in R:
                sel = [r for r in recs if r['chi'] == chi]
                prev = [r for r in recs if r['chi'] < chi]
                Es.append((sel[-1]['E'], sel[-1]['E'] - sel[-2]['E'] if len(sel) > 1 else np.nan,
                           sel[-1]['E'] - prev[-1]['E'], recs[0]['cfg']['N'], sel[-1]['trunc']))
            N = Es[0][3]
            d = Es[2][0] + Es[1][0] - 2 * Es[0][0]
            res.setdefault(nom, []).append(dict(L=L, N=N, D=d, dlast=max(abs(e[1]) for e in Es),
                                                dchi=max(abs(e[2]) for e in Es), chi=chi, src='C176'))
    return res


def main(dirs):
    st = {}
    for d in dirs:
        if os.path.isdir(d):
            for k, v in load_dir(d).items():
                key = (v['nom'], v['L'], v['N'])
                if key not in st or v['chi'] > st[key]['chi']:
                    st[key] = v
    table = {}
    for (nom, L, N), v in st.items():
        if v['N'] % 2:
            continue
        m, p = st.get((nom, L, N - 1)), st.get((nom, L, N + 1))
        if not (m and p) or not (m['chi'] == v['chi'] == p['chi']):
            continue
        d = p['E'] + m['E'] - 2 * v['E']
        table.setdefault(nom, []).append(dict(L=L, N=N, D=d, dlast=max(abs(x['dlast']) for x in (v, m, p)),
                                              dchi=max(abs(x['dchi']) for x in (v, m, p)) if v['dchi'] is not None else None,
                                              chi=v['chi'], src='C194'))
    for nom, rows in c176().items():
        table.setdefault(nom, []).extend(rows)
    out = {}
    for nom in sorted(table, key=lambda s: eval(s)):
        rows = sorted(table[nom], key=lambda r: r['L'])
        n = eval(nom)
        print(f"=== n = {nom} ({n:.4f})")
        for r in rows:
            print(f"  L={r['L']:4d} N={r['N']:3d} chi={r['chi']:4d} Dpb={r['D']:.6f} L*Dpb={r['L']*r['D']:.4f} "
                  f"max|last-sweep dE|={r['dlast']:.1e} max|dE(chi)|={r['dchi'] if r['dchi'] is None else format(r['dchi'], '.1e')} [{r['src']}]")
        if len(rows) >= 3:
            f = fit_gap([r['L'] for r in rows], [r['D'] for r in rows])
            print(f"  fit D + a/L + b/L^2: D_inf = {f['D']:+.5f} +- {f['err']:.5f} ({f['D']/f['err']:+.1f} err; se {f['se']:.1e}, "
                  f"spread {f['spread']:.1e}), a = {f['a']:.3f}, b = {f['b']:.2f}")
            # largest three sizes only (sub-fit check, not registered)
            if len(rows) >= 4:
                g = fit_gap([r['L'] for r in rows[-3:]], [r['D'] for r in rows[-3:]])
                print(f"  (check, largest three L) D_inf = {g['D']:+.5f} +- {g['err']:.5f}")
            out[nom] = dict(n=n, rows=rows, fit=f)
        else:
            out[nom] = dict(n=n, rows=rows)
    json.dump(out, open('gap_table.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main(sys.argv[1:] or ['out_main', 'typhon_out/out_fin'])
