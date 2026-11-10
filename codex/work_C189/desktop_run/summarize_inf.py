"""C189 infinite-chain summary (registered estimators, RESULT.md Part 2 and its amendment A2).
Usage: python3 summarize_inf.py <dir with inf_*.jsonl> [tag prefix, default inf_; use loc_ for local/out]
S_bar = mean of the two cell-boundary entropies of consecutive iterations (the converged state can be 2q periodic).
xi_X = per-site correlation length of the leading eigenvalue in sector X (T0: dQ = 0, second eigenvalue; T1b, T1f:
dQ = 1 without and with the fermion string; T2: dQ = 2), xi = max over the sectors.
c_fit = 6 x least-squares slope of S_bar against ln xi over the largest four chi.
"""
import json, glob, os, sys
import numpy as np

EXACT = {'inf_c1': -(2 / np.pi) * 0.75 * np.sin(np.pi / 3)}


def ff_exact(n):
    k = np.linspace(-np.pi, np.pi, 400000, endpoint=False)
    e = np.sort(2 * np.cos(k) + 2 * np.cos(2 * k))
    return e[:int(round(len(k) * n))].sum() / len(k)


def load(d, prefix='inf_'):
    out = {}
    for f in sorted(glob.glob(os.path.join(d, prefix + '*.jsonl'))):
        if f.endswith('.iter.jsonl'):
            continue
        recs = [json.loads(l) for l in open(f) if l.strip()]
        if recs:
            out[recs[0]['tag']] = recs
    return out


def xis(rec):
    m = rec['meas']
    o = {}
    if 'T0' in m and len(m['T0']) > 1:
        o['xi0'] = m['T0'][1][2]
        o['k0'] = m['T0'][1][3]
    for nm in ('T1b', 'T1f', 'T2', 'T0f'):
        if nm in m and m[nm]:
            o['xi' + nm[1:]] = m[nm][0][2]
            o['k' + nm[1:]] = m[nm][0][3]
    vals = [v for k, v in o.items() if k.startswith('xi') and k != 'xi0f' and np.isfinite(v)]
    o['xi'] = max(vals) if vals else np.nan
    return o


def cfit(Ls, Ss):
    A = np.vstack([np.log(Ls), np.ones(len(Ls))]).T
    return 6 * np.linalg.lstsq(A, np.asarray(Ss), rcond=None)[0][0]


def main(d, prefix='inf_'):
    np.seterr(all='ignore')
    runs = load(d, prefix)
    for tag, recs in runs.items():
        c = recs[0]['cfg']
        n = c['p'] / c['q']
        ex = EXACT.get('inf_' + tag[len(prefix):], ff_exact(n) if c.get('fermionic') else np.nan)
        print(f"=== {tag}  n = {c['p']}/{c['q']}  exact e = {ex:.10f}")
        rows = []
        for r in sorted(recs, key=lambda r: r['chi']):
            x = xis(r)
            Sb = 0.5 * (r['S'] + r['S_prev'])
            rows.append((r['chi'], Sb, x, r))
            dens = r['meas'].get('dens', [])
            print(f"chi={r['chi']:5d} it={r['conv']['iters']:3d} e={r['e_site']:.10f} (e-ex={r['e_site'] - ex:.1e}) "
                  f"S={r['S']:.4f}/{r['S_prev']:.4f} Sbar={Sb:.4f} trunc={r['trunc']:.1e} "
                  f"xi0={x.get('xi0', np.nan):.1f} xi1b={x.get('xi1b', np.nan):.1f} xi1f={x.get('xi1f', np.nan):.1f} "
                  f"xi2={x.get('xi2', np.nan):.1f} lam0={r['meas']['T0'][0][0] if 'T0' in r['meas'] else np.nan:.5f} "
                  f"dens[min,max]={min(dens) if dens else np.nan:.4f},{max(dens) if dens else np.nan:.4f} t={r['secs']:.0f}s")
        if len(rows) >= 2:
            print("  local c between consecutive chi (xi = max sector):",
                  ' '.join(f"{rows[i][0]}->{rows[i + 1][0]}:{6 * (rows[i + 1][1] - rows[i][1]) / np.log(rows[i + 1][2]['xi'] / rows[i][2]['xi']):.3f}"
                           for i in range(len(rows) - 1)))
            print("  local c with xi0:",
                  ' '.join(f"{rows[i][0]}->{rows[i + 1][0]}:{6 * (rows[i + 1][1] - rows[i][1]) / np.log(rows[i + 1][2]['xi0'] / rows[i][2]['xi0']):.3f}"
                           for i in range(len(rows) - 1)))
        if len(rows) >= 4:
            last = rows[-4:]
            cf = cfit([r[2]['xi'] for r in last], [r[1] for r in last])
            print(f"  c_fit (largest four chi, xi = max) = {cf:.3f}")
            # registered readings
            byc = {r[0]: r for r in rows}
            top = rows[-1][0]
            half = top // 2
            if half in byc:
                a, b = byc[half][2], rows[-1][2]
                # numpy floats, so a frozen run with zero correlation length prints nan instead of raising
                a = {k: np.float64(v) for k, v in a.items()}
                b = {k: np.float64(v) for k, v in b.items()}
                g1 = b.get('xi1f', np.nan) / a.get('xi1f', np.nan)
                g2 = b.get('xi2', np.nan) / a.get('xi2', np.nan)
                g0 = b.get('xi0', np.nan) / a.get('xi0', np.nan)
                rat = (b['xi1f'] / b['xi2']) / (a['xi1f'] / a['xi2'])
                print(f"  last doubling {half}->{top}: growth xi1f x{g1:.3f}, xi2 x{g2:.3f}, xi0 x{g0:.3f}; "
                      f"(xi1f/xi2) changes by x{rat:.3f}")
                liquid = cf >= 1.6 and (1 / 1.5 < rat < 1.5)
                gap = (g1 < 1.10 and (g2 > 1.30 or g0 > 1.30)) or cf < 1.4
                print(f"  registered reading: c = 2 liquid: {liquid}; pseudospin gap: {gap}; "
                      f"{'inconclusive' if liquid == gap else ''}")


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'out', sys.argv[2] if len(sys.argv) > 2 else 'inf_')
