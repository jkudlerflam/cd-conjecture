#!/usr/bin/env python3
"""Fits for contrib/d2_staggered.

Open chains (data/obc/*.json), per (point, L):
  c_obc    S(l) = (c/6) ln[(2L/pi) sin(pi l/L)] + a, fitted separately on even and odd bonds with
           L/8 <= l <= 7L/8; reported as the mean, with half the even/odd spread as a spread.
  K_F      charge variance of the block, F(l) = (K/(2 pi^2)) ln[(2L/pi) sin(pi l/L)] + a (same fit).
  gaps     D1 = E0(N0+1)-E0(N0), D2 = E0(N0+2)-E0(N0), Dn1, Dn2 = first two neutral gaps (sector N0),
           Dc1 = first gap inside the N0+1 sector.
  v_c      = 2 K_F L D1 / pi   (charge tower E(q) = (pi v/L) q^2/(2K) of one open-boundary boson)
  x_n1     = Dn1/(2 K_F D1)    neutral gap in units of the charged oscillator pi v_c/L
           (= 1 for a single boson; an extra neutral sector with speed v' and weight h gives h v'/v_c).
  c_size   6 dS(L/2)/d ln L between consecutive sizes.
FES (data/fes/*.json): c from S = (c/6) ln xi over consecutive chi pairs and a joint fit of the
largest four chi.
usage: analyze.py [--prefix stag] [--json summary.json]
Written for contrib/d2_staggered (Claude Opus 5.5 subagent, 2026-10-28).
"""
import os, sys, json, glob, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def chord_fit(y, L, lo=None, hi=None, pbc=False):
    y = np.asarray(y)
    l = np.arange(len(y))
    lo = lo or max(2, L // 8); hi = hi or L - max(2, L // 8)
    out = []
    for par in (0, 1):
        sel = (l >= lo) & (l <= hi) & (l % 2 == par)
        x = np.log(((1 if pbc else 2) * L / np.pi) * np.sin(np.pi * l[sel] / L))
        p = np.polyfit(x, y[sel], 1)
        out.append(p[0])
    return float(np.mean(out)), float(abs(out[0] - out[1]) / 2)


def obc_row(r):
    """Open chain: S = (c/6) ln chord, F = (K/2pi^2) ln chord, E(q) = (pi v/L) q^2/(2K), neutral unit pi v/L.
    Ring (bc = pbc; the MPS bond l is a two-endpoint interval of the ring): S = (c/3) ln chord,
    F = (K/pi^2) ln chord, E(q) = (2 pi v/L) q^2/(4K), neutral unit 2 pi v/L.  In both cases
    v = 2 K L D1/pi and x_n = Dn/(2 K D1).  Field names keep the c_obc label for both."""
    L = r['L']; s = r['sectors']; pbc = r.get('bc') == 'pbc'
    g0 = s['0'][0]
    E0 = g0['E']
    c, dc = chord_fit(g0['S'], L, pbc=pbc)
    k, dk = chord_fit(g0['F'], L, pbc=pbc)
    fc = 3 if pbc else 6; fk = np.pi ** 2 if pbc else 2 * np.pi ** 2
    K = fk * k; dK = fk * dk
    row = dict(L=L, E0=E0, c_obc=fc * c, c_obc_spread=fc * dc, K_F=K, K_F_spread=dK, S_mid=g0['S'][L // 2],
               trunc=max(h['trunc'] for x in s.values() for y in x for h in y['hist'][-1:]),
               chi=max(h['maxchi'] for x in s.values() for y in x for h in y['hist'][-1:]))
    ex = [x['E'] - E0 for x in s['0'][1:]]
    row['Dn'] = ex
    if '1' in s:
        row['D1'] = s['1'][0]['E'] - E0
        row['Dc'] = [x['E'] - s['1'][0]['E'] for x in s['1'][1:]]
    if '2' in s:
        row['D2'] = s['2'][0]['E'] - E0
    if 'D1' in row and row['D1'] > 0:
        row['v_c'] = 2 * K * L * row['D1'] / np.pi
        if 'D2' in row:
            row['D2_over_D1'] = row['D2'] / row['D1']
        if ex:
            unit = (4 if pbc else 2) * K * row['D1']
            row['x_n'] = [e / unit for e in ex]
            row['x_c'] = [e / unit for e in row.get('Dc', [])]
            row['x_n1'] = ex[0] / unit
            row['v_n'] = L * ex[0] / np.pi / (2 if pbc else 1)
            if len(ex) > 1:
                row['x_n2'] = ex[1] / unit
    if ex:
        row['L_Dn1'] = L * ex[0]
    return row


def load_obc(pattern, sub='obc'):
    groups = {}
    for f in sorted(glob.glob(os.path.join(HERE, 'data', sub, pattern))):
        r = json.load(open(f))
        if '0' not in r['sectors'] or 'S' not in r['sectors']['0'][0]:
            continue
        key = (r['name'], round(r['u'], 6), round(r['b'], 6))
        groups.setdefault(key, {})[r['L']] = obc_row(r)
    out = []
    for key, rows in sorted(groups.items()):
        Ls = sorted(rows)
        for a, b in zip(Ls[:-1], Ls[1:]):
            fc = 3 if sub == 'pbc' else 6
            rows[b]['c_size'] = fc * (rows[b]['S_mid'] - rows[a]['S_mid']) / np.log(b / a)
        out.append(dict(name=key[0], u=key[1], b=key[2], rows=[rows[L] for L in Ls]))
    return out


def fes_fit(path):
    d = json.load(open(path))
    R = [x for x in d['results'] if x['err'] < 1e-5]
    S = np.array([np.mean(x['S']) for x in R]); xi = np.array([x['xi'] for x in R])
    chi = [x['chi'] for x in R]
    pairs = [(chi[i], chi[i + 1], float(6 * (S[i + 1] - S[i]) / np.log(xi[i + 1] / xi[i])))
             for i in range(len(R) - 1)]
    joint = float(6 * np.polyfit(np.log(xi[-4:]), S[-4:], 1)[0]) if len(R) >= 4 else None
    return dict(name=d['name'], u=d['u'], b=d['b'], chi=chi, S=S.tolist(), xi=xi.tolist(),
                e=[x['e'] for x in R], c_pairs=pairs, c_joint4=joint)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--pattern', default='[gsti]*.json')
    ap.add_argument('--json', default=os.path.join(HERE, 'summary.json'))
    ap.add_argument('--sub', default='obc')
    a = ap.parse_args()
    obc = load_obc(a.pattern, a.sub)
    fes = [fes_fit(f) for f in sorted(glob.glob(os.path.join(HERE, 'data', 'fes', '*.json')))]
    json.dump(dict(obc=obc, fes=fes), open(a.json, 'w'), indent=1)
    for g in obc:
        print(f"{g['name']} u={g['u']} b={g['b']}")
        for r in g['rows']:
            print('  L=%3d c_obc=%.4f(%.4f) c_size=%s K_F=%.3f D2/D1=%s x_n1=%s x_n2=%s v_c=%s v_n=%s LDn1=%s chi=%d tr=%.1e' % (
                r['L'], r['c_obc'], r['c_obc_spread'], '%.4f' % r['c_size'] if 'c_size' in r else '-', r['K_F'],
                '%.3f' % r['D2_over_D1'] if 'D2_over_D1' in r else '-', '%.3f' % r['x_n1'] if 'x_n1' in r else '-',
                '%.3f' % r['x_n2'] if 'x_n2' in r else '-', '%.3f' % r['v_c'] if 'v_c' in r else '-',
                '%.3f' % r['v_n'] if 'v_n' in r else '-', '%.3f' % r['L_Dn1'] if 'L_Dn1' in r else '-', r['chi'], r['trunc']))
    for f in fes:
        print(f"FES {f['name']} u={f['u']} b={f['b']} c_joint4={f['c_joint4']} pairs=" +
              ' '.join(f'{p[2]:.3f}' for p in f['c_pairs']))
