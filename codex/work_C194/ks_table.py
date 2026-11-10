"""C194 K_s(n) on the gapless side, frozen estimator (RESULT.md "Frozen after the controls").
K_s = 1/(4 eta_b - eta_pair); eta_b = boson exponent at the valley momentum (offline.py's fit, both valleys
averaged), eta_pair = mean of the d = 2 and d = 3 inter-valley pair exponents at momentum pi (open-chain forms).
Central value window 4, error = spread of K_s over windows 2 to 4.
Inputs: final records with 'meas' in out_main/ (and typhon_out/out_fin/ if present), and C176's sym runs at 1/8,
5/32, 3/16. Writes ks_table.json.  The fitting code is offline.py's, imported as functions here (same numbers)."""
import json, glob, os, re, sys
import numpy as np

d = lambda x, L: np.sin(np.pi * x / (2 * L))


def sm(h, n):
    for _ in range(n):
        h = np.convolve(h, [.25, .5, .25], mode='same')
    return h


def peaks(rs, g, npk, lo=-np.pi, hi=np.pi):
    ks = np.linspace(lo, hi, 8001)[:-1]; w = np.cos(np.pi * (rs - rs[0]) / (2 * (rs[-1] - rs[0] + 1))) ** 2
    F = np.abs(np.exp(-1j * np.outer(ks, rs)) @ (g * w)); out = []
    for _ in range(npk):
        i = np.argmax(F); out.append((ks[i], F[i])); dk = np.abs(np.angle(np.exp(1j * (ks - ks[i])))); F = np.where(dk < 0.6, 0, F)
    return out


def fit(x, y):
    p = np.polyfit(x, y, 1); return -p[0], np.std(y - np.polyval(p, x))


def exponents(m, L):
    """offline.py's boson and pair fits (windows 1..4); returns eta_b[4] (valley mean), eta_pair[4] (d=2,3 mean)."""
    G = np.array(m['Gb_re']) + 1j * np.array(m['Gb_im']); W = G.shape[0]; c0 = W // 2; nc = 6
    rs = np.arange(1, W - 2 * nc - 2)
    g = np.array([np.mean([G[c - r // 2, c - r // 2 + r] for c in range(c0 - nc, c0 + nc + 1)]) for r in rs])
    pk = peaks(rs, g, 2)
    etas = []
    for (k, amp) in pk:
        hh = sm(g * np.exp(-1j * k * rs), 4); sl = (rs >= 8) & (rs <= W // 2)
        k = k + np.polyfit(rs[sl], np.unwrap(np.angle(hh[sl])), 1)[0]
        h = np.abs(sm(g * np.exp(-1j * k * rs), 4)); x = np.log(np.tan(np.pi * rs / (2 * L))); y = np.log(h)
        etas.append([fit(x[(rs >= lo) & (rs <= hi)], y[(rs >= lo) & (rs <= hi)])[0]
                     for lo, hi in ((6, W // 4), (10, W // 3), (16, W // 2), (24, int(.7 * W)))])
    eb = np.mean(etas, 0)
    pe = []
    for sep in ('2', '3'):
        dat = []
        for key, v in m['Dpair'].items():
            if not key.startswith(sep + '_'):
                continue
            x0 = int(key.split('_')[1]) + 0.5 + int(sep) / 2; a = np.array(v); D = a[:, 0] + 1j * a[:, 1]
            r = np.arange(len(D)) + int(sep) + 1
            dat.append((x0, r, D))
        n0 = min(len(D) for _, _, D in dat)
        Q = peaks(dat[0][1][:n0], np.mean([D[:n0] for _, _, D in dat], 0), 1, lo=np.pi - 1, hi=np.pi + 1)[0][0]
        res = []
        for lo, hi in ((6, L // 8), (10, L // 5), (16, L // 3), (24, int(.4 * L))):
            X = []; Y = []
            for x0, r, D in dat:
                h = np.abs(sm(D * np.exp(-1j * Q * r), 2)); F = d(r, L) * d(2 * x0 + r, L) / np.sqrt(d(2 * x0, L) * d(2 * (x0 + r), L))
                s = (r >= lo) & (r <= hi); X += list(np.log(F[s])); Y += list(np.log(h[s]))
            res.append(fit(np.array(X), np.array(Y))[0])
        pe.append(res)
    return eb, np.mean(pe, 0), [p[0][0] for p in [pk]][0]


def last_meas(f):
    recs = [json.loads(l) for l in open(f) if l.strip()]
    for r in reversed(recs):
        if 'meas' in r and r['meas'] and 'Gb_re' in r['meas']:
            return r
    return None


def main():
    files = []
    for dd in ('out_main', 'typhon_out/out_fin'):
        files += [(f, None) for f in glob.glob(os.path.join(dd, '*.jsonl'))]
    for key, nom in (('sym18', '1/8'), ('sym532', '5/32'), ('sym316', '3/16')):
        files += [(f, nom) for f in glob.glob(f'../work_C176/typhon_out/out/{key}_L*.jsonl')]
    rows = {}
    for f, nom in files:
        r = last_meas(f)
        if r is None:
            continue
        c = r['cfg']; L = c['L']; nom = nom or c['nom']
        if L < 120:
            continue
        eb, ep, _ = exponents(r['meas'], L)
        ks = 1 / (4 * eb - ep)
        row = dict(L=L, N=c['N'], eta_b=float(eb[3]), eta_pair=float(ep[3]), Ks=float(ks[3]),
                   Ks_err=float(np.max(ks[1:]) - np.min(ks[1:])) / 2, Ks_win=[float(x) for x in ks],
                   src='C176' if 'work_C176' in f else os.path.dirname(f))
        old = rows.setdefault(nom, {}).get(L)
        if old is None or 'out_fin' in row['src']:
            rows[nom][L] = row
    out = {}
    for nom in sorted(rows, key=lambda s: eval(s)):
        print(f"=== n = {nom} ({eval(nom):.4f})")
        for L in sorted(rows[nom]):
            r = rows[nom][L]
            print(f"  L={L:4d} N={r['N']:3d} eta_b={r['eta_b']:.4f} eta_pair={r['eta_pair']:.4f} K_s={r['Ks']:.3f} +- {r['Ks_err']:.3f} "
                  f"(windows 1-4: {' '.join(f'{x:.3f}' for x in r['Ks_win'])}) [{r['src']}]")
        Ls = sorted(rows[nom])
        if len(Ls) >= 2:
            a, b = rows[nom][Ls[0]], rows[nom][Ls[-1]]
            print(f"  drift K_s({Ls[-1]}) - K_s({Ls[0]}) = {b['Ks'] - a['Ks']:+.3f}")
        out[nom] = dict(n=eval(nom), rows=[rows[nom][L] for L in Ls])
    json.dump(out, open('ks_table.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
