"""C193 estimators (Claude Opus 5.5 agent, 2026-11-10).  Run from codex/work_C193.

E1  one-body exponent from a stored window matrix G (Gb or Gf): pairs symmetric about the window centre,
    averaged over 2*NC+1 centres, demodulated at a peak momentum found from the data (refined by the
    phase slope), low-passed with NSM binomial passes, |h| fit to F^(-eta) with the open-chain form
    F(x, y) = f(x - y) f(x + y) / sqrt(f(2x) f(2y)), f(u) = sin(pi u / (2 (L + 1))), x = site + 1.
E2  the C_P correlator from site i0: dominant momentum from the tapered Fourier transform, demodulated,
    low-passed, fit to the same form.
Usage: python3 c193.py controls | python3 c193.py candidates
"""
import json
import sys
import glob
import os
import numpy as np

BASE = '../work_C190'
GPU = BASE + '/gpu/typhon_out/outG/'
NC = 4
NSM = 4


def f(u, L):
    return np.sin(np.pi * u / (2 * (L + 1)))


def Fform(x, y, L):
    return f(np.abs(x - y), L) * f(x + y, L) / np.sqrt(f(2 * x, L) * f(2 * y, L))


def sm(h, n=NSM):
    for _ in range(n):
        h = np.convolve(h, [.25, .5, .25], mode='same')
    return h


def peaks(rs, g, npk, lo=-np.pi, hi=np.pi, excl=0.6):
    ks = np.linspace(lo, hi, 8001)[:-1]
    w = np.cos(np.pi * (rs - rs[0]) / (2 * (rs[-1] - rs[0] + 1))) ** 2
    F = np.abs(np.exp(-1j * np.outer(ks, rs)) @ (g * w))
    out = []
    for _ in range(npk):
        i = np.argmax(F)
        out.append((ks[i], F[i]))
        dk = np.abs(np.angle(np.exp(1j * (ks - ks[i]))))
        F = np.where(dk < excl, 0, F)
    return out


def fit(X, Y):
    p = np.polyfit(X, Y, 1)
    return -p[0], float(np.std(Y - np.polyval(p, X)))


def load(fn):
    d = json.load(open(fn))
    assert d['status'] == 'done', fn
    return d


def energy(d):
    e = d.get('E_phys')
    return e if e is not None else d['hist'][-1]['E_phys']


def onebody_seq(G, L, i_lo):
    W = G.shape[0]
    c0 = W // 2
    rmax = 2 * (W // 2 - 1 - NC)
    rs = np.arange(1, rmax + 1)
    g = np.array([np.mean([G[c - r // 2, c - r // 2 + r] for c in range(c0 - NC, c0 + NC + 1)]) for r in rs])
    xs = i_lo + c0 - rs // 2 + 1.0
    ys = xs + rs
    return rs, g, Fform(xs, ys, L), rmax


def E1(G, L, i_lo, ranges, npk=2, rel=0.3):
    rs, g, F, rmax = onebody_seq(G, L, i_lo)
    pk = peaks(rs, g, npk)
    pk = [p for p in pk if p[1] >= rel * pk[0][1]]
    res = []
    for k, amp in pk:
        hh = sm(g * np.exp(-1j * k * rs))
        s = (rs >= 6) & (rs <= rmax - NSM)
        k = k + np.polyfit(rs[s], np.unwrap(np.angle(hh[s])), 1)[0]
        h = np.abs(sm(g * np.exp(-1j * k * rs)))
        ets = [fit(np.log(F[(rs >= lo) & (rs <= hi)]), np.log(h[(rs >= lo) & (rs <= hi)])) for lo, hi in ranges]
        res.append(dict(k=float(np.angle(np.exp(1j * k))), amp=float(amp), eta=[e for e, _ in ets],
                        rms=[s_ for _, s_ in ets]))
    return res, rmax


def E2(m, L, ranges, qwin=None):
    cp = np.array(m['CP'])
    r = cp[:, 0]
    v = cp[:, 1]
    i0 = m['CP_i0']
    s0 = r >= 3
    if qwin is None:
        q = peaks(r[s0], v[s0], 1, lo=0.0, hi=np.pi, excl=0.3)[0][0]
    else:
        q = peaks(r[s0], v[s0], 1, lo=qwin[0], hi=qwin[1], excl=0.3)[0][0]
    if 0.05 < q < np.pi - 0.05:
        hh = sm(v * np.exp(-1j * q * r))
        s = (r >= 8) & (r <= r.max() - NSM)
        q = q + np.polyfit(r[s], np.unwrap(np.angle(hh[s])), 1)[0]
        h = np.abs(sm(v * np.exp(-1j * q * r)))
    else:
        qq = 0.0 if q < 1 else np.pi
        h = np.abs(sm(v * np.cos(qq * r)))
        q = qq
    F = Fform(i0 + 1.0, i0 + 1.0 + r, L)
    out = []
    for lo, hi in ranges:
        s = (r >= lo) & (r <= hi) & (h > 0)
        out.append(fit(np.log(F[s]), np.log(h[s])))
    return dict(q=float(q), eta=[e for e, _ in out], rms=[s_ for _, s_ in out], sign=float(np.sign(np.mean(v[(r >= 8) & (r <= 30)]))))


# frozen ranges (see RESULT.md calibration section)
def ranges_E1(W):
    rmax = 2 * (W // 2 - 1 - NC) - NSM
    if W == 48:
        return [(4, 12), (6, 20), (8, 28), (12, rmax)]
    return [(4, 8), (5, 12), (6, 16), (8, rmax)]


def ranges_E2(L):
    return [(6, L // 8), (10, L // 5), (16, L // 3), (24, int(.45 * L))]


def summarize(vals):
    a = np.array(vals)
    return float(a.mean()), float((a.max() - a.min()) / 2)


def run_file(fn, gf=True, qwin=None):
    d = load(fn)
    m = d['meas']
    L = d['cfg']['L']
    i_lo, W = m['gwin']
    out = dict(L=L, P=d['cfg']['P'], E=energy(d), nmid=float(np.mean(m['dens'][L // 4:3 * L // 4])))
    Gb = np.array(m['Gb_re']) + 1j * np.array(m['Gb_im'])
    out['b'], out['rmax'] = E1(Gb, L, i_lo, ranges_E1(W))
    if gf:
        Gf = np.array(m['Gf_re']) + 1j * np.array(m['Gf_im'])
        out['f'], _ = E1(Gf, L, i_lo, ranges_E1(W))
    out['cp'] = E2(m, L, ranges_E2(L), qwin)
    return out


def fmt(o):
    s = f"L={o['L']} P={o['P']} nmid={o['nmid']:.4f}\n"
    for lab in ('b', 'f'):
        for p in o.get(lab, []):
            s += f"  G{lab} k={p['k']:+.4f} amp={p['amp']:.3g} eta=" + ' '.join(f'{e:.4f}' for e in p['eta']) + \
                 '  rms=' + ' '.join(f'{x:.3f}' for x in p['rms']) + '\n'
    c = o['cp']
    s += f"  CP q={c['q']:.4f} sign={c['sign']:+.0f} eta=" + ' '.join(f'{e:.4f}' for e in c['eta']) + \
         '  rms=' + ' '.join(f'{x:.3f}' for x in c['rms']) + '\n'
    return s


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'controls':
        files = [GPU + 'xx2_L64_P0_rand1.json', GPU + 'xx2_L128_P0_rand1.json', GPU + 'xx2_L256_P0_rand1.json',
                 GPU + 'xx1_L64_P1_rand1.json', GPU + 'xx1_L128_P1_rand1.json', BASE + '/typhon_out/outB/xx1_L256_P0_rand1.json']
        res = {}
        for fn in files:
            o = run_file(fn)
            res[os.path.basename(fn)] = o
            print(os.path.basename(fn), fmt(o))
        json.dump(res, open('controls.json', 'w'), indent=1)

    if mode == 'candidates':
        IE1, IE2 = (2, 3), (1, 2, 3)          # frozen
        INVKC, BAND = 0.62, (0.55, 0.64)
        keys = ['p153c', 'sym40_n316', 'rat_n316', 'sym40_n18', 'p153_n316', 'sym80_n316']
        allres = {}
        for key in keys:
            for L in (64, 128, 256):
                sec = {}
                for P in (0, 1):
                    fns = sorted(glob.glob(GPU + f'{key}_L{L}_P{P}_*.json'))
                    ds = [(energy(load(fn)), fn) for fn in fns]
                    E, fn = min(ds)
                    o = run_file(fn)
                    eb = [p['eta'][i] for p in o['b'] for i in IE1]
                    o['eta_b'], o['deta_b'] = summarize(eb)
                    o['deta_b'] = float(np.hypot(o['deta_b'], 0.005))
                    o['eta_cp'], o['deta_cp'] = summarize([o['cp']['eta'][i] for i in IE2])
                    o['eta_f'] = summarize([p['eta'][i] for p in o['f'] for i in IE1])[0]
                    x = 4 * o['eta_b'] - INVKC
                    o['Ks'] = 1 / x
                    o['dKs'] = float(4 * o['deta_b'] / x ** 2)
                    o['Ks_band'] = [1 / (4 * o['eta_b'] - b) for b in BAND]
                    o['file'] = os.path.basename(fn)
                    o['dE_inits'] = float(max(e for e, _ in ds) - E)
                    sec[P] = o
                gb = min(sec, key=lambda p: sec[p]['E'])
                allres[f'{key}_L{L}'] = dict(sec={str(p): sec[p] for p in sec}, ground=gb)
                for P in (0, 1):
                    o = sec[P]
                    ks = sorted(abs(p['k']) for p in o['b'])
                    dq = abs(np.angle(np.exp(1j * (o['b'][0]['k'] - o['b'][1]['k'])))) if len(o['b']) > 1 else float('nan')
                    print(f"{key:11s} L={L:3d} P={P}{'*' if P == gb else ' '} nmid={o['nmid']:.4f} "
                          f"k={','.join(f'{k:.3f}' for k in ks)} |dk|={dq:.3f} "
                          f"eta_b={o['eta_b']:.4f}({o['deta_b']:.4f}) Ks={o['Ks']:.3f}({o['dKs']:.3f}) "
                          f"band={o['Ks_band'][0]:.2f}-{o['Ks_band'][1]:.2f} eta_f={o['eta_f']:.3f} "
                          f"CP q={o['cp']['q']:.3f} eta={o['eta_cp']:.3f}({o['deta_cp']:.3f}) rms={max(o['cp']['rms'][1:]):.2f} "
                          f"| R3,R4 per valley " + ' '.join(f"{p['eta'][2]:.4f},{p['eta'][3]:.4f}" for p in o['b']))
        json.dump(allres, open('candidates.json', 'w'), indent=1)
