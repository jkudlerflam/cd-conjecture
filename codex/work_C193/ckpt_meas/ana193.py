#!/usr/bin/env python3
"""C193b analysis (Claude Opus 5.5 agent, 2026-11-10).  Run from codex/work_C193/ckpt_meas.

Estimators (PROOF.md §3, §4), all at fixed r/L:
  EB  boson exponent eta_b.  g(r) = mean of <sm_i sp_j> over pairs with j - i = r and |i + j - (L - 1)| <= 8
      (centres within 4 sites of the chain centre), demodulated at the two largest Fourier peaks (refined by the
      phase slope), NSM binomial low-pass passes, |h| fit to F^(-eta), F the open-chain form at the centred pair.
  EP  pair exponent.  Connected D(i0, j) = <sm_i0 sm_i0+d sp_j sp_j+d> - <sm sm><sp sp>, rows from the four
      reference sites, separation r = j - i0, demodulated at the dominant momentum Q (found from the data),
      low-passed, all rows fit jointly to F(x1, x2)^(-eta) with x1, x2 the pair centres (site + 1 convention).
      Primary d = 3 (the separation at which the two-leg control's pair is the inter-leg pair), d = 2 reported.
  ES  structure factor.  C(r) = mean over the reference rows of the connected <n n>, S(q) = 2 sum_{r>=1} w C cos qr
      with w = cos^2(pi r / 2R), r <= R, fit a + b q + c q^2 on q in [lo, hi] pi/R; K_c = pi b (C192 sq2.py).
Usage
  python3 -E -s ana193.py controls [<meas dir>]          trial ranges on the two controls
  python3 -E -s ana193.py all [<meas dir>] [--eff]       frozen ranges on everything, plus the C176 effective model
  python3 -E -s ana193.py states <file_meas.json> ...    frozen estimators on any measured states (new runs)
"""
import os
import sys
import json
import glob
import numpy as np

NSM = 4
DC = 8                       # |i + j - (L - 1)| <= DC
EFF_DIR = '../../work_C176/typhon_out/out/'


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


def linfit(X, Y):
    p = np.polyfit(X, Y, 1)
    return -p[0], float(np.std(Y - np.polyval(p, X)))


# ---------------------------------------------------------------------------------------------------------------
# data access: a common dict form for the qubit outputs and the C176 effective-model records
# ---------------------------------------------------------------------------------------------------------------
def from_qubit(fn):
    d = json.load(open(fn))
    m = d['meas']
    L = m['L']
    G = {}
    for i_s, (re, im) in m['Gb_rows'].items():
        i = int(i_s)
        for t, (a, b) in enumerate(zip(re, im)):
            G[(i, i + 1 + t)] = a + 1j * b
    D = {}
    for key, row in m['Dpair'].items():
        dd, i0 = map(int, key.split('_'))
        mmv = complex(*m['mm'][key])
        pp = m['pp'][str(dd)]
        vals = []
        for t, (a, b) in enumerate(row):
            j = i0 + dd + 1 + t
            vals.append(a + 1j * b - mmv * complex(*pp[str(j)]))
        D[(dd, i0)] = np.array(vals)
    C = {int(k): np.array(v) for k, v in m['Cnn'].items()}
    hist = d.get('cfg', {})
    E = d['check']['E_stored']
    return dict(L=L, G=G, D=D, C=C, dens=np.array(m['dens']), E=E, P=d['cfg']['P'], tag=d['tag'], check=d['check'],
                Dmax=d.get('Dmax'))


def from_eff(tag):
    rec = json.loads(open(EFF_DIR + f'{tag}.jsonl').read().strip().split('\n')[-1])
    m = rec['meas']
    L = rec['cfg']['L']
    Gm = np.array(m['Gb_re']) + 1j * np.array(m['Gb_im'])
    a0 = m['gwin'][0]
    G = {}
    W = Gm.shape[0]
    for a in range(W):
        for b in range(a + 1, W):
            G[(a0 + a, a0 + b)] = Gm[a, b]
    D = {}
    for key, row in m['Dpair'].items():
        dd, i0 = map(int, key.split('_'))
        D[(dd, i0)] = np.array([a + 1j * b for a, b in row])      # N conserved: <bb> = 0, already connected
    C = {int(k): np.array(v) for k, v in m['Cnn'].items()}
    return dict(L=L, G=G, D=D, C=C, dens=np.array(m['dens']), tag=tag)


# ---------------------------------------------------------------------------------------------------------------
# estimators
# ---------------------------------------------------------------------------------------------------------------
def onebody_seq(G, L):
    rmax = max(j - i for (i, j) in G)
    rs, g, F = [], [], []
    for r in range(1, rmax + 1):
        vals = []
        for dl in range(-DC, DC + 1):
            s = L - 1 + dl
            if (s - r) % 2:
                continue
            i = (s - r) // 2
            if (i, i + r) in G:
                vals.append(G[(i, i + r)])
        if len(vals) < 3:
            break
        rs.append(r)
        g.append(np.mean(vals))
        x = (L + 1 - r) / 2.0
        F.append(Fform(x, x + r, L))
    return np.array(rs), np.array(g), np.array(F)


def EB(dat, ranges, npk=2, rel=0.3):
    L = dat['L']
    rs, g, F = onebody_seq(dat['G'], L)
    rmax = rs[-1]
    pk = peaks(rs, g, npk)
    pk = [p for p in pk if p[1] >= rel * pk[0][1]]
    res = []
    for k, amp in pk:
        hh = sm(g * np.exp(-1j * k * rs))
        s = (rs >= 6) & (rs <= rmax - NSM)
        k = k + np.polyfit(rs[s], np.unwrap(np.angle(hh[s])), 1)[0]
        h = np.abs(sm(g * np.exp(-1j * k * rs)))
        ets = []
        for lo, hi in ranges:
            s = (rs >= lo * L) & (rs <= hi * L) & (rs <= rmax - NSM)
            if s.sum() < 4 or hi * L > rmax - NSM + 1:
                ets.append((np.nan, np.nan))
            else:
                ets.append(linfit(np.log(F[s]), np.log(h[s])))
        res.append(dict(k=float(np.angle(np.exp(1j * k))), amp=float(amp), eta=[e for e, _ in ets],
                        rms=[x for _, x in ets]))
    return res, int(rmax)


def EP(dat, d, ranges):
    L = dat['L']
    rows = sorted((i0, v) for (dd, i0), v in dat['D'].items() if dd == d and len(v) > 8)
    # momentum from the mean row
    n0 = min(len(v) for _, v in rows)
    rs = np.arange(n0) + d + 1.0
    mean = np.mean([v[:n0] for _, v in rows], axis=0)
    s0 = rs >= 3
    Q = peaks(rs[s0], mean[s0], 1, lo=-np.pi, hi=np.pi, excl=0.3)[0][0]
    hh = sm(mean * np.exp(-1j * Q * rs))
    s = (rs >= 6) & (rs <= rs.max() - NSM)
    Q = Q + np.polyfit(rs[s], np.unwrap(np.angle(hh[s])), 1)[0]
    data = []
    for i0, v in rows:
        r = np.arange(len(v)) + d + 1.0
        h = np.abs(sm(v * np.exp(-1j * Q * r)))
        x1 = i0 + 1 + d / 2.0
        Fr = Fform(x1, x1 + r, L)
        data.append((r, Fr, h))
    out = []
    rmx = min(r.max() for r, _, _ in data) - NSM
    for lo, hi in ranges:
        if hi * L > rmx + 1:
            out.append((np.nan, np.nan))
            continue
        X, Y = [], []
        for r, Fr, h in data:
            s = (r >= lo * L) & (r <= hi * L) & (h > 0)
            X += list(np.log(Fr[s]))
            Y += list(np.log(h[s]))
        out.append(linfit(np.array(X), np.array(Y)) if len(X) > 4 else (np.nan, np.nan))
    sign = float(np.sign(np.mean(np.real(mean * np.exp(-1j * Q * rs))[(rs >= 6) & (rs <= rs.max() / 2)])))
    return dict(Q=float(np.angle(np.exp(1j * Q))), eta=[e for e, _ in out], rms=[x for _, x in out], sign=sign)


def ES(dat, settings):
    L = dat['L']
    vs = list(dat['C'].values())
    n0 = min(len(v) for v in vs)
    Cs = np.mean([v[:n0] for v in vs], axis=0)
    r = np.arange(n0) + 1.0
    out = []
    for rho, lo, hi in settings:
        R = rho * L
        if R > r.max():
            out.append(np.nan)
            continue
        s = r <= R
        w = np.cos(np.pi * r[s] / (2 * R)) ** 2
        qs = np.linspace(lo * np.pi / R, hi * np.pi / R, 16)
        y = np.array([2 * np.sum(w * Cs[s] * np.cos(q * r[s])) for q in qs])
        A = np.stack([np.ones_like(qs), qs, qs ** 2], 1)
        c, a, b = np.linalg.lstsq(A, y, rcond=None)[0]
        out.append(float(a * np.pi))
    return out


# ---------------------------------------------------------------------------------------------------------------
# trial and frozen ranges (fractions of L)
# ---------------------------------------------------------------------------------------------------------------
TRIAL_B = [(0.04, 0.12), (0.06, 0.20), (0.10, 0.30), (0.15, 0.40), (0.20, 0.60), (0.30, 0.80)]
TRIAL_P = [(0.04, 0.12), (0.06, 0.20), (0.10, 0.30), (0.15, 0.40)]
TRIAL_S = [(0.25, 3, 8), (0.25, 4, 10), (0.35, 3, 8), (0.35, 4, 10), (0.45, 3, 8), (0.45, 4, 10)]
FROZEN = None                # filled from frozen.json once the controls are done
SYS = None


def load_frozen():
    global FROZEN
    fz = json.load(open('frozen.json'))
    FROZEN = fz
    return fz


def summarize(vals):
    a = np.array([v for v in vals if np.isfinite(v)])
    if len(a) == 0:
        return np.nan, np.nan
    return float(a.mean()), float((a.max() - a.min()) / 2)


def frozen_estimates(dat, fz):
    o = {}
    b, rmax = EB(dat, fz['B'])
    o['eta_b'], sp = summarize([e for p in b for e in p['eta']])
    o['deta_b'] = float(np.hypot(sp, fz['sys_b']))
    o['b'] = b
    o['k'] = [p['k'] for p in b]
    for d in (3, 2):
        p = EP(dat, d, fz['P'])
        p['eta_mean'], spp = summarize(p['eta'])
        p['deta'] = float(np.hypot(spp, fz['sys_p']))
        o[f'pair{d}'] = p
    o['eta_p'], o['deta_p'] = o['pair3']['eta_mean'], o['pair3']['deta']
    s = ES(dat, fz['S'])
    o['Kc'], spk = summarize(s)
    o['dKc'] = float(np.hypot(spk, fz['sys_s'] * o['Kc']))
    o['Kc_list'] = s
    for lab, inv, dinv in (('pair', o['eta_p'], o['deta_p']), ('sf', 1 / o['Kc'], o['dKc'] / o['Kc'] ** 2)):
        x = 4 * o['eta_b'] - inv
        o[f'Ks_{lab}'] = float(1 / x)
        o[f'dKs_{lab}'] = float(np.hypot(4 * o['deta_b'], dinv) / x ** 2)
    o['nmid'] = float(np.mean(dat['dens'][dat['L'] // 4: 3 * dat['L'] // 4]))
    return o


def bound(K1, K2, dK1, dK2, dl=np.log(2)):
    """drift.py: |y1| from -dK over dl (central and with dK - sigma), y_par = 2(K2 - 1) (central, lowered)."""
    dK = K2 - K1
    sig = float(np.hypot(dK1, dK2))
    Kb = 0.5 * (K1 + K2)
    y2c = max(0.0, -2 * dK / (Kb ** 2 * dl))
    y2u = max(0.0, -2 * (dK - sig) / (Kb ** 2 * dl))
    return dict(dK=dK, sig=sig, y1c=float(np.sqrt(y2c)), y1u=float(np.sqrt(y2u)), ypc=2 * (K2 - 1),
                ypl=2 * (K2 - dK2 - 1))


KEYS = ['p153c', 'sym40_n316', 'rat_n316', 'sym40_n18', 'p153_n316', 'sym80_n316']


def main():
    mode = sys.argv[1]
    args = [a for a in sys.argv[2:] if not a.startswith('--')]
    mdir = args[0] if args else 'typhon_out/out'
    if mode == 'controls':
        files = ['xx2_L64_P0_rand1', 'xx2_L128_P0_rand1', 'xx2_L256_P0_rand1',
                 'xx1_L64_P1_rand1', 'xx1_L128_P1_rand1', 'xx1_L256_P0_rand1']
        res = {}
        for t in files:
            fn = os.path.join(mdir, t + '_meas.json')
            if not os.path.exists(fn):
                print('missing', t)
                continue
            dat = from_qubit(fn)
            b, rmax = EB(dat, TRIAL_B)
            p3 = EP(dat, 3, TRIAL_P)
            p2 = EP(dat, 2, TRIAL_P)
            s = ES(dat, TRIAL_S)
            res[t] = dict(b=b, rmax=rmax, p3=p3, p2=p2, s=s, check=dat['check'])
            print(f"{t}: rmax={rmax} check dE={dat['check']['dE']:.2e} dS6={dat['check']['dSbar6']:.2e} "
                  f"dGb={dat['check']['dGb_window']:.2e}")
            for p in b:
                print(f"   EB k={p['k']:+.4f} eta=" + ' '.join(f'{e:.4f}' for e in p['eta']) +
                      '  rms=' + ' '.join(f'{x:.3f}' for x in p['rms']))
            for lab, p in (('d3', p3), ('d2', p2)):
                print(f"   EP {lab} Q={p['Q']:+.4f} sign={p['sign']:+.0f} eta=" + ' '.join(f'{e:.4f}' for e in p['eta']) +
                      '  rms=' + ' '.join(f'{x:.3f}' for x in p['rms']))
            print('   ES pi*slope=' + ' '.join(f'{x:.4f}' for x in s))
        if '--eff' in sys.argv:
            for t in ['c2_L64', 'c2_L128', 'c2_L256']:
                dat = from_eff(t)
                b, rmax = EB(dat, TRIAL_B)
                p3 = EP(dat, 3, TRIAL_P)
                s = ES(dat, TRIAL_S)
                print(f'eff {t}: rmax={rmax} EB ' + ' | '.join(' '.join(f'{e:.4f}' for e in p['eta']) for p in b) +
                      '  EP3 ' + ' '.join(f'{e:.4f}' for e in p3['eta']) + '  ES ' + ' '.join(f'{x:.4f}' for x in s))
        json.dump(res, open('controls_trial.json', 'w'), indent=1, default=float)
        return

    fz = load_frozen()
    if mode == 'states':
        for fn in args:
            dat = from_qubit(fn)
            o = frozen_estimates(dat, fz)
            print(f"{dat['tag']}: E={dat['E']:.9f} nmid={o['nmid']:.4f} eta_b={o['eta_b']:.4f}({o['deta_b']:.4f}) "
                  f"pair3={o['eta_p']:.4f}({o['deta_p']:.4f}) pair2={o['pair2']['eta_mean']:.4f} Kc={o['Kc']:.3f}({o['dKc']:.3f}) "
                  f"Ks_pair={o['Ks_pair']:.3f}({o['dKs_pair']:.3f}) Ks_sf={o['Ks_sf']:.3f}({o['dKs_sf']:.3f}) "
                  f"check dE={dat['check']['dE']:.1e} dS6={dat['check']['dSbar6']:.1e} dGb={dat['check']['dGb_window']:.1e}")
        return

    # mode all
    allres = {}
    print('== controls (frozen)')
    for t in ['xx2_L64_P0_rand1', 'xx2_L128_P0_rand1', 'xx2_L256_P0_rand1', 'xx1_L64_P1_rand1', 'xx1_L128_P1_rand1',
              'xx1_L256_P0_rand1']:
        fn = os.path.join(mdir, t + '_meas.json')
        dat = from_qubit(fn)
        o = frozen_estimates(dat, fz)
        allres[t] = o
        print(f"{t:18s} eta_b={o['eta_b']:.4f}({o['deta_b']:.4f}) pair3={o['eta_p']:.4f}({o['deta_p']:.4f}) "
              f"Kc={o['Kc']:.4f}({o['dKc']:.4f})")
    print('== candidates')
    cand = {}
    for key in KEYS:
        for L in (64, 128, 256):
            sec = {}
            for P in (0, 1):
                fns = sorted(glob.glob(os.path.join(mdir, f'{key}_L{L}_P{P}_*_meas.json')))
                if not fns:
                    print('missing', key, L, P)
                    continue
                dat = from_qubit(fns[0])
                o = frozen_estimates(dat, fz)
                o['E'] = dat['E']
                o['check'] = dat['check']
                o['file'] = os.path.basename(fns[0])
                sec[P] = o
            if not sec:
                continue
            gb = min(sec, key=lambda p: sec[p]['E'])
            cand[f'{key}_L{L}'] = dict(sec={str(p): sec[p] for p in sec}, ground=gb)
            for P in sorted(sec):
                o = sec[P]
                c = o['check']
                print(f"{key:11s} L={L:3d} P={P}{'*' if P == gb else ' '} n={o['nmid']:.4f} k={','.join(f'{k:+.3f}' for k in o['k'])} "
                      f"eta_b={o['eta_b']:.4f}({o['deta_b']:.4f}) pair3={o['eta_p']:.4f}({o['deta_p']:.4f}) Q3={o['pair3']['Q']:+.3f} "
                      f"pair2={o['pair2']['eta_mean']:.4f} Kc={o['Kc']:.3f}({o['dKc']:.3f}) "
                      f"Ks_pair={o['Ks_pair']:.3f}({o['dKs_pair']:.3f}) Ks_sf={o['Ks_sf']:.3f}({o['dKs_sf']:.3f}) "
                      f"| dE={c['dE']:.1e} dS6={c['dSbar6']:.1e} dGb={c['dGb_window']:.1e}")
    allres['cand'] = cand
    print('== drift and KT bound (Delta l = ln 2 per doubling)')
    drift = {}
    for key in KEYS:
        for lab in ('pair', 'sf'):
            for br in ('ground', '0', '1'):
                try:
                    o = [cand[f'{key}_L{L}']['sec'][str(cand[f'{key}_L{L}']['ground']) if br == 'ground' else br]
                         for L in (64, 128, 256)]
                except KeyError:
                    continue
                K = [x[f'Ks_{lab}'] for x in o]
                dK = [x[f'dKs_{lab}'] for x in o]
                b1 = bound(K[0], K[1], dK[0], dK[1])
                b2 = bound(K[1], K[2], dK[1], dK[2])
                acc = b2['dK'] - b1['dK']
                sacc = float(np.hypot(b1['sig'], b2['sig']))
                drift[f'{key}_{lab}_{br}'] = dict(K=K, dK=dK, d1=b1['dK'], d2=b2['dK'], s1=b1['sig'], s2=b2['sig'],
                                                   y1u=b2['y1u'], y1c=b2['y1c'], ypl=b2['ypl'], ypc=b2['ypc'],
                                                   acc=acc, sacc=sacc)
                print(f"{key:11s} {lab:4s} {br:6s} Ks 64/128/256 = {K[0]:.3f}/{K[1]:.3f}/{K[2]:.3f} "
                      f"(err {dK[0]:.3f}/{dK[1]:.3f}/{dK[2]:.3f})  d(64->128)={b1['dK']:+.3f}({b1['sig']:.3f}) "
                      f"d(128->256)={b2['dK']:+.3f}({b2['sig']:.3f})  |y1| {b2['y1c']:.3f}/{b2['y1u']:.3f}  "
                      f"y_par {b2['ypc']:.3f}/{b2['ypl']:.3f}  d2-d1={acc:+.3f}({sacc:.3f})")
    allres['drift'] = drift
    if '--eff' in sys.argv:
        print('== effective model (C176), same frozen estimators')
        eff = {}
        for tag in ['c2_L64', 'c2_L128', 'c2_L256', 'sym316_L64', 'sym316_L128', 'sym316_L256',
                    'sym18_L64', 'sym18_L128', 'sym18_L256']:
            dat = from_eff(tag)
            o = frozen_estimates(dat, fz)
            eff[tag] = o
            print(f"{tag:12s} n={o['nmid']:.4f} eta_b={o['eta_b']:.4f}({o['deta_b']:.4f}) [ranges " +
                  ' '.join(f'{e:.4f}' for p in o['b'] for e in p['eta']) + f"] pair3={o['eta_p']:.4f}({o['deta_p']:.4f}) "
                  f"pair2={o['pair2']['eta_mean']:.4f} Kc={o['Kc']:.3f}({o['dKc']:.3f}) "
                  f"Ks_pair={o['Ks_pair']:.3f}({o['dKs_pair']:.3f}) Ks_sf={o['Ks_sf']:.3f}({o['dKs_sf']:.3f})")
        for t in ('sym316', 'sym18'):
            for lab in ('pair', 'sf'):
                K = [eff[f'{t}_L{L}'][f'Ks_{lab}'] for L in (64, 128, 256)]
                dK = [eff[f'{t}_L{L}'][f'dKs_{lab}'] for L in (64, 128, 256)]
                b1 = bound(K[0], K[1], dK[0], dK[1])
                b2 = bound(K[1], K[2], dK[1], dK[2])
                print(f"eff {t:7s} {lab:4s} Ks {K[0]:.3f}/{K[1]:.3f}/{K[2]:.3f} d1={b1['dK']:+.3f}({b1['sig']:.3f}) "
                      f"d2={b2['dK']:+.3f}({b2['sig']:.3f}) |y1| {b2['y1c']:.3f}/{b2['y1u']:.3f} y_par {b2['ypc']:.3f}/{b2['ypl']:.3f}")
        allres['eff'] = eff
    json.dump(allres, open('analysis.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
