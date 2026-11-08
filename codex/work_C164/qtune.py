"""C164: continuation of a two-valley band to its onset (Claude Opus agent, 2026-11-08).

Given a bond with two isolated valleys A, B (qopt.evaluate at L1 = 12, L2 = 14), take minimal-norm Newton steps in
the 9 bond coordinates that move (E_A, E_B) toward (g_t, g_t), lowering the common gap g_t geometrically while keeping
E_A = E_B. Valleys are tracked by nearest momentum. Each accepted step must keep the qopt quality conditions
(uniform background, smooth isolated band, barrier, incommensurate valleys). Stops at the smallest gap reached.
The derivatives use central differences with step 1e-3 (in the unnormalized v; the bond is normalized after).
usage: qtune.py in.jsonl out.jsonl [Othr] [index range a:b, or stride i/n]
"""
import sys, json, time
import numpy as np
import qopt

TWO_PI = 2 * np.pi


def track(v, kA0, kB0, L1=12, L2=14):
    O, inf = qopt.evaluate(v, L1, L2)
    if 'kB' not in inf:
        return None, O, inf
    kA, kB = inf['kA'], inf['kB']
    if qopt.cdist(kA, kB0) + qopt.cdist(kB, kA0) < qopt.cdist(kA, kA0) + qopt.cdist(kB, kB0):
        inf = dict(inf); inf['kA'], inf['kB'] = kB, kA; inf['EA'], inf['EB'] = inf['EB'], inf['EA']
    moved = qopt.cdist(inf['kA'], kA0) + qopt.cdist(inf['kB'], kB0)
    return (inf['EA'], inf['EB']), O, inf if moved < 0.6 else dict(inf, lost=True)


def quality_ok(inf, bmin=0.01, tmax=0.05, margin=0.3):
    if 'kB' not in inf or inf.get('lost'):
        return False
    if inf['barrier'] < bmin or max(inf['tail']) > tmax or not inf['iso']:
        return False
    for k in (inf['kA'], inf['kB']):
        if qopt.bad_momentum(k, margin):
            return False
    if qopt.cdist(inf['kA'] + inf['kB'], 0.0) < margin:
        return False
    return inf['gap14'] / max(inf['gap12'], 1e-12) >= 0.8


def tune(v, steps=14, shrink=0.6, h=1e-3, log=None):
    v = np.asarray(v, float); v /= np.linalg.norm(v)
    O, inf = qopt.evaluate(v)
    kA0, kB0 = inf['kA'], inf['kB']
    E, O, inf = track(v, kA0, kB0)
    hist = [dict(v=v.tolist(), O=O, info=inf)]
    gt = 0.5 * (E[0] + E[1])
    for it in range(steps):
        gt *= shrink
        Jm = np.zeros((2, 9)); ok = True
        for i in range(9):
            dv = np.zeros(9); dv[i] = h
            Ep, _, ip = track(v + dv, kA0, kB0); Em, _, im = track(v - dv, kA0, kB0)
            if Ep is None or Em is None or ip.get('lost') or im.get('lost'):
                ok = False; break
            Jm[:, i] = (np.array(Ep) - np.array(Em)) / (2 * h)
        if not ok:
            break
        r = np.array([gt, gt]) - np.array(E)
        dv = np.linalg.pinv(Jm) @ r
        nrm = np.linalg.norm(dv)
        if nrm > 0.05:
            dv *= 0.05 / nrm
        acc = False
        for lam in (1.0, 0.5, 0.25):
            vn = v + lam * dv; vn /= np.linalg.norm(vn)
            En, On, inn = track(vn, kA0, kB0)
            if En is None or inn.get('lost'):
                continue
            # quality must stay acceptable (penalties of qopt other than gap/degeneracy below 0.15)
            if quality_ok(inn) and min(En) > -0.02:
                v, E, O, inf, acc = vn, En, On, inn, True
                kA0, kB0 = inf['kA'], inf['kB']
                break
        if not acc:
            break
        hist.append(dict(v=v.tolist(), O=O, info=inf))
        if log:
            log(f"  it {it} gap target {gt:.4f}: E = ({E[0]:.4f}, {E[1]:.4f}) kA={inf['kA']:.3f} kB={inf['kB']:.3f} O={O:.3f}")
    return hist


if __name__ == '__main__':
    fin, fout = sys.argv[1], sys.argv[2]
    thr = float(sys.argv[3]) if len(sys.argv) > 3 else 0.3
    R = [json.loads(l) for l in open(fin)]
    R = [r for r in R if r['O'] < thr and 'kB' in r['info']]
    R.sort(key=lambda r: r['O'])
    if len(sys.argv) > 4:
        if '/' in sys.argv[4]:
            i, n = map(int, sys.argv[4].split('/')); R = R[i::n]
        else:
            a, b = map(int, sys.argv[4].split(':')); R = R[a:b]
    with open(fout, 'a') as fo:
        for r in R:
            t0 = time.time()
            hist = tune(r['v'], log=lambda s: print(s, flush=True))
            fo.write(json.dumps(dict(src=dict(task=r['task'], restart=r['restart']), hist=hist, t=time.time() - t0)) + '\n')
            fo.flush()
            last = hist[-1]['info']
            print('done', r['task'], r['restart'], 'final', {k: last.get(k) for k in ('kA', 'kB', 'EA', 'EB', 'bw')}, f'{time.time() - t0:.0f}s', flush=True)
