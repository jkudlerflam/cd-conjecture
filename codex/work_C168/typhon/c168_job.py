"""C168 Typhon job (Claude Opus 5.5, 2026-11-08): the sampling bound on open c > 1 regions of B_2.

Stage P  proxy calibration of the ring screen on qubit chains with known c (range 3, outside B_2):
         Jordan-Wigner fermions with complex NN and NNN hopping (one or two Fermi seas), a NN density
         interaction V Z_j Z_{j+1}, and a random generic NN perturbation W of norm epsw (no on-site symmetry).
         Families: 'free' (V = 0, epsw = 0.1), 'int' (V uniform in [-0.6, 0.6], epsw = 0.25).
         Rings (12, 16) for all samples, rings (16, 20) for those with c_S(12, 16) > 1.5.
Stage S1 every C155 stage-A uniform point with c_S(12, 16) > 1.5 (and the 28 C155 gapless flags):
         rings (16, 20).
Stage S2 every stage-A point with c_S(16, 20) > 1.5 or c_obj(16, 20) > 1.1, and every C155 gapless flag:
         open-chain DMRG at L = 32, 64, 128.
Rules fixed before any output: screen F1 = c_S(12,16) > 1.5, F2 = c_S(16,20) > 1.5; a c > 1 candidate
needs c_two(64 -> 128) and c_prof(128) above 1.2 with truncation below 1e-7 (then FES, C163 rule).
Usage: python3 c168_job.py STAGE NPROC
"""
import json, os, sys, time, glob
import numpy as np
from multiprocessing import Pool

KEYS = ('c_S', 'c_fit', 'c_cas', 'c_obj', 'crit', 'ndeg', 'm0', 'Shalf')
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], complex); Y = np.array([[0, -1j], [1j, 0]], complex); Z = np.diag([1., -1.]).astype(complex)
SP = np.array([[0, 1], [0, 0]], complex); SM = SP.T.copy()
P4 = [I2, X, Y, Z]


def k3(a, b, c):
    return np.kron(np.kron(a, b), c)


def seas(t1, a1, t2, a2, mu, n=20000):
    k = np.linspace(-np.pi, np.pi, n, endpoint=False)
    e = -2 * t1 * np.cos(k - a1) - 2 * t2 * np.cos(2 * k - a2) - mu
    occ = e < 0
    if occ.all() or not occ.any():
        return 0, []
    d = np.diff(np.concatenate([occ, occ[:1]]).astype(int))
    starts = np.nonzero(d == 1)[0]; ends = np.nonzero(d == -1)[0]
    widths = []
    for s in starts:
        ea = ends[ends > s]
        en = ea[0] if len(ea) else ends[0] + n
        widths.append(float((en - s) * 2 * np.pi / n))
    return len(starts), widths


def proxy_bond(t1, a1, t2, a2, mu, V, W):
    h = -t1 * np.exp(1j * a1) * k3(SP, SM, I2) - t2 * np.exp(1j * a2) * k3(SP, -Z, SM)
    h = h + h.conj().T
    h = h - (mu / 2) * k3(Z, I2, I2) + V * k3(Z, Z, I2) + np.kron(W, I2)
    return h


def proxy_task(args):
    fam, i = args
    import est_c155 as est
    rng = np.random.default_rng((1 if fam == 'free' else 2) * 1000003 + i)
    t1 = 1.0; t2 = rng.uniform(0.5, 1.5); a1 = rng.uniform(0, 2 * np.pi); a2 = rng.uniform(0, 2 * np.pi)
    mu = rng.uniform(-2 * (t1 + t2), 2 * (t1 + t2))
    V = 0.0 if fam == 'free' else rng.uniform(-0.6, 0.6)
    epsw = 0.1 if fam == 'free' else 0.25
    nr, w = seas(t1, a1, t2, a2, mu)
    c = rng.normal(size=(4, 4)); c[0, 0] = 0
    W = sum(c[a, b] * np.kron(P4[a], P4[b]) for a in range(4) for b in range(4)) * epsw / np.linalg.norm(c)
    h = proxy_bond(t1, a1, t2, a2, mu, V, W)
    out = dict(fam=fam, i=i, t2=t2, a1=a1, a2=a2, mu=mu, V=V, epsw=epsw, nseas=nr, widths=w)
    try:
        e = est.estimate(h, Ls=(12, 16), r=3, nev=2); out['e1216'] = {k: e[k] for k in KEYS}
        if e['c_S'] > 1.5:
            e = est.estimate(h, Ls=(16, 20), r=3, nev=2); out['e1620'] = {k: e[k] for k in KEYS}
    except Exception as ex:
        out['error'] = str(ex)
    return out


def s1_task(item):
    import est_c155 as est
    out = dict(item)
    try:
        e = est.estimate(est.bond_v(item['v']), Ls=(16, 20), nev=2); out['e1620'] = {k: e[k] for k in KEYS}
    except Exception as ex:
        out['error'] = str(ex)
    return out


def s2_task(item):
    import dmrg_c_c155 as dc
    out = dict(item)
    try:
        o = dc.run(item['v'], Ls=(32, 64, 128), chis=(16, 32, 64, 128, 128))
        for L in o['L']:
            o['L'][L].pop('S', None)
        o['c_two_64_128'] = None
        out['dmrg'] = o
    except Exception as ex:
        out['error'] = str(ex)
    return out


def run_pool(fn, items, nproc, outname):
    done = set()
    if os.path.exists(outname):
        for l in open(outname):
            try:
                done.add(json.loads(l)['_key'])
            except Exception:
                pass
    todo = [it for it in items if it['_key'] not in done]
    with Pool(nproc) as p, open(outname, 'a') as f:
        for r in p.imap_unordered(fn, todo):
            f.write(json.dumps(r) + '\n'); f.flush()


def ptask(it):
    r = proxy_task((it['fam'], it['i'])); r['_key'] = it['_key']; return r


if __name__ == '__main__':
    stage, nproc = sys.argv[1], int(sys.argv[2])
    os.makedirs('out', exist_ok=True)
    if stage == 'P':
        items = [dict(fam=f, i=i, _key=f'{f}_{i}') for f in ('free', 'int') for i in range(800)]
        run_pool(ptask, items, nproc, 'out/proxy.jsonl')
    elif stage == 'S1':
        pts = json.load(open('stageA_cS15.json'))
        items = [dict(p, _key=str(p['idx'])) for p in pts]
        run_pool(s1_task, items, nproc, 'out/s1.jsonl')
    elif stage == 'S2':
        sel = []
        flags = set(x['idx'] for x in json.load(open('stageA_flags.json')))   # the 28 C155 gapless flags
        for l in open('out/s1.jsonl'):
            r = json.loads(l); e = r.get('e1620')
            if e and (e['c_S'] > 1.5 or e['c_obj'] > 1.1 or r['idx'] in flags):
                sel.append(dict(idx=r['idx'], v=r['v'], e1620=e, _key=r['_key']))
        json.dump(sel, open('out/s2_selection.json', 'w'))
        run_pool(s2_task, sel, nproc, 'out/s2.jsonl')
