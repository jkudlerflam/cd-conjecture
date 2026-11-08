"""C155 stage E: directed scan of the XZ neighbourhood (Claude Opus agent, 2026-11-07).

Base point h = X (x) Z (J_xz = 1). Directions in J-space (12 coordinates, search.py): YX, ZY, -YX, -ZY, YX + ZY,
YX - ZY, 12 random directions inside C144's tuned hyperplane T (dJ_zx = (2/pi) dJ_yy) and 8 random directions
in the full normal space. Distances delta in {0.1, 0.2, 0.3, 0.45}. For each point: labels, ring ED at
(12, 16) and (16, 20), open-chain DMRG at L = 64 and 128 (chi up to 128), and L = 256 when any L = 128 estimate
exceeds 1.05. Also the Kitaev chain in a [110] field (C156): XZ + B (x + z)/sqrt 2, B = 0.05 to 1.0, and 4 random
kicks of size 0.1 around B = 0.2, 0.3, 0.4. Output xz/pt_<i>.json (resumable).
"""
import json, os, sys, time
import numpy as np
from multiprocessing import Pool
from search import W0, E_YX, E_ZY, B_XZ, B_XZT, F_BIS, w_to_v


def directions():
    rng = np.random.default_rng(1557)
    D = [('YX', E_YX), ('ZY', E_ZY), ('-YX', -E_YX), ('-ZY', -E_ZY),
         ('YX+ZY', (E_YX + E_ZY) / np.sqrt(2)), ('YX-ZY', (E_YX - E_ZY) / np.sqrt(2))]
    for i in range(12):
        d = B_XZT @ rng.normal(size=B_XZT.shape[1]); D.append((f'T{i}', d / np.linalg.norm(d)))
    for i in range(8):
        d = B_XZ @ rng.normal(size=B_XZ.shape[1]); D.append((f'N{i}', d / np.linalg.norm(d)))
    return D


def points():
    P = []
    for name, d in directions():
        for delta in (0.1, 0.2, 0.3, 0.45):
            P.append((name, delta, (W0 + delta * d).tolist()))
    # Kitaev chain in a [110] field (C156, arXiv:2209.06221): XZ + B (x + z)/sqrt 2 on every site
    rng = np.random.default_rng(1558)
    for B in (0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0):
        P.append(('KitaevB', B, (W0 + B * F_BIS).tolist()))
    for B in (0.2, 0.3, 0.4):
        for j in range(4):
            k = B_XZ @ rng.normal(size=B_XZ.shape[1]); k *= 0.1 / np.linalg.norm(k)
            P.append((f'KitaevB+kick{j}', B, (W0 + B * F_BIS + k).tolist()))
    return P


def work(item):
    i, (name, delta, w), chi = item
    fn = os.path.join(os.environ.get('XZ_OUT', 'xz'), f'pt_{i}.json')
    if os.path.exists(fn):
        return fn
    import est, dmrg_c
    from crit import labels
    w = np.array(w); v = w_to_v(w)
    out = dict(i=i, dir=name, delta=delta, w=w.tolist(), v=v.tolist(), labels=labels(v))
    keys = ('c_S', 'c_fit', 'c_cas', 'c_obj', 'z1', 'z2', 'crit', 'ndeg', 'v', 'm0', 'gap1', 'gap2')
    for Ls in ((12, 16), (16, 20)):
        try:
            e = est.estimate(est.bond_v(v), Ls=Ls, nev=3)
            out[f'ed_{Ls[0]}_{Ls[1]}'] = {k: e[k] for k in keys}
        except Exception as ex:
            out[f'ed_{Ls[0]}_{Ls[1]}'] = dict(error=str(ex))
    try:
        o = dmrg_c.run(v, Ls=(64, 128), chis=(16, 32, 64, chi, chi, chi))
        big = max(max(x['c_prof'], x['c_sm']) for x in o['L'].values())
        if big > 1.05:
            o2 = dmrg_c.run(v, Ls=(256,), chis=(16, 32, 64, chi, chi, chi))
            o['L'][256] = o2['L'][256]
            mids = [o['L'][L]['Smid'] for L in (128, 256)]
            o['c_two_128_256'] = float(6 * (mids[1] - mids[0]) / np.log(2))
        out['dmrg'] = o
    except Exception as ex:
        out['dmrg'] = dict(error=str(ex))
    json.dump(out, open(fn, 'w'))
    return fn


if __name__ == '__main__':
    procs = int(sys.argv[1]) if len(sys.argv) > 1 else 48
    chi = int(sys.argv[2]) if len(sys.argv) > 2 else 128
    os.makedirs('xz', exist_ok=True)
    P = points()
    # largest delta first (fastest DMRG convergence), so partial runs already cover every direction
    order = sorted(range(len(P)), key=lambda i: -P[i][1])
    with Pool(procs) as p:
        for fn in p.imap_unordered(work, [(i, P[i], chi) for i in order]):
            print('done', fn, flush=True)
