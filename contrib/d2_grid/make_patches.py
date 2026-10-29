#!/usr/bin/env python3
"""Junction candidates and local patches.

On every sphere a point is typed G (c_S(16) < 0.25), I (0.3 < c_S < 0.75 at L = 12 and 16, Ising-like),
B (0.8 < c_S < 1.3 at both sizes, boson-like) or O (other).  A junction candidate is an I point whose
eight nearest Fibonacci neighbours contain both a B and a G point, i.e. an Ising-like line touching a
c = 1 region and a gapped region.  Around each candidate a 7 x 7 patch in the tangent plane of the
sphere (half-width --w, about 1.4 neighbour spacings) is written to patch_points.json for scan.py.
"""
import json, glob, os, sys, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import strata

ap = argparse.ArgumentParser()
ap.add_argument('--w', type=float, default=0.2)
ap.add_argument('--n', type=int, default=7)
a = ap.parse_args()
NN = {'NS': 600, 'Z2': 400, 'TR': 500, 'RA': 500}
pts, cands = [], []
for f in sorted(glob.glob(os.path.join(HERE, 'data', 'scan', '*_*.jsonl'))):
    k = os.path.basename(f)[:-6]
    if k.split('_')[0] not in NN or len(k.split('_')) != 2:
        continue
    st, seed = k.split('_'); seed = int(seed)
    X = strata.fibonacci(NN[st]); U = strata.sphere_basis(st, seed)
    rows = {int(json.loads(l)['tag'].split('_')[-1]): json.loads(l) for l in open(f)}
    idx = np.array(sorted(rows)); Xi = X[idx]; nb = np.argsort(-(Xi @ Xi.T), 1)[:, :9]
    typ = []
    for i in idx:
        c12, c16 = rows[i]['L12']['cS'], rows[i]['L16']['cS']
        typ.append('G' if c16 < 0.25 else 'I' if (0.3 < c12 < 0.75 and 0.3 < c16 < 0.75)
                   else 'B' if (0.8 < c12 < 1.3 and 0.8 < c16 < 1.3) else 'O')
    typ = np.array(typ)
    for k2, i in enumerate(idx):
        s = set(typ[nb[k2]])
        if not (typ[k2] == 'I' and 'B' in s and 'G' in s):
            continue
        x0 = X[i]
        e1 = np.cross(x0, [0.3, 0.5, 0.8]); e1 /= np.linalg.norm(e1); e2 = np.cross(x0, e1)
        tag = f'{st}_{seed}_{i}'
        cands.append(tag)
        g = np.linspace(-a.w, a.w, a.n)
        for p, s1 in enumerate(g):
            for q, s2 in enumerate(g):
                y = x0 + s1 * e1 + s2 * e2; y /= np.linalg.norm(y)
                pts.append(dict(name=f'{tag}_p{p}{q}', v=(y @ U).tolist()))
json.dump(pts, open(os.path.join(HERE, 'patch_points.json'), 'w'))
json.dump(cands, open(os.path.join(HERE, 'data', 'junction_candidates.json'), 'w'))
print(len(cands), 'candidates,', len(pts), 'patch points')
