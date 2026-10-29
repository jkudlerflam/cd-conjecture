#!/usr/bin/env python3
"""Compact copies of the raw ED scans (the raw files, with all sector energies and entropy profiles,
stay on Typhon).  Per point: tag, the 9-vector v, and for each L: c_S, ground-state momentum index m0,
S(L/2), and the four lowest scaled excitation energies L (E_n - E_0) over all momenta.
Writes data/scan_compact.jsonl.gz (sphere and staggered scans) and data/patch_compact.jsonl.gz."""
import os, glob, json, gzip
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))


def row(r):
    o = dict(tag=r['tag'], v=[round(x, 8) for x in r['v']])
    for k in r:
        if k.startswith('L'):
            L = int(k[1:])
            e = np.sort(np.concatenate([np.asarray(x) for x in r[k]['E']]))
            o[k] = dict(cS=round(r[k]['cS'], 5), m0=r[k]['m0'], Shalf=round(r[k]['S'][-1], 5),
                        gaps=[round(float(L * (x - e[0])), 5) for x in e[1:5]])
    return o


for out, pat in (('scan_compact.jsonl.gz', lambda f: not os.path.basename(f).startswith('patch')),
                 ('patch_compact.jsonl.gz', lambda f: os.path.basename(f).startswith('patch_points'))):
    files = [f for f in sorted(glob.glob(os.path.join(HERE, 'data', 'scan', '*.jsonl'))) if pat(f)]
    n = 0
    with gzip.open(os.path.join(HERE, 'data', out), 'wt') as fo:
        for f in files:
            for l in open(f):
                if l.strip():
                    fo.write(json.dumps(row(json.loads(l))) + '\n'); n += 1
    print(out, n)
