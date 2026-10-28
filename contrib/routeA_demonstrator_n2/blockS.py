#!/usr/bin/env python3
"""Block entanglement entropies S(l) of l whole cells in a saved uMPS (second central-charge estimator).

For a critical chain approximated at finite chi, S(l) = (c/3) ln l + s0 for 1 << l << xi.  The estimator fits
S(l) on l in [lmin, lmax] with lmax well below xi, and is calibrated by applying the same window (in units of
xi) to the logical chains.
usage: blockS.py STATE.npz --cell 3 --lmax 12 --out f.json   (cell = sites per cell: 3 for cell3, 1 for d = 4)
Claude (Opus 5.5), 2026-10-28.
"""
import os, sys, json, argparse
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import vumps2 as V


def load(f):
    z = np.load(f)
    n = len([k for k in z.files if k.startswith('AL')])
    return V.UMPS([z[f'AL{k}'] for k in range(n)], [z[f'AR{k}'] for k in range(n)], [z[f'C{k}'] for k in range(n)])


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('state')
    ap.add_argument('--cell', type=int, default=3)
    ap.add_argument('--lmax', type=int, default=12)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    psi = load(a.state)
    xi = -psi.n / np.log(abs(V.transfer_spectrum(psi)[1]))
    S = V.block_entropies(psi, 0, a.lmax * a.cell, every=a.cell)
    res = dict(state=os.path.basename(a.state), xi_cells=float(xi / a.cell), chi=int(psi.C[0].shape[0]),
               S=[(l // a.cell, s) for l, s in S])
    json.dump(res, open(a.out, 'w'), indent=1)
    print(res)
