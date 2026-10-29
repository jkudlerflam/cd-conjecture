#!/usr/bin/env python3
"""Charge- and momentum-resolved ring spectra through the ED engine of
../counterexample_search/ced.py (Sector class), for the bonds that conserve additive
site charges.  Returns E_0 over all requested charge targets, the ground momentum k_0,
the lowest level at k_0 +- 2 pi/L in the ground sector (the "current" state, whose gap
times L/(2 pi) estimates the velocity when its scaling dimension is 1), the overall gap
and the half-ring entropy.  numpy and scipy only.
"""
import os, sys, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(HERE, '..', 'counterexample_search'))
from ced import Sector  # noqa: E402


def all_targets(charges, L):
    """Every reachable total charge of L sites (charges: d x nq integer array)."""
    charges = np.asarray(charges, dtype=np.int64)
    if charges.ndim == 1:
        charges = charges[:, None]
    tot = {tuple(np.zeros(charges.shape[1], dtype=np.int64))}
    for _ in range(L):
        tot = {tuple(np.array(t) + c) for t in tot for c in charges}
    return sorted(tot)


def symring(h, d, L, charges, targets=None, ks=None, nlev=3, tol=1e-10, entropy=True,
            target_filter=None):
    charges = np.asarray(charges, dtype=np.int64)
    if charges.ndim == 1:
        charges = charges[:, None]
    if targets is None:
        targets = all_targets(charges, L)
        if target_filter is not None:
            targets = [t for t in targets if target_filter(t, L)]
    out = {'L': L, 'spec': {}}
    best = (np.inf, None, None, None, None)
    for t in targets:
        sec = Sector(d, L, [h], charges=charges, target=t, cache=False)
        if len(sec.basis) == 0:
            continue
        kk = range(L) if ks is None else [k % L for k in ks]
        for k in kk:
            e, v, res = sec.solve([1.0], k, nlev, tol=tol)
            if len(e) == 0:
                continue
            out['spec'][(tuple(int(x) for x in t), int(k))] = [float(x) for x in e]
            if e[0] < best[0] - 1e-12:
                best = (float(e[0]), int(k), tuple(int(x) for x in t), v[:, 0], sec)
    out['E0'], out['k0'], out['target'] = best[0], best[1], best[2]
    # gaps
    allE = sorted(x for v in out['spec'].values() for x in v)
    out['gap'] = allE[1] - allE[0] if len(allE) > 1 else np.nan
    k0 = out['k0']
    cur = [out['spec'][(out['target'], k)][0] for k in ((k0 + 1) % L, (k0 - 1) % L)
           if (out['target'], k) in out['spec']]
    out['gap_current'] = (min(cur) - out['E0']) if cur else np.nan
    out['v_est'] = out['gap_current'] * L / (2 * np.pi)
    if entropy:
        sec = best[4]
        psi = sec.expand(best[1], best[3])
        out['S_half'] = float(sec.entropies(psi, cuts=[L // 2])[0])
    return out
