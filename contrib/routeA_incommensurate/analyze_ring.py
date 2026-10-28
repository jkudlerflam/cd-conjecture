#!/usr/bin/env python3
"""Compare the full-chain ring spectra (data/bigring_T_LAM.json, whole logical band per parity
sector) with (a) the second-order logical chain and (b) the cluster-ED logical chain
H = sum h2(lam) + sum h3(lam) on the same ring (data/heff_T.npz).  Differences in units of J.
On a ring of L = 3 cells there are also loop processes around the ring at order lam^3 that no open
cluster contains, so (b) is not expected to be exact there.
usage: analyze_ring.py INC 0.05 0.1 0.2
Written by Claude (Opus 5.5) for contrib/routeA_incommensurate, 2026-10-28.
"""
import sys, json
import numpy as np
from scipy import linalg
import cell as CE


def ring_dense(terms, N, L):
    """terms: list of (r, op) with op on r consecutive sites; ring of L sites."""
    dim = N ** L
    H = np.zeros((dim, dim), complex)
    I = np.eye(dim).reshape([N] * L + [dim])
    for r, op in terms:
        T = op.reshape([N] * (2 * r))
        for j in range(L):
            sites = [(j + s) % L for s in range(r)]
            t = np.tensordot(T, I, (list(range(r, 2 * r)), sites))
            H += np.moveaxis(t, list(range(r)), sites).reshape(dim, dim)
    return 0.5 * (H + H.conj().T)


def main():
    T = sys.argv[1]
    z = np.load(f'data/heff_{T}.npz')
    J0 = float(z['J0'])
    out = []
    for lam in map(float, sys.argv[2:]):
        d = json.load(open(f'data/bigring_{T}_{lam}.json'))
        J = J0 * lam ** 2
        for r in d['results']:
            L = r['L']
            N = 4
            full = np.sort(np.concatenate([s['E'] for s in r['sectors']]))
            sec2 = np.sort(np.concatenate([s['E_log2'] for s in r['sectors']]))
            terms = [(2, z[f'h2_{lam}'])]
            if f'h3_{lam}' in z.files:
                terms.append((3, z[f'h3_{lam}']))
            Hc = ring_dense(terms, N, L)
            clus = np.sort(linalg.eigvalsh(Hc))
            m = len(full)
            rec = dict(lam=lam, L=L, levels=m,
                       full_vs_2nd_abs=float(np.max(np.abs(full - sec2[:m])) / J),
                       full_vs_2nd_rel=float(np.max(np.abs((full - full[0]) - (sec2[:m] - sec2[0]))) / J),
                       full_vs_cluster_abs=float(np.max(np.abs(full - clus[:m])) / J),
                       full_vs_cluster_rel=float(np.max(np.abs((full - full[0]) - (clus[:m] - clus[0]))) / J),
                       bandwidth=float((full[-1] - full[0]) / J))
            out.append(rec)
            print(json.dumps(rec))
    json.dump(out, open(f'data/ring_compare_{T}.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
