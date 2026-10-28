"""Driver for the exact-diagonalization checks reported in README.md.

Writes checks_results.json.  Runtime: a few minutes on one core.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import json
import numpy as np
from checks import run
from model import Chain

PARAMS = dict(eT=1.0, eH=1.0, eD=1.0, wT=0.3, wH=0.3, wD=-0.3, phiT=0.3,
              tauA=0.3, tauB=0.3, kappa=1.0, muc=0.0)
LAMS = [0.08, 0.04, 0.02]
CASES = [(1, 2, 2), (1, 2, 4), (1, 3, 2), (2, 2, 2), (2, 3, 2)]


def conserves_number(n, M, K):
    ch = Chain(n, M, K, PARAMS)
    H0, V = ch.hamiltonians()
    worst = 0
    for A in (H0, V):
        r, c = A.nonzero()
        worst = max(worst, int(np.max(np.abs(ch.number[r] - ch.number[c]))) if len(r) else 0)
    return worst


if __name__ == '__main__':
    results = []
    for (n, M, K) in CASES:
        out, ch, H2, h = run(n, M, K, LAMS, PARAMS, verbose=True)
        out['number_changing_matrix_elements'] = conserves_number(n, M, K) if M == 2 else 'not rerun'
        results.append(out)
        with open('checks_results.json', 'w') as f:
            json.dump(results, f, indent=1, default=str)
