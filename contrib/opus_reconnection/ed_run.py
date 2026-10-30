"""Interacting ED checks of the exact decomposition, the clustering formula
and the candidate inequality G <= R_L.

Usage: python3 ed_run.py [local|typhon|typhon_qubit|typhon_su3]
  local  : qubit L = 4, 6 (doubled 8, 12), SU(3) L = 3 (doubled 6)
  typhon : qubit L = 8 (doubled 16), SU(3) L = 6 (doubled 12)
Output: ed_<mode>.json in this folder.
"""
import sys, json, time
from pathlib import Path
import numpy as np
from ed_core import blocks, decomposition

ROOT = Path(__file__).resolve().parent

MODELS = [  # name, model, Delta, hz, M (gauge-optimized half-spread)
    ('XX', 'xxz', 0.0, 0.0, 0.5),
    ('XXZ_m05', 'xxz', -0.5, 0.0, 0.5),
    ('XXZ_p05', 'xxz', 0.5, 0.0, 0.5),
    ('Heisenberg', 'xxz', 1.0, 0.0, 0.5),
    ('XXZ_p05_field', 'xxz', 0.5, 0.4, None),  # M computed below
    ('SU3', 'su3', 1.0, 0.0, 1.0),
]


def bond_M(model, delta, hz):
    """Gauge-optimized half-spread by direct minimization over o (2x2 or 3x3
    Hermitian) of the eigenvalue spread of h + o x 1 - 1 x o."""
    from scipy.optimize import minimize
    if model == 'su3':
        return 1.0
    sx = np.array([[0, 1], [1, 0]]) / 2; sy = np.array([[0, -1j], [1j, 0]]) / 2
    sz = np.diag([0.5, -0.5]); I = np.eye(2)
    h = (np.kron(sx, sx) + np.kron(sy, sy) + delta * np.kron(sz, sz)
         + 0.5 * hz * (np.kron(sz, I) + np.kron(I, sz)))
    basis = [sx, sy, sz]

    def spread(c):
        o = sum(ci * b for ci, b in zip(c, basis))
        ev = np.linalg.eigvalsh(h + np.kron(o, I) - np.kron(I, o))
        return (ev[-1] - ev[0]) / 2
    best = min((minimize(spread, x0, method='Nelder-Mead',
                         options=dict(xatol=1e-10, fatol=1e-12, maxiter=4000))
                for x0 in [np.zeros(3), np.array([.3, -.2, .1])]), key=lambda r: r.fun)
    return float(best.fun)


def run(mode):
    start = time.time()
    rows = []
    for name, model, delta, hz, M in MODELS:
        if mode == 'typhon_qubit' and model == 'su3': continue
        if mode == 'typhon_su3' and model != 'su3': continue
        if M is None:
            M = bond_M(model, delta, hz)
        if mode == 'local':
            Ls = [3] if model == 'su3' else [4, 6]
            grid = [(a, N) for a in [0.25, 1.0, 2.0] for N in [1, 2, 4]]
        else:
            Ls = [6] if model == 'su3' else [8]
            grid = [(a, N) for a in [0.25, 1.0, 2.0] for N in [1, 2, 4]]
        for L in Ls:
            sym = (hz == 0.0)
            small = blocks(L, model, delta, hz, reduce=sym)
            small_full = blocks(L, model, delta, hz, reduce=False)
            big = blocks(2 * L, model, delta, hz, reduce=sym)
            maxblock = max(len(b['ev']) for b in big['blocks'])
            for a, N in grid:
                beta = N * a * L / M
                row = decomposition(small, big, beta, N, small_full)
                row.update(model=name, a=a, M=M, maxblock=maxblock)
                d = row['d']
                # proved sandwich (C121 (33a)) and pinching bounds
                assert -1e-9 <= row['epsL'] <= N * np.log(d) + 1e-9
                assert -1e-9 <= row['epsD'] <= 2 * N * np.log(d) + 1e-9
                assert -1e-9 <= row['H2'] <= N * np.log(d) + 1e-9
                assert abs(row['identity_residual']) < 1e-8
                row['G_le_R'] = bool(row['G'] <= row['R'] + 1e-10)
                rows.append(row)
                print(f"{name:14s} L={L} a={a} N={N} R={row['R']:.5f} G={row['G']:.5f} "
                      f"G-R={row['GmR']:.5f} pred={row['cluster_pred_GmR']:.5f} "
                      f"epsL={row['epsL']:.4f}/{row['cluster_pred_epsL']:.4f} "
                      f"t={time.time()-start:.0f}s", flush=True)
            (ROOT / f'ed_{mode}.json').write_text(json.dumps(rows, indent=1) + '\n')
    print('done', len(rows), 'cases', time.time() - start)


if __name__ == '__main__':
    run(sys.argv[1] if len(sys.argv) > 1 else 'local')
