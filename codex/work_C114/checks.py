"""C114 ring-only checks. Full mode is restricted to a SLURM allocation.

Local mode uses Hilbert spaces of dimension at most 729. No shared files
are imported or written. NumPy and SciPy implement every numerical step.
"""
import argparse
import json
import os
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh
from scipy.integrate import quad
from scipy.special import digamma


def ring_matrix(h, d, length):
    size = d**length
    ids = np.arange(size, dtype=np.int64)
    powers = d**np.arange(length, dtype=np.int64)
    rows, cols, vals = [], [], []
    for i in range(length):
        j = (i+1) % length
        a, b = (ids//powers[i]) % d, (ids//powers[j]) % d
        source = d*a+b
        for target, initial in zip(*np.nonzero(h)):
            chosen = ids[source == initial]
            rows.append(chosen+(target//d-initial//d)*powers[i]
                        +(target % d-initial % d)*powers[j])
            cols.append(chosen)
            vals.append(np.full(chosen.size, h[target, initial]))
    return sparse.coo_matrix((np.concatenate(vals),
                              (np.concatenate(rows), np.concatenate(cols))),
                             shape=(size, size)).tocsr()


def ground(H):
    rng = np.random.default_rng(114+H.shape[0])
    e, v = eigsh(H, k=1, which='SA', tol=2e-12,
                 v0=rng.normal(size=H.shape[0]))
    psi = v[:, 0]
    residual = float(np.linalg.norm(H@psi-e[0]*psi))
    assert residual < 2e-9, residual
    return float(e[0]), psi, residual


def translate(v, d, n):
    # The amplitudes are indexed with site zero as the least significant digit.
    return np.moveaxis(v.reshape([d]*n), 0, -1).reshape(-1)


def swap_rung(v, d, L, j):
    ids = np.arange(v.size, dtype=np.int64)
    pa, pb = d**j, d**(j+L)
    a, b = (ids//pa) % d, (ids//pb) % d
    target = ids+(b-a)*pa+(a-b)*pb
    return v[target]


def purity(psi, d, length, r):
    mat = psi.reshape(d**(length-r), d**r)
    rho = mat.conj().T@mat
    return float(np.trace(rho@rho).real)


def coherent_trial(H, psi, d, L, e):
    phi = np.kron(psi, psi)
    states = [phi]
    for j in range(1, L):
        states.append(swap_rung(states[-1], d, L, j-1))
    X = np.column_stack(states)
    G = X.conj().T@X
    K = X.conj().T@(H@X)-2*e*G
    g, u = eigh(G)
    keep = g > 1e-11
    W = u[:, keep]/np.sqrt(g[keep])
    k = W.conj().T@K@W
    best, coeff = eigh((k+k.conj().T)/2)
    trial = X@(W@coeff[:, 0])
    best_direct = float(np.vdot(trial, H@trial).real-2*e)
    assert abs(best_direct-best[0]) < 2e-10
    phase = np.vdot(psi, translate(psi, d, L))
    assert abs(abs(phase)-1) < 1e-9
    # Prefix swaps give interval purities without any momentum convention.
    predicted = np.array([purity(psi, d, L, r) for r in range(L)])
    gram_error = float(np.max(np.abs(G[0]-predicted)))
    assert gram_error < 1e-9, (phase, G[0], predicted)
    q = float(K[0, 0].real)
    real_kernel_error = float(np.max(np.abs(K.real-q*np.eye(L))))
    assert real_kernel_error < 2e-9
    rownorm = float(np.sum(G[0]).real)
    uniform = float(np.sum(K).real/np.sum(G).real)
    assert abs(uniform-q/rownorm) < 2e-9
    if np.max(np.abs(K.imag)) < 1e-10 and q > 0:
        assert abs(best[0]-q/rownorm) < 2e-9
    return dict(bare=float(K[0, 0].real), coherent=float(best[0]),
                uniform=uniform, gram_min=float(g.min()),
                gram_purity_error=gram_error, norm_per_position=rownorm,
                real_kernel_error=real_kernel_error,
                energy_kernel_row=[float(x.real) for x in K[0]],
                gram_row=[float(x.real) for x in G[0]])


def models():
    sx = np.array([[0., 1.], [1., 0.]])/2
    sy = np.array([[0., -1j], [1j, 0.]])/2
    sz = np.diag([.5, -.5])
    xx = (np.kron(sx, sx)+np.kron(sy, sy)).real
    zz = np.kron(sz, sz)
    data = []
    for name, delta in [('XX', 0.), ('XXZ_minus', -.5),
                        ('XXZ_plus', .5), ('Heisenberg', 1.)]:
        A = np.pi/2 if delta == 1 else np.pi*np.sqrt(1-delta**2)/(2*np.arccos(delta))
        data.append((name, xx+delta*zz, 2, .5, A))
    p = np.zeros((9, 9))
    for a in range(3):
        for b in range(3):
            p[3*b+a, 3*a+b] = 1.
    data.append(('SU3', p, 3, 1., 4*np.pi/3))
    return data


def inherited_bulk(name):
    if name == 'SU3':
        return 1+2*(digamma(1/3)+np.euler_gamma)/3
    delta = {'XX': 0., 'XXZ_minus': -.5, 'XXZ_plus': .5, 'Heisenberg': 1.}[name]
    if delta == 1.:
        return .25-np.log(2.)
    gamma = np.arccos(delta)
    def integrand(x):
        if x == 0.:
            return (np.pi-gamma)/np.pi
        return (2*np.exp(-2*gamma*x)*(-np.expm1(-2*(np.pi-gamma)*x))
                /((-np.expm1(-2*np.pi*x))*(1+np.exp(-2*gamma*x))))
    integral, error = quad(integrand, 0, np.inf, epsabs=1e-12, epsrel=1e-12)
    assert error < 1e-9
    return float(delta/4-np.sin(gamma)*integral)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--full', action='store_true')
    args = parser.parse_args()
    if args.full and not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Full ED requires a SLURM allocation. Use small mode locally.')
    report = ['# C114 ring checks', '',
              '[Numerical] Full spin-ring diagonalization with sparse matrices. '
              'Continuum A values are inherited calibration inputs, not fitted results.', '',
              '| model | L | E(L) | E(2L) | D_L | L D_L/M | L coherent/M | bare cost |',
              '|:--|--:|--:|--:|--:|--:|--:|--:|']
    results = []
    max_residual = 0.
    max_dimension = 0
    for name, h, d, M, A in models():
        lengths = ([4, 6, 8] if d == 2 else [3, 6]) if args.full else ([4] if d == 2 else [3])
        for L in lengths:
            H = ring_matrix(h, d, L)
            H2 = ring_matrix(h, d, 2*L)
            max_dimension = max(max_dimension, H2.shape[0])
            if not args.full:
                assert H2.shape[0] <= 729
            e, psi, residual = ground(H)
            e2, psi2, residual2 = ground(H2)
            max_residual = max(max_residual, residual, residual2)
            D = e2-2*e
            deficit = L*inherited_bulk(name)-e
            assert deficit > 1e-8
            variational = coherent_trial(H2, psi, d, L, e)
            assert D > 1e-8, (name, L, D)
            assert D <= variational['coherent']+1e-9
            assert variational['coherent'] <= variational['bare']+1e-9
            assert abs(D) <= 4*M+1e-9
            convex_gap = e2/(2*L)-e/L
            assert abs(convex_gap-D/(2*L)) < 1e-12
            if name == 'XX':
                assert abs(e+1/np.sin(np.pi/L)) < 1e-10
                assert abs(e2+1/np.sin(np.pi/(2*L))) < 1e-10
            if name == 'Heisenberg' and L == 4:
                assert D > np.pi**2*M/(4*L)
            row = dict(model=name, d=d, M=M, A=A, L=L,
                       energy=e, doubled_energy=e2, D=D,
                       deficit=deficit,
                       convex_gap=convex_gap, residual=max(residual, residual2),
                       **variational)
            results.append(row)
            report.append(f'| {name} | {L} | {e:.12f} | {e2:.12f} | {D:.9f} | '
                          f'{L*D/M:.9f} | {L*variational["coherent"]/M:.9f} | '
                          f'{variational["bare"]:.9f} |')
            print(name, L, 'D', D, 'coherent', variational['coherent'], flush=True)
    report += ['', '[Numerical] Every positive D_L refutes the zero-cost lift and '
               'the zero-error comparison E(L)/L >= E(2L)/(2L). The inequalities '
               'D_L <= optimized coherent cost <= bare cost <= 4M pass. '
               'Using the bulk energies inherited from C98, every tested E(L) '
               'is also strictly below L e0, excluding an exact finite-range '
               'positive certificate that wraps on that size. No approximate '
               'SOS certificate has been constructed or numerically tested.', '',
               f'[Numerical] Maximum eigenpair residual {max_residual:.3g}. '
               f'Largest Hilbert-space dimension {max_dimension}.', '',
               '[Numerical] All assertions passed. The Gram matrix agrees with '
               'single-ring interval purities. No asymptotic bound follows from '
               'these finite sizes.']
    suffix = '_FULL' if args.full else ''
    Path('CHECKS'+suffix+'.md').write_text('\n'.join(report)+'\n')
    Path('checks'+suffix+'.json').write_text(json.dumps(results, indent=2)+'\n')
    print('\n'.join(report))


if __name__ == '__main__':
    main()
