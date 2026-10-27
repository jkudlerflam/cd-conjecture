"""Small checks for C7_result.md. No input files, network, or large computations.

Run with python3 c7.py. These checks supplement the proofs in the note;
they do not establish a conformal counterexample.
"""

import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh


# Alphabet: marker, L, X, R. The permitted pairs concatenate arrows
# marker=(R,L), L=(L,L), X=(L,R), R=(R,R).
MARKER, LEFT, HEAD, RIGHT = range(4)
ALLOWED = {(0, 1), (0, 2), (1, 1), (1, 2),
           (2, 0), (2, 3), (3, 0), (3, 3)}
W = np.array([[int((a, b) not in ALLOWED) for b in range(4)]
              for a in range(4)], dtype=float)


def a_cell(b, t=1.0):
    # Stable evaluation of 2t(1-cos(pi/b)).
    return 4.0 * t * np.sin(np.pi / (2.0 * np.asarray(b))) ** 2


def delta_boundary(b, t=1.0):
    return (b + 1) * a_cell(b, t) - b * a_cell(b + 1, t)


def selected_delta(b, t=1.0):
    return (delta_boundary(b, t) + delta_boundary(b - 1, t)) / 2


def algebra_checks():
    hop = np.zeros((16, 16))
    hop[4 * LEFT + HEAD, 4 * HEAD + RIGHT] = -1
    hop[4 * HEAD + RIGHT, 4 * LEFT + HEAD] = -1
    q = np.diag(W.reshape(-1))
    marker = np.diag([1., 0., 0., 0.])
    eye = np.eye(4)
    pairs = [(q, hop),
             (np.kron(q, eye), np.kron(eye, hop)),
             (np.kron(eye, q), np.kron(hop, eye)),
             (np.kron(marker, eye), hop),
             (np.kron(eye, marker), hop)]
    errors = [np.max(np.abs(x @ y - y @ x)) for x, y in pairs]
    assert max(errors) == 0
    print('Local marker and forbidden-bond commutators:', errors)


def period_checks():
    print('Period selection: b, delta, energy/site, neighbor cell cost')
    for b in (5, 8, 16, 32, 64):
        delta = selected_delta(b)
        lengths = np.arange(2, max(1024, 8 * b) + 1)
        energies = (a_cell(lengths) - delta) / lengths
        best = int(lengths[np.argmin(energies)])
        assert best == b
        e = (a_cell(b) - delta) / b
        neighbor_cost = min((a_cell(k) - delta) - k * e
                            for k in (b - 1, b + 1))
        assert e < 0 and neighbor_cost > 0 and 0 < delta < 2
        print(f'  {b:3d}  {delta:.12g}  {e:.12g}  {neighbor_cost:.12g}')


def ring_hamiltonian(length, delta, penalty=8.0):
    size = 4 ** length
    indices = np.arange(size, dtype=np.int64)
    powers = 4 ** np.arange(length, dtype=np.int64)
    digits = (indices[:, None] // powers[None, :]) % 4
    diagonal = (2.0 - delta) * np.sum(digits == MARKER, axis=1)
    rows, cols, vals = [], [], []
    for i in range(length):
        j = (i + 1) % length
        diagonal += penalty * W[digits[:, i], digits[:, j]]
        forward = indices[(digits[:, i] == LEFT) & (digits[:, j] == HEAD)]
        backward = indices[(digits[:, i] == HEAD) & (digits[:, j] == RIGHT)]
        rows.extend([forward, backward])
        cols.extend([forward + powers[i] + powers[j],
                     backward - powers[i] - powers[j]])
        vals.extend([-np.ones(len(forward)), -np.ones(len(backward))])
    matrix = sparse.coo_matrix((np.concatenate(vals),
                                (np.concatenate(rows), np.concatenate(cols))),
                               shape=(size, size)).tocsr()
    return matrix + sparse.diags(diagonal, format='csr')


def ring_checks():
    rng = np.random.default_rng(7)
    for b in (5, 6):
        delta = selected_delta(b)
        matrix = ring_hamiltonian(b, delta)
        # An unshifted Lanczos run can omit the isolated all-L/all-R
        # zero-energy states. Shift the matrix and explicitly check them.
        shift = 0.371
        shifted = matrix + shift * sparse.eye(matrix.shape[0], format='csr')
        values = np.sort(eigsh(shifted, k=24, which='SA', tol=1e-11,
                               v0=rng.normal(size=matrix.shape[0]),
                               return_eigenvectors=False)) - shift
        for label in (LEFT, RIGHT):
            product_index = sum(label * 4**j for j in range(b))
            assert matrix.getrow(product_index).nnz == 0
        assert np.sum(np.abs(values) < 1e-8) == 2
        predicted = a_cell(b) - delta
        ground_count = int(np.sum(np.abs(values - predicted) < 1e-8))
        assert abs(values[0] - predicted) < 1e-8
        assert ground_count == b
        gap = values[ground_count] - values[0]
        assert gap > 0
        print(f'Ring L=b={b}: E0={values[0]:.12g}, '
              f'prediction={predicted:.12g}, degeneracy={ground_count}, '
              f'next gap={gap:.12g}')


def modular_functional_checks():
    bound = 1.0 / (2.0 * np.pi * np.e)
    for c in (1.0, 2.0, 10.0, 100.0):
        C = c / 12.0
        coefficient = C / 2.0 * np.exp(-2 * np.pi * C) / (-np.expm1(-3*np.pi*C))

        def functional(E):
            return (E * np.exp(-2*np.pi*E)
                    + coefficient * (np.exp(-4*np.pi*E) - np.exp(-np.pi*E)))

        vacuum = functional(-C)
        target = -C / 2.0 * np.exp(2*np.pi*C)
        assert np.isclose(vacuum, target, rtol=1e-12)
        assert np.max(functional(np.linspace(-C, 0, 10001))) <= 1e-12
        assert np.max(functional(np.linspace(0, 2*C+2, 10001))) <= bound + 1e-12
        assert np.max(functional(np.linspace(2*C+2, 2*C+20, 10001))) <= 1e-12
        print(f'Modular functional c={c:g}: signs and maximum checked; '
              f'Delta cutoff={3*C+2:g}')


def swap_rank_checks():
    for n in (2, 3, 4):
        swap = np.zeros((n, n, n, n))
        for i in range(n):
            for j in range(n):
                swap[i, j, j, i] = 1
        identity = np.eye(n*n).reshape(n, n, n, n)
        connected = swap - identity / n
        reshuffled = connected.transpose(0, 2, 1, 3).reshape(n*n, n*n)
        rank = np.linalg.matrix_rank(reshuffled)
        assert rank == n*n - 1
        print(f'Connected swap N={n}: operator Schmidt rank {rank}')


if __name__ == '__main__':
    start = time.monotonic()
    algebra_checks()
    period_checks()
    ring_checks()
    modular_functional_checks()
    swap_rank_checks()
    print(f'All checks passed in {time.monotonic() - start:.3f} seconds.')
