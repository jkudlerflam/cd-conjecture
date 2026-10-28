"""Closed-book verification of the explicit Clifford/positive-weight recipe.
All tensor orders are p_1,...,p_m,s_1,...,s_m,c.
Run with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 verify.py.
"""
import itertools
import json
import time
import numpy as np

I = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)
H = (X + Z) / np.sqrt(2)
S = np.diag([1, 1j])
PAULI = [I, X, Y, Z]
NAMES = 'IXYZ'

def kron_all(matrices):
    out = np.array([[1]], dtype=complex)
    for a in matrices:
        out = np.kron(out, a)
    return out

def pauli(labels):
    return kron_all([PAULI[k] for k in labels])

def eta(labels):
    return (-1) ** labels.count(2)

def character(a, b):
    return (-1) ** sum(x != 0 and y != 0 and x != y for x, y in zip(a, b))

# Derived directly from 2-by-2 matrix products, not hard-coded phases.
PRODUCT = {}
for a in range(4):
    for b in range(4):
        prod = PAULI[a] @ PAULI[b]
        coeffs = [np.trace(p.conj().T @ prod) / 2 for p in PAULI]
        k = int(np.argmax(np.abs(coeffs)))
        PRODUCT[a, b] = (k, coeffs[k])

def multiply(a, b):
    result, phase = [], 1 + 0j
    for x, y in zip(a, b):
        k, v = PRODUCT[x, y]
        result.append(k)
        phase *= v
    return tuple(result), phase

def gate_matrix(gate, n):
    name, *sites = gate
    if name == 'SWAP':
        a, b = sites
        return gate_matrix(('CX', a, b), n) @ gate_matrix(('CX', b, a), n) @ gate_matrix(('CX', a, b), n)
    if name == 'CX':
        control, target = sites
        d = 2 ** n
        out = np.zeros((d, d), complex)
        for col in range(d):
            row = col ^ ((1 << (n - 1 - target)) if (col >> (n - 1 - control)) & 1 else 0)
            out[row, col] = 1
        return out
    mats = [I] * n
    mats[sites[0]] = {'H': H, 'S': S, 'Sd': S.conj().T,
                      'X': X, 'Y': Y, 'Z': Z}[name]
    return kron_all(mats)

def inverse_gate(gate):
    if gate[0] == 'S':
        return ('Sd', gate[1])
    if gate[0] == 'Sd':
        return ('S', gate[1])
    return gate

def normalizer(labels):
    """Chronological gates for L D L^dagger = Z_c, no sign."""
    n, gates = len(labels), []
    c = n - 1
    support = [j for j, x in enumerate(labels) if x]
    assert support
    for j, x in enumerate(labels):
        if x == 1:
            gates.append(('H', j))
        elif x == 2:
            gates.extend([('Sd', j), ('H', j)])
    if c not in support:
        j = support[0]
        gates.extend([('CX', j, c), ('CX', c, j), ('CX', j, c)])
        support.remove(j)
        support.append(c)
    for j in support:
        if j != c:
            gates.append(('CX', j, c))
    return gates

def circuit_matrix(gates, n, cache):
    out = np.eye(2 ** n, dtype=complex)
    for gate in gates:
        if gate not in cache:
            cache[gate] = gate_matrix(gate, n)
        out = cache[gate] @ out
    return out

def target_terms(m, qlabels):
    terms = []
    for a in (1, 2):
        f = (a,) + tuple(qlabels)
        for p in itertools.product(range(4), repeat=m):
            if character(f, p) == -1:
                continue
            r, phase = multiply(f, p)
            assert abs(phase.imag) < 1e-14
            if not any(r):
                continue
            sigma = int(round(phase.real)) * eta(p)
            terms.append((r, p + (a,), sigma))
    assert len(terms) == 4 ** m - 2
    return terms

def verify_one(m, qlabels):
    q, n = 2 ** m, m + 1
    d = 2 * q * q
    omega = np.eye(q).reshape(q * q) / np.sqrt(q)
    pi0 = np.outer(omega, omega)
    pi = np.kron(pi0, I)
    target = np.zeros((d, d), complex)
    local_sc = np.zeros((2 * q, 2 * q), complex)
    for a in (1, 2):
        f = pauli((a,) + tuple(qlabels))
        target += np.kron(np.kron(f, np.eye(q)) @ pi0 + pi0 @ np.kron(f, np.eye(q)), PAULI[a])
        local_sc += 2 / q ** 2 * np.kron(f.T, PAULI[a])
    local = np.kron(np.eye(q), local_sc)
    terms = target_terms(m, qlabels)
    expansion = sum((2 / q ** 2 * sigma * np.kron(pauli(r), pauli(b))
                     for r, b, sigma in terms), np.zeros_like(target))
    accumulated = np.zeros_like(target)
    cache = {}
    frames = list(itertools.product(range(4), repeat=m))
    frame_matrices = {u: np.kron(pauli(u), I) for u in frames}
    frame_count, max_gate_count = 0, 0
    max_mapping_error, max_normalizer_error = 0., 0.
    for r, b, sigma in terms:
        lr, lb = normalizer(r + (0,)), normalizer(b)
        gates_c = lb + [inverse_gate(g) for g in reversed(lr)]
        C = circuit_matrix(gates_c, n, cache)
        source = np.kron(pauli(r).T, I)
        max_mapping_error = max(max_mapping_error, float(np.max(np.abs(C.conj().T @ source @ C - eta(r) * pauli(b)))))
        for labels, gates_l in ((r + (0,), lr), (b, lb)):
            L = circuit_matrix(gates_l, n, cache)
            max_normalizer_error = max(max_normalizer_error, float(np.max(np.abs(L @ pauli(labels) @ L.conj().T - np.kron(np.eye(q), Z)))))
        t = sigma * eta(r)
        for u in frames:
            if character(r, u) != t:
                continue
            W = frame_matrices[u] @ C
            # Columns are (I_p tensor W^dagger)(Omega tensor |a>_c).
            V = W.conj().T.reshape(2 * q, q, 2).transpose(1, 0, 2).reshape(d, 2) / np.sqrt(q)
            accumulated += 4 / q ** 2 * (V @ V.conj().T)
            frame_count += 1
            max_gate_count = max(max_gate_count, len(gates_c) + sum(x != 0 for x in u))
    alpha = 2 - 4 / q ** 2
    residual = accumulated - alpha * np.eye(d) + local - target
    # This partial trace detects the impossible s-only Q conjugation.
    base = np.zeros_like(target)
    for a in (1, 2):
        f = pauli((a,) + (0,) * (m - 1))
        base += np.kron(np.kron(f, np.eye(q)) @ pi0 + pi0 @ np.kron(f, np.eye(q)), PAULI[a])
    qs = pauli((0,) + tuple(qlabels)).T
    transform = np.kron(np.eye(q), np.kron(qs, I))
    difference = transform.conj().T @ base @ transform - target
    difference_trace_s = np.einsum('ascbsd->acbd', difference.reshape(q, q, 2, q, q, 2)).reshape(2*q, 2*q)
    return {
        'm': m, 'Q': ''.join(NAMES[x] for x in qlabels) or 'empty',
        'frames': frame_count, 'weight': 4 / q ** 2, 'alpha': alpha,
        'g': 1, 'max_gates': max_gate_count,
        'residual_max_abs': float(np.max(np.abs(residual))),
        'residual_fro': float(np.linalg.norm(residual)),
        'relative_fro': float(np.linalg.norm(residual) / np.linalg.norm(target)),
        'pauli_expansion_max_abs': float(np.max(np.abs(expansion - (target - local)))),
        'mapping_max_abs': max_mapping_error,
        'normalizer_max_abs': max_normalizer_error,
        'uniformity_trace_s_fro': float(np.linalg.norm(difference_trace_s)),
        'trace_error': float(abs(np.trace(accumulated) - 2 * frame_count * 4 / q ** 2)),
    }

def check_m1_bell_frames():
    omega = np.eye(2).reshape(4) / np.sqrt(2)
    pi = np.kron(np.outer(omega, omega), I)
    swap = gate_matrix(('CX', 0, 1), 2) @ gate_matrix(('CX', 1, 0), 2) @ gate_matrix(('CX', 0, 1), 2)
    result = np.zeros((8, 8), complex)
    wrong_order_errors = []
    # f1=1, f2=3, f=2, with order I,Z,X,Y.
    for label, weight in ((0, 2), (3, 2), (1, 3), (2, 1)):
        pc = np.kron(I, PAULI[label])
        W = swap @ pc
        big = np.kron(I, W)
        result += weight * big.conj().T @ pi @ big
        wrong = np.kron(I, pc @ swap)
        swapped = np.kron(I, swap)
        wrong_order_errors.append(float(np.max(np.abs(wrong.conj().T @ pi @ wrong - swapped.conj().T @ pi @ swapped))))
    port = (kron_all([X, I, X]) + kron_all([Y, I, Y])) / 2
    return {'residual_max_abs': float(np.max(np.abs(result - 2*np.eye(8) - port))),
            'wrong_order_frame_variation': max(wrong_order_errors)}

if __name__ == '__main__':
    start = time.monotonic()
    rows = []
    for m in (1, 2, 3):
        for labels in itertools.product(range(4), repeat=m-1):
            row = verify_one(m, labels)
            rows.append(row)
            print(json.dumps(row), flush=True)
    output = {'rows': rows, 'm1_four_bell_frames': check_m1_bell_frames(),
              'elapsed_seconds': time.monotonic() - start}
    with open('verification.json', 'w') as handle:
        json.dump(output, handle, indent=2)
    print(json.dumps({'m1_four_bell_frames': output['m1_four_bell_frames'],
                      'elapsed_seconds': output['elapsed_seconds']}), flush=True)
