"""Finite checks of the graph-charge target; no thermodynamic inference."""
from itertools import combinations, combinations_with_replacement, product
import json
from math import pi
from pathlib import Path
import numpy as np


def projective_representatives(p):
    reps = []
    for x in product(range(p), repeat=3):
        if not any(x):
            continue
        first = next(a for a in x if a)
        if first == 1:
            reps.append(x)
    return np.asarray(reps, dtype=int)


def projective_graph(p):
    x = projective_representatives(p)
    incidence = (x @ x.T) % p == 0
    m = len(x)
    edges = [(i, m + j) for i, j in zip(*np.nonzero(incidence))]
    return incidence, edges


def incidence_matrix(q, edges):
    w = np.zeros((q, len(edges)), dtype=int)
    for a, (left, right) in enumerate(edges):
        w[left, a] = -1
        w[right, a] = 1
    return w


def short_kernel_vectors(w, max_l1=4):
    # Each canonical signed multiset enumerates every nonzero integer
    # vector up to this l1 radius; opposite signs on the same index are
    # excluded so the multiset size is its actual l1 norm.
    n = w.shape[1]
    checked = 0
    witnesses = []
    for size in range(1, max_l1 + 1):
        for indices in combinations_with_replacement(range(n), size):
            unique = sorted(set(indices))
            for signs in product((-1, 1), repeat=len(unique)):
                z = np.zeros(n, dtype=int)
                for a, sign in zip(unique, signs):
                    z[a] = sign * indices.count(a)
                checked += 1
                if not np.any(w @ z):
                    witnesses.append(z.tolist())
    return {"vectors_checked": checked, "nonzero_neutral_vectors": witnesses}


def fermion_creation(n, a):
    size = 1 << n
    op = np.zeros((size, size))
    for state in range(size):
        if state & (1 << a):
            continue
        parity = (state & ((1 << a) - 1)).bit_count() % 2
        op[state | (1 << a), state] = (-1) ** parity
    return op


def choi_check():
    # Three vertices on each side, the complete bipartite graph minus a
    # perfect matching, hence one six-cycle and degree two.
    q, n = 6, 6
    edges = [(a, 3 + b) for a in range(3) for b in range(3) if a != b]
    dim = 1 << n
    parity = np.diag([(-1) ** s.bit_count() for s in range(dim)])
    norms = []
    for right in (False, True):
        choi = np.zeros((q * dim, q * dim))
        for a, (left, head) in enumerate(edges):
            fdag = fermion_creation(n, a)
            val = fdag @ parity if right else fdag
            choi[head*dim:(head+1)*dim, left*dim:(left+1)*dim] = val
            choi[left*dim:(left+1)*dim, head*dim:(head+1)*dim] = val.T
        norms.append(float(np.max(np.abs(np.linalg.eigvalsh(choi)))))
    assert max(norms) <= 2 + 1e-12
    return {"graph": "six-cycle", "q": q, "N": n,
            "degree_bound": 2, "left_right_choi_norms": norms,
            "epsilon_sufficient": 1 / (4 * q * 2)}


def momentum_check(n):
    k = pi / 4 + np.arange(1, n + 1) * pi / (12 * (n + 1))
    dist = lambda x: np.abs((x + pi) % (2 * pi) - pi)
    minimum = float(np.min(dist(2 * k)))
    for a, b in combinations(range(n), 2):
        minimum = min(minimum, float(dist(2 * (k[a] - k[b]))),
                      float(dist(2 * (k[a] + k[b]))))
    claimed = pi / (6 * (n + 1))
    assert abs(minimum - claimed) < 5e-14
    length_k = pi * np.sqrt(2) / 12
    mixed_min = float(min(np.min(dist(2*(k+length_k))),
                          np.min(dist(2*(k-length_k)))))
    assert mixed_min > pi / 4
    assert float(dist(2 * length_k)) > pi / 5
    return {"N": n, "minimum_low_dimension_mismatch": minimum,
            "claimed_bound": claimed,
            "mobile_mixed_minimum": mixed_min,
            "pure_length_mismatch": float(dist(2 * length_k)),
            "velocity_range": [float(np.tan(k).min()), float(np.tan(k).max())]}


def main():
    result = {"projective_planes": [], "momenta": []}
    for p in (2, 3, 5, 7):
        incidence, edges = projective_graph(p)
        m = p*p + p + 1
        assert incidence.shape == (m, m)
        assert np.all(incidence.sum(axis=0) == p+1)
        assert np.all(incidence.sum(axis=1) == p+1)
        intersections = incidence.astype(int) @ incidence.astype(int).T
        offdiag = intersections[~np.eye(m, dtype=bool)]
        assert np.all(offdiag == 1)
        q, n = 2*m, len(edges)
        w = incidence_matrix(q, edges)
        rank = int(np.linalg.matrix_rank(w.astype(float)))
        assert rank == q - 1
        row = {"prime": p, "port_dimension": q, "channels": n,
               "degree": p+1, "four_cycles": 0,
               "incidence_rank": rank, "kernel_rank": n - rank}
        if p == 2:
            row["short_integer_kernel_check"] = short_kernel_vectors(w)
            assert not row["short_integer_kernel_check"]["nonzero_neutral_vectors"]
        result["projective_planes"].append(row)
        result["momenta"].append(momentum_check(n))
    result["compression"] = choi_check()
    path = Path(__file__).with_name("checks_results.json")
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
