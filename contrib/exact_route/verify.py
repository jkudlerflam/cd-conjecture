"""Independent dense checks of the sparse covariant two-port isometry.

Run with Python 3 and numpy. At most four fermion modes are used; no
many-body diagonalization or infrared inference is performed.
"""

from pathlib import Path
import json
import time
from collections import Counter

import numpy as np


def annihilators(n):
    dim = 1 << n
    ops = []
    for edge in range(n):
        op = np.zeros((dim, dim))
        for state in range(dim):
            if (state >> edge) & 1:
                sign = (-1) ** ((state & ((1 << edge) - 1)).bit_count())
                op[state ^ (1 << edge), state] = sign
        ops.append(op)
    return ops


def kraus(u, v, edges, annihilation, a, b):
    dim = annihilation[0].shape[0]
    count = len(edges)
    out = np.zeros((u + v, u + count, dim, dim))
    for alpha in range(u):
        out[alpha, alpha] = a * np.eye(dim)
    for index, (alpha, beta) in enumerate(edges):
        out[u + beta, alpha] = b * annihilation[index].T
        out[u + beta, u + index] = b * annihilation[index]
    return out


def channel(k, op):
    dim = op.shape[0]
    result = np.zeros((dim, dim))
    for entry in k.reshape((-1, dim, dim)):
        result += entry.T @ op @ entry
    return result


def max_abs(matrix):
    return float(np.max(np.abs(matrix)))


def stable_route(source, target):
    working = list(source)
    swaps = []
    for position, wanted in enumerate(target):
        found = next(j for j in range(position, len(working)) if working[j] == wanted)
        for j in range(found - 1, position - 1, -1):
            working[j], working[j + 1] = working[j + 1], working[j]
            swaps.append(j)
    assert working == list(target)
    return swaps


def basis_transducer_check(u, v, edges, left, right):
    """Exhaust every basis branch, including its physical word and inverse."""
    n, q = len(edges), u + v
    species = np.eye(q, dtype=int)
    reference = tuple(s for s in range(q) for _ in range(2))
    flags = (
        [2 * np.ones(q, dtype=int) - species[alpha] for alpha in range(u)]
        + [
            2 * np.ones(q, dtype=int) + species[alpha] - 2 * species[u + beta]
            for alpha, beta in edges
        ]
    )
    checked = max_swaps = 0
    for dressed, kraus_map in ((False, left), (True, right)):
        seen = set()
        for source in range(1 << n):
            input_word = tuple(
                alpha if (source >> e) & 1 else u + beta
                for e, (alpha, beta) in enumerate(edges)
            ) + reference
            for branch in range(u + n):
                if branch < u:
                    port, flag, out, sign = branch, branch, source, 1
                else:
                    e = branch - u
                    alpha, beta = edges[e]
                    occupied = (source >> e) & 1
                    port = u + beta
                    flag = u + e if occupied else alpha
                    out = source ^ (1 << e)
                    sign = (-1) ** ((source & ((1 << e) - 1)).bit_count())
                    if dressed:
                        sign *= (-1) ** (source.bit_count() - occupied)
                amplitude = kraus_map[port, flag, out, source]
                assert np.isclose(amplitude, sign / np.sqrt(2 * (u if branch < u else n)))
                output_word = (port,) + tuple(
                    s for s in range(q) for _ in range(int(flags[flag][s]))
                ) + tuple(
                    alpha if (out >> e) & 1 else u + beta
                    for e, (alpha, beta) in enumerate(edges)
                )
                assert Counter(input_word) == Counter(output_word)
                assert output_word not in seen
                seen.add(output_word)
                swaps = stable_route(input_word, output_word)
                max_swaps = max(max_swaps, len(swaps))

                # The inverse recovers the branch and the original occupation.
                if port < u:
                    recovered_branch, recovered_source = port, out
                elif flag < u:
                    e = edges.index((flag, port - u))
                    assert (out >> e) & 1
                    recovered_branch, recovered_source = u + e, out ^ (1 << e)
                else:
                    e = flag - u
                    assert not ((out >> e) & 1)
                    assert edges[e][1] == port - u
                    recovered_branch, recovered_source = u + e, out ^ (1 << e)
                assert (recovered_branch, recovered_source) == (branch, source)
                checked += 1
    return {"basis_branches_checked": checked, "largest_stable_routing_swap_count": max_swaps}


def head_gate_check(q):
    """Check the complete local permutation gates, and their 3-gate action."""
    ell = (q - 1).bit_length()
    types = (
        [("D", s) for s in range(q)]
        + [("blank",), ("bit", 0), ("bit", 1)]
        + [("H", j, b) for j in range(ell) for b in range(2)]
        + [("T", j, b, s) for j in range(ell) for b in range(2) for s in range(q)]
    )

    def charge(x):
        return Counter({x[1] if x[0] == "D" else x[3]: 1}) if x[0] in ("D", "T") else Counter()

    def load(pair):
        x, y = pair
        if x[0] == "D" and y[0] == "H":
            return ("blank",), ("T", y[1], y[2], x[1])
        if x == ("blank",) and y[0] == "T":
            return ("D", y[3]), ("H", y[1], y[2])
        return pair

    def exchange(pair):
        x, y = pair
        if x[0] == "T" and y[0] == "D" and x[2] == 1:
            return ("T", x[1], 1, y[1]), ("D", x[3])
        return pair

    def query(pair):
        x, y = pair
        if x[0] == "D" and y[0] == "H":
            return x, ("H", y[1], y[2] ^ ((x[1] >> y[1]) & 1))
        return pair

    for gate in (load, exchange, query):
        images = set()
        for x in types:
            for y in types:
                result = gate((x, y))
                images.add(result)
                assert gate(result) == (x, y)
                assert charge(x) + charge(y) == charge(result[0]) + charge(result[1])
        assert len(images) == len(types) ** 2
    checked = 0
    for j in range(ell):
        for b in range(2):
            for s in range(q):
                for t in range(q):
                    x, head = load((("D", s), ("H", j, b)))
                    head, y = exchange((head, ("D", t)))
                    x, head = load((x, head))
                    assert head == ("H", j, b)
                    assert (x, y) == ((("D", t), ("D", s)) if b else (("D", s), ("D", t)))
                    checked += 1
    return {"q": q, "local_alphabet": len(types), "controlled_swaps_checked": checked}


def check(u, v, edges):
    n = len(edges)
    dim = 1 << n
    q = u + v
    flags = u + n
    a = 1 / np.sqrt(2 * u)
    b = 1 / np.sqrt(2 * n)
    eps = a * b
    f = annihilators(n)
    parity = np.diag([(-1) ** state.bit_count() for state in range(dim)])
    g = [parity @ op for op in f]
    left = kraus(u, v, edges, f, a, b)
    right = kraus(u, v, edges, g, a, b)

    # Compose explicit operator coefficients in output order
    # left port, left flag, right port, right flag, residual Fock, input.
    joint = np.einsum("rsij,abjk->abrsik", right, left)
    matrix = joint.reshape((-1, dim))
    errors = {"normalization": max_abs(matrix.T @ matrix - np.eye(dim))}
    errors["one_port_normalization"] = max_abs(
        channel(left, np.eye(dim)) - np.eye(dim)
    )

    edge_index = {edge: index for index, edge in enumerate(edges)}
    left_error = right_error = missing_error = 0.0
    for alpha in range(u):
        for beta in range(v):
            lo = joint[alpha].reshape((-1, dim))
            li = joint[u + beta].reshape((-1, dim))
            ro = joint[:, :, alpha].reshape((-1, dim))
            ri = joint[:, :, u + beta].reshape((-1, dim))
            compressed_left = lo.T @ li
            compressed_right = ro.T @ ri
            edge = (alpha, beta)
            if edge in edge_index:
                creation = f[edge_index[edge]].T
                left_error = max(left_error, max_abs(compressed_left - eps * creation))
                right_error = max(
                    right_error,
                    max_abs(compressed_right - eps * (1 - b * b) * creation @ parity),
                )
            else:
                missing_error = max(
                    missing_error, max_abs(compressed_left), max_abs(compressed_right)
                )
    errors.update(left_port=left_error, right_port=right_error, missing_edge=missing_error)
    errors["dressed_channel"] = max(
        max_abs(channel(left, op.T @ parity) - (1 - b * b) * op.T @ parity)
        for op in f
    )
    errors["undressed_channel"] = max(
        max_abs(channel(left, op.T) - (1 - (2 * n - 1) * b * b) * op.T)
        for op in f
    )
    errors["dressed_car"] = max(
        max_abs(g[i] @ g[j].T + g[j].T @ g[i] - (np.eye(dim) if i == j else 0))
        for i in range(n) for j in range(n)
    )

    # Charge covariance is checked on all nonzero coefficients of the
    # full composed isometry, not inferred from the individual stages.
    species = np.eye(q, dtype=int)
    weights = np.array([species[alpha] - species[u + beta] for alpha, beta in edges])
    logical_charge = np.array([
        sum((weights[e] for e in range(n) if (state >> e) & 1), np.zeros(q, dtype=int))
        for state in range(dim)
    ])
    offset = 2 * np.ones(q, dtype=int)
    flag_charge = np.array(
        [offset - species[alpha] for alpha in range(u)]
        + [offset + species[alpha] - 2 * species[u + beta] for alpha, beta in edges]
    )
    charge_error = parity_error = 0
    for lp, la, rp, ra, output, source in np.argwhere(joint != 0):
        out_charge = (
            species[lp] + flag_charge[la] + species[rp] + flag_charge[ra]
            + logical_charge[output]
        )
        expected = 2 * offset + logical_charge[source]
        charge_error = max(charge_error, int(np.max(np.abs(out_charge - expected))))
        out_parity = int(out_charge[:u].sum()) % 2
        parity_error = max(parity_error, int(out_parity != source.bit_count() % 2))
    errors["species_covariance"] = charge_error
    errors["parity_covariance"] = parity_error
    errors["flag_length"] = int(np.max(np.abs(flag_charge.sum(axis=1) - (2 * q - 1))))
    errors["negative_flag_charge"] = int(max(0, -flag_charge.min()))

    assert max(errors.values()) < 1e-12, errors
    classical_checks = basis_transducer_check(u, v, edges, left, right)
    # The physical exchange blocks are built independently as matrix units.
    link = np.zeros((q * q, q * q), dtype=complex)
    hopping = [complex(index + 1, 0.3 * index) / (n + 1) for index in range(n)]
    for (alpha, beta), coefficient in zip(edges, hopping):
        source = (u + beta) * q + alpha
        target = alpha * q + u + beta
        link[target, source] = coefficient
        link[source, target] = coefficient.conjugate()
    errors["link_norm"] = abs(float(np.max(np.abs(np.linalg.eigvalsh(link)))) - max(map(abs, hopping)))
    assert max(errors.values()) < 1e-12, errors
    return {
        "u": u, "v": v, "edges": edges, "logical_modes": n,
        "port_dimension": q, "flags_per_stage": flags,
        "joint_matrix_shape": list(matrix.shape),
        "left_visibility": eps, "right_visibility": eps * (1 - b * b),
        "logical_hopping_factor": eps * eps * (1 - b * b),
        "max_errors": errors,
        "basis_transducer": classical_checks,
    }


def main():
    started = time.perf_counter()
    cases = [
        (1, 1, [(0, 0)]),
        (2, 2, [(0, 0), (0, 1), (1, 1)]),
        (2, 2, [(0, 0), (0, 1), (1, 0), (1, 1)]),
        (3, 2, [(0, 0), (0, 1), (1, 1)]),
    ]
    results = [check(*case) for case in cases]
    heads = [head_gate_check(q) for q in (2, 3, 4, 5)]
    output = {"cases": results, "head_gates": heads, "runtime_seconds": time.perf_counter() - started}
    path = Path(__file__).with_name("verification.json")
    path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({
        "cases": len(results),
        "maximum_error": max(max(case["max_errors"].values()) for case in results),
        "runtime_seconds": output["runtime_seconds"],
        "output": str(path),
    }, indent=2))


if __name__ == "__main__":
    main()
