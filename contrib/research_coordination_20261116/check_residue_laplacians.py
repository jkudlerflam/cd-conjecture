"""Independent complete-space check of the two sixth-order residue models.

This checks the published effective operators, not their microscopic
Schrieffer--Wolff derivation or any finite-coupling conformal limit.
All matrix-identity and graph-kernel checks use exact integer arithmetic.
"""
import json
from pathlib import Path
import time

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.csgraph import connected_components


def rotate(state, distance, length):
    distance %= length
    return ((state >> distance) | (state << (length - distance))) & ((1 << length) - 1)


def pair_count(state, distance, length):
    return (state & rotate(state, distance, length)).bit_count()


def residue_word(state, period, length):
    word = tuple(j % period for j in range(length) if state >> j & 1)
    return min((word[j:] + word[:j] for j in range(len(word))), default=())


def swap_laplacian(length):
    answer = np.zeros((1 << length, 1 << length), dtype=np.int64)
    for state in range(1 << length):
        for j in range(length):
            k = (j + 1) % length
            if (state >> j & 1) != (state >> k & 1):
                answer[state, state] += 1
                answer[state ^ (1 << j) ^ (1 << k), state] -= 1
    return answer


def check(period, length):
    forbidden = (1, 2, 3) if period == 4 else (1, 2, 4)
    states = [s for s in range(1 << length)
              if all(pair_count(s, r, length) == 0 for r in forbidden)]
    index = {s: i for i, s in enumerate(states)}
    numbers = np.array([s.bit_count() for s in states], dtype=np.int64)
    pairs = {r: np.array([pair_count(s, r, length) for s in states], dtype=np.int64)
             for r in range(3, 8)}
    rows, cols = [], []
    words = [residue_word(s, period, length) for s in states]
    for col, state in enumerate(states):
        for j in range(length):
            if not (state >> j & 1):
                continue
            for step in (-period, period):
                k = (j + step) % length
                if state >> k & 1:
                    continue
                target = state ^ (1 << j) ^ (1 << k)
                row = index.get(target)
                if row is not None:
                    assert words[row] == words[col]
                    rows.append(row)
                    cols.append(col)
    adjacency = coo_matrix((np.ones(len(rows), dtype=np.int64), (rows, cols)),
                           shape=(len(states), len(states))).tocsr()
    assert (adjacency - adjacency.T).nnz == 0
    degree = np.asarray(adjacency.sum(axis=0)).ravel()
    if period == 4:
        hopping_times_three = 3
        diagonal = 6 * numbers - 6 * pairs[4] + 10 * pairs[5] + 14 * pairs[6]
        potential = 16 * pairs[5] + 20 * pairs[6] + 6 * pairs[7]
    else:
        hopping_times_three = 128
        gap8 = np.array([
            sum((s >> j & 1) * (s >> ((j + 8) % length) & 1)
                * (1 - (s >> ((j + 3) % length) & 1))
                * (1 - (s >> ((j + 5) % length) & 1))
                for j in range(length)) for s in states], dtype=np.int64)
        diagonal = 256 * numbers - 256 * pairs[3] + 24 * gap8
        potential = 256 * pairs[5] + 256 * pairs[7] + 24 * gap8
    assert np.array_equal(diagonal, hopping_times_three * degree + potential)
    assert np.all(potential >= 0)
    h_times_three = diags(diagonal, dtype=np.int64) - hopping_times_three * adjacency
    components, labels = connected_components(adjacency, directed=False)
    zero_components = []
    for label in range(components):
        members = np.flatnonzero(labels == label)
        if np.all(potential[members] == 0):
            zero_components.append(label)
            assert all(len(set(words[i])) <= 1 for i in members)
    assert len(zero_components) == length + 1

    coarse_length = length // period
    expected = hopping_times_three * swap_laplacian(coarse_length)
    for residue in range(period):
        mapping = [index[sum(((s >> j) & 1) << (residue + period * j)
                             for j in range(coarse_length))]
                   for s in range(1 << coarse_length)]
        actual = h_times_three.tocsr()[mapping][:, mapping].toarray()
        assert np.array_equal(actual, expected)
    return dict(period=period, length=length, full_constrained_dimension=len(states),
                oriented_hops=len(rows), connected_components=components,
                exact_zero_components=len(zero_components),
                integer_laplacian_identity=True, ordered_residue_conservation=True,
                all_single_residue_swap_identities=True)


def main():
    start = time.perf_counter()
    results = [check(m, length) for m, lengths in ((4, (12, 16, 20)),
                                                  (3, (9, 12, 15, 18)))
               for length in lengths]
    output = dict(status="Exact integer identities for the stated finite effective spaces",
                  microscopic_derivation_checked=False, conformal_claim=False,
                  results=results, runtime_seconds=time.perf_counter() - start)
    Path(__file__).with_name("residue_laplacian_checks.json").write_text(
        json.dumps(output, indent=2) + "\n")
    for row in results:
        print(row, flush=True)


if __name__ == "__main__":
    main()
