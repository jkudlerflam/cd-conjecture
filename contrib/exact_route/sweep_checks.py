"""Check the explicit fresh/updated single-head sweep interface."""

from itertools import product
from pathlib import Path
import json
import time

import numpy as np


def alphabet(d, k):
    return [(a,) for a in range(d)] + [
        (phase, direction, opcode, a)
        for phase in (0, 1)
        for direction in ("R", "L")
        for opcode in range(k)
        for a in range(d)
    ]


def local_column(left, right, gates):
    d = int(round(np.sqrt(gates[0].shape[0])))
    head_left, head_right = len(left) == 4, len(right) == 4
    mode = None
    if head_left and not head_right:
        phase, direction, opcode, a = left
        b = right[0]
        if (phase, direction) == (0, "R"):
            mode = "forward_R"
        elif (phase, direction) == (1, "L"):
            mode = "reverse_L"
    elif not head_left and head_right:
        a = left[0]
        phase, direction, opcode, b = right
        if (phase, direction) == (0, "L"):
            mode = "forward_L"
        elif (phase, direction) == (1, "R"):
            mode = "reverse_R"
    if mode is None:
        return {(left, right): 1.0}
    gate = gates[opcode]
    if mode.startswith("reverse"):
        gate = gate.conj().T
    output = {}
    for index in np.flatnonzero(gate[:, a * d + b]):
        aa, bb = divmod(int(index), d)
        if mode == "forward_R":
            pair = ((aa,), (1, "R", opcode, bb))
        elif mode == "forward_L":
            pair = ((1, "L", opcode, aa), (bb,))
        elif mode == "reverse_R":
            pair = ((0, "R", opcode, aa), (bb,))
        else:
            pair = ((aa,), (0, "L", opcode, bb))
        output[pair] = gate[index, a * d + b]
    return output


def apply_bond(state, bond, gates):
    out = {}
    for word, amplitude in state.items():
        for pair, coefficient in local_column(word[bond], word[bond + 1], gates).items():
            changed = word[:bond] + pair + word[bond + 2:]
            out[changed] = out.get(changed, 0) + amplitude * coefficient
    return {word: value for word, value in out.items() if abs(value) > 1e-15}


def phase_flip(word):
    return tuple((1 - x[0], *x[1:]) if len(x) == 4 else x for x in word)


def vector_error(a, b):
    return max((abs(a.get(word, 0) - b.get(word, 0)) for word in a.keys() | b.keys()), default=0)


def sweep_check(d, gates, lengths):
    count = 0
    worst = 0.0
    for n in lengths:
        for payload in product(range(d), repeat=n):
            for direction in ("R", "L"):
                for opcode, gate in enumerate(gates):
                    for position in range(n):
                        new_position = position + (1 if direction == "R" else -1)
                        if new_position not in range(n):
                            continue
                        word = tuple(
                            (0, direction, opcode, a) if j == position else (a,)
                            for j, a in enumerate(payload)
                        )
                        actual = {word: 1.0}
                        for bond in range(n - 1):
                            actual = apply_bond(actual, bond, gates)
                        actual = {phase_flip(w): z for w, z in actual.items()}
                        bond = min(position, new_position)
                        column = payload[bond] * d + payload[bond + 1]
                        expected = {}
                        for row in np.flatnonzero(gate[:, column]):
                            aa, bb = divmod(int(row), d)
                            target = list(payload)
                            target[bond:bond + 2] = aa, bb
                            target = tuple(
                                (0, direction, opcode, a) if j == new_position else (a,)
                                for j, a in enumerate(target)
                            )
                            expected[target] = gate[row, column]
                        worst = max(worst, vector_error(actual, expected))
                        count += 1
    assert worst < 1e-12
    return {"basis_inputs": count, "maximum_sweep_error": float(worst)}


def quantum_test():
    d = 2
    swap = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], complex)
    theta = 0.37
    gate = np.diag([np.exp(0.13j), 1, 1, np.exp(-0.21j)]).astype(complex)
    gate[1:3, 1:3] = [[np.cos(theta), 1j * np.sin(theta)], [1j * np.sin(theta), np.cos(theta)]]
    gates = [swap, gate]
    states = alphabet(d, len(gates))
    pairs = list(product(states, repeat=2))
    index = {pair: i for i, pair in enumerate(pairs)}
    matrix = np.zeros((len(pairs), len(pairs)), complex)
    charge_error = 0
    for pair, column in index.items():
        for target, amplitude in local_column(*pair, gates).items():
            matrix[index[target], column] = amplitude
            charge_error = max(charge_error, abs(sum(x[-1] for x in pair) - sum(x[-1] for x in target)))
    error = float(np.max(np.abs(matrix.conj().T @ matrix - np.eye(len(pairs)))))
    involution_error = float(np.max(np.abs(matrix @ matrix - np.eye(len(pairs)))))
    assert max(error, involution_error, charge_error) < 1e-12
    return {
        "local_alphabet": len(states), "two_site_dimension": len(pairs),
        "unitarity_error": error, "involution_error": involution_error,
        "species_conservation_error": charge_error,
        "sweep": sweep_check(d, gates, (2, 3, 4)),
    }


def loaded_exchange_test(q):
    # Ordinary D_s are 0,...,q-1; loaded T_s are q,...,2q-1.
    # Both have species charge e_s. T_s is a TAPE symbol, not head memory.
    d = 2 * q
    gate = np.eye(d * d)
    for s in range(q):
        for t in range(q):
            source = (q + s) * d + t
            target = (q + t) * d + s
            gate[source, source] = 0
            gate[target, source] = 1
    states = alphabet(d, 1)
    images = set()
    for pair in product(states, repeat=2):
        result = local_column(*pair, [gate])
        assert len(result) == 1
        target, amplitude = next(iter(result.items()))
        assert amplitude == 1
        assert sorted(x[-1] % q for x in pair) == sorted(x[-1] % q for x in target)
        assert local_column(*target, [gate]) == {pair: 1.0}
        images.add(target)
    assert len(images) == len(states) ** 2
    return {
        "species": q, "tape_alphabet": d, "head_alphabet": len(states),
        "two_site_basis_pairs": len(images), "exact_permutation_and_charge_checks": True,
        "sweep": sweep_check(d, [gate], (2, 3)),
    }


def main():
    start = time.perf_counter()
    output = {
        "quantum_gate_case": quantum_test(),
        "loaded_exchange_cases": [loaded_exchange_test(q) for q in (2, 3, 5)],
        "runtime_seconds": time.perf_counter() - start,
    }
    path = Path(__file__).with_name("sweep_verification.json")
    path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
