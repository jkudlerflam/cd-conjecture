#!/usr/bin/env python3
"""Exact finite-matrix checks for small-d closure and vacancy lifting.

Requires Python 3 and numpy. No central charge is fitted or inferred here.
All rings have at least three sites, and all bonds are oriented (j,j+1).
"""

from __future__ import annotations

import argparse
import itertools
import json
import platform
import time
from pathlib import Path

import numpy as np


def words(q, length):
    return list(itertools.product(range(q), repeat=length))


def ring_matrix(h, q, length):
    basis = words(q, length)
    lookup = {w: i for i, w in enumerate(basis)}
    out = np.zeros((len(basis), len(basis)), dtype=complex)
    for col, w in enumerate(basis):
        for j in range(length):
            k = (j + 1) % length
            for a, b in itertools.product(range(q), repeat=2):
                value = h[a * q + b, w[j] * q + w[k]]
                if value:
                    z = list(w)
                    z[j], z[k] = a, b
                    out[lookup[tuple(z)], col] += value
    return out


def seed_bond():
    eye = np.eye(2)
    x = np.array([[0, 1], [1, 0]], dtype=complex)
    y = np.array([[0, -1j], [1j, 0]])
    z = np.diag([1, -1])
    return (
        0.5 * (np.kron(x, x) - np.kron(y, y))
        + (np.kron(x, y) - np.kron(y, x)) / 160
        + 10 * np.kron(z, z)
        - (19.976047 / 2) * (np.kron(z, eye) + np.kron(eye, z))
    )


def normalize_bond(h):
    ev = np.linalg.eigvalsh(h)
    shift = float((ev[-1] + ev[0]) / 2)
    scale = float((ev[-1] - ev[0]) / 2)
    return (h - shift * np.eye(len(h))) / scale, shift, scale


def stacked_bond(ha, a, hb, b, gamma):
    """Local d=ab ordering is |alpha,beta>; this is not spatial blocking."""
    d = a * b
    out = np.zeros((d * d, d * d), dtype=complex)
    for ai, bi, aj, bj, ak, bk, al, bl in itertools.product(
        range(a), range(b), range(a), range(b),
        range(a), range(b), range(a), range(b)
    ):
        row = (ai * b + bi) * d + aj * b + bj
        col = (ak * b + bk) * d + al * b + bl
        if bi == bk and bj == bl:
            out[row, col] += ha[ai * a + aj, ak * a + al]
        if ai == ak and aj == al:
            out[row, col] += gamma * hb[bi * b + bj, bk * b + bl]
    return out


def lifted_bond(h, q, t=1.0, coupling=1.0, mu=0.0):
    """State 0 is vacant; states 1,...,q carry the original payload."""
    d = q + 1
    out = np.zeros((d * d, d * d), dtype=complex)
    for a in range(1, d):
        out[a, a * d] = out[a * d, a] = -t
    for a, b, c, e in itertools.product(range(q), repeat=4):
        out[(a + 1) * d + b + 1, (c + 1) * d + e + 1] += (
            coupling * h[a * q + b, c * q + e]
        )
    for a, b in itertools.product(range(d), repeat=2):
        out[a * d + b, a * d + b] -= mu * ((a != 0) + (b != 0)) / 2
    return out


def cyclic_zero_basis(q, number):
    basis = words(q, number)
    lookup = {w: i for i, w in enumerate(basis)}
    unused = set(basis)
    columns = []
    while unused:
        w = min(unused)
        orbit = {w[k:] + w[:k] for k in range(number)}
        unused.difference_update(orbit)
        col = np.zeros(len(basis))
        for z in orbit:
            col[lookup[z]] = 1 / np.sqrt(len(orbit))
        columns.append(col)
    return np.column_stack(columns)


def projection_check(h, q, length, number, name):
    """Retain the periodic boundary hop, including its color-word rotation."""
    positions = list(itertools.combinations(range(length), number))
    colors = words(q, number)
    basis = []
    for x in positions:
        for spin in colors:
            w = [0] * length
            for j, a in zip(x, spin):
                w[j] = a + 1
            basis.append(tuple(w))
    lookup = {w: i for i, w in enumerate(basis)}
    pos_lookup = {x: i for i, x in enumerate(positions)}
    dim = len(basis)
    kinetic = np.zeros((dim, dim), dtype=complex)
    interaction = np.zeros_like(kinetic)
    charge = np.zeros((len(positions), len(positions)))

    for col, x in enumerate(positions):
        occupied = set(x)
        for j in range(length):
            k = (j + 1) % length
            if (j in occupied) != (k in occupied):
                z = set(occupied)
                z.symmetric_difference_update((j, k))
                charge[pos_lookup[tuple(sorted(z))], col] -= 1
    for col, w in enumerate(basis):
        for j in range(length):
            k = (j + 1) % length
            if bool(w[j]) != bool(w[k]):
                z = list(w)
                z[j], z[k] = z[k], z[j]
                kinetic[lookup[tuple(z)], col] -= 1
            elif w[j] and w[k]:
                for a, b in itertools.product(range(q), repeat=2):
                    value = h[a * q + b, (w[j] - 1) * q + w[k] - 1]
                    if value:
                        z = list(w)
                        z[j], z[k] = a + 1, b + 1
                        interaction[lookup[tuple(z)], col] += value

    ev, vectors = np.linalg.eigh(charge)
    phi = vectors[:, 0]
    spin_zero = cyclic_zero_basis(q, number)
    isometry = np.kron(phi[:, None], spin_zero)
    squeezed = ring_matrix(h, q, number)
    contacts = np.array([
        sum((j in x) and ((j + 1) % length in x) for j in range(length))
        for x in positions
    ])
    g_ed = float(np.dot(abs(phi) ** 2, contacts) / number)
    rho = number / length
    g_formula = (length / number) * (
        rho ** 2 - (np.sin(np.pi * rho) / (length * np.sin(np.pi / length))) ** 2
    )
    projected = isometry.conj().T @ interaction @ isometry
    target = g_formula * spin_zero.conj().T @ squeezed @ spin_zero
    return {
        "name": name, "payload_dimension": q, "length": length,
        "particles": number, "full_sector_dimension": dim,
        "payload_zero_momentum_dimension": spin_zero.shape[1],
        "charge_energy": float(ev[0]),
        "charge_energy_formula_error": float(abs(
            ev[0] + 2 * np.sin(np.pi * rho) / np.sin(np.pi / length)
        )),
        "isometry_error": float(np.max(abs(isometry.conj().T @ isometry
                                            - np.eye(isometry.shape[1])))),
        "kinetic_eigenspace_error": float(np.max(abs(
            kinetic @ isometry - ev[0] * isometry
        ))),
        "contact_coefficient_ed": g_ed,
        "contact_coefficient_formula": float(g_formula),
        "contact_error": float(abs(g_ed - g_formula)),
        "first_order_projection_error": float(np.max(abs(projected - target))),
        "physical_interaction_hermiticity_error": float(np.max(abs(
            interaction - interaction.conj().T
        ))),
    }


def embedding_check(h, length):
    """Check the defect-gap bound in every extra-level occupancy sector."""
    q, d = 2, 3
    norm = float(np.max(abs(np.linalg.eigvalsh(h))))
    penalty = 2 * norm + 0.75
    bond = lifted_bond(h, q, t=0, coupling=1)
    # lifted_bond calls level 0 the vacancy; here it is a static penalized state.
    for a, b in itertools.product(range(d), repeat=2):
        bond[a * d + b, a * d + b] += penalty * ((a == 0) + (b == 0)) / 2
    physical = ring_matrix(bond, d, length)
    original_e0 = float(np.linalg.eigvalsh(ring_matrix(h, q, length))[0])
    basis = words(d, length)
    margins = []
    count = 0
    for mask in itertools.product((False, True), repeat=length):
        defects = sum(mask)
        if not defects:
            continue
        indices = [i for i, w in enumerate(basis)
                   if tuple(a == 0 for a in w) == mask]
        e0 = float(np.linalg.eigvalsh(physical[np.ix_(indices, indices)])[0])
        bound = original_e0 + (penalty - 2 * norm) * defects
        margins.append(e0 - bound)
        count += 1
    keep = [i for i, w in enumerate(basis) if all(w)]
    original = ring_matrix(h, q, length)
    return {
        "length": length, "penalty": penalty, "bond_norm": norm,
        "defect_patterns_checked": count,
        "smallest_margin_above_proved_bound": float(min(margins)),
        "payload_block_error": float(np.max(abs(
            physical[np.ix_(keep, keep)] - original
        ))),
    }


def entropy(vector, q, length, cut):
    if cut == 0 or cut == length:
        return 0.0
    weights = np.linalg.svd(vector.reshape(q ** cut, q ** (length - cut)),
                            compute_uv=False) ** 2
    weights = weights[weights > 1e-15]
    return float(-np.dot(weights, np.log(weights)))


def entropy_identity_check(q, length, number, cut, rng):
    """Any factorized charge / ordered-payload state satisfies this identity.

    Random complex states deliberately test the kinematics independently
    of Gaussian physics or a ground-state ansatz.
    """
    positions = list(itertools.combinations(range(length), number))
    phi = rng.normal(size=len(positions)) + 1j * rng.normal(size=len(positions))
    phi /= np.linalg.norm(phi)
    chi = rng.normal(size=q ** number) + 1j * rng.normal(size=q ** number)
    chi /= np.linalg.norm(chi)
    physical = np.zeros((q + 1) ** length, dtype=complex)
    charge = np.zeros(2 ** length, dtype=complex)
    probabilities = np.zeros(number + 1)
    for i, x in enumerate(positions):
        charge_index = sum(2 ** (length - 1 - j) for j in x)
        charge[charge_index] = phi[i]
        probabilities[sum(j < cut for j in x)] += abs(phi[i]) ** 2
        for ci, spin in enumerate(words(q, number)):
            w = [0] * length
            for j, a in zip(x, spin):
                w[j] = a + 1
            index = sum(a * (q + 1) ** (length - 1 - j) for j, a in enumerate(w))
            physical[index] = phi[i] * chi[ci]
    s_charge = entropy(charge, 2, length, cut)
    s_payload_average = float(sum(
        probabilities[k] * entropy(chi, q, number, k) for k in range(number + 1)
    ))
    s_physical = entropy(physical, q + 1, length, cut)
    return {
        "payload_dimension": q, "length": length, "particles": number,
        "cut": cut, "physical_entropy": s_physical, "charge_entropy": s_charge,
        "mean_payload_entropy": s_payload_average,
        "entropy_identity_error": abs(s_physical - s_charge - s_payload_average),
    }


def closure_envelope(limit, conditional_seed):
    # Only benchmark SU(d)_1 seeds, static embeddings and onsite stacking.
    best = [0] * (limit + 1)
    for d in range(2, limit + 1):
        best[d] = max(d - 1, best[d - 1], 2 if d == 2 and conditional_seed else 0)
        for a in range(2, d):
            for b in range(2, d // a + 1):
                best[d] = max(best[d], best[a] + best[b])
    return best[2:]


def serialize_matrix(h):
    return {"real": h.real.tolist(), "imaginary": h.imag.tolist()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("checks.json"))
    args = parser.parse_args()
    start = time.perf_counter()
    h, shift, scale = normalize_bond(seed_bond())
    gamma = 1.37
    h4 = stacked_bond(h, 2, h, 2, gamma)
    stack_results = []
    for length in (3, 4):
        e = np.linalg.eigvalsh(ring_matrix(h, 2, length))
        expected = np.sort((e[:, None] + gamma * e[None, :]).ravel())
        actual = np.linalg.eigvalsh(ring_matrix(h4, 4, length))
        stack_results.append({"length": length, "dimension": 4 ** length,
                              "spectrum_error": float(np.max(abs(actual - expected)))})

    rng = np.random.default_rng(20261116)
    projection_results = []
    for length, number in ((6, 3), (7, 3), (8, 4)):
        projection_results.append(projection_check(h, 2, length, number, "qubit_seed"))
    for q, length, number in ((2, 7, 4), (3, 6, 3), (3, 7, 3)):
        raw = rng.normal(size=(q * q, q * q)) + 1j * rng.normal(size=(q * q, q * q))
        random_bond, _, _ = normalize_bond(raw + raw.conj().T)
        projection_results.append(projection_check(
            random_bond, q, length, number, "generic_complex_hermitian"
        ))

    eta = 1e-6
    h3 = lifted_bond(h, 2, coupling=eta)
    n = np.diag([0, 1, 1])
    number_bond = np.kron(n, np.eye(3)) + np.kron(np.eye(3), n)
    results = {
        "description": "Finite algebra checks only; no thermodynamic c estimate.",
        "runtime": {"python": platform.python_version(), "numpy": np.__version__},
        "seed_shift": shift, "seed_scale": scale,
        "normalized_seed_bond": serialize_matrix(h),
        "qutrit_example": {
            "t": 1.0, "eta": eta, "mu": 0.0,
            "status": "explicit bond, not a verified c=3 point",
            "bond": serialize_matrix(h3),
            "hermiticity_error": float(np.max(abs(h3 - h3.conj().T))),
            "number_commutator_error": float(np.max(abs(h3 @ number_bond - number_bond @ h3))),
        },
        "stack_gamma": gamma, "stack_spectra": stack_results,
        "vacancy_first_order_projection": projection_results,
        "embedding": [embedding_check(h, length) for length in (3, 4, 5)],
        "factorized_state_entropy": [
            entropy_identity_check(q, length, number, cut, rng)
            for q, length, number, cut in ((2, 6, 3, 2), (2, 6, 3, 3),
                                           (2, 8, 4, 4), (3, 6, 3, 3))
        ],
        "benchmark_closure_d_2_through_32": closure_envelope(32, False),
        "conditional_c2_seed_closure_d_2_through_32": closure_envelope(32, True),
    }
    errors = [x["spectrum_error"] for x in stack_results]
    for row in projection_results:
        errors += [v for k, v in row.items() if k.endswith("error")]
    errors += [results["qutrit_example"]["hermiticity_error"],
               results["qutrit_example"]["number_commutator_error"]]
    for row in results["embedding"]:
        errors.append(row["payload_block_error"])
        assert row["smallest_margin_above_proved_bound"] >= -1e-10
    errors += [row["entropy_identity_error"] for row in results["factorized_state_entropy"]]
    assert max(errors) < 1e-10, max(errors)
    results["maximum_identity_error"] = float(max(errors))
    results["elapsed_seconds"] = time.perf_counter() - start
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps({"output": str(args.output), "maximum_identity_error": max(errors),
                      "elapsed_seconds": results["elapsed_seconds"],
                      "projection_cases": len(projection_results)}, indent=2))


if __name__ == "__main__":
    main()
