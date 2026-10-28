#!/usr/bin/env python3
"""Check the single composed table, including endpoint operators at every vertex."""

from __future__ import annotations

import json
from pathlib import Path
import random

from bridge import AutonomousBridge
from checks import apply_pair, matrices
from integrated_readout import IntegratedReadout


def micro_sweep(machine, config, callback=None):
    """Yield every actual G/F vertex, retaining fresh/updated head phases."""
    phase = 0
    moves = 0
    for left in range(len(config.tape) - 1):
        control = machine.control_by_name[config.control]
        if control.direction == 1:
            forward = phase == 0 and config.position == left
            backward = phase == 1 and config.position == left + 1
            target = left + 1 if forward else left
        else:
            forward = phase == 0 and config.position == left + 1
            backward = phase == 1 and config.position == left
            target = left if forward else left + 1
        if forward or backward:
            if control.oracle is not None and callback is not None:
                callback(control.oracle, left, control.adjoint ^ backward)
            config.position, phase = target, 1 - phase
            moves += 1
        yield "G", left, phase
    assert moves == 1 and phase == 1
    for site in reversed(range(len(config.tape))):
        if config.position == site:
            source = machine.encode(config.control, config.tape[site])
            assert phase == 1 and source in machine.partial
            table = machine.collision if phase else machine.inverse
            config.control, config.tape[site] = machine.decode(table[source])
            phase = 1 - phase
        yield "F", site, phase
    assert phase == 0


def active(config, endpoint):
    return bool(config.tape[endpoint].cursor and config.position != endpoint)


def inner(a, b):
    return sum(x.conjugate() * y for x, y in zip(a, b))


def raised(vector, n, site):
    mask = 1 << (n - 1 - site)
    result = [0j] * len(vector)
    for index, amplitude in enumerate(vector):
        if not index & mask:
            result[index | mask] = amplitude
    return result


def pullback(columns, n, site):
    operated = [raised(v, n, site) for v in columns]
    return [[inner(u, v) for v in operated] for u in columns]


def matrix_error(a, b):
    return max(abs(x - y) for row, other in zip(a, b) for x, y in zip(row, other))


def main():
    machine = IntegratedReadout()
    assert (machine.B, machine.K, machine.local_factor) == (80, 40, 6480)
    assert len(machine.partial) == 497
    assert json.loads(Path(__file__).with_name("integrated_transition_table.json").read_text()) == machine.manifest()
    for source in range(machine.B * machine.K):
        target = machine.collision[source]
        assert machine.inverse[target] == source
        assert machine.collision[machine.inverse[source]] == source
        assert machine.decode(source)[1].static == machine.decode(target)[1].static

    trajectories = [machine.run_cycle(n) for n in range(4, 129)]
    # Separately test every F input, including malformed dynamic records.
    size = machine.B * machine.K
    def f(label):
        if label < machine.B:
            return label
        phase, x = divmod(label - machine.B, size)
        return machine.B + (1 - phase) * size + (machine.inverse[x] if phase == 0 else machine.collision[x])
    for label in range(machine.local_factor):
        assert f(f(label)) == label

    # Coherent four-dimensional input code with no product-state or
    # computational-basis restriction.  Pure Python keeps this check free
    # of optional numerical dependencies.
    n, code_dimension = 6, 4
    rng = random.Random(2026102801)
    initial = []
    for _ in range(code_dimension):
        column = [complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(2**n)]
        for previous in initial:
            overlap = inner(previous, column)
            column = [x - overlap * y for x, y in zip(column, previous)]
        norm = inner(column, column).real ** .5
        initial.append([x / norm for x in column])
    gates = matrices()
    gates["SWAP"] = [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]
    prepared = [v[:] for v in initial]
    for opcode, left, adjoint in AutonomousBridge.expected_forward_events(n):
        prepared = [apply_pair(v, n, left, gates[opcode], adjoint) for v in prepared]
    targets = [pullback(prepared, n, site) for site in (1, n - 2)]
    columns = [v[:] for v in initial]
    def callback(opcode, left, adjoint):
        nonlocal columns
        columns = [apply_pair(v, n, left, gates[opcode], adjoint) for v in columns]

    config, reference = machine.initialize(n), machine.initialize(n)
    active_vertices = [0, 0]
    sums = [[[0j] * code_dimension for _ in range(code_dimension)] for _ in range(2)]
    unflagged = [[[0j] * code_dimension for _ in range(code_dimension)] for _ in range(2)]
    errors = [0.0, 0.0]
    total_vertices = 0
    simultaneous_active = 0
    for _ in range(machine.cycle_steps(n)):
        machine.step(reference)
        for operation, site, phase in micro_sweep(machine, config, callback):
            total_vertices += 1
            simultaneous_active += int(active(config, 0) and active(config, n - 1))
            for which, endpoint in enumerate((0, n - 1)):
                value = pullback(columns, n, endpoint)
                if config.position != endpoint:
                    for j in range(code_dimension):
                        for k in range(code_dimension):
                            unflagged[which][j][k] += value[j][k]
                if active(config, endpoint):
                    errors[which] = max(errors[which], matrix_error(value, targets[which]))
                    active_vertices[which] += 1
                    for j in range(code_dimension):
                        for k in range(code_dimension):
                            sums[which][j][k] += value[j][k]
        assert config.snapshot() == reference.snapshot()
    assert config.snapshot() == machine.initialize(n).snapshot()
    assert max(errors) < 1e-13 and min(active_vertices) > 0
    assert simultaneous_active > 0
    cleanup_error = max(abs(x-y) for v, w in zip(columns, initial) for x, y in zip(v, w))
    assert cleanup_error < 1e-13
    weights = [count / total_vertices for count in active_vertices]
    weighted_errors, negative_errors = [], []
    for j in range(2):
        actual = [[x / total_vertices for x in row] for row in sums[j]]
        expected = [[weights[j] * x for x in row] for row in targets[j]]
        without_ready = [[x / total_vertices for x in row] for row in unflagged[j]]
        weighted_errors.append(matrix_error(actual, expected))
        negative_errors.append(matrix_error(without_ready, expected))
    assert max(weighted_errors) < 1e-13
    assert min(negative_errors) > 1e-3

    # Check endpoint READY safety by tracking record identities at every
    # local gate, for a range of lengths.  Record identities refer to the
    # prepared state at the start of readout; computational gates don't
    # permute them, and no READY is active during those phases.
    safe_vertices = 0
    for width in range(4, 25):
        config = machine.initialize(width)
        records = list(range(width))
        def swap_record(opcode, left, adjoint):
            if opcode == "SWAP":
                records[left], records[left+1] = records[left+1], records[left]
        for _ in range(machine.cycle_steps(width)):
            for operation, site, phase in micro_sweep(machine, config, swap_record):
                for endpoint, output in ((0, 1), (width-1, width-2)):
                    if active(config, endpoint):
                        assert records[endpoint] == output
                        assert config.control.startswith("RO:")
                safe_vertices += 1
        assert records == list(range(width))
        assert config.snapshot() == machine.initialize(width).snapshot()

    result = {
        "all_passed": True,
        "B_K_factor_qubit_dimension": [machine.B, machine.K, machine.local_factor, 2 * machine.local_factor],
        "complete_collision_pairs_and_inverses": size,
        "prescribed_injective_entries": len(machine.partial),
        "full_classical_F_involution_basis_checked": machine.local_factor,
        "blank_input_lengths_checked": [4, 128],
        "complete_cycles": len(trajectories),
        "total_elementary_steps": sum(run["cycle_steps"] for run in trajectories),
        "readout_safety_lengths_checked": [4, 24],
        "readout_safety_local_G_F_vertices": safe_vertices,
        "coherent_operator_test": {
            "physical_sites": n, "arbitrary_coherent_code_dimension": code_dimension,
            "complete_cycle_vertices": total_vertices,
            "active_endpoint_vertices": active_vertices,
            "simultaneously_active_vertices": simultaneous_active,
            "maximum_active_vertex_operator_errors": errors,
            "uniform_vertex_sample_weights": weights,
            "weighted_endpoint_operator_errors": weighted_errors,
            "missing_READY_negative_control_errors": negative_errors,
            "final_quantum_cleanup_error": cleanup_error,
        },
        "weight_scope": "Uniform samples after every listed G/F gate; actual initialized-clock weights still require that clock",
        "limitations": [
            "Static output-mark grammar and initial fresh head remain assumed",
            "Binary-digit circuit demonstration, not full sparse-encoding generator",
            "No full Hamiltonian legal-sector or infrared theorem established by these tests",
        ],
    }
    Path(__file__).with_name("integrated_results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
