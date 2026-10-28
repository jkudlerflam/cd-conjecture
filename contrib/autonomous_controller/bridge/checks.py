#!/usr/bin/env python3
"""Independent checks of the composed local table and coherent legal orbits."""

from __future__ import annotations

import json
from pathlib import Path
import random

from bridge import AutonomousBridge, Configuration, Symbol


def apply_pair(vector, width, left, matrix, adjoint=False):
    out = [0j] * len(vector)
    if adjoint:
        matrix = [[matrix[b][a].conjugate() for b in range(4)] for a in range(4)]
    for word, amplitude in enumerate(vector):
        shift = width - left - 2
        pair = (word >> shift) & 3
        erased = word & ~(3 << shift)
        for new_pair in range(4):
            out[erased | (new_pair << shift)] += matrix[new_pair][pair] * amplitude
    return out


def matrices():
    a, b = (1 + 1j) / 2, (1 - 1j) / 2
    return {
        "CNOT": [[1+0j, 0j, 0j, 0j], [0j, 1+0j, 0j, 0j],
                 [0j, 0j, 0j, 1+0j], [0j, 0j, 1+0j, 0j]],
        "CV": [[1+0j, 0j, 0j, 0j], [0j, 1+0j, 0j, 0j],
               [0j, 0j, a, b], [0j, 0j, b, a]],
    }


def ordered_g(machine, control_name, position, width, callback=None):
    """Apply every bond's complete one-head G blocks, in spatial order.

    Fresh -> updated applies U_control; the reverse block applies its
    adjoint.  We execute both possibilities, so a second unintended move
    would be visible.  Neutral records remain attached to physical sites.
    """
    control = machine.control_by_name[control_name]
    fresh, moves = True, []
    for left in range(width - 1):
        if control.direction == 1:
            forward = fresh and position == left
            backward = not fresh and position == left + 1
            target = left + 1 if forward else left
        else:
            forward = fresh and position == left + 1
            backward = not fresh and position == left
            target = left if forward else left + 1
        if forward or backward:
            position, fresh = target, not fresh
            moves.append((left, 1 if forward else -1))
            if control.oracle is not None and callback is not None:
                callback(control.oracle, left, control.adjoint ^ backward)
    return position, fresh, moves


def sweep_step(machine, config, callback=None):
    """Actual spatial G sweep followed by F sweep, with no phase repairs."""
    position, fresh, moves = ordered_g(machine, config.control, config.position,
                                      len(config.tape), callback)
    assert len(moves) == 1 and moves[0][1] == 1
    assert not fresh  # Otherwise F would apply C^-1, an invalid legal step.
    config.position = position
    # F is the identity away from the head.  At an updated head it applies
    # C and turns the head fresh; no host-supplied control or tape change.
    encoded = machine.encode(config.control, config.tape[position])
    assert encoded in machine.partial
    config.control, config.tape[position] = machine.decode(machine.collision[encoded])


def main():
    machine = AutonomousBridge()
    size = machine.B * machine.K
    assert (machine.B, machine.K, machine.local_factor) == (48, 30, 2928)
    assert len(machine.partial) == 112
    for x in range(size):
        assert machine.inverse[machine.collision[x]] == x
        assert machine.collision[machine.inverse[x]] == x
        assert machine.decode(x)[1].static == machine.decode(machine.collision[x])[1].static
    assert json.loads(Path(__file__).with_name("transition_table.json").read_text()) == machine.manifest()

    # F has 48 no-head symbols and 2*1440 headed symbols.  Check the entire
    # classical basis, including malformed tapes and controls; the quantum
    # record is the identity in F.  Fresh (phase0) uses C^-1; updated uses C.
    def f(label):
        if label < machine.B:
            return label
        phase, x = divmod(label - machine.B, size)
        y = machine.inverse[x] if phase == 0 else machine.collision[x]
        return machine.B + (1 - phase) * size + y
    for label in range(machine.local_factor):
        assert f(f(label)) == label
    sweeps = 0
    for width in range(2, 15):
        for control in machine.controls:
            for position in range(width):
                target = position + control.direction
                if 0 <= target < width:
                    dest, fresh, moves = ordered_g(machine, control.name, position, width)
                    assert dest == target and not fresh
                    assert moves == [(min(position, target), 1)]
                    sweeps += 1

    trajectory_count, total_steps = 0, 0
    examples = []
    for n in range(4, 129):
        run = machine.run_cycle(n)
        t_counter = n * n + 3 * n - 4 - 2 * (n - 1).bit_count()
        t_program = t_counter + 2 * (n - 1)
        assert run["forward_steps"] == n * n + 8 * n - 11 - 2 * (n - 1).bit_count()
        for time, expected_control in ((t_counter, "F:COPY"), (t_program, "F:DISPATCH"),
                                       (run["forward_steps"], "R:DISPATCH")):
            position, control, tape = run["checkpoints"][time]
            assert control == expected_control
            assert position == (n - 1 if time == run["forward_steps"] else 0)
            assert sum(a.bit << j for j, a in enumerate(tape[1:-1])) == n - 1
            assert [j for j, a in enumerate(tape) if a.frontier] == [n - 1]
            assert all(a.cursor == 0 for a in tape)
            assert sum(a.program << j for j, a in enumerate(tape[1:-1])) == (0 if time == t_counter else n - 1)
        if n in (4, 5, 8, 16, 128):
            examples.append({"sites": n, "forward_steps": run["forward_steps"],
                             "cycle_steps": run["cycle_steps"],
                             "binary_program_little_endian": [(n - 1) >> j & 1 for j in range(n - 2)],
                             "forward_gates": len(run["events"]) // 2})
        trajectory_count += 1
        total_steps += run["cycle_steps"]

    # Full G/F sweeps for each time step, with coherent quantum data.  The
    # direct expected circuit is calculated from n, independently of the
    # generated tape.  Different inputs must yield identical classical
    # trajectories, and CV ensures the inverse must use actual adjoints.
    gates = matrices()
    rng = random.Random(240317242)
    coherent_trials, sweep_steps = 0, 0
    maximum_forward_error, maximum_return_error = 0.0, 0.0
    noncommuting_gate_error = 0.0
    for n in range(4, 11):
        baseline_trajectory = None
        for trial in range(4):
            initial = [complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) for _ in range(2 ** n)]
            norm = sum(abs(x) ** 2 for x in initial) ** .5
            initial = [x / norm for x in initial]
            vector = initial[:]
            expected = initial[:]
            for j in range(n - 2):
                opcode = "CV" if (n - 1) >> j & 1 else "CNOT"
                expected = apply_pair(expected, n, j + 1, gates[opcode])
            reference = machine.initialize(n)
            config = machine.initialize(n)
            trajectory = []
            def callback(opcode, position, adjoint):
                nonlocal vector
                vector = apply_pair(vector, n, position, gates[opcode], adjoint)
            for tick in range(1, 2 * machine.forward_steps(n) + 1):
                machine.step(reference)
                sweep_step(machine, config, callback)
                assert config.snapshot() == reference.snapshot()
                trajectory.append(config.snapshot())
                if tick == machine.forward_steps(n):
                    maximum_forward_error = max(maximum_forward_error,
                                                max(abs(x - y) for x, y in zip(vector, expected)))
                sweep_steps += 1
            maximum_return_error = max(maximum_return_error,
                                      max(abs(x - y) for x, y in zip(vector, initial)))
            assert config.snapshot() == machine.initialize(n).snapshot()
            if baseline_trajectory is None:
                baseline_trajectory = trajectory
            else:
                assert trajectory == baseline_trajectory
            coherent_trials += 1
        # A negative control: adjacent controlled gates cannot in general
        # be reordered.  Thus matching just a multiset would be inadequate.
        xy = apply_pair(apply_pair(initial, n, 1, gates["CV"]), n, 2, gates["CNOT"])
        yx = apply_pair(apply_pair(initial, n, 2, gates["CNOT"]), n, 1, gates["CV"])
        noncommuting_gate_error = max(noncommuting_gate_error, max(abs(x-y) for x,y in zip(xy,yx)))
    assert maximum_forward_error < 1e-13
    assert maximum_return_error < 1e-13
    assert noncommuting_gate_error > 1e-3

    result = {
        "all_passed": True,
        "B_K_factor_qubit_dimension": [machine.B, machine.K, machine.local_factor, 2 * machine.local_factor],
        "complete_collision_pairs_and_inverses": size,
        "prescribed_injective_entries": len(machine.partial),
        "static_marks_conserved_on_every_local_input": True,
        "full_classical_F_involution_basis_checked": machine.local_factor,
        "ordered_G_sweep_movement_cases": sweeps,
        "blank_input_lengths_checked": [4, 128],
        "complete_cycles": trajectory_count,
        "total_elementary_steps": total_steps,
        "coherent_random_input_trials": coherent_trials,
        "full_G_then_F_sweep_steps_with_coherent_data": sweep_steps,
        "maximum_forward_circuit_amplitude_error": maximum_forward_error,
        "maximum_inverse_cleanup_amplitude_error": maximum_return_error,
        "classical_trajectory_independent_of_quantum_input": True,
        "noncommuting_gate_negative_control_error": noncommuting_gate_error,
        "all_handoffs_are_collision_entries": True,
        "program_initially_blank": True,
        "examples": examples,
        "limitations": [
            "Generates a binary-digit circuit, not the sparse-encoding circuit",
            "Static boundary marks and the unique initial head remain an initialization contract",
            "No Hamiltonian legality penalties, READY readout, or port transport are instantiated here",
            "Dimension is for this controller plus one qubit, not the proposed counterexample",
        ],
    }
    Path(__file__).with_name("results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
