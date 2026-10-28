#!/usr/bin/env python3
"""Independent finite checks for the preloaded-program dispatcher."""

from __future__ import annotations

import itertools
import json
from pathlib import Path
import random

from dispatcher import Dispatcher, MOVES, Symbol, lower_gate_stream


def apply_pair(vector, width, left, matrix, adjoint=False):
    if adjoint:
        matrix = [[matrix[b][a].conjugate() for b in range(4)] for a in range(4)]
    out = [0j] * len(vector)
    for word, amp in enumerate(vector):
        right_bit = width - left - 2
        pair = (word >> right_bit) & 3
        erased = word & ~(3 << right_bit)
        for new_pair in range(4):
            out[erased | (new_pair << right_bit)] += matrix[new_pair][pair] * amp
    return out


def quantum_matrices():
    s = 2 ** -0.5
    a, b = (1 + 1j) / 2, (1 - 1j) / 2
    return {
        "X": [[0j, 0j, 1+0j, 0j], [0j, 0j, 0j, 1+0j],
              [1+0j, 0j, 0j, 0j], [0j, 1+0j, 0j, 0j]],
        "H": [[s+0j, 0j, s+0j, 0j], [0j, s+0j, 0j, s+0j],
              [s+0j, 0j, -s+0j, 0j], [0j, s+0j, 0j, -s+0j]],
        "CNOT": [[1+0j, 0j, 0j, 0j], [0j, 1+0j, 0j, 0j],
                 [0j, 0j, 0j, 1+0j], [0j, 0j, 1+0j, 0j]],
        "CV": [[1+0j, 0j, 0j, 0j], [0j, 1+0j, 0j, 0j],
               [0j, 0j, a, b], [0j, 0j, b, a]],
    }


def ordered_g_one_head(control, position, width):
    """Apply all bonds' G blocks explicitly, including reverse blocks.

    The event sign is +1 on the forward U block and -1 on its adjoint.
    Starting fresh must trigger exactly one forward block.  This test does
    not replace a quantum matrix check; the block form proves unitarity.
    """
    fresh = True
    events = []
    for left in range(width - 1):
        if control.direction == 1:
            if fresh and position == left:
                position, fresh = left + 1, False
                events.append((left, 1))
            elif not fresh and position == left + 1:
                position, fresh = left, True
                events.append((left, -1))
        else:
            if fresh and position == left + 1:
                position, fresh = left, False
                events.append((left, 1))
            elif not fresh and position == left:
                position, fresh = left + 1, True
                events.append((left, -1))
    return position, fresh, events


def is_legal(program, width):
    cursor = 0
    for command in program:
        if command == "MOVE_R":
            cursor += 1
        elif command == "MOVE_L":
            cursor -= 1
        elif cursor == width - 1:
            return False
        if not 0 <= cursor < width:
            return False
    return True


def main():
    results = {}
    dispatcher = Dispatcher()
    count = dispatcher.B * dispatcher.K
    assert count == 13200
    assert dispatcher.B == 88 and dispatcher.K == 150
    assert dispatcher.local_factor == 26488
    # Every malformed as well as valid local pair has exactly one image and
    # one preimage.  The fresh/updated F is therefore exactly an involution.
    assert len(set(dispatcher.collision)) == count
    for x in range(count):
        assert dispatcher.inverse[dispatcher.collision[x]] == x
        assert dispatcher.collision[dispatcher.inverse[x]] == x
        assert dispatcher.decode(x)[1].program == dispatcher.decode(dispatcher.collision[x])[1].program
    assert all(dispatcher.collision[s] == t for s, t in dispatcher.partial.items())
    results["all_collision_pairs_and_inverses"] = count
    results["prescribed_injective_collision_entries"] = len(dispatcher.partial)
    results["program_symbol_conserved_on_all_local_inputs"] = True
    results["B_K_and_factor"] = [dispatcher.B, dispatcher.K, dispatcher.local_factor]
    saved = json.loads(Path(__file__).with_name("transition_table.json").read_text())
    assert saved == dispatcher.manifest()
    results["saved_table_exactly_reproduced"] = True

    sweeps = 0
    for width in range(2, 9):
        for control in dispatcher.controls:
            for position in range(width):
                target = position + control.direction
                if not 0 <= target < width:
                    continue
                dest, fresh, events = ordered_g_one_head(control, position, width)
                assert dest == target and not fresh
                assert events == [(min(position, target), 1)]
                sweeps += 1
    results["ordered_sweep_one_move_cases"] = sweeps

    # Exhaust all short programs for a nontrivial four-oracle menu.  The
    # validator is independent of the dispatcher: it tracks only the data
    # cursor and bounds.  Every accepted program must run entirely through
    # prescribed rules, reproduce its oracle order, and restore all neutral
    # records and the head after inverse execution.
    small = Dispatcher(("X", "H", "CNOT", "CV"))
    exhaustive = 0
    maximum_steps = 0
    for width in (2, 3, 4):
        for length in range(5):
            for program in itertools.product(small.commands, repeat=length):
                if not is_legal(program, width):
                    continue
                run = small.run_cycle(program, width)
                offset = len(program) + 2
                cursor = 0
                expected = []
                for command in program:
                    if command == "MOVE_R":
                        cursor += 1
                    elif command == "MOVE_L":
                        cursor -= 1
                    else:
                        expected.append((command, offset + cursor, False))
                expected += [(op, at, True) for op, at, _ in reversed(expected)]
                assert run["events"] == expected
                assert run["final"].snapshot() == small.make_configuration(program, width).snapshot()
                maximum_steps = max(maximum_steps, run["steps"])
                exhaustive += 1
    results["exhaustive_short_legal_program_cycles"] = exhaustive
    results["maximum_short_program_steps"] = maximum_steps

    # Random addressed scripts stress both routing directions and labels.
    rng = random.Random(17242)
    for trial in range(100):
        width = rng.randrange(2, 12)
        stream = [(rng.choice(dispatcher.oracles), rng.randrange(width - 1))
                  for _ in range(rng.randrange(1, 25))]
        program = lower_gate_stream(stream, width)
        run = dispatcher.run_cycle(program, width)
        offset = len(program) + 2
        expected = [(op, offset + left, False) for op, left in stream]
        expected += [(op, left, True) for op, left, _ in reversed(expected)]
        assert run["events"] == expected
    results["random_full_menu_addressed_scripts"] = 100

    # Actually evolve coherent quantum states, including a non-self-adjoint
    # gate.  Compare the forward half to direct circuit evaluation and the
    # full cycle to identity.  Repeat with distinct inputs and compare the
    # entire classical trajectory, not merely the instruction list.
    matrices = quantum_matrices()
    width = 4
    stream = [("H", 0), ("CNOT", 0), ("CV", 2), ("H", 1),
              ("X", 2), ("CV", 0), ("CNOT", 1), ("H", 2)]
    program = lower_gate_stream(stream, width, small.oracles)
    offset = len(program) + 2
    previous_trajectory = None
    max_error = 0.0
    for trial in range(20):
        initial = [complex(rng.uniform(-1, 1), rng.uniform(-1, 1))
                   for _ in range(2 ** width)]
        norm = sum(abs(x) ** 2 for x in initial) ** .5
        initial = [x / norm for x in initial]
        state = initial[:]
        direct = initial[:]
        for op, left in stream:
            direct = apply_pair(direct, width, left, matrices[op])
        gate_count = 0
        midpoint_state = None

        def callback(op, physical_left, adjoint):
            nonlocal state, gate_count, midpoint_state
            state = apply_pair(state, width, physical_left - offset, matrices[op], adjoint)
            gate_count += 1
            if gate_count == len(stream):
                midpoint_state = state[:]

        run = small.run_cycle(program, width, callback)
        if previous_trajectory is not None:
            assert run["snapshots"] == previous_trajectory
        previous_trajectory = run["snapshots"]
        err = max(abs(a - b) for a, b in zip(state, initial))
        err = max(err, max(abs(a - b) for a, b in zip(midpoint_state, direct)))
        assert err < 1e-12
        max_error = max(max_error, err)
    results["coherent_input_independent_trajectory_trials"] = 20
    results["max_coherent_forward_and_cleanup_error"] = max_error
    results["all_passed"] = True
    destination = Path(__file__).with_name("results.json")
    destination.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
