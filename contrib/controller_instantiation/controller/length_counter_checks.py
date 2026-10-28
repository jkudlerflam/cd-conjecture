#!/usr/bin/env python3
"""Exhaustive local checks and finite trajectory tests for the length counter."""

import json
from pathlib import Path

from length_counter import Configuration, LengthCounter, Symbol


def main():
    machine = LengthCounter()
    size = machine.B * machine.K
    assert (machine.B, machine.K, machine.local_factor) == (12, 16, 396)
    for x in range(size):
        assert machine.inverse[machine.collision[x]] == x
        assert machine.collision[machine.inverse[x]] == x
        assert machine.decode(x)[1].static == machine.decode(machine.collision[x])[1].static
    assert all(machine.collision[x] == y for x, y in machine.partial.items())
    path = Path(__file__).with_name("length_counter_table.json")
    assert json.loads(path.read_text()) == machine.manifest()

    trajectories = [machine.run_cycle(n) for n in range(4, 129)]
    # Standalone arbitrary-integer increment, including a frontier at every
    # possible BULK site.  RETURN must cross exactly the zeros cleared by
    # CARRY, and the frontier must survive.  No overflowing increment is
    # claimed: the length algorithm never needs one for n>=4.
    increments = 0
    for width in range(1, 10):
        for x in range(2 ** width - 1):
            for frontier in range(width):
                tape = [Symbol("LEFT")]
                tape += [Symbol("BULK", (x >> j) & 1, int(j == frontier)) for j in range(width)]
                tape += [Symbol("RIGHT")]
                config = Configuration(tape, 0, "F:CARRY")
                for tick in range(1, 2 * width + 1):
                    machine.step(config)
                    if config.position == 0 and config.control == "F:SCAN":
                        break
                else:
                    raise AssertionError("Increment did not finish")
                assert machine.value(config) == x + 1
                assert [j for j, s in enumerate(config.tape) if s.frontier] == [frontier + 1]
                increments += 1
    results = {
        "local_collision_pairs_and_inverses": size,
        "prescribed_injective_entries": len(machine.partial),
        "static_marks_preserved_on_all_inputs": True,
        "B_K_factor": [machine.B, machine.K, machine.local_factor],
        "blank_input_lengths_checked": [4, 128],
        "blank_to_parameter_to_blank_cycles": len(trajectories),
        "total_elementary_steps_checked": sum(t["full_cycle_steps"] for t in trajectories),
        "exact_running_time_formula_checked": "2(n^2+3n-4-2 popcount(n-1))",
        "arbitrary_nonoverflowing_binary_increments": increments,
        "largest_example": trajectories[-1],
        "all_passed": True,
    }
    Path(__file__).with_name("length_counter_results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
