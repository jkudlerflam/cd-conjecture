#!/usr/bin/env python3
"""Connect the concrete sparse circuit stream to the dispatcher alphabet.

Only the N=1 fixture is expanded.  Its machine has a long quadratic-time
classical travel schedule, so this checks and fingerprints the exact program
without pretending to have simulated every interpreter step.
"""

from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "circuit"))
from sparse_circuit import Graph, LOCAL_OPCODES, Localizer, make_circuit
from dispatcher import DEFAULT_ORACLES, Dispatcher, iter_lower_gate_stream


def main():
    assert tuple(LOCAL_OPCODES) == DEFAULT_ORACLES
    circuit = make_circuit(Graph(1, 1, ((0, 0),)))
    localizer = Localizer(circuit)
    counts = Counter()
    digest = sha256()
    cursor = 0
    gate_count = 0
    # Use two independent stream traversals to check address reconstruction
    # rather than trusting the lowerer to count its own emitted instructions.
    original = iter(localizer.stream())
    for command in iter_lower_gate_stream(localizer.stream(), localizer.size):
        counts[command] += 1
        digest.update((command + "\n").encode())
        if command == "MOVE_R":
            cursor += 1
        elif command == "MOVE_L":
            cursor -= 1
        else:
            opcode, address = next(original)
            assert command == opcode and cursor == address
            assert 0 <= cursor < localizer.size - 1
            gate_count += 1
    assert next(original, None) is None
    m = sum(counts.values())
    n = m + 2 + localizer.size
    dispatcher = Dispatcher()
    results = {
        "scope": "Exact finite gate-stream/program interface; the full long dispatcher trajectory was not simulated",
        "fixture": "u=v=N=1, one-edge sparse isometry",
        "data_width": localizer.size,
        "quantum_gate_calls": gate_count,
        "program_length": m,
        "program_command_counts": dict(sorted(counts.items())),
        "program_sha256_newline_separated_utf8": digest.hexdigest(),
        "physical_tape_sites": n,
        "full_forward_inverse_step_upper_bound": 2 * (m + 1) * (2 * n + 6),
        "neutral_alphabet_B": dispatcher.B,
        "direction_fixed_control_count_K": dispatcher.K,
        "dispatcher_single_site_factor": dispatcher.local_factor,
        "all_reconstructed_operand_addresses_match": True,
        "all_passed": True,
    }
    destination = Path(__file__).with_name("integration_results.json")
    destination.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
