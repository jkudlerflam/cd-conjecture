#!/usr/bin/env python3
"""Independent integer oracle tests of the finite-state digit tables.

Host integers occur only in fixture generation and expected-value calculation.
The implementation under test never receives those integers, only framed bits.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json
import operator
from pathlib import Path
import random
import time

from machine import ALPHABET, OPERATIONS, OUT, TABLE, Runner, run_binary, write_table

ORACLES = {
    "ADD": operator.add, "SUB": operator.sub, "MUL": operator.mul,
    "DIV": operator.floordiv, "MOD": operator.mod,
    "SHL": operator.lshift, "SHR": operator.rshift,
    "AND": operator.and_, "XOR": operator.xor,
    "EQ": operator.eq, "NE": operator.ne, "LT": operator.lt,
    "LE": operator.le, "GT": operator.gt, "GE": operator.ge,
}


def bits(value):
    return tuple(map(int, bin(abs(value))[:1:-1]))


def integer(result):
    value = sum(bit << index for index, bit in enumerate(result["bits_le"]))
    return -value if result["sign"] else value


def check_table():
    states = set(TABLE["active_tape"])
    halts = set(TABLE["halting_states"])
    assert len(TABLE["transitions"]) == len(states) * len(ALPHABET)
    assert set(TABLE["entries"]) == set(OPERATIONS)
    events = {}
    for event_id, ((state, read), (nxt, write, move)) in enumerate(TABLE["transitions"].items(), 1):
        assert state in states and read in ALPHABET
        assert nxt in states | halts and write in ALPHABET and move in (-1, 0, 1)
        assert 0 <= TABLE["active_tape"][state] < TABLE["tape_count"]
        # Adding a finite transition ID makes every local write injective.
        # This tests the record embedding, not a spatial history-tape compiler.
        target = (nxt, write, event_id)
        assert target not in events
        events[target] = (state, read, move)
        # Check the actual interpreter's inverse for every possible local rule,
        # including malformed-symbol rejection. Other tapes are spectators.
        runner = Runner("ADD", (0,), (0,), journal=True)
        c = runner.config
        c.state = state
        tape = TABLE["active_tape"][state]
        c.heads[tape] = 0
        c.tapes[tape] = {} if read == "_" else {0: read}
        before = c.snapshot()
        runner.step()
        assert c.state == nxt and c.heads[tape] == move
        runner.undo()
        assert c.snapshot() == before and not runner.history
    return len(events)


def check_case(op, left, right, stats, journal=False):
    runner = Runner(op, bits(left), bits(right), left < 0, right < 0, journal=journal)
    initial = runner.config.snapshot()
    input_tapes = [t.copy() for t in runner.config.tapes[:2]]
    result = runner.run()
    assert integer(result) == ORACLES[op](left, right), (op, left, right, result)
    assert runner.config.tapes[:2] == input_tapes
    b = result["bits_le"]
    assert len(b) == 1 or b[-1] == 1
    assert b != (0,) or result["sign"] == 0
    end = len(b) + 1
    assert set(runner.config.tapes[OUT]) == set(range(end + 1))
    assert all(position >= 0 for tape in runner.config.tapes for position in tape)
    stats["cases"] += 1
    stats["microsteps"] += result["steps"]
    stats["max_microsteps"] = max(stats["max_microsteps"], result["steps"])
    stats["max_input_bits"] = max(stats["max_input_bits"], len(bits(left)), len(bits(right)))
    if journal:
        count = len(runner.history)
        assert count == result["steps"]
        while runner.history:
            runner.undo()
        assert runner.config.snapshot() == initial
        stats["rewound_cases"] += 1
        stats["rewound_steps"] += count


def main():
    start = time.monotonic()
    local_rules = check_table()
    by_operation = {}
    rng = random.Random(731114)
    for op in OPERATIONS:
        stats = Counter()
        # Exhaust all signed magnitudes through 16. Shift counts are bounded
        # because their numerical value controls the output length.
        rights = range(13) if op in ("SHL", "SHR") else range(-16, 17)
        for left in range(-16, 17):
            for right in rights:
                if op in ("DIV", "MOD") and right == 0:
                    continue
                check_case(op, left, right, stats,
                           journal=abs(left) <= 4 and abs(right) <= 4)
        # Variable-width and carry/borrow boundaries, including long negative
        # infinite-two's-complement tails, use an independent host oracle.
        for width in (1, 2, 4, 8, 16, 32, 64, 128, 256):
            edge = (1 << width) - 1
            pairs = [(edge, 1), (-edge, 1), (-edge, -1), (1 << width, -edge)]
            pairs += [(rng.randrange(-edge, edge + 1), rng.randrange(-edge, edge + 1)) for _ in range(2)]
            for left, right in pairs:
                if op in ("DIV", "MOD") and right == 0:
                    right = 1
                if op in ("SHL", "SHR"):
                    right = abs(right) % 35
                check_case(op, left, right, stats, journal=True)
        by_operation[op] = dict(stats)

    error_cases = 0
    for op, left, right, expected in (
        ("DIV", 23, 0, "ERROR_DIV_ZERO"), ("MOD", -23, 0, "ERROR_DIV_ZERO"),
        ("SHL", 7, -1, "ERROR_NEGATIVE_SHIFT"), ("SHR", -7, -3, "ERROR_NEGATIVE_SHIFT"),
    ):
        runner = Runner(op, bits(left), bits(right), left < 0, right < 0, journal=True)
        initial = runner.config.snapshot()
        try:
            runner.run()
        except ArithmeticError as exc:
            assert str(exc) == expected
        else:
            raise AssertionError("Expected arithmetic exception")
        while runner.history:
            runner.undo()
        assert runner.config.snapshot() == initial
        error_cases += 1

    table_path = write_table()
    payload = json.loads(table_path.read_text())
    restored = {(row[0], row[1]): tuple(row[2:]) for row in payload["transitions"]}
    assert restored == TABLE["transitions"]
    summary = {
        "status": "pass", "scope": "Fixed-control digit TM; history records are not yet a local Hamiltonian",
        "operations": list(OPERATIONS), "states": len(TABLE["active_tape"]),
        "halting_states": len(TABLE["halting_states"]), "tapes": TABLE["tape_count"],
        "alphabet_size": len(ALPHABET), "transition_rules": len(TABLE["transitions"]),
        "every_local_rule_forward_inverse_checked": local_rules,
        "error_roundtrips": error_cases, "by_operation": by_operation,
        "cases": sum(s["cases"] for s in by_operation.values()),
        "microsteps": sum(s["microsteps"] for s in by_operation.values()),
        "rewound_cases": sum(s["rewound_cases"] for s in by_operation.values()),
        "rewound_steps": sum(s["rewound_steps"] for s in by_operation.values()),
        "table_sha256": hashlib.sha256(table_path.read_bytes()).hexdigest(),
        "elapsed_seconds": round(time.monotonic() - start, 3),
    }
    Path(__file__).with_name("results.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({key: value for key, value in summary.items() if key != "by_operation"}, indent=2))


if __name__ == "__main__":
    main()
