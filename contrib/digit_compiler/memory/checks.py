#!/usr/bin/env python3
"""Reference semantics appear here only; machine.py never indexes an array."""
from __future__ import annotations

from itertools import product
import json
from pathlib import Path
import random

from machine import ALPHABET, COPY_TABLE, TABLE, Machine, bits_frame, clean_run, memory_frame


def bits(n):
    return tuple(map(int, bin(n)[2:][::-1]))


def main():
    count = steps = max_steps = 0
    coverage = set()
    rollback_prefixes = 0

    def check(op, words, address=(0,), value=(0,)):
        nonlocal count, steps, max_steps
        result = clean_run(op, words, address, value)
        index = sum(bit << i for i, bit in enumerate(address))
        accepted = op in ("APPEND", "LEN") or index < len(words)
        assert result["accepted"] == accepted, (op, words, address, result)
        if accepted:
            if op == "GET":
                expected = bits_frame(words[index])
            elif op == "SET":
                changed = list(words)
                changed[index] = value
                expected = memory_frame(changed)
            elif op == "APPEND":
                expected = memory_frame(list(words) + [value])
            else:
                expected = bits_frame(bits(len(words)))
            assert result["result"] == expected, (op, words, address, result, expected)
        assert result["workspace_clean"]
        count += 1
        steps += result["forward_steps"] + result["copy_steps"] + result["forward_steps"]
        max_steps = max(max_steps, result["forward_steps"])

    # Exhaust all arrays up to four words from a three-word, variable-width set.
    values = [(0,), (1,), (1, 0)]
    for n in range(5):
        for words in product(values, repeat=n):
            for address in range(n + 3):
                for op in ("GET", "SET"):
                    for value in values if op == "SET" else [(0,)]:
                        check(op, words, bits(address), value)
            for value in values:
                check("APPEND", words, value=value)
            check("LEN", words)

    rng = random.Random(17242)
    for _ in range(120):
        words = [tuple(rng.randrange(2) for _ in range(rng.randrange(1, 65))) for _ in range(rng.randrange(33))]
        address = bits(rng.randrange(len(words) + 8)) + (0,) * rng.randrange(4)
        value = tuple(rng.randrange(2) for _ in range(rng.randrange(1, 130)))
        for op in ("GET", "SET", "APPEND", "LEN"):
            check(op, words, address, value)
    for width in (1, 2, 8, 32, 128, 513):
        for op in ("GET", "SET"):
            check(op, [(0,), (1,)], (1,) * width, (1,) * width)
            check(op, [(0,), (1,)], (0,) * width, (1,) * width)
    for n in (0, 1, 2, 3, 7, 8, 15, 16, 31, 32, 127, 128, 255, 256):
        check("LEN", [(0,)] * n)

    # Walk and reverse every prefix of selected operations, including failures.
    for op, words, address, value in [
        ("GET", values, (0,), (0,)), ("GET", values, (0, 1), (0,)),
        ("GET", values, (1, 1), (0,)), ("SET", values, (1,), (1, 1, 0, 1)),
        ("SET", values, (1, 1), (1,)), ("APPEND", [], (0,), (0, 1)),
        ("APPEND", values, (0,), (1, 0)), ("LEN", [(0,)] * 16, (0,), (0,)),
    ]:
        m = Machine(op, memory_frame(words), bits_frame(address), bits_frame(value))
        initial = m.snapshot()
        snapshots = [initial]
        while m.state != "ACCEPT":
            tape = m.tapes[TABLE.active[m.state]]
            key = m.state, tape.read()
            if key not in TABLE.rules:
                break
            coverage.add(key)
            assert m.step()
            snapshots.append(m.snapshot())
        for previous in reversed(snapshots[:-1]):
            assert m.reverse_step()
            assert m.snapshot() == previous
            rollback_prefixes += 1
        assert m.snapshot() == initial

    malformed = [
        ("GET", "^1#$", "^$", "^0$"),
        ("GET", "^1#$", "^1#0$", "^0$"),
        ("GET", "^#$", "^0$", "^0$"),
        ("GET", "^$", "^0$", "^0$"),
        ("SET", "^1#$", "^0$", "^$"),
        ("SET", "^1#$", "^0$", "^1#0$"),
        ("APPEND", "^1#0$", "^0$", "^1$"),
        ("APPEND", "^1##$", "^0$", "^1$"),
        ("LEN", "^1#0$", "^0$", "^1$"),
        ("LEN", "^#$", "^0$", "^1$"),
    ]
    for op, memory, address, value in malformed:
        m = Machine(op, memory, address, value)
        initial = m.snapshot()
        assert not m.run(), (op, memory, address, value)
        m.rollback()
        assert m.snapshot() == initial
    assert all(r.move in (-1, 0, 1) and r.read in ALPHABET and r.write in ALPHABET for r in TABLE.rules.values())
    # Determinism and one active tape per state are constructor invariants.
    assert json.loads(Path(__file__).with_name("transition_table.json").read_text()) == json.loads(json.dumps(TABLE.payload()))
    copy_payload = COPY_TABLE.payload()
    copy_payload["initial_state"] = "COPY"
    assert json.loads(Path(__file__).with_name("copy_table.json").read_text()) == json.loads(json.dumps(copy_payload))
    result = {
        "status": "pass", "semantic_cases": count,
        "base_forward_plus_reverse_and_copy_steps": steps,
        "maximum_base_forward_steps": max_steps,
        "every_prefix_exact_rollback_steps": rollback_prefixes,
        "malformed_input_rejections": len(malformed),
        "sampled_execution_rule_coverage": len(coverage),
        "base_control_states": len(TABLE.active), "base_rules": len(TABLE.rules),
        "base_alphabet_size": len(ALPHABET), "base_tapes": 5,
        "copy_control_states": len(COPY_TABLE.active), "copy_rules": len(COPY_TABLE.rules),
        "history_alphabet_size_including_blank": len(TABLE.rules) + 1,
        "maximum_tested_address_width": 513,
        "claims": ["serial digit access", "GET", "functional SET", "APPEND", "LEN", "exact initialized rollback", "clean scratch and history after output copy"],
        "not_claimed": ["heap allocator", "VM integration", "spatial nearest-neighbor single-head compiler", "full-space Hamiltonian legality", "c-d violation"],
    }
    Path(__file__).with_name("results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
