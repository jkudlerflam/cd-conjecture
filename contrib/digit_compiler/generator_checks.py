#!/usr/bin/env python3
"""Replay real frozen-generator arithmetic through the fixed digit tables.

This is an integration test, NOT a compiled VM. Python still controls the
word-machine's registers, heap, call stack, and opcode dispatch. Repeated
arithmetic input pairs are replayed once and counted separately as events.
"""
from collections import Counter
import importlib.util
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from alu.machine import run_binary, OPERATIONS, TABLE

spec = importlib.util.spec_from_file_location("frozen_word_vm", HERE.parent / "autonomous_controller/generator/machine.py")
vm = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = vm
spec.loader.exec_module(vm)


def bits(n):
    return tuple(int(x) for x in bin(abs(n))[2:][::-1])


def main():
    started = time.perf_counter()
    program = json.loads((HERE.parent / "autonomous_controller/generator/bytecode.json").read_text())
    assert program == vm.compile_program(), "Frozen program has changed"
    cache = {}
    events = Counter()
    digit_steps = 0
    runs = []
    for entry, arguments in [("plane", [p]) for p in (2, 3, 5)] + [("width", [n]) for n in (1, 7, 65)] + [("tiny", [1])]:
        machine = vm.Machine(program, entry, arguments)
        before_events = sum(events.values())
        while not machine.halted:
            op = program["code"][machine.pc][0]
            if op in OPERATIONS:
                a, b = machine.stack[-2:]
                events[op] += 1
                key = op, a, b
                if key not in cache:
                    result = run_binary(op, bits(a), bits(b), a < 0, b < 0)
                    magnitude = sum(bit << j for j, bit in enumerate(result["bits_le"]))
                    cache[key] = -magnitude if result["sign"] else magnitude
                    digit_steps += result["steps"]
                expected = cache[key]
                machine.step()
                assert machine.stack[-1] == expected, key
            else:
                machine.step()
        row = dict(entry=entry, arguments=arguments, word_steps=machine.steps,
                   arithmetic_events=sum(events.values()) - before_events)
        if entry == "plane":
            p, = arguments
            graph = machine.heap[machine.stack[-1]]
            points = p*p+p+1
            edges = machine.heap[graph[2]]
            assert graph[:2] == [points, points]
            assert len(edges) == (p+1)*points
            degree_left, degree_right = Counter(), Counter()
            for ptr in edges:
                a, b = machine.heap[ptr]
                degree_left[a] += 1
                degree_right[b] += 1
            assert set(degree_left.values()) == set(degree_right.values()) == {p+1}
            row.update(species=2*points, fermion_channels=len(edges))
        elif entry == "tiny":
            graph = machine.heap[machine.stack[-1]]
            row["abstract_gates"] = len(machine.heap[graph[2]])
            assert row["abstract_gates"] == 906
        runs.append(row)
    result = dict(status="PASS", frozen_program=vm.program_manifest(program), runs=runs,
        arithmetic_events=sum(events.values()), unique_operand_replays=len(cache),
        digit_steps_for_unique_replays=digit_steps, operations=dict(sorted(events.items())),
        table_states=len(TABLE["active_tape"]), table_rules=len(TABLE["transitions"]),
        runtime_seconds=round(time.perf_counter()-started, 3),
        scope="Host-selected actual generator arithmetic; VM control and memory remain host operations")
    (HERE / "generator_results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
