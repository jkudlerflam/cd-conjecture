#!/usr/bin/env python3
"""Replay actual frozen-generator memory events through the finite tape table.

The host VM still selects heap blocks and schedules operations. This is exact
operation-trace validation, not an autonomous spatial compiler or a heap VM.
"""
from __future__ import annotations

from collections import Counter
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
CONTRIB = HERE.parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


serial = load("serial_memory_trace_component", HERE / "machine.py")
vm = load("word_ram_trace_component", CONTRIB / "autonomous_controller/generator/machine.py")
circuit = load("sparse_circuit_trace_reference", CONTRIB / "controller_instantiation/circuit/sparse_circuit.py")


def unsigned_bits(value):
    assert value >= 0
    return tuple(map(int, bin(value)[2:][::-1]))


def signed_word(value):
    # Sign is an opaque first bit to memory; remaining bits are magnitude LE.
    return (int(value < 0),) + unsigned_bits(abs(value))


def expected_frame(op, output):
    if op == "GET":
        return serial.bits_frame(signed_word(output))
    if op == "LEN":
        return serial.bits_frame(unsigned_bits(output))
    return serial.memory_frame(tuple(signed_word(value) for value in output))


def run_case(program, entry, arguments, label):
    machine = vm.Machine(program, entry, arguments)
    events, replayed = Counter(), Counter()
    cache = {}
    checked_events = total_transitions = maximum_transitions = 0
    max_array_words = max_word_bits = max_address_bits = 0
    trace_digest = sha256()
    while not machine.halted:
        op = program["code"][machine.pc][0]
        if op not in ("GET", "SET", "APPEND", "LEN"):
            machine.step()
            continue
        pc = machine.pc
        events[op] += 1
        if op == "GET":
            pointer, index = machine.stack[-2:]
            value = 0
        elif op == "SET":
            pointer, index, value = machine.stack[-3:]
        elif op == "APPEND":
            pointer, value = machine.stack[-2:]
            index = 0
        else:
            pointer = machine.stack[-1]
            index = value = 0
        # This selection/serialization is the explicitly UNCOMPILED host adapter.
        words = tuple(machine.heap[pointer])
        key = op, words, index, value
        max_array_words = max(max_array_words, len(words))
        max_word_bits = max(max_word_bits, *(len(signed_word(x)) for x in words), len(signed_word(value)))
        max_address_bits = max(max_address_bits, len(unsigned_bits(index)))
        machine.step()
        expected = machine.stack[-1] if op in ("GET", "LEN") else tuple(machine.heap[pointer])
        if key not in cache:
            result = serial.clean_run(op, tuple(signed_word(x) for x in words), unsigned_bits(index), signed_word(value))
            assert result["accepted"] and result["workspace_clean"], (label, pc, key)
            assert result["result"] == expected_frame(op, expected), (label, pc, op, expected, result)
            cache[key] = expected
            transitions = 2 * result["forward_steps"] + result["copy_steps"]
            total_transitions += transitions
            maximum_transitions = max(maximum_transitions, transitions)
            replayed[op] += 1
        else:
            assert cache[key] == expected
        trace_digest.update((json.dumps([pc, op, words, index, value, expected], separators=(",", ":")) + "\n").encode())
        checked_events += 1
    result_pointer = machine.stack[-1]
    if entry == "plane":
        actual = machine.heap[result_pointer]
        edges = tuple(tuple(machine.heap[p]) for p in machine.heap[actual[2]])
        reference = circuit.projective_plane(arguments[0])
        assert (actual[0], actual[1], edges) == (reference.u, reference.v, reference.edges)
        output_check = dict(
            exact_projective_graph_and_edge_order=True, p=arguments[0],
            left_species=actual[0], right_species=actual[1], fermion_modes=len(edges),
            edge_sha256=sha256(json.dumps(edges, separators=(",", ":")).encode()).hexdigest(),
        )
    else:
        wire_count, charged_count, operations_pointer = machine.heap[result_pointer]
        op_names = ("", "X", "CX", "CCX", "H", "CH", "R", "QUERY", "CSWAPD", "P")
        actual = []
        for pointer in machine.heap[operations_pointer]:
            word = machine.heap[pointer]
            actual.append((op_names[word[0]], *word[1:]))
        reference = circuit.make_circuit(circuit.Graph(1, 1, ((0, 0),)))
        assert (wire_count, charged_count) == (reference.wires, reference.data_sites)
        assert actual == list(reference.abstract_ops())
        output_check = dict(
            exact_complete_abstract_circuit=True, abstract_gates=len(actual),
            abstract_sha256=sha256("".join(json.dumps(op, separators=(",", ":")) + "\n" for op in actual).encode()).hexdigest(),
        )
    return dict(
        case=label, frozen_vm_steps=machine.steps,
        checked_memory_events=checked_events, event_counts=dict(sorted(events.items())),
        unique_finite_table_replays=len(cache), replay_counts=dict(sorted(replayed.items())),
        finite_table_forward_reverse_copy_transitions=total_transitions,
        maximum_transitions_per_replay=maximum_transitions,
        maximum_array_words=max_array_words, maximum_signed_word_bits=max_word_bits,
        maximum_address_bits=max_address_bits, event_trace_sha256=trace_digest.hexdigest(),
        uncompiled_NEW_events=machine.op_counts["NEW"],
        every_replayed_operation_has_clean_workspace=True,
        **output_check,
    )


def main():
    frozen_path = CONTRIB / "autonomous_controller/generator/bytecode.json"
    program = json.loads(frozen_path.read_text())
    assert program == vm.compile_program()
    checks = []
    for entry, arguments, label in [("plane", [2], "plane_p2"), ("plane", [3], "plane_p3"), ("tiny", [1], "complete_one_edge_abstract_circuit")]:
        result = run_case(program, entry, arguments, label)
        checks.append(result)
        print(json.dumps(result), flush=True)
    payload = dict(
        status="pass", program_sha256=vm.program_manifest(program)["sha256"],
        total_checked_memory_events=sum(c["checked_memory_events"] for c in checks),
        total_unique_replays=sum(c["unique_finite_table_replays"] for c in checks),
        total_finite_table_transitions=sum(c["finite_table_forward_reverse_copy_transitions"] for c in checks),
        checks=checks,
        scope="Actual generator memory events replayed through fixed finite tape tables, with exact cleanup and final graph/circuit comparison",
        uncompiled=["VM instruction scheduling", "host heap-block selection", "NEW allocation", "alias-preserving heap-block replacement", "serialization/sign adapter", "spatial nearest-neighbor head routing"],
        deduplication="Within each case, identical operation, full input array, index and replacement values replay once; every occurrence is checked against the cached exact answer",
    )
    (HERE / "vm_results.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({k: v for k, v in payload.items() if k != "checks"}, indent=2))


if __name__ == "__main__":
    main()
