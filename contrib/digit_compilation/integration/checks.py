#!/usr/bin/env python3
"""Run real generator workloads with local digit kernels replacing word opcodes.

This is a hybrid compiler check. Stack, heap, branch dispatch, unsupported
arithmetic and operand placement remain word-machine operations. Supported
arithmetic is actually computed by the finite local rewrite table; no native
arithmetic result is substituted into those instructions.
"""
from __future__ import annotations

from collections import Counter
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import time


HERE = Path(__file__).resolve().parent
CONTRIB = HERE.parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result


vm = module("digit_integration_word_machine",
            CONTRIB / "autonomous_controller/generator/machine.py")
kernel = module("digit_integration_arithmetic",
                HERE.parent / "arithmetic/arithmetic.py")
sys.path.insert(0, str(CONTRIB / "controller_instantiation/circuit"))
from sparse_circuit import Graph, make_circuit, projective_plane


LOCAL_OPS = frozenset(("ADD", "SUB", "AND", "XOR", "EQ", "NE", "LT", "GT"))
ALL_ARITHMETIC = frozenset(vm.BINOPS.values())
ABSTRACT = ("", "X", "CX", "CCX", "H", "CH", "R", "QUERY", "CSWAPD", "P")


class DigitMachine(vm.Machine):
    """Preserve VM semantics and journal schema, replacing eight opcodes."""
    def __init__(self, *args, width=16, **kwargs):
        super().__init__(*args, **kwargs)
        self.width = width
        self.local_counts = Counter()
        self.remaining_arithmetic = Counter()

    def step(self):
        op = self.program["code"][self.pc][0]
        if op not in LOCAL_OPS:
            if op in ALL_ARITHMETIC:
                self.remaining_arithmetic[op] += 1
            return super().step()
        assert not self.halted
        old_pc = self.pc
        a, b = self.stack[-2:]
        # apply() initializes aligned digit tracks, follows the fixed table
        # through compute/copy/uncompute, checks cleanup, and decodes output.
        value = kernel.apply(op, a, b, self.width)
        del self.stack[-2:]
        self.stack.append(value)
        self.pc += 1
        if self.journal is not None:
            self.journal.append((old_pc, [a, b], 1, None))
        self.steps += 1
        self.op_counts[op] += 1
        self.local_counts[op] += 1
        self.max_stack = max(self.max_stack, len(self.stack))
        self.max_frames = max(self.max_frames, len(self.frames))
        self.max_integer_bits = max(self.max_integer_bits, abs(value).bit_length())


def check_case(program, entry, args, label, width=16, journal=True):
    actual = DigitMachine(program, entry, args, journal=journal, width=width)
    expected = vm.Machine(program, entry, args, journal=journal)
    initial = actual.state()
    start = time.monotonic()
    while not expected.halted:
        # The native VM is an independent result oracle, not the executed
        # arithmetic path used by actual to decide its later branches.
        assert (actual.pc, actual.stack) == (expected.pc, expected.stack)
        actual.step()
        expected.step()
    assert actual.state() == expected.state()
    assert actual.journal == expected.journal
    result = actual.stack[-1]
    if entry == "plane":
        u, v, edges_at = actual.array(result)[:3]
        edges = tuple(tuple(actual.array(x)) for x in actual.array(edges_at))
        reference = projective_plane(args[0])
        assert (u, v, edges) == (reference.u, reference.v, reference.edges)
        output_count = len(edges)
    else:
        wires, charged, ops_at = actual.array(result)
        ops = [actual.array(x) for x in actual.array(ops_at)]
        decoded = [(ABSTRACT[x[0]], *x[1:]) for x in ops]
        graph = projective_plane(args[0]) if entry == "main" else Graph(1, 1, ((0, 0),))
        reference = make_circuit(graph)
        assert (wires, charged) == (reference.wires, reference.data_sites)
        assert decoded == list(reference.abstract_ops())
        output_count = len(decoded)
    stats = actual.statistics()
    if journal:
        while actual.journal:
            actual.undo()
        assert actual.state() == initial
    return dict(case=label, signed_width=actual.width, output_count=output_count,
                word_steps=stats["steps"],
                local_arithmetic_calls=sum(actual.local_counts.values()),
                local_opcode_counts=dict(sorted(actual.local_counts.items())),
                remaining_host_arithmetic=dict(sorted(actual.remaining_arithmetic.items())),
                exact_stepwise_stack_match=True, exact_complete_state_match=True,
                exact_reference_output_match=True,
                exact_journal_and_rollback=True if journal else "not retained in this larger run",
                runtime_seconds=round(time.monotonic()-start, 3))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full", action="store_true", help="Also execute the entire p=2 abstract circuit with 128-bit local arithmetic")
    args = parser.parse_args()
    path = CONTRIB / "autonomous_controller/generator/bytecode.json"
    program = json.loads(path.read_text())
    assert program == vm.compile_program()
    results = [check_case(program, "plane", [2], "p=2 incidence graph"),
               check_case(program, "plane", [3], "p=3 incidence graph"),
               check_case(program, "tiny", [1], "complete one-edge abstract circuit")]
    if args.full:
        results.append(check_case(program, "main", [2, 1], "complete p=2 21-mode abstract circuit",
                                  width=128, journal=False))
    report = dict(scope="Actual hybrid generator execution using finite local digit kernels; memory and dispatch are still abstract",
                  all_passed=True, local_opcode_families=sorted(LOCAL_OPS),
                  total_local_arithmetic_calls=sum(r["local_arithmetic_calls"] for r in results),
                  program_sha256=vm.program_manifest(program)["sha256"], checks=results)
    (HERE / ("full_results.json" if args.full else "results.json")).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
