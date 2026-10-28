#!/usr/bin/env python3
"""Independent comparisons of a fixed bytecode against the existing circuit."""
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "controller_instantiation" / "circuit"))
sys.path.insert(0, str(HERE.parents[1] / "controller_instantiation" / "controller"))
from sparse_circuit import Graph, LOCAL_OPCODES, Localizer, make_circuit, projective_plane
from dispatcher import iter_lower_gate_stream
from machine import Machine, compile_program, program_manifest
from resource_analysis import analyze, evaluate

ABSTRACT = ("", "X", "CX", "CCX", "H", "CH", "R", "QUERY", "CSWAPD", "P")
COMMANDS = ("",) + LOCAL_OPCODES + ("MOVE_R", "MOVE_L")


def abstract_output(machine, result):
    wires, charged, pointer = machine.array(result)
    ops = [machine.array(p) for p in machine.array(pointer)]
    return wires, charged, [(ABSTRACT[op[0]], *op[1:]) for op in ops]


def check_abstract(program, p=None, journal=False):
    machine = Machine(program, "tiny" if p is None else "main", [1] if p is None else [p, 1], journal=journal)
    initial = machine.state()
    result = machine.run()
    wires, charged, actual = abstract_output(machine, result)
    graph = Graph(1, 1, ((0, 0),)) if p is None else projective_plane(p)
    reference = make_circuit(graph)
    expected = list(reference.abstract_ops())
    assert (wires, charged) == (reference.wires, reference.data_sites)
    assert actual == expected
    stats = machine.statistics()
    stats.update(case="one_edge" if p is None else f"p={p}", abstract_gates=len(actual),
                 exact_abstract_stream_match=True,
                 abstract_sha256=sha256("".join(json.dumps(op, separators=(",", ":")) + "\n" for op in actual).encode()).hexdigest())
    if journal:
        while machine.journal:
            machine.undo()
        assert machine.state() == initial
        stats["complete_forward_reverse_restores_input_and_empty_heap"] = True
    return stats


def check_expanded(program, mode):
    circuit = make_circuit(Graph(1, 1, ((0, 0),)))
    localizer = Localizer(circuit)
    if mode == 2:
        reference = iter(localizer.stream())
    else:
        reference = iter(iter_lower_gate_stream(localizer.stream(), localizer.size))
    digest = sha256()
    count = 0
    def sink(op, address):
        nonlocal count
        actual = (COMMANDS[op], address) if mode == 2 else COMMANDS[op]
        assert actual == next(reference), (count, actual)
        digest.update((json.dumps(actual, separators=(",", ":")) + "\n").encode())
        count += 1
    machine = Machine(program, "tiny", [mode], sink=sink)
    machine.run()
    assert next(reference, None) is None
    stats = machine.statistics()
    stats.update(case="one_edge_local" if mode == 2 else "one_edge_dispatcher",
                 exact_stream_match=True, output_count=count, output_sha256=digest.hexdigest())
    return stats


def check_graph_and_reversal(program):
    machine = Machine(program, "plane", [2], journal=True)
    initial = machine.state()
    result = machine.run()
    values = machine.array(result)
    edges = tuple(tuple(machine.array(p)) for p in machine.array(values[2]))
    reference = projective_plane(2)
    assert (values[0], values[1], edges) == (reference.u, reference.v, reference.edges)
    stats = machine.statistics()
    while machine.journal:
        machine.undo()
    assert machine.state() == initial
    stats.update(case="p=2_graph", exact_incidence_order_match=True, complete_reversal=True)
    return stats


def check_prefix_reversal(program):
    # Cross into the output-lowering stage and include actual EMIT operations.
    machine = Machine(program, "tiny", [3], journal=True)
    initial = machine.state()
    while machine.emissions < 500:
        machine.step()
    assert machine.output and not machine.halted
    stats = machine.statistics()
    while machine.journal:
        machine.undo()
    assert machine.state() == initial
    stats.update(case="one_edge_command_prefix", emitted_symbols=500,
                 arbitrary_stop_reversal_restores_input_and_empty_output=True)
    return stats


def main():
    start = time.monotonic()
    program = compile_program()
    frozen = json.loads((HERE / "bytecode.json").read_text())
    assert program == frozen
    results = [check_graph_and_reversal(program), check_abstract(program, journal=True),
               check_abstract(program, 2), check_abstract(program, 3),
               check_expanded(program, 2), check_expanded(program, 3), check_prefix_reversal(program)]
    bounds = analyze()
    assert bounds == json.loads((HERE / "resource_bounds.json").read_text())
    for result in results:
        p = 3 if result["case"] == "p=3" else 2
        graph = (projective_plane(p) if result["case"].startswith("p=")
                 else Graph(1, 1, ((0, 0),)))
        L = graph.n + 4 * graph.q
        w = max(1, graph.species_width, graph.label_width)
        A = L * L * w
        bound = evaluate(bounds["T_polynomial_coefficients_in_ascending_powers_of_B"], 128 * A)
        assert result["steps"] <= bound
        assert result["heap_words"] <= 100 * result["steps"] + 40000
        assert result["max_integer_bits"] <= 2 * A + bound.bit_length() + 64
        result["fits_source_derived_resource_bound"] = True
    payload = dict(scope="Exact fixed word-RAM generator and journal inverse, not digit-level local TM or initialized Hamiltonian",
                   program=program_manifest(program), checks=results,
                   resource_bound_degree=bounds["polynomial_degree"],
                   elapsed_seconds=round(time.monotonic()-start, 3), all_passed=True)
    (HERE / "results.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
