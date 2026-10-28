#!/usr/bin/env python3
"""Independent exhaustive small-instance audit of physical-head lowering."""
from itertools import product
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from tm import AsyncTable, Configuration, step
from single_head import Compiler, Runner, REJECT, Cell


def checks():
    macro_cases = local_cases = moves = 0
    manifests = []
    alphabet = ("_", "0", "1")
    for k in (1, 2, 3):
        # Include every active tape and direction in a single finite machine.
        active = {f"q{i}:{d}": i for i in range(k) for d in (-1, 0, 1)}
        transitions = {(q, a): ("halt", "0" if a == "1" else "1", int(q.split(":")[1]))
                       for q in active for a in alphabet}
        table = AsyncTable(alphabet, "_", k, active, transitions,
                           frozenset({"halt"}), "q0:0")
        table.validate()
        compiler = Compiler(table)
        physical_alphabet = compiler.alphabet()
        physical_alphabet_set = set(physical_alphabet)
        for control, cell in product(compiler.controls, physical_alphabet):
            target, write, direction = compiler.local_rule(control, cell)
            assert target in compiler.controls
            assert write in physical_alphabet_set
            assert direction in (-1, 0, 1)
            local_cases += 1
        manifests.append(compiler.manifest())
        for q, heads, scanned_symbol in product(active, product((-1, 0, 1), repeat=k), alphabet):
            tapes = [{p: str((p+i) % 2) for p in range(-2, 3)} for i in range(k)]
            i = active[q]
            if scanned_symbol == "_":
                tapes[i].pop(heads[i], None)
            else:
                tapes[i][heads[i]] = scanned_symbol
            reference = Configuration(q, tapes, heads)
            physical = Runner(compiler, reference, -2, 2)
            step(table, reference)
            count = physical.macrostep()
            assert physical.decode().snapshot() == reference.snapshot()
            assert count <= 2*len(physical.cells)+k+4
            macro_cases += 1
            moves += count

        # Missing the selected marker must reject after its scan. Duplicate
        # markers need not reject: the initialized-format promise is separate.
        reference = Configuration("q0:0", [{} for _ in range(k)], [0]*k)
        physical = Runner(compiler, reference, -2, 2)
        physical.cells = [Cell(c.tape, c.symbol, False) if isinstance(c, Cell) and c.tape == 0
                          else c for c in physical.cells]
        while physical.step():
            pass
        assert physical.control == REJECT
    # A concrete merge shows why deterministically routing a reversible
    # source step is not itself a reversible local simulator.
    identity = AsyncTable(("_", "0"), "_", 1, {"q": 0},
                          {("q", a): ("halt", a, 0) for a in ("_", "0")},
                          frozenset({"halt"}), "q")
    compiled = Compiler(identity)
    initial = Configuration("q", [{0: "0"}], [0])
    first = Runner(compiled, initial, 0, 1)
    second = Runner(compiled, initial, 0, 1)
    first.position, first.control = 1, ("SCAN", "q")
    second.position, second.control = 2, ("RETURN", "halt")
    assert (first.position, first.control) != (second.position, second.control)
    first.step()
    second.step()
    assert (first.position, first.control, first.cells) == (second.position, second.control, second.cells)
    return {"status": "pass", "source_tape_counts": [1, 2, 3],
            "complete_local_control_symbol_pairs": local_cases,
            "initialized_source_macrosteps_compared": macro_cases,
            "physical_microsteps_executed": moves,
            "missing_selected_head_rejections": 3,
            "explicit_noninjective_microstep_witness": True,
            "fixture_manifests": manifests,
            "scope": "Exact deterministic spatial simulation on initialized tagged tapes; no reversibility or general malformed-format theorem"}


if __name__ == "__main__":
    result = checks()
    Path(__file__).with_name("single_head_results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
