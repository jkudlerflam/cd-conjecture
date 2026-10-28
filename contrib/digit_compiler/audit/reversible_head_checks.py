"""Exact local injectivity, autonomous execution and complete cleanup tests."""
from itertools import product
import json
from pathlib import Path
import random
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from reversible_head import Source, Compiler, Runner, Work
from tm import AsyncTable, Configuration
from single_head import Compiler as SpatialCompiler, Runner as SpatialRunner
from alu.machine import Runner as ArithmeticRunner, machine_for


def table_checks():
    transitions = {(q, a): ("halt", "_", d)
                   for q, d in (("left", -1), ("stay", 0), ("right", 1))
                   for a in ("_", "0", "1")}
    compiler = Compiler(Source.from_dictionary(("_", "0", "1"), transitions, {"halt"}))
    controls, alphabet, rules = compiler.materialize()
    inverse, direction = {}, {}
    for (q, a), (r, b, d) in rules.items():
        assert r in controls and b in alphabet
        assert (r, b) not in inverse, "Merged incoming local transition"
        inverse[r, b] = q, a, d
        assert direction.setdefault(r, d) == d, "Target has inconsistent incoming direction"
        assert compiler.incoming_direction(r) == d
        assert compiler.inverse_rule(r, b) == (q, a, d)
    for r, b in product(controls, alphabet):
        assert compiler.inverse_rule(r, b) == inverse.get((r, b))
    rng = random.Random(39221)
    entries = list(rules.items())
    for _ in range(10000):
        (q, a), _ = rng.choice(entries)
        runner = Runner(compiler, "stay", ["0"], 0, 2)
        runner.cells = [rng.choice(alphabet) for _ in range(7)]
        runner.position = 3
        runner.control = q
        runner.cells[3] = a
        initial = runner.snapshot()
        assert runner.step() and runner.backward()
        assert runner.snapshot() == initial
    return dict(controls=len(controls), alphabet=len(alphabet), rules=len(rules),
                inverse_pairs_exhausted=len(controls)*len(alphabet),
                full_configuration_local_roundtrips=10000,
                fixed_incoming_direction=True, injective_target_written_pairs=True)


def destructive_runs():
    transitions = {("erase", a): ("erase", "_", 1) for a in "01"}
    transitions["erase", "#"] = "halt", "#", 0
    compiler = Compiler(Source.from_dictionary(("_", "0", "1", "#"), transitions, {"halt"}))
    cases = moves = 0
    for length in range(7):
        for word in product("01", repeat=length):
            runner = Runner(compiler, "erase", [*word, "#"], 0, length+2)
            initial = runner.snapshot()
            while runner.step():
                pass
            assert runner.control == ("F", "halt")
            assert runner.decode_work() == ["_"]*length+["#"]
            count = runner.steps
            for _ in range(count):
                assert runner.backward()
            assert runner.snapshot() == initial
            cases += 1
            moves += count
    return dict(all_binary_strings_through_length=6, complete_runs=cases,
                forward_microsteps=moves, inverse_microsteps=moves,
                all_work_and_history_restored=True)


def lowered_arithmetic():
    cases = []
    for op, a, b, expected in (("EQ", 0, 1, 0), ("ADD", 1, 1, 2)):
        arithmetic = ArithmeticRunner(op, (a,), (b,))
        table = AsyncTable.from_export(machine_for(op))
        spatial_compiler = SpatialCompiler(table)
        config = Configuration(arithmetic.config.state, arithmetic.config.tapes, arithmetic.config.heads)
        spatial = SpatialRunner(spatial_compiler, config, 0, 5)
        initial_work = list(spatial.cells)
        while spatial.step():
            pass
        source = Source.from_spatial_compiler(spatial_compiler)
        compiler = Compiler(source)
        runner = Runner(compiler, ("SCAN", config.state), initial_work, 0, spatial.steps+1)
        initial = runner.snapshot()
        completed = 0
        while runner.control != ("F", spatial.control):
            runner.source_step()
            completed += 1
            assert completed <= spatial.steps
        assert completed == spatial.steps
        assert runner.decode_work() == spatial.cells
        assert runner.position == spatial.position+1
        result = spatial.decode()
        out = result.tapes[2]
        value, j = 0, 1
        while out[j] != "R":
            value += int(out[j]) << (j-1)
            j += 1
        if out[0] == "M":
            value = -value
        assert value == expected
        count = runner.steps
        for _ in range(count):
            assert runner.backward()
        assert runner.snapshot() == initial
        cases.append(dict(operation=op, operands=[a, b], answer=value,
                          deterministic_source_microsteps=spatial.steps,
                          reversible_physical_microsteps=count,
                          inverse_physical_microsteps=count,
                          exact_full_cleanup=True, manifest=compiler.manifest()))
    return cases


if __name__ == "__main__":
    start = time.monotonic()
    result = dict(local_table=table_checks(), destructive_source=destructive_runs(),
                  actual_lowered_alu=lowered_arithmetic(), status="pass",
                  runtime_seconds=round(time.monotonic()-start, 3))
    Path(__file__).with_name("reversible_head_results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))
