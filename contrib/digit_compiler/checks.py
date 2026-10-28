#!/usr/bin/env python3
"""Independent source/table and one-physical-head simulation checks."""
import json
from pathlib import Path
import random
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from tm import AsyncTable, Configuration, step
from single_head import Compiler, Runner, REJECT
from alu.machine import Runner as AluRunner, machine_for, OPERATIONS


def bits(n):
    return tuple(int(x) for x in bin(abs(n))[2:][::-1])


def compare(table, initial, lo, hi, source_limit=100000):
    source = Configuration(initial.state, initial.tapes, initial.heads)
    physical = Runner(Compiler(table), source, lo, hi)
    assert physical.decode().snapshot() == source.snapshot()
    while source.state not in table.halting_states and source.steps < source_limit:
        step(table, source)
        physical.macrostep()
        assert physical.decode().snapshot() == source.snapshot()
    assert source.state in table.halting_states or source.steps == source_limit
    return source.steps, physical.steps


def main():
    started = time.perf_counter()
    total_source = total_physical = cases = 0
    manifests = {}
    # Every arithmetic opcode, both signs and unequal word lengths. The ALU
    # suite separately supplies broader semantic tests; here each individual
    # source step is compared to the spatially lowered machine's full state.
    for op in OPERATIONS:
        for a, b in ((5, 2), (-3, 1), (0, -2)):
            original = AluRunner(op, bits(a), bits(b), a < 0, b < 0).config
            table = AsyncTable.from_export(machine_for(op))
            counts = compare(table, original, -1, 14)
            total_source += counts[0]
            total_physical += counts[1]
            cases += 1
    manifests["arithmetic"] = Compiler(table).manifest()

    data = json.loads((HERE / "memory/transition_table.json").read_text())
    names = data["tapes"]
    for op in ("GET", "SET", "APPEND", "LEN"):
        for address in ("^0$", "^1$", "^01$"):
            data["initial_state"] = op
            table = AsyncTable.from_export(data)
            inputs = {"M": "^1#001#01#$", "A": address, "V": "^1101$"}
            original = Configuration(op, [dict(enumerate(inputs.get(name, ""))) for name in names],
                                     [0] * len(names))
            counts = compare(table, original, -1, 22)
            total_source += counts[0]
            total_physical += counts[1]
            cases += 1
    manifests["memory"] = Compiler(table).manifest()

    # Independently generated arbitrary tables, not arithmetic-specific
    # examples. Source heads may visit negative coordinates.
    rng = random.Random(20261028)
    random_cases = 0
    for k in (1, 2, 4):
        for _ in range(8):
            active = {str(q): rng.randrange(k) for q in range(4)}
            rules = {(q, a): (rng.choice(list(active) + ["HALT"]), rng.choice(("_", "0", "1")),
                              rng.choice((-1, 0, 1))) for q in active for a in ("_", "0", "1")}
            table = AsyncTable.from_export(dict(alphabet=["_", "0", "1"], blank="_",
                tape_count=k, active_tape=active, transitions=rules, halting_states=["HALT"], initial_state="0"))
            original = Configuration("0", [{0: rng.choice(("0", "1"))} for _ in range(k)], [0] * k)
            counts = compare(table, original, -35, 35, source_limit=30)
            total_source += counts[0]
            total_physical += counts[1]
            random_cases += 1

    # Capacity is a precondition, never an implicit growing tape or a host
    # operation. A move outside the reserved interval visibly rejects.
    tiny = AsyncTable.from_export(dict(alphabet=["_"], blank="_", tape_count=1,
        active_tape={"q": 0}, transitions={("q", "_"): ("q", "_", 1)},
        halting_states=[], initial_state="q"))
    physical = Runner(Compiler(tiny), Configuration("q", [{}], [0]), 0, 0)
    try:
        physical.macrostep()
    except ValueError:
        assert physical.control == REJECT
    else:
        raise AssertionError("Insufficient capacity accepted")

    result = dict(status="PASS", arithmetic_and_memory_executions=cases,
        random_table_executions=random_cases, source_steps_compared=total_source,
        single_physical_head_microsteps=total_physical, capacity_negative_control="PASS",
        manifests=manifests, runtime_seconds=round(time.perf_counter() - started, 3),
        scope="Exact classical spatial lowering; no reversible single-head or Hamiltonian claim")
    (HERE / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
