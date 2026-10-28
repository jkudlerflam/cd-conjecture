#!/usr/bin/env python3
"""Autonomous cycles made from fixed full-space local permutation gates."""
from itertools import product
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "audit"))
from reversible_head import Source, Compiler, Runner, Work, History
from reversible_gates import UnitMoves, Gates, Site


def main():
    started = time.perf_counter()
    transitions = {("erase", a): ("erase", "_", 1) for a in "01"}
    transitions["erase", "#"] = "halt", "#", 0
    reversible = Compiler(Source.from_dictionary(("_", "0", "1", "#"), transitions, {"halt"}))
    normalized = UnitMoves(reversible)
    gates = Gates(normalized)
    cases = []
    for length in range(6):
        for word in product("01", repeat=length):
            tape = Runner(reversible, "erase", [*word, "#"], 0, length+2)
            initial = [Site(a, ("BASE", tape.control) if i == tape.position else None)
                       for i, a in enumerate(tape.cells)]
            sites = initial
            saw_answer = False
            max_records = 0
            for sweep in range(1, 20001):
                old = sites
                sites = gates.forward_sweep(sites)
                assert gates.backward_sweep(sites) == old
                heads = [s for s in sites if s.control is not None]
                assert len(heads) == 1
                head, = heads
                max_records = max(max_records, sum(isinstance(s.symbol, History) for s in sites))
                if head.control == ("BASE", ("F", "halt")) and head.phase == 0:
                    output = [s.symbol.symbol for s in sites[1:len(word)+2]]
                    assert output == ["_"]*length + ["#"]
                    saw_answer = True
                if sites == initial:
                    assert saw_answer and max_records == length+1
                    cases.append(dict(word="".join(word), physical_sites=len(sites),
                                      fixed_sweeps_per_cycle=sweep, peak_history_records=max_records))
                    break
            else:
                raise AssertionError(("No autonomous return", word))
    result = dict(status="PASS", autonomous_cycles=len(cases),
        complete_fixed_sweeps=sum(c["fixed_sweeps_per_cycle"] for c in cases),
        all_forward_answers_and_exact_initial_returns=True,
        scheduling="Repeat exactly the same onsite collision layer and left-to-right bond movement sweep",
        scope="Classical erasure fixture, initialized tape. No output retained after reset or full Hamiltonian claim",
        cases=cases, runtime_seconds=round(time.perf_counter()-started, 3))
    (HERE / "reversible_results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k:v for k,v in result.items() if k != "cases"}, indent=2))


if __name__ == "__main__":
    main()
