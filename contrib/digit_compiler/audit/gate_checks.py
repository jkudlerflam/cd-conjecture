"""Independent full-space audit of unit moves, collision and movement gates."""
from itertools import product
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from reversible_head import Source, Compiler, Runner, Work, History
from reversible_gates import UnitMoves, Gates, Site


def checks():
    source = Source.from_dictionary(("0", "1"),
            {("q", a): ("halt", "0", 0) for a in "01"}, {"halt"})
    compiler = Compiler(source)
    controls, alphabet, rules = compiler.materialize()
    normalized = UnitMoves(compiler)
    unit_controls = [("BASE", q) for q in controls]
    unit_controls += [("BOUNCE", q) for q in controls if compiler.incoming_direction(q) == 0]
    inverse, directions = {}, {}
    for q, a in product(unit_controls, alphabet):
        out = normalized.local_rule(q, a)
        if out is None:
            continue
        r, b, d = out
        assert r in unit_controls and b in alphabet and d in (-1, 1)
        assert (r, b) not in inverse
        inverse[r, b] = q, a, d
        assert directions.setdefault(r, d) == d
        assert normalized.incoming_direction(r) == d
        assert normalized.inverse_rule(r, b) == (q, a, d)
    for r, b in product(unit_controls, alphabet):
        assert normalized.inverse_rule(r, b) == inverse.get((r, b))

    # Every source rule, every possible symbol in the bounce's adjacent cell.
    normalized_cases = 0
    for (q, a), (r, b, d) in rules.items():
        for neighbor in alphabet:
            tape = [neighbor, a, neighbor]
            control, position = ("BASE", q), 1
            for _ in range(2 if d == 0 else 1):
                control, tape[position], move = normalized.local_rule(control, tape[position])
                position += move
            assert control == ("BASE", r) and position == 1+d
            assert tape == [neighbor, b, neighbor]
            normalized_cases += 1

    gates = Gates(normalized)
    sites = [Site(a) for a in alphabet]
    sites += [Site(a, q, phase) for a, q, phase in product(alphabet, unit_controls, (0, 1))]
    collision_outputs = set()
    for s in sites:
        out = gates.collision(s)
        assert gates.collision(out) == s
        collision_outputs.add(out)
    assert len(collision_outputs) == len(sites)
    sectors = {0: 0, 1: 0, 2: 0}
    for left, right in product(sites, repeat=2):
        out = gates.movement(left, right)
        assert gates.movement(*out) == (left, right)
        assert out[0].symbol == left.symbol and out[1].symbol == right.symbol
        heads = int(left.control is not None)+int(right.control is not None)
        assert sum(s.control is not None for s in out) == heads
        sectors[heads] += 1

    # The same gate word reflects at both ends of its extended computation
    # orbit. These checks intentionally use no inverse call during evolution.
    erase_source = Source.from_dictionary(("_", "0", "1", "#"),
            {("erase", "0"): ("erase", "_", 1),
             ("erase", "1"): ("erase", "_", 1),
             ("erase", "#"): ("halt", "#", 0)}, {"halt"})
    erase_compiler = Compiler(erase_source)
    erase_gates = Gates(UnitMoves(erase_compiler))
    cycles = []
    for word in ("", "0", "1", "010", "10110"):
        runner = Runner(erase_compiler, "erase", [*word, "#"], 0, len(word)+2)
        start = [Site(a, ("BASE", runner.control) if j == runner.position else None)
                 for j, a in enumerate(runner.cells)]
        evolved, answer = start, False
        for count in range(1, 20000):
            evolved = erase_gates.forward_sweep(evolved)
            heads = [(j, s) for j, s in enumerate(evolved) if s.control is not None]
            assert len(heads) == 1
            _, head = heads[0]
            if head.control == ("BASE", ("F", "halt")) and head.phase == 0:
                assert [s.symbol.symbol for s in evolved[1:len(word)+2]] == ["_"]*len(word)+["#"]
                answer = True
            if evolved == start:
                assert answer
                cycles.append(dict(word=word, sweeps=count))
                break
        else:
            raise AssertionError("No exact return")
    return {"status": "pass", "normalized_controls": len(unit_controls),
            "normalized_rules": len(inverse), "source_rule_bounce_contexts": normalized_cases,
            "collision_full_space_inputs": len(sites),
            "movement_full_space_pairs": len(sites)**2,
            "movement_pair_counts_by_head_number": sectors,
            "independent_forward_only_autonomous_cycles": cycles,
            "scope": "Full finite fixture gate spaces and initialized finite erasure orbits; no general malformed-sector Hamiltonian claim"}


if __name__ == "__main__":
    result = checks()
    Path(__file__).with_name("gate_results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))
