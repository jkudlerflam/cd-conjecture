#!/usr/bin/env python3
"""Independent exact checks of finite-table logging and implementation adapters."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import random
import sys
import time

from logged_lift import LoggedLift, SerialLoggedLift, SerialRule, Rule, Tape


ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def snapshot(control, tapes, history):
    return control, tuple(t.snapshot() for t in tapes), history.snapshot()


def full_space_fixture():
    alphabet = ("_", "0", "1", "#")
    rules = [Rule("scan", (bit, x), "scan", ("_", "_"), (1, -1))
             for bit in "01" for x in alphabet]
    rules += [Rule("scan", ("#", x), "halt", ("#", x), (0, 0)) for x in alphabet]
    lift = LoggedLift(rules, [alphabet, alphabet])
    collision, erasure = lift.finite_tables()
    for x, y in collision:
        assert lift.collision(*y) == x
    for x, y in erasure:
        assert (lift.erase_control(*y), y[1]) == x
    assert len({y for _, y in collision}) == len(collision)
    assert len({y for _, y in erasure}) == len(erasure)

    # Global malformed configurations exercise arbitrary controls, old history,
    # displaced heads and spectator cells, rather than only initialized runs.
    rng = random.Random(821104)
    for _ in range(10000):
        control = rng.choice(lift.controls)
        tapes = [Tape(cells={i: rng.choice(alphabet) for i in range(-2, 3)},
                      head=rng.randrange(-2, 3)) for _ in range(2)]
        history = Tape(cells={i: rng.choice(lift.history_alphabet) for i in range(-2, 3)},
                       head=rng.randrange(-2, 3))
        before = snapshot(control, tapes, history)
        control = lift.forward(control, tapes, history)
        control = lift.backward(control, tapes, history)
        assert snapshot(control, tapes, history) == before
    return dict(**lift.manifest(), exhaustive_collision_inputs=len(collision),
                exhaustive_erasure_inputs=len(erasure),
                malformed_full_configuration_roundtrips=10000)


def check_serial_table(lift):
    checked = 0
    for i, r in enumerate(lift.rules):
        # Every supported transposition, every same-symbol spectator filling,
        # and malformed history examples. Other spectator combinations are
        # identities by construction, as proved in README.
        for fill in lift.alphabets[0]:
            symbols = [fill] * lift.k
            symbols[r.tape] = r.read
            x = ("Q", r.source), tuple(symbols), "_"
            y = lift.collision(*x)
            assert y[0] == ("P", i) and y[2] == i
            assert y[1][r.tape] == r.write
            assert lift.collision(*y) == x
            assert lift.erase_control(("P", i), i) == ("Q", r.target)
            assert lift.erase_control(("Q", r.target), i) == ("P", i)
            malformed = ("Q", r.source), tuple(symbols), i
            assert lift.collision(*malformed) == malformed
            checked += 1
    return checked


def bitword(value):
    value = abs(value)
    return tuple(map(int, reversed(bin(value)[2:])))


def alu_checks():
    alu = load("independent_digit_alu", ROOT / "alu" / "machine.py")
    lift = SerialLoggedLift.from_specification(alu.TABLE)
    transposition_checks = check_serial_table(lift)
    examples = [(-13, 5), (17, 9), (0, 3), (7, -4)]
    cases = steps = 0
    for op in alu.OPERATIONS:
        for a, b in examples:
            if op in ("SHL", "SHR"):
                b = abs(b)
            runner = alu.Runner(op, bitword(a), bitword(b), a < 0, b < 0)
            tapes = [Tape(cells=t, head=h) for t, h in zip(runner.config.tapes, runner.config.heads)]
            history = Tape()
            control = ("Q", runner.config.state)
            initial = snapshot(control, tapes, history)
            count = 0
            while runner.config.state not in alu.TABLE["halting_states"]:
                runner.step()
                control = lift.forward(control, tapes, history)
                count += 1
                assert count < 1000000
            assert runner.config.state == "HALT"
            assert control == ("Q", runner.config.state)
            assert all(t.cells == source and t.head == head for t, source, head in
                       zip(tapes, runner.config.tapes, runner.config.heads))
            expected = {
                "ADD": lambda: a+b, "SUB": lambda: a-b, "MUL": lambda: a*b,
                "DIV": lambda: a//b, "MOD": lambda: a%b, "SHL": lambda: a<<b,
                "SHR": lambda: a>>b, "AND": lambda: a&b, "XOR": lambda: a^b,
                "EQ": lambda: int(a==b), "NE": lambda: int(a!=b),
                "LT": lambda: int(a<b), "LE": lambda: int(a<=b),
                "GT": lambda: int(a>b), "GE": lambda: int(a>=b),
            }[op]()
            out = tapes[alu.OUT]
            value, j = 0, 1
            while out.cells[j] != "R":
                value += int(out.cells[j]) << (j-1)
                j += 1
            if out.cells[0] == "M":
                value = -value
            assert value == expected, (op, a, b, value, expected)
            for _ in range(count):
                control = lift.backward(control, tapes, history)
            assert snapshot(control, tapes, history) == initial
            cases += 1
            steps += count
    return dict(**lift.manifest(), supported_collision_checks=transposition_checks,
                exact_source_run_comparisons=cases, source_steps_logged_and_reversed=steps,
                arithmetic_results_checked_against_host_reference=cases)


def memory_checks():
    memory = load("independent_digit_memory", ROOT / "memory" / "machine.py")
    rules = [SerialRule(q, memory.TAPES.index(r.tape), symbol, r.target, r.write, r.move)
             for (q, symbol), r in sorted(memory.TABLE.rules.items())]
    lift = SerialLoggedLift(rules, [memory.ALPHABET] * len(memory.TAPES), ["ACCEPT"])
    transposition_checks = check_serial_table(lift)
    cases = steps = 0
    for words in ([], [(0,)], [(1,), (0, 1), (1, 1, 0, 1)]):
        for op in ("GET", "SET", "APPEND", "LEN"):
            for address in range(len(words) if op in ("GET", "SET") else 1):
                value = (1, 0, 1, 1, 0, 1)
                source = memory.Machine(op, memory.memory_frame(words),
                                        memory.bits_frame(bitword(address)), memory.bits_frame(value))
                tapes = [Tape(cells=source.tapes[name].cells, head=source.tapes[name].head)
                         for name in memory.TAPES]
                history = Tape()
                control = ("Q", source.state)
                initial = snapshot(control, tapes, history)
                count = 0
                while source.step():
                    control = lift.forward(control, tapes, history)
                    count += 1
                    assert count < 1000000
                assert source.state == "ACCEPT"
                assert control == ("Q", source.state)
                assert all(t.cells == source.tapes[name].cells and t.head == source.tapes[name].head
                           for name, t in zip(memory.TAPES, tapes))
                for _ in range(count):
                    control = lift.backward(control, tapes, history)
                assert snapshot(control, tapes, history) == initial
                cases += 1
                steps += count
    return dict(**lift.manifest(), supported_collision_checks=transposition_checks,
                exact_source_run_comparisons=cases, source_steps_logged_and_reversed=steps)


if __name__ == "__main__":
    started = time.monotonic()
    result = {"full_space_fixture": full_space_fixture(),
              "actual_arithmetic_table": alu_checks(),
              "actual_serial_memory_table": memory_checks(),
              "status": "pass", "runtime_seconds": round(time.monotonic()-started, 3)}
    Path(__file__).with_name("results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
