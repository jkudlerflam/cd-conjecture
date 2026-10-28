#!/usr/bin/env python3
"""One finite reversible controller: blank -> length -> program -> circuit -> blank.

This is an autonomous, length-dependent demonstration program, not the full
sparse-encoding circuit generator.  Every phase handoff is a collision entry.
There is no quantum-data-dependent transition and no host-side tape mutation.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import importlib.util
import json
from pathlib import Path
import sys


COUNTER_PATH = (Path(__file__).resolve().parents[2] / "controller_instantiation"
                / "controller" / "length_counter.py")
_spec = importlib.util.spec_from_file_location("bridge_length_counter", COUNTER_PATH)
_counter = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _counter
_spec.loader.exec_module(_counter)


@dataclass(frozen=True)
class Symbol:
    static: str
    bit: int = 0
    frontier: int = 0
    program: int = 0
    cursor: int = 0


@dataclass(frozen=True)
class Control:
    name: str
    direction: int
    oracle: str | None = None
    adjoint: bool = False


@dataclass
class Configuration:
    tape: list[Symbol]
    position: int = 0
    control: str = "F:ADVANCE"

    def snapshot(self):
        return self.position, self.control, tuple(self.tape)


EXTRA_CONTROLS = (
    ("COPY", 1, None), ("COPY_BACK", -1, None), ("DISPATCH", 1, None),
    ("EXEC_0", 1, "CNOT"), ("EXEC_1", 1, "CV"),
    ("GATE_BACK_0", -1, None), ("GATE_BACK_1", -1, None),
)


class AutonomousBridge:
    """Generate a finite collision permutation with an initialized legal orbit.

    A step first moves the direction-fixed head and applies its two-site
    oracle, then collides with the arrived-at neutral symbol.  Reverse
    controls have opposite directions and adjoint oracle matrices.  The
    symbols retain their static LEFT/BULK/RIGHT component on every input,
    including malformed sectors filled by permutation completion.
    """

    def __init__(self):
        self.symbols = tuple(
            Symbol(mark, bit, frontier, program, cursor)
            for mark in ("LEFT", "BULK", "RIGHT")
            for bit in range(2) for frontier in range(2)
            for program in range(2) for cursor in range(2)
        )
        self.symbol_index = {a: i for i, a in enumerate(self.symbols)}
        forward_controls = tuple((q, d, None) for q, d in _counter.FORWARD) + EXTRA_CONTROLS
        self.controls = tuple(
            Control(f"{sign}:{q}", direction * (1 if sign == "F" else -1),
                    oracle, sign == "R")
            for sign in ("F", "R") for q, direction, oracle in forward_controls
        )
        self.control_index = {q.name: i for i, q in enumerate(self.controls)}
        self.control_by_name = {q.name: q for q in self.controls}
        self.B, self.K = len(self.symbols), len(self.controls)
        self.local_factor = (1 + 2 * self.K) * self.B
        self.partial, self.ranges = {}, {}
        self.forward_rules = []
        self._build_forward()
        for q, a, qp, b in self.forward_rules:
            self._add("F:" + q, a, "F:" + qp, b)
            self._add("R:" + qp, b, "R:" + q, a)
        # All handoffs, including both turnarounds, are literal local rules.
        right = Symbol("RIGHT", frontier=1)
        self._add("F:DISPATCH", right, "R:DISPATCH", right)
        left = Symbol("LEFT")
        self._add("R:ADVANCE", left, "F:ADVANCE", left)
        self.collision = self._complete()
        self.inverse = [0] * (self.B * self.K)
        for src, dst in enumerate(self.collision):
            self.inverse[dst] = src

    def encode(self, q, a):
        return self.control_index[q] * self.B + self.symbol_index[a]

    def decode(self, value):
        q, a = divmod(value, self.B)
        return self.controls[q].name, self.symbols[a]

    def _add(self, q, a, qp, b):
        src, dst = self.encode(q, a), self.encode(qp, b)
        if src in self.partial:
            raise AssertionError(f"Repeated local domain: {(q, a)}")
        if dst in self.ranges:
            raise AssertionError(f"Repeated local range: {(qp, b)}")
        self.partial[src], self.ranges[dst] = dst, src

    def _build_forward(self):
        def add(q, a, qp, b):
            self.forward_rules.append((q, a, qp, b))

        # Reuse the published, explicit counter rules, restricted to blank
        # program/cursor tracks.  Its stand-alone turnaround is NOT reused.
        counter = _counter.LengthCounter()
        for q, a, qp, b in counter.forward_rules:
            add(q, Symbol(a.static, a.bit, a.frontier),
                qp, Symbol(b.static, b.bit, b.frontier))

        left = Symbol("LEFT")
        right = Symbol("RIGHT", frontier=1)
        add("FINAL_RETURN", left, "COPY", left)
        for bit in range(2):
            for program in range(2):
                a = Symbol("BULK", bit=bit, program=program)
                # Reversible copy into the program track: defined as XOR on
                # BOTH target values, not a destructive assignment.
                add("COPY", a, "COPY", Symbol("BULK", bit=bit, program=program ^ bit))
                add("COPY_BACK", a, "COPY_BACK", a)
                marked = Symbol("BULK", bit=bit, program=program, cursor=1)
                add("DISPATCH", a, f"EXEC_{program}", marked)
                add(f"GATE_BACK_{program}", marked, "DISPATCH", a)
                for opcode in range(2):
                    # The gate was applied on the move into this site.
                    add(f"EXEC_{opcode}", a, f"GATE_BACK_{opcode}", a)
        add("COPY", right, "COPY_BACK", right)
        add("COPY_BACK", left, "DISPATCH", left)
        for opcode in range(2):
            add(f"EXEC_{opcode}", right, f"GATE_BACK_{opcode}", right)

    def _complete(self):
        size = self.B * self.K
        result = [-1] * size
        for src, dst in self.partial.items():
            result[src] = dst
        domain = set(range(size)) - self.partial.keys()
        image = set(range(size)) - self.ranges.keys()
        for mark in ("LEFT", "BULK", "RIGHT"):
            sources = {x for x in domain if self.decode(x)[1].static == mark}
            targets = {x for x in image if self.decode(x)[1].static == mark}
            assert len(sources) == len(targets)
            common = sources & targets
            for x in common:
                result[x] = x
            for src, dst in zip(sorted(sources - common), sorted(targets - common)):
                result[src] = dst
        assert sorted(result) == list(range(size))
        return result

    @staticmethod
    def initialize(n):
        if n < 4:
            raise ValueError("At least four physical sites are required")
        return Configuration([Symbol("LEFT")] + [Symbol("BULK")] * (n - 2)
                             + [Symbol("RIGHT")])

    def step(self, config, callback=None, require_prescribed=True):
        control = self.control_by_name[config.control]
        target = config.position + control.direction
        if not 0 <= target < len(config.tape):
            raise ValueError("Head left the physical tape")
        event = None
        if control.oracle is not None:
            event = (control.oracle, min(config.position, target), control.adjoint)
            if callback is not None:
                callback(*event)
        config.position = target
        src = self.encode(config.control, config.tape[target])
        if require_prescribed and src not in self.partial:
            raise ValueError("Trajectory entered a non-prescribed sector")
        config.control, config.tape[target] = self.decode(self.collision[src])
        return event

    @staticmethod
    def counter_value(config):
        return sum(a.bit << j for j, a in enumerate(config.tape[1:-1]))

    @staticmethod
    def program_value(config):
        return sum(a.program << j for j, a in enumerate(config.tape[1:-1]))

    @staticmethod
    def forward_steps(n):
        return _counter.LengthCounter.forward_steps(n) + 5 * n - 7

    @staticmethod
    def expected_forward_events(n):
        return [("CV" if (n - 1) >> j & 1 else "CNOT", j + 1, False)
                for j in range(n - 2)]

    def run_cycle(self, n, callback=None, retain_trajectory=False):
        config = self.initialize(n)
        initial = config.snapshot()
        forward = self.forward_steps(n)
        counter_time = _counter.LengthCounter.forward_steps(n)
        program_time = counter_time + 2 * (n - 1)
        events, trajectory, checkpoints = [], [], {}
        for tick in range(1, 2 * forward + 1):
            event = self.step(config, callback)
            if event:
                events.append(event)
            if retain_trajectory:
                trajectory.append(config.snapshot())
            if tick in (counter_time, program_time, forward):
                checkpoints[tick] = config.snapshot()
        assert config.snapshot() == initial
        expected = self.expected_forward_events(n)
        expected += [(op, at, True) for op, at, _ in reversed(expected)]
        assert events == expected
        return {"sites": n, "forward_steps": forward, "cycle_steps": 2 * forward,
                "events": events, "checkpoints": checkpoints, "trajectory": trajectory}

    def manifest(self):
        return {
            "schema": "cd-autonomous-length-program-bridge-v1",
            "scope": "Blank dynamic tracks to n-1 binary instruction tape and local circuit, then inverse",
            "neutral_alphabet_B": self.B,
            "direction_fixed_controls_K": self.K,
            "single_site_factor_with_fresh_updated_heads": self.local_factor,
            "single_site_dimension_with_one_quantum_qubit": 2 * self.local_factor,
            "prescribed_collision_entries": len(self.partial),
            "forward_time_formula": "n^2+8n-11-2 popcount(n-1)",
            "controls": [{"name": c.name, "direction": c.direction,
                          "oracle": c.oracle, "adjoint": c.adjoint} for c in self.controls],
            "neutral_symbols": [[a.static, a.bit, a.frontier, a.program, a.cursor]
                                for a in self.symbols],
            "collision_permutation": self.collision,
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--sites", type=int, default=16)
    args = parser.parse_args()
    machine = AutonomousBridge()
    if args.output:
        args.output.write_text(json.dumps(machine.manifest(), indent=2) + "\n")
    run = machine.run_cycle(args.sites)
    print(json.dumps({k: run[k] for k in ("sites", "forward_steps", "cycle_steps", "events")}, indent=2))


if __name__ == "__main__":
    main()
