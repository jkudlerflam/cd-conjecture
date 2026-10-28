#!/usr/bin/env python3
"""Compose the autonomous bridge and endpoint readout in one collision table.

The first/last BULK sites have static OUT_L/OUT_R marks.  The bridge's
cursor track doubles as the endpoint READY track, with disjoint phase
controls.  Both transport legs, READY intervals, all handoffs, and inverse
computation are performed by this one finite machine.
"""

from __future__ import annotations

from dataclasses import replace
import importlib.util
import json
from pathlib import Path
import sys

from bridge import AutonomousBridge, Configuration, Control, Symbol


_path = Path(__file__).resolve().parents[1] / "readout" / "transport.py"
_spec = importlib.util.spec_from_file_location("bridge_readout", _path)
_readout = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _readout
_spec.loader.exec_module(_readout)
MARKS = ("LEFT", "BULK", "RIGHT", "OUT_L", "OUT_R")


class IntegratedReadout(AutonomousBridge):
    def __init__(self):
        bridge = AutonomousBridge()
        readout = _readout.Readout()
        self.symbols = tuple(
            Symbol(mark, bit, frontier, program, cursor)
            for mark in MARKS for bit in range(2) for frontier in range(2)
            for program in range(2) for cursor in range(2)
        )
        self.symbol_index = {a: j for j, a in enumerate(self.symbols)}
        self.controls = bridge.controls + tuple(
            Control("RO:" + q, direction, None if oracle == "I" else oracle)
            for q, direction, oracle in _readout.CONTROLS
        ) + (Control("TO_READOUT", -1), Control("FROM_READOUT", 1))
        self.control_index = {q.name: j for j, q in enumerate(self.controls)}
        self.control_by_name = {q.name: q for q in self.controls}
        self.B, self.K = len(self.symbols), len(self.controls)
        self.local_factor = (1 + 2 * self.K) * self.B
        self.partial, self.ranges = {}, {}

        # OUT_L and OUT_R act exactly as BULK throughout the computation.
        # Remove only the bridge's forward-to-adjoint turnaround.
        for src, dst in bridge.partial.items():
            q, a = bridge.decode(src)
            qp, b = bridge.decode(dst)
            if q == "F:DISPATCH" and a.static == "RIGHT":
                continue
            marks = ("BULK", "OUT_L", "OUT_R") if a.static == "BULK" else (a.static,)
            for mark in marks:
                self._add(q, replace(a, static=mark), qp, replace(b, static=mark))

        # Readout changes only READY (=cursor) and leaves the other three
        # dynamic tracks intact, on all payload values.  Replace its final
        # cycle edge with a handoff back to reverse computation.
        for src, dst in readout.partial.items():
            q, (mark, ready) = readout.decode(src)
            qp, (target_mark, target_ready) = readout.decode(dst)
            assert mark == target_mark
            target_control = ("FROM_READOUT" if q == "FIND_START" and mark == "OUT_L"
                              else "RO:" + qp)
            for bit in range(2):
                for frontier in range(2):
                    for program in range(2):
                        a = Symbol(mark, bit, frontier, program, ready)
                        b = replace(a, cursor=target_ready)
                        self._add("RO:" + q, a, target_control, b)

        right = Symbol("RIGHT", frontier=1)
        self._add("F:DISPATCH", right, "TO_READOUT", right)
        self._add("FROM_READOUT", right, "R:DISPATCH", right)
        for a in self.symbols:
            if a.cursor:
                continue
            if a.static in ("BULK", "OUT_R"):
                self._add("TO_READOUT", a, "TO_READOUT", a)
                self._add("FROM_READOUT", a, "FROM_READOUT", a)
            elif a.static == "OUT_L":
                self._add("TO_READOUT", a, "RO:LEFT_OUT", a)

        self.collision = self._complete()
        self.inverse = [0] * (self.B * self.K)
        for src, dst in enumerate(self.collision):
            self.inverse[dst] = src

    def _complete(self):
        size = self.B * self.K
        result = [-1] * size
        for src, dst in self.partial.items():
            result[src] = dst
        domain = set(range(size)) - self.partial.keys()
        image = set(range(size)) - self.ranges.keys()
        for mark in MARKS:
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
        config = AutonomousBridge.initialize(n)
        config.tape[1] = replace(config.tape[1], static="OUT_L")
        config.tape[-2] = replace(config.tape[-2], static="OUT_R")
        return config

    @staticmethod
    def readout_extra_steps(n):
        return 2 * (n - 2) + 4 * (n - 1)

    @staticmethod
    def cycle_steps(n):
        return 2 * AutonomousBridge.forward_steps(n) + IntegratedReadout.readout_extra_steps(n)

    @staticmethod
    def expected_transport_events(n):
        # OUT_L and OUT_R are adjacent to the corresponding endpoints.
        return [("SWAP", 0, False), ("SWAP", n - 2, False),
                ("SWAP", 0, False), ("SWAP", n - 2, False)]

    def run_cycle(self, n, callback=None, retain_trajectory=False):
        config = self.initialize(n)
        initial = config.snapshot()
        events, trajectory = [], []
        for _ in range(self.cycle_steps(n)):
            event = self.step(config, callback)
            if event:
                events.append(event)
            if retain_trajectory:
                trajectory.append(config.snapshot())
        assert config.snapshot() == initial
        forward = self.expected_forward_events(n)
        expected = forward + self.expected_transport_events(n)
        expected += [(op, at, True) for op, at, _ in reversed(forward)]
        assert events == expected
        return {"sites": n, "cycle_steps": self.cycle_steps(n),
                "events": events, "trajectory": trajectory}

    def manifest(self):
        return {
            "schema": "cd-autonomous-length-program-with-readout-v1",
            "scope": "Blank dynamic tracks to binary program, quantum circuit, endpoint READY, inverse cleanup",
            "static_grammar": "LEFT OUT_L BULK^(n-4) OUT_R RIGHT",
            "neutral_alphabet_B": self.B,
            "direction_fixed_controls_K": self.K,
            "single_site_factor_with_fresh_updated_heads": self.local_factor,
            "single_site_dimension_with_one_quantum_qubit": 2 * self.local_factor,
            "prescribed_collision_entries": len(self.partial),
            "cycle_time_formula": "2(n^2+8n-11-2 popcount(n-1))+6n-8",
            "controls": [{"name": c.name, "direction": c.direction,
                          "oracle": c.oracle, "adjoint": c.adjoint} for c in self.controls],
            "neutral_symbols": [[a.static, a.bit, a.frontier, a.program, a.cursor]
                                for a in self.symbols],
            "collision_permutation": self.collision,
        }


if __name__ == "__main__":
    machine = IntegratedReadout()
    Path(__file__).with_name("integrated_transition_table.json").write_text(
        json.dumps(machine.manifest(), indent=2) + "\n")
    run = machine.run_cycle(8)
    print(json.dumps({k: run[k] for k in ("sites", "cycle_steps", "events")}, indent=2))
