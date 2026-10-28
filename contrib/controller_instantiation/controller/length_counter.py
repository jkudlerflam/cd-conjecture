#!/usr/bin/env python3
"""Explicit reversible blank-to-length computation and exact uncomputation.

On static LEFT BULK^(n-2) RIGHT, n>=4, both dynamic bits initially zero,
the unique head computes n-1 in little-endian BULK bits and places a frontier
at RIGHT.  It then reverses its complete trajectory and restores blank bits.
There is no program tape or quantum-data oracle in this primitive.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path


@dataclass(frozen=True)
class Symbol:
    static: str
    bit: int = 0
    frontier: int = 0


@dataclass
class Configuration:
    tape: list[Symbol]
    position: int = 0
    control: str = "F:ADVANCE"

    def snapshot(self):
        return self.position, self.control, tuple(self.tape)


FORWARD = (
    ("ADVANCE", 1), ("BACK", -1), ("CARRY", 1), ("RETURN", -1),
    ("SCAN", 1), ("FINAL_BACK", -1), ("FINAL_CARRY", 1), ("FINAL_RETURN", -1),
)


class LengthCounter:
    def __init__(self):
        self.symbols = tuple(Symbol(mark, bit, frontier)
                             for mark in ("LEFT", "BULK", "RIGHT")
                             for bit in range(2) for frontier in range(2))
        self.symbol_index = {s: i for i, s in enumerate(self.symbols)}
        self.controls = tuple(f"{sign}:{name}" for sign in ("F", "R") for name, _ in FORWARD)
        self.control_index = {q: i for i, q in enumerate(self.controls)}
        self.direction = {f"{sign}:{name}": direction * (1 if sign == "F" else -1)
                          for sign in ("F", "R") for name, direction in FORWARD}
        self.B = len(self.symbols)
        self.K = len(self.controls)
        self.local_factor = (1 + 2 * self.K) * self.B
        self.partial = {}
        self.ranges = {}
        self.forward_rules = []
        self._build_forward()
        for q, a, qp, b in self.forward_rules:
            self._add("F:" + q, a, "F:" + qp, b)
            # Reversed labels carry opposite movement directions.  After a
            # reversed move the collision is the inverse forward collision.
            self._add("R:" + qp, b, "R:" + q, a)
        left = Symbol("LEFT")
        # The terminal collision and the initial reverse collision do not
        # occur in forward_rules. Their otherwise-unused ranges let us close
        # the reversible orbit without an absorbing halt or a local merge.
        self._add("F:FINAL_RETURN", left, "R:FINAL_RETURN", left)
        self._add("R:ADVANCE", left, "F:ADVANCE", left)
        self.collision = self._complete()
        self.inverse = [0] * (self.K * self.B)
        for src, dst in enumerate(self.collision):
            self.inverse[dst] = src

    def encode(self, q, symbol):
        return self.control_index[q] * self.B + self.symbol_index[symbol]

    def decode(self, x):
        q, a = divmod(x, self.B)
        return self.controls[q], self.symbols[a]

    def _add(self, q, a, qp, b):
        src, dst = self.encode(q, a), self.encode(qp, b)
        if src in self.partial:
            raise AssertionError(f"Repeated local domain: {(q,a)}")
        if dst in self.ranges:
            raise AssertionError(f"Repeated local range: {(qp,b)}")
        self.partial[src] = dst
        self.ranges[dst] = src

    def _build_forward(self):
        def add(q, a, qp, b):
            self.forward_rules.append((q, a, qp, b))
        left = Symbol("LEFT")
        for bit in range(2):
            add("ADVANCE", Symbol("BULK", bit), "BACK", Symbol("BULK", bit, 1))
            for prefix in ("", "FINAL_"):
                add(prefix + "BACK", Symbol("BULK", bit), prefix + "BACK", Symbol("BULK", bit))
            add("SCAN", Symbol("BULK", bit), "SCAN", Symbol("BULK", bit))
            add("SCAN", Symbol("BULK", bit, 1), "ADVANCE", Symbol("BULK", bit))
        add("ADVANCE", Symbol("RIGHT"), "FINAL_BACK", Symbol("RIGHT", 0, 1))
        for prefix in ("", "FINAL_"):
            add(prefix + "BACK", left, prefix + "CARRY", left)
            for frontier in range(2):
                add(prefix + "CARRY", Symbol("BULK", 1, frontier),
                    prefix + "CARRY", Symbol("BULK", 0, frontier))
                add(prefix + "CARRY", Symbol("BULK", 0, frontier),
                    prefix + "RETURN", Symbol("BULK", 1, frontier))
                add(prefix + "RETURN", Symbol("BULK", 0, frontier),
                    prefix + "RETURN", Symbol("BULK", 0, frontier))
        add("RETURN", left, "SCAN", left)

    def _complete(self):
        count = self.B * self.K
        table = [-1] * count
        for src, dst in self.partial.items():
            table[src] = dst
        unused_domain = set(range(count)) - self.partial.keys()
        unused_range = set(range(count)) - self.ranges.keys()
        for mark in ("LEFT", "BULK", "RIGHT"):
            domain = {x for x in unused_domain if self.decode(x)[1].static == mark}
            image = {x for x in unused_range if self.decode(x)[1].static == mark}
            assert len(domain) == len(image)
            common = domain & image
            for x in common:
                table[x] = x
            for src, dst in zip(sorted(domain - common), sorted(image - common)):
                table[src] = dst
        assert sorted(table) == list(range(count))
        return table

    def initialize(self, n):
        if n < 4:
            raise ValueError("At least four sites are required for the counter to fit")
        return Configuration([Symbol("LEFT")] + [Symbol("BULK")] * (n - 2) + [Symbol("RIGHT")])

    def step(self, config, prescribed=True):
        config.position += self.direction[config.control]
        if not 0 <= config.position < len(config.tape):
            raise ValueError("Head left tape")
        src = self.encode(config.control, config.tape[config.position])
        if prescribed and src not in self.partial:
            raise ValueError("Computation entered an unspecified sector or counter overflowed")
        config.control, config.tape[config.position] = self.decode(self.collision[src])

    @staticmethod
    def value(config):
        return sum(symbol.bit << i for i, symbol in enumerate(config.tape[1:-1]))

    @staticmethod
    def forward_steps(n):
        return n * n + 3 * n - 4 - 2 * (n - 1).bit_count()

    def run_cycle(self, n):
        config = self.initialize(n)
        initial = config.snapshot()
        midpoint = None
        total = 2 * self.forward_steps(n)
        for tick in range(1, total + 1):
            self.step(config)
            if config.control == "R:FINAL_RETURN" and config.position == 0:
                assert midpoint is None
                midpoint = {
                    "tick": tick, "counter": self.value(config),
                    "frontiers": [i for i, a in enumerate(config.tape) if a.frontier],
                }
        assert config.snapshot() == initial
        assert midpoint == {"tick": total // 2, "counter": n - 1, "frontiers": [n - 1]}
        return {"sites": n, "midpoint": midpoint, "full_cycle_steps": total}

    def manifest(self):
        return {
            "schema": "cd-blank-length-counter-v1",
            "scope": "Blank bits to n-1 at a reversible midpoint; no program generation",
            "neutral_alphabet_B": self.B,
            "direction_fixed_controls_K": self.K,
            "single_site_factor_with_fresh_updated_heads": self.local_factor,
            "prescribed_collision_entries": len(self.partial),
            "controls": [{"name": q, "direction": self.direction[q]} for q in self.controls],
            "neutral_symbols": [[s.static, s.bit, s.frontier] for s in self.symbols],
            "collision_permutation": self.collision,
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--sites", type=int, default=16)
    args = parser.parse_args()
    machine = LengthCounter()
    if args.output:
        args.output.write_text(json.dumps(machine.manifest(), indent=2) + "\n")
    print(json.dumps(machine.run_cycle(args.sites), indent=2))


if __name__ == "__main__":
    main()
