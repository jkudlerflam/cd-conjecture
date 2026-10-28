#!/usr/bin/env python3
"""Execute a frozen adjacent-cell rewrite table; no runtime word arithmetic.

Integers enter/leave only through encode/decode adapters. Rewriting itself
does not parse digits, calculate results, or receive the operand width.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
FROZEN = json.loads((HERE / "rules.json").read_text())
SUPPORTED = frozenset(FROZEN["operations"])
FORWARD = {tuple(rule["left"]): tuple(rule["right"]) for rule in FROZEN["rules"]}
BACKWARD = {right: left for left, right in FORWARD.items()}
assert len(FORWARD) == len(BACKWARD) == len(FROZEN["rules"])


class Rewriter:
    def __init__(self, tape, head_position):
        self.tape = list(tape)
        self.head_position = head_position
        self.steps = 0

    def step(self, inverse=False):
        relation = BACKWARD if inverse else FORWARD
        matches = []
        # A single head can occur in at most the two adjacent pairs below.
        # The finite rules determine both direction and phase, not the driver.
        for position in (self.head_position - 1, self.head_position):
            if 0 <= position < len(self.tape) - 1:
                pair = (self.tape[position], self.tape[position + 1])
                replacement = relation.get(pair)
                if replacement is not None:
                    matches.append((position, replacement))
        if not matches:
            return False
        assert len(matches) == 1, "Nondeterministic local transition"
        position, replacement = matches[0]
        self.tape[position:position + 2] = replacement
        self.head_position = position if replacement[0].startswith("H:") else position + 1
        self.steps += 1
        return True

    def run(self, inverse=False):
        while self.step(inverse=inverse):
            pass


def encode(op, a, b, width):
    """Host-side fixture adapter, not part of the local machine."""
    if op not in SUPPORTED:
        raise ValueError(f"Unsupported digit operation {op}")
    if width < 1 or not -(1 << (width - 1)) <= a < (1 << (width - 1)) or not -(1 << (width - 1)) <= b < (1 << (width - 1)):
        raise ValueError("Operands must fit the requested signed width")
    tape = [f"H:{op}:enter:_:L"]
    for i in range(width):
        kind = "F" if i == 0 else "D"
        tape.append(f"{kind}{(a >> i) & 1}{(b >> i) & 1}00_")
    tape.append("R")
    return tuple(tape)


def decode(tape, signed=True):
    """Host-side output adapter, never called from a local transition."""
    value = sum(int(cell[4]) << i for i, cell in enumerate(tape[1:-1]))
    width = len(tape) - 2
    if signed and tape[-2][4] == "1":
        value -= 1 << width
    return value


@dataclass
class Run:
    operation: str
    output: int
    steps: int
    initial_tape: tuple[str, ...]
    final_tape: tuple[str, ...]

    def verify_inverse(self):
        machine = Rewriter(self.final_tape, 0)
        machine.run(inverse=True)
        assert tuple(machine.tape) == self.initial_tape
        assert machine.steps == self.steps
        return True


def execute(op, a, b, width):
    initial = encode(op, a, b, width)
    machine = Rewriter(initial, 0)
    machine.run()
    assert machine.tape[0] == f"H:{op}:done:_:L", "Dirty/malformed tape stalled before completion"
    assert machine.head_position == 0
    for before, after in zip(initial[1:-1], machine.tape[1:-1]):
        assert before[:3] == after[:3] and after[3] == "0" and after[5:] == "_"
    assert machine.steps == 2 * width + 3
    result = decode(machine.tape, signed=op not in ("EQ", "NE", "LT", "LE", "GT", "GE"))
    return Run(op, result, machine.steps, initial, tuple(machine.tape))


def apply(op, a, b, width):
    """Return signed b-bit arithmetic (booleans unsigned); ADD/SUB wrap."""
    return execute(op, a, b, width).output
