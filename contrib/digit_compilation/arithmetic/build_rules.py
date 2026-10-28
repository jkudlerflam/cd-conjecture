#!/usr/bin/env python3
"""Freeze a finite, width-independent two-cell arithmetic rewrite relation.

This is a rule *compiler*. All arithmetic here is on at most two input bits
and one member of a fixed, three-element comparison state. Runtime execution
loads the emitted JSON and performs dictionary lookups only.
"""
from __future__ import annotations

from itertools import product
import hashlib
import json
from pathlib import Path


ARITHMETIC = ("ADD", "SUB", "AND", "XOR", "COPY")
COMPARISONS = ("EQ", "NE", "LT", "LE", "GT", "GE")
OPERATIONS = ARITHMETIC + COMPARISONS
HERE = Path(__file__).resolve().parent


def digit(kind, a, b, scratch=0, output=0, history="_"):
    return f"{kind}{a}{b}{scratch}{output}{history}"


def head(op, phase, state, cell):
    return f"H:{op}:{phase}:{state}:{cell}"


def states(op):
    return ("L", "E", "G") if op in COMPARISONS else (("0", "1") if op in ("ADD", "SUB") else ("0",))


def initial(op):
    return "E" if op in COMPARISONS else ("1" if op == "SUB" else "0")


def truth(op, q):
    return int({"EQ": q == "E", "NE": q != "E", "LT": q == "L",
                "LE": q != "G", "GT": q == "G", "GE": q != "L"}[op])


def bit_transition(op, a, b, q, last):
    if op in ("ADD", "SUB"):
        effective_b = b if op == "ADD" else 1 - b
        total = a + effective_b + int(q)
        return str(total // 2), total % 2
    if op == "AND":
        return q, a & b
    if op == "XOR":
        return q, a ^ b
    if op == "COPY":
        return q, a
    # More significant unequal bits supersede earlier comparison state.
    # The last (sign) bit reverses the unsigned ordering.
    if a != b:
        return ("L" if a < b else "G") if not last else ("G" if a < b else "L"), 0
    return q, 0


def build():
    rules = {}

    def add(left, right):
        lhs, rhs = tuple(left), tuple(right)
        assert lhs not in rules or rules[lhs] == rhs, (lhs, rhs, rules.get(lhs))
        rules[lhs] = rhs

    clean_next = [digit("D", a, b) for a, b in product((0, 1), repeat=2)] + ["R"]
    returned_next = [digit("D", a, b, output=o) for a, b, o in product((0, 1), repeat=3)] + ["R"]
    for op in OPERATIONS:
        for a, b in product((0, 1), repeat=2):
            first = digit("F", a, b)
            add([head(op, "enter", "_", "L"), first],
                ["L", head(op, "scan", initial(op), first)])
        for kind, a, b, q, nxt in product(("F", "D"), (0, 1), (0, 1), states(op), clean_next):
            new_q, scratch = bit_transition(op, a, b, q, nxt == "R")
            clean = digit(kind, a, b)
            used = digit(kind, a, b, scratch=scratch, history=q)
            add([head(op, "scan", q, clean), nxt],
                [used, head(op, "scan", new_q, nxt)])
            if nxt == "R":
                back_state = new_q + str(truth(op, new_q)) if op in COMPARISONS else new_q
                add([used, head(op, "scan", new_q, "R")],
                    [used, head(op, "back", back_state, "R")])
        for kind, a, b, q, nxt in product(("F", "D"), (0, 1), (0, 1), states(op), returned_next):
            new_q, scratch = bit_transition(op, a, b, q, nxt == "R")
            used = digit(kind, a, b, scratch=scratch, history=q)
            for verdict in ((0, 1) if op in COMPARISONS else (None,)):
                before = new_q + str(verdict) if verdict is not None else new_q
                after = q + str(verdict) if verdict is not None else q
                out = (verdict if kind == "F" else 0) if verdict is not None else scratch
                clean = digit(kind, a, b, output=out)
                add([used, head(op, "back", before, nxt)],
                    [head(op, "back", after, clean), nxt])
        for a, b, out in product((0, 1), repeat=3):
            first = digit("F", a, b, output=out)
            back_state = initial(op) + str(out) if op in COMPARISONS else initial(op)
            add(["L", head(op, "back", back_state, first)],
                [head(op, "done", "_", "L"), first])
    serial = [{"left": list(left), "right": list(right)} for left, right in sorted(rules.items())]
    inverse = [tuple(rule["right"]) for rule in serial]
    assert len(inverse) == len(set(inverse)), "Local relation must be injective"
    return {"format": "two-adjacent-cell-rewrite-v1", "operations": list(OPERATIONS), "rules": serial}


def manifest(payload):
    symbols = {s for rule in payload["rules"] for side in ("left", "right") for s in rule[side]}
    control_states = {tuple(s.split(":")[1:4]) for s in symbols if s.startswith("H:")}
    tape_symbols = {s.split(":")[-1] for s in symbols}
    return {"rules": len(payload["rules"]), "symbols_used": len(symbols),
            "tape_symbols_used": len(tape_symbols), "control_states_used": len(control_states),
            "sha256": hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
            "width_dependent_rules_or_states": False,
            "scope": "Aligned signed digit kernels with autonomous compute/copy/uncompute; no addressed memory or full VM"}


if __name__ == "__main__":
    payload = build()
    (HERE / "rules.json").write_text(json.dumps(payload, separators=(",", ":")) + "\n")
    info = manifest(payload)
    (HERE / "manifest.json").write_text(json.dumps(info, indent=2) + "\n")
    print(json.dumps(info, indent=2))
