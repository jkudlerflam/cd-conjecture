#!/usr/bin/env python3
"""Structural proofs-by-enumeration and independent semantic kernel tests."""
from collections import Counter, defaultdict
from itertools import product
import json
from pathlib import Path
import random
import time

from arithmetic import BACKWARD, FORWARD, FROZEN, Rewriter, SUPPORTED, encode, execute
from build_rules import build, manifest


HERE = Path(__file__).resolve().parent


def oracle(op, a, b, width):
    if op == "ADD":
        value = a + b
    elif op == "SUB":
        value = a - b
    elif op == "AND":
        value = a & b
    elif op == "XOR":
        value = a ^ b
    elif op == "COPY":
        value = a
    else:
        return int({"EQ": a == b, "NE": a != b, "LT": a < b,
                    "LE": a <= b, "GT": a > b, "GE": a >= b}[op])
    unsigned = value % (1 << width)
    return unsigned - (1 << width) if unsigned >= 1 << (width - 1) else unsigned


def structure():
    assert FROZEN == build(), "The reviewed finite generator must reproduce the frozen rules"
    assert manifest(FROZEN) == json.loads((HERE / "manifest.json").read_text())
    assert len(FORWARD) == len(BACKWARD) == len(FROZEN["rules"])
    sides = [defaultdict(set), defaultdict(set)]
    movements = Counter()
    for lhs, rhs in FORWARD.items():
        heads = []
        for relation_side, pair in enumerate((lhs, rhs)):
            assert len(pair) == 2
            head_indices = [i for i, symbol in enumerate(pair) if symbol.startswith("H:")]
            assert len(head_indices) == 1
            at = head_indices[0]
            heads.append(at)
            sides[relation_side][pair[at]].add(at)
        movements[heads[1] - heads[0]] += 1
    # For every fully decorated head symbol, its applicable neighbor is
    # unique, in both directions. Together with unique pair keys this proves
    # global successor/predecessor uniqueness on EVERY one-head tape, not
    # only on successful test trajectories.
    assert all(len(positions) == 1 for side in sides for positions in side.values())
    return {**manifest(FROZEN), "unique_forward_and_backward_pair_keys": True,
            "unique_neighbor_side_per_head_symbol_both_directions": True,
            "head_moves": dict(sorted(movements.items()))}


def check_case(op, a, b, width):
    run = execute(op, a, b, width)
    assert run.output == oracle(op, a, b, width), (op, a, b, width, run.output)
    assert run.verify_inverse()
    return run.steps


def exhaustive():
    cases = steps = 0
    for width in range(1, 7):
        numbers = range(-(1 << (width - 1)), 1 << (width - 1))
        for op, a, b in product(sorted(SUPPORTED), numbers, numbers):
            steps += check_case(op, a, b, width)
            cases += 1
    return {"widths": list(range(1, 7)), "operand_pairs_times_operations": cases,
            "forward_steps": steps, "semantic_match_and_full_inverse_every_case": True}


def large_width():
    rng = random.Random(20261028)
    cases = steps = 0
    widths = [7, 31, 64, 257, 1024, 4096]
    for width in widths:
        low, high = -(1 << (width - 1)), (1 << (width - 1)) - 1
        # Sign boundaries, full carry/borrow chains, exact equality and
        # independent random strings. No width is passed to Rewriter.
        operands = [(low, low), (high, high), (low, high), (-1, 1),
                    (high, 1), (low, 1), (0, 0), (-1, -1)]
        operands += [(rng.randrange(low, high + 1), rng.randrange(low, high + 1)) for _ in range(8)]
        for op, (a, b) in product(sorted(SUPPORTED), operands):
            steps += check_case(op, a, b, width)
            cases += 1
    return {"widths": widths, "cases": cases, "forward_steps": steps,
            "semantic_match_and_full_inverse_every_case": True}


def every_prefix():
    count = 0
    for op in sorted(SUPPORTED):
        initial = encode(op, -73, 31, 9)
        forward = Rewriter(initial, 0)
        prefix_length = 0
        while True:
            inverse = Rewriter(forward.tape, forward.head_position)
            inverse.run(inverse=True)
            assert tuple(inverse.tape) == initial
            assert inverse.steps == prefix_length
            count += 1
            if not forward.step():
                break
            prefix_length += 1
    return {"prefix_states_checked": count, "every_prefix_has_exact_inverse_to_start": True}


def initialization():
    stalled = 0
    initial = encode("ADD", 3, 5, 7)
    for at in range(1, len(initial) - 1):
        for track, bad in ((3, "1"), (4, "1"), (5, "G")):
            dirty = list(initial)
            old = dirty[at]
            dirty[at] = old[:track] + bad + old[track + 1:]
            machine = Rewriter(dirty, 0)
            machine.run()
            assert machine.tape[0] != "H:ADD:done:_:L"
            stalled += 1
    return {"dirty_track_cases_stall": stalled,
            "scope": "A domain test, not a Hamiltonian penalty or illegal-sector gap proof"}


def main():
    start = time.monotonic()
    results = {"scope": "Fixed aligned digit kernels; neither a complete VM nor a counterexample",
               "structural": structure(), "exhaustive": exhaustive(),
               "large_width": large_width(), "prefix_inverse": every_prefix(),
               "initialization": initialization(), "all_passed": True}
    results["elapsed_seconds"] = round(time.monotonic() - start, 3)
    (HERE / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
