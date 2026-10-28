#!/usr/bin/env python3
"""Regression-check the explicit loop-accounted bounds, not a proof by sampling."""
from collections import Counter
import json
from pathlib import Path
import random
import time

from machine import OPERATIONS, run_binary


def magnitude_bits(value):
    return tuple(map(int, bin(abs(value))[:1:-1]))


def operation_bound(op, width, shift_count=0):
    assert width >= 1
    if op in ("SHL", "SHR"):
        assert shift_count >= 0
        return 64 * (shift_count + 1) * (width + shift_count + 5)
    return 128 * (width + 1) * (2 * width + 5)


def conditional_generator_bound(resource_a, word_steps):
    assert resource_a >= 1 and word_steps >= 1
    # ceil(log2(T+1)) = bit_length(T) for positive integers T.
    width = 2 * resource_a + word_steps.bit_length() + 64
    per_op = 128 * (width + 1) * (2 * width + 5)
    return {"A": resource_a, "word_steps_T": word_steps, "operand_bit_bound_B": width,
            "digit_steps_per_arithmetic_opcode": per_op,
            "ALU_digit_steps": word_steps * per_op}


def main():
    start = time.monotonic()
    rng = random.Random(77031)
    count = 0
    total_moves = 0
    worst_ratio = 0
    tightest = None
    per_operation = Counter()
    for op in OPERATIONS:
        for width in (1, 2, 3, 4, 8, 16, 32, 64, 128, 256):
            cap = (1 << width) - 1
            pairs = [(cap, cap), (-cap, 1), (1 << width, -cap)]
            pairs += [(rng.randrange(-cap, cap + 1), rng.randrange(-cap, cap + 1)) for _ in range(2)]
            for left, right in pairs:
                if op in ("SHL", "SHR"):
                    right = abs(right) % 97
                if op in ("DIV", "MOD") and right == 0:
                    right = 1
                left_bits, right_bits = magnitude_bits(left), magnitude_bits(right)
                b = max(len(left_bits), len(right_bits))
                result = run_binary(op, left_bits, right_bits, left < 0, right < 0)
                upper = operation_bound(op, b, right if op in ("SHL", "SHR") else 0)
                assert result["steps"] <= upper, (op, left, right, result["steps"], upper)
                count += 1
                total_moves += result["steps"]
                per_operation[op] += 1
                ratio = result["steps"] / upper
                if ratio > worst_ratio:
                    worst_ratio = ratio
                    tightest = {"op": op, "magnitude_width": b,
                                "shift_count": right if op in ("SHL", "SHR") else None,
                                "actual_steps": result["steps"], "bound": upper}
    grid = []
    for a in (1, 2, 7, 100, 10**6):
        for steps in (1, 3, 10, 10**4, 10**12):
            row = conditional_generator_bound(a, steps)
            b = row["operand_bit_bound_B"]
            assert operation_bound("SHL", b, 2*a) <= row["digit_steps_per_arithmetic_opcode"]
            assert row["ALU_digit_steps"] <= 256 * steps * (b + 3)**2
            grid.append(row)
    polynomial_checks = []
    for a in (1, 2, 3, 7, 100, 10**6, 10**30):
        word_bound = 18463 * (128*a)**6
        row = conditional_generator_bound(a, word_bound)
        assert word_bound < 2**57 * a**6
        assert row["operand_bit_bound_B"] <= 123*a <= 128*a
        assert row["ALU_digit_steps"] < 2**80 * a**8
        polynomial_checks.append({"A": a, "operand_bit_bound_B": row["operand_bit_bound_B"],
                                  "coarse_ALU_digit_step_bound": 2**80 * a**8})
    result = {"status": "pass", "scope": "Regression evidence for the documented loop-accounted bounds",
              "cases": count, "microsteps": total_moves, "by_operation": dict(per_operation),
              "largest_actual_to_bound_ratio": worst_ratio, "tightest_test_case": tightest,
              "conditional_generator_parameter_grid": grid,
              "conditional_polynomial_bound_checks": polynomial_checks,
              "elapsed_seconds": round(time.monotonic() - start, 3)}
    Path(__file__).with_name("bounds_results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "conditional_generator_parameter_grid"}, indent=2))


if __name__ == "__main__":
    main()
