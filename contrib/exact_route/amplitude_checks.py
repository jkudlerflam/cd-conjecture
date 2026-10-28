"""Numerical checks of exact finite-menu uniform interval preparation."""

import json
import math
from pathlib import Path

import numpy as np


def prepare_uniform(m, width):
    state = np.zeros(1 << width)
    state[0] = 1
    exceptional = 0
    for depth in range(width):
        block = 1 << (width - depth)
        half = block >> 1
        partial_at_depth = 0
        for prefix in range(1 << depth):
            start = prefix * block
            count = min(block, max(0, m - start))
            if count == 0:
                continue
            if count == block:
                cosine = sine = 1 / math.sqrt(2)
            else:
                partial_at_depth += 1
                exceptional += 1
                left = min(half, count)
                cosine = math.sqrt(left / count)
                sine = math.sqrt((count - left) / count)
            # Apply the unitary to every suffix, not just its populated input.
            a = state[start : start + half].copy()
            b = state[start + half : start + block].copy()
            state[start : start + half] = cosine * a - sine * b
            state[start + half : start + block] = sine * a + cosine * b
        assert partial_at_depth <= 1
    assert exceptional <= width
    return state, exceptional


def main():
    error = 0.0
    exceptional_max = 0
    interval_cases = 0
    branch_cases = 0
    for m in range(1, 257):
        width = (m - 1).bit_length()
        state, exceptional = prepare_uniform(m, width)
        expected = np.zeros_like(state)
        expected[:m] = 1 / math.sqrt(m)
        error = max(error, float(np.max(np.abs(state - expected))))
        error = max(error, abs(float(np.dot(state, state)) - 1))
        exceptional_max = max(exceptional_max, exceptional)
        interval_cases += 1
    for upper in range(1, 17):
        for edges in range(1, 65):
            width = (max(upper, edges) - 1).bit_length()
            a, count_a = prepare_uniform(upper, width)
            b, count_b = prepare_uniform(edges, width)
            joint = np.concatenate([a, b]) / math.sqrt(2)
            expected = np.zeros_like(joint)
            expected[:upper] = 1 / math.sqrt(2 * upper)
            expected[(1 << width) : (1 << width) + edges] = 1 / math.sqrt(2 * edges)
            error = max(error, float(np.max(np.abs(joint - expected))))
            error = max(error, abs(float(np.dot(joint, joint)) - 1))
            assert count_a + count_b <= 2 * width
            branch_cases += 1
    assert error < 2e-14
    report = {
        "interval_cases": interval_cases,
        "two_branch_cases": branch_cases,
        "maximum_exceptional_rotations_one_interval": exceptional_max,
        "maximum_error": error,
        "status": "passed",
    }
    Path(__file__).with_name("amplitude_verification.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
