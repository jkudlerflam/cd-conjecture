#!/usr/bin/env python3
"""An explicit, initialized autonomous-qubit demonstration Hamiltonian.

One n-independent two-site bond and two one-site boundary terms are given
as sparse columns on their ENTIRE Hilbert spaces. Local dimension 194400.
This is not the high-channel sparse code and is not a c-d counterexample.
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bridge"))
from local_gates import WORK, WORK_DIM, decode as work_decode
from local_gates import f_column, g_column, initial_allowed, seed_column


# These are exactly the full-sweep rules of codex/work_C8/sweeping_clock_checks.py.
# Keeping the small combinatorial definition here avoids a NumPy/SciPy runtime
# requirement for querying this Hamiltonian's local matrix elements.
ALPHABET = ("L", "0", "1", "C0", "C1", "R0", "R1", "I", "Z",
            "A0", "A1", "D0", "D1", "E0", "E1")
CLOCK_INDEX = {letter: index for index, letter in enumerate(ALPHABET)}
COUNTER = ("C0", "C1", "R0", "R1")
SWEEP = ("D0", "D1", "E0", "E1")
HEADS = set(COUNTER + SWEEP + ("I",))
BITS = ("0", "1")
PAIRS = {("L", "L"), ("I", "Z"), ("Z", "Z")}
PAIRS |= {("L", h) for h in COUNTER + ("I",)}
PAIRS |= {(h, b) for h in COUNTER + SWEEP for b in BITS}
PAIRS |= {(a, b) for a in BITS for b in BITS}
PAIRS |= {("A" + a, "A" + b) for a in BITS for b in BITS}
PAIRS |= {("A" + a, h) for a in BITS for h in SWEEP}
LEFT = {"L", "A0", "A1"} | HEADS
RIGHT = set(BITS) | HEADS | {"Z"}
RULES = []
for b in BITS:
    RULES.append((("C1", b), ("L", "C" + b), "carry"))
    RULES.append((("L", "R" + b), ("R0", b), "return"))
    for a in BITS:
        RULES.append((("D" + a, b), ("A" + a, "D" + b), "out"))
        RULES.append((("A" + a, "E" + b), ("E" + a, b), "back"))
RULES.append((("I", "Z"), ("L", "I"), "initialize"))

STATIC_PAIRS = {("LEFT", "OUT_L"), ("OUT_L", "BULK"),
                ("OUT_L", "OUT_R"), ("BULK", "BULK"),
                ("BULK", "OUT_R"), ("OUT_R", "RIGHT")}
LOCAL_DIM = len(ALPHABET) * WORK_DIM


def encode(clock, work):
    if clock not in CLOCK_INDEX or not 0 <= work < WORK_DIM:
        raise ValueError("Invalid clock/work basis label")
    return CLOCK_INDEX[clock] * WORK_DIM + work


def decode(index):
    if not 0 <= index < LOCAL_DIM:
        raise ValueError("Core basis label out of range")
    clock, work = divmod(index, WORK_DIM)
    return ALPHABET[clock], work


def cleaned(column):
    return {row: value for row, value in column.items() if value != 0}


def onsite_column(index):
    """Clock turn projector plus initialization pin, before bond splitting."""
    clock, work = decode(index)
    result = defaultdict(complex)
    if clock in ("C0", "R1"):
        target = "R1" if clock == "C0" else "C0"
        result[index] += 0.5
        result[encode(target, work)] -= 0.5
    if clock == "I" and not initial_allowed(work):
        result[index] += 1
    return cleaned(result)


def pair_unitary(label, left, right):
    """The work unitary on a directed bulk clock edge; all are involutions."""
    if label == "out":
        return g_column(left, right)
    if label == "back":
        return {(left, target): coefficient for target, coefficient in f_column(right).items()}
    if label == "initialize":
        return {(target, right): coefficient for target, coefficient in seed_column(left).items()}
    return {(left, right): 1+0j}


def bond_column(left, right):
    """Column of h_c on C^194400 tensor C^194400, for arbitrary basis inputs.

    Each propagation term is 1/2 [[I,-U†],[-U,I]] in its clock-edge block.
    Clock/static defects have weight one. Onsite terms have half their weight
    on each adjacent bond; boundary_column supplies the missing endpoint half.
    """
    ca, wa = decode(left)
    cb, wb = decode(right)
    a, _, _, _ = work_decode(wa)
    b, _, _, _ = work_decode(wb)
    source = (left, right)
    result = defaultdict(complex)
    result[source] += int((ca, cb) not in PAIRS)
    result[source] += int((a.static, b.static) not in STATIC_PAIRS)
    for target, coefficient in onsite_column(left).items():
        result[(target, right)] += 0.5 * coefficient
    for target, coefficient in onsite_column(right).items():
        result[(left, target)] += 0.5 * coefficient
    for before, after, label in RULES:
        if (ca, cb) == before:
            target_clock = after
        elif (ca, cb) == after:
            target_clock = before
        else:
            continue
        result[source] += 0.5
        # U†=U for the complete G, F and seed matrices, on every work input.
        for (ta, tb), coefficient in pair_unitary(label, wa, wb).items():
            target = (encode(target_clock[0], ta), encode(target_clock[1], tb))
            result[target] -= 0.5 * coefficient
    return cleaned(result)


def boundary_column(side, index):
    """The left (a_c) or right (b_c) one-site term on its full Hilbert space."""
    if side not in ("left", "right"):
        raise ValueError(side)
    clock, work = decode(index)
    symbol, _, _, _ = work_decode(work)
    result = defaultdict(complex)
    allowed_clock = LEFT if side == "left" else RIGHT
    required_static = "LEFT" if side == "left" else "RIGHT"
    result[index] += int(clock not in allowed_clock)
    result[index] += int(symbol.static != required_static)
    for target, coefficient in onsite_column(index).items():
        result[target] += 0.5 * coefficient
    if side == "left":
        rules = [("R0", "D0", False), ("R1", "D1", False),
                 ("E0", "C0", True), ("E1", "C1", True)]
    else:
        rules = [("I", "R0", False), ("D0", "E0", False), ("D1", "E1", False)]
    for before, after, apply_f in rules:
        if clock == before:
            target_clock = after
        elif clock == after:
            target_clock = before
        else:
            continue
        result[index] += 0.5
        column = f_column(work) if apply_f else {work: 1+0j}
        for target, coefficient in column.items():
            result[encode(target_clock, target)] -= 0.5 * coefficient
    return cleaned(result)


def left_column(index):
    return boundary_column("left", index)


def right_column(index):
    return boundary_column("right", index)


def clock_defects(word):
    return ((word[0] not in LEFT),
            *((a, b) not in PAIRS for a, b in zip(word, word[1:])),
            (word[-1] not in RIGHT))


def clock_forward(word):
    result = []
    for j in range(len(word) - 1):
        for before, after, label in RULES:
            if word[j:j+2] == before:
                result.append((word[:j] + after + word[j+2:], label, j))
    for j, clock in enumerate(word):
        if clock == "C0":
            result.append((word[:j] + ("R1",) + word[j+1:], "turn", j))
    if word[0] in ("R0", "R1"):
        result.append((("D" + word[0][1],) + word[1:], "begin", 0))
    if word[0] in ("E0", "E1"):
        result.append((("C" + word[0][1],) + word[1:], "finish", 0))
    if word[-1] == "I":
        result.append((word[:-1] + ("R0",), "activate", len(word) - 1))
    if word[-1] in ("D0", "D1"):
        result.append((word[:-1] + ("E" + word[-1][1],), "bounce", len(word) - 1))
    return result


def clock_path(n):
    word = ("I",) + ("Z",) * (n - 1)
    words, steps = [word], []
    seen = {word}
    while True:
        candidates = clock_forward(word)
        assert len(candidates) <= 1
        if not candidates:
            break
        word, label, site = candidates[0]
        assert word not in seen and not any(clock_defects(word))
        words.append(word)
        steps.append((label, site))
        seen.add(word)
    assert len(words) == (2*n + 4)*2**n + n - 4
    return words, steps


def clock_work_column(label, site, work_word):
    """Unitary column on the full work tensor product at one legal clock edge."""
    if label in ("out", "back", "initialize"):
        result = {}
        for (a, b), coefficient in pair_unitary(label, work_word[site], work_word[site+1]).items():
            result[work_word[:site] + (a, b) + work_word[site+2:]] = coefficient
        return result
    if label == "finish":
        return { (target,) + work_word[1:]: coefficient
                 for target, coefficient in f_column(work_word[0]).items() }
    return {work_word: 1+0j}


def initial_work_word(n, quantum_bits):
    """One accepted t=0 input; quantum_bits is an arbitrary basis bit string."""
    if n < 4 or len(quantum_bits) != n:
        raise ValueError("The demonstration grammar requires n>=4")
    config = WORK.machine.initialize(n)
    return tuple(WORK.encode(symbol, bit) for symbol, bit in zip(config.tape, quantum_bits))
