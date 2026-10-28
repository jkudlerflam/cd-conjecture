"""Adversarial check of the initialized input-format query interface.

This checks finite local permutations and the actual helper walk, not a
fully instantiated reversible scheduler or the macroscopic conjecture.
Run with ordinary Python 3; no dependencies are required.
"""

from collections import defaultdict
from itertools import product
import json
from pathlib import Path


def alphabet(q, m):
    return (
        [("D", s) for s in range(q)]
        + [("blank",), ("bit", 0), ("bit", 1)]
        + [("H", j, c) for j in range(m) for c in range(2)]
        + [("T", j, c, s) for j in range(m) for c in range(2) for s in range(q)]
    )


def charge(symbol, q):
    out = [0] * q
    if symbol[0] == "D":
        out[symbol[1]] = 1
    elif symbol[0] == "T":
        out[symbol[3]] = 1
    return tuple(out)


def query_value(symbol, kind, bit):
    if kind == "species":
        return (symbol[1] >> bit) & 1 if symbol[0] == "D" else 0
    types = {"D": 0, "blank": 1, "bit": 2, "H": 3, "T": 4}
    return (types[symbol[0]] >> bit) & 1


def signature(symbol, q):
    return tuple(query_value(symbol, "type", j) for j in range(3)) + tuple(
        query_value(symbol, "species", j) for j in range((q - 1).bit_length())
    )


def local_query(pair, kind, type_bit=None):
    left, right = pair
    if right[0] != "H":
        return pair
    _, index, value = right
    bit = index if kind == "species" else type_bit
    return left, ("H", index, value ^ query_value(left, kind, bit))


def shift(symbol, delta, m):
    if symbol[0] == "H":
        return ("H", (symbol[1] + delta) % m, symbol[2])
    if symbol[0] == "T":
        return ("T", (symbol[1] + delta) % m, *symbol[2:])
    return symbol


def blank_bit_conversion(symbol):
    if symbol == ("blank",):
        return ("bit", 0)
    if symbol == ("bit", 0):
        return ("blank",)
    return symbol


def local_checks(q, m):
    symbols = alphabet(q, m)
    classes = defaultdict(list)
    for value in symbols:
        classes[signature(value, q)].append(value)
    # Concrete failure of an unqualified "arbitrary prescribed state" format.
    assert signature(("bit", 0), q) == signature(("bit", 1), q)
    assert signature(("H", 0, 0), q) == signature(("H", m - 1, 1), q)
    # Repair: every permitted untrusted single-symbol reference is distinguishable.
    for value in [("blank",)] + [("D", s) for s in range(q)]:
        assert classes[signature(value, q)] == [value]
    allowed_occupation = {("D", 0), ("D", q - 1)}
    accept_signatures = {signature(value, q) for value in allowed_occupation}
    assert {x for x in symbols if signature(x, q) in accept_signatures} == allowed_occupation
    pair_count = 0
    for kind, bit in [("species", None)] + [("type", j) for j in range(3)]:
        images = set()
        for pair in product(symbols, repeat=2):
            target = local_query(pair, kind, bit)
            assert local_query(target, kind, bit) == pair
            assert tuple(a + b for a, b in zip(charge(pair[0], q), charge(pair[1], q))) == tuple(
                a + b for a, b in zip(charge(target[0], q), charge(target[1], q))
            )
            images.add(target)
            pair_count += 1
        assert len(images) == len(symbols) ** 2
    for value in symbols:
        assert shift(shift(value, 1, m), -1, m) == value
        assert blank_bit_conversion(blank_bit_conversion(value)) == value
        assert charge(value, q) == charge(blank_bit_conversion(value), q)
    return {
        "q": q,
        "m": m,
        "alphabet_size": len(symbols),
        "whole_space_query_columns_checked": pair_count,
        "bit_zero_one_alias_demonstrated": True,
        "all_blank_and_species_references_separated": True,
        "blank_bit_zero_conversion_is_charge_preserving_involution": True,
    }


def helper_round_trip(word, target, kind, bit, q, m, scratch):
    initial = tuple(word)
    work = list(word)
    helper = len(work) - 1
    assert work[helper] == ("H", 0, 0)
    swaps = []
    while helper > target + 1:
        work[helper - 1], work[helper] = work[helper], work[helper - 1]
        swaps.append(helper - 1)
        helper -= 1
    index = bit if kind == "species" else 0
    for _ in range(index):
        work[helper] = shift(work[helper], 1, m)
    work[target], work[helper] = local_query((work[target], work[helper]), kind, bit)
    scratch ^= work[helper][2]
    work[target], work[helper] = local_query((work[target], work[helper]), kind, bit)
    for _ in range(index):
        work[helper] = shift(work[helper], -1, m)
    for left in reversed(swaps):
        work[left], work[left + 1] = work[left + 1], work[left]
    assert tuple(work) == initial
    return scratch


def walk_checks(q, m, n):
    symbols = alphabet(q, m)
    queries = [("type", j) for j in range(3)] + [
        ("species", j) for j in range((q - 1).bit_length())
    ]
    count = 0
    for values in product(symbols, repeat=n - 1):
        word = values + (("H", 0, 0),)
        for target in range(n - 1):
            for kind, bit in queries:
                for initial_scratch in (0, 1):
                    result = helper_round_trip(word, target, kind, bit, q, m, initial_scratch)
                    assert result == initial_scratch ^ query_value(word[target], kind, bit)
                    count += 1
    return {"q": q, "m": m, "n": n, "basis_round_trips": count}


def product_format_checks(q, m):
    symbols = alphabet(q, m)
    # Two independent logical occupations, one charged reference and one blank.
    prescribed = [
        {("D", 0), ("D", 1)},
        {("D", 1), ("D", q - 1)},
        {("D", 0)},
        {("blank",)},
    ]
    counts = []
    for allowed in prescribed:
        sigs = {signature(x, q) for x in allowed}
        accepted = [x for x in symbols if signature(x, q) in sigs]
        assert set(accepted) == allowed
        counts.append(len(accepted))
    dimension = 1
    for count in counts:
        dimension *= count
    assert dimension == 4
    return {"accepted_slot_dimensions": counts, "accepted_space_dimension": dimension}


def main():
    report = {
        "status": "passed",
        "scope": "finite gate and input-format audit, not a full controller or CFT verification",
        "local_gates": [local_checks(q, m) for q, m in [(2, 2), (3, 4), (5, 6)]],
        "helper_walks": [walk_checks(2, 2, 3), walk_checks(3, 4, 3)],
        "repaired_product_format": product_format_checks(5, 6),
    }
    path = Path(__file__).with_name("format_query_results.json")
    path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
