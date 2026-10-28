"""Exhaust the trusted-helper round trip, including malformed H/T records."""

import itertools
import json
from pathlib import Path


def alphabet(q):
    width = (q - 1).bit_length()
    values = [("D", s) for s in range(q)]
    values += [("blank",), ("bit", 0), ("bit", 1)]
    values += [("H", j, c) for j in range(width) for c in range(2)]
    values += [
        ("T", j, c, s)
        for j in range(width)
        for c in range(2)
        for s in range(q)
    ]
    return values, width


def charges(word, q):
    result = [0] * q
    for symbol in word:
        if symbol[0] == "D":
            result[symbol[1]] += 1
        elif symbol[0] == "T":
            result[symbol[3]] += 1
    return result


def query_value(symbol, kind, index):
    if kind == "species":
        return ((symbol[1] >> index) & 1) if symbol[0] == "D" else 0
    type_code = {"D": 0, "blank": 1, "bit": 2, "H": 3, "T": 4}
    return (type_code[symbol[0]] >> index) & 1


def round_trip(source, target, kind, index, q):
    word = list(source)
    original_charge = charges(word, q)
    helper_position = len(word) - 1
    assert word[helper_position] == ("H", 0, 0)
    path = []
    while helper_position > target + 1:
        left = helper_position - 1
        word[left], word[helper_position] = word[helper_position], word[left]
        path.append((left, helper_position))
        helper_position = left
    assert charges(word, q) == original_charge
    helper_index = index if kind == "species" else 0
    word[helper_position] = ("H", helper_index, 0)
    value = query_value(word[target], kind, index)
    word[helper_position] = ("H", helper_index, value)
    scratch = word[helper_position][2]  # One exact copy into trusted zero bit.
    word[helper_position] = (
        "H",
        helper_index,
        value ^ query_value(word[target], kind, index),
    )
    assert word[helper_position][2] == 0
    word[helper_position] = ("H", 0, 0)
    for left, right in reversed(path):
        word[left], word[right] = word[right], word[left]
    assert tuple(word) == source
    assert charges(word, q) == original_charge
    assert scratch == query_value(source[target], kind, index)


def main():
    cases = 0
    malformed_cases = 0
    for q in (2, 3, 4):
        values, width = alphabet(q)
        for length in ((2, 3, 4) if q == 2 else (2, 3)):
            for untrusted in itertools.product(values, repeat=length - 1):
                source = untrusted + (("H", 0, 0),)
                for target in range(length - 1):
                    queries = [("species", bit) for bit in range(width)]
                    queries += [("type", bit) for bit in range(3)]
                    for kind, bit in queries:
                        round_trip(source, target, kind, bit, q)
                        cases += 1
                        if source[target][0] in ("H", "T"):
                            malformed_cases += 1
    report = {
        "basis_query_round_trips": cases,
        "round_trips_targeting_H_or_T_records": malformed_cases,
        "charge_checks": "exact",
        "restoration_checks": "exact",
        "status": "passed",
    }
    Path(__file__).with_name("helper_verification.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
