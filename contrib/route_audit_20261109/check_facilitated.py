"""Exact periodic second-order identities for the proposed transverse-DM repair.

Standard library only; all amplitudes use fractions. Does not assert c=2.
"""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import json


def occ(s, j, L):
    return (s >> (j % L)) & 1


def adjacent(s, L):
    return sum(occ(s, j, L) * occ(s, j + 1, L) for j in range(L))


def transitions(s, L, kind):
    out = defaultdict(int)
    if kind == "pair":
        for j in range(L):
            k = (j + 1) % L
            if occ(s, j, L) == occ(s, k, L):
                out[s ^ (1 << j) ^ (1 << k)] += 1
    elif kind == "facilitation":
        for j in range(L):
            out[s ^ (1 << j)] += occ(s, j - 1, L) - occ(s, j + 1, L)
    return {q: a for q, a in out.items() if a}


def second_order(s, L, first, second):
    out = defaultdict(F)
    for q, a in transitions(s, L, first).items():
        penalty = adjacent(q, L)
        assert penalty > 0
        for f, b in transitions(q, L, second).items():
            if adjacent(f, L) == 0:
                out[f] -= F(a * b, penalty)
    return {f: a for f, a in out.items() if a}


def predicted_mixed(s, L):
    out = defaultdict(F)
    for j in range(L):
        f = s ^ (1 << j)
        if adjacent(f, L) == 0:
            out[f] -= occ(s, j - 2, L) - occ(s, j + 2, L)
    return {f: a for f, a in out.items() if a}


def predicted_squared(s, L):
    out = defaultdict(F)
    out[s] += -2 * s.bit_count() + 2 * sum(
        occ(s, j, L) * occ(s, j + 2, L) for j in range(L)
    )
    for j in range(L):
        k = (j + 1) % L
        if occ(s, j, L) != occ(s, k, L):
            f = s ^ (1 << j) ^ (1 << k)
            if adjacent(f, L) == 0:
                out[f] += 1
    return {f: a for f, a in out.items() if a}


def main():
    cases = []
    for L in (7, 8, 9, 10):
        basis = [s for s in range(1 << L) if adjacent(s, L) == 0]
        entries = 0
        for s in basis:
            a = second_order(s, L, "pair", "facilitation")
            b = second_order(s, L, "facilitation", "pair")
            mixed = defaultdict(F, a)
            for f, amp in b.items():
                mixed[f] += amp
            mixed = {f: amp for f, amp in mixed.items() if amp}
            assert mixed == predicted_mixed(s, L)
            squared = second_order(s, L, "facilitation", "facilitation")
            assert squared == predicted_squared(s, L)
            entries += len(mixed) + len(squared)
        cases.append(dict(L=L, legal_states=len(basis), nonzero_verified_entries=entries))

    out = dict(
        status="all exact identities passed",
        conventions="H0=V sum nn, pair amplitude +Delta, facilitation lambda sum(nX-Xn)",
        mixed_coefficient="-Delta*lambda/V times projected distance-two antisymmetric facilitation",
        squared_coefficient="lambda^2/V times projected [NN hop+h.c.-2N+2 sum n_j n_(j+2)]",
        cases=cases,
        scope="second-order legal effective operator, including charge-changing terms; no infrared claim",
    )
    text = json.dumps(out, indent=2) + "\n"
    (Path(__file__).parent / "facilitated_checks.json").write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
