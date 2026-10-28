"""Enumerate ALL clock/work states of a two-site initialized toy core.

The quantum payload is a charged two-valued record. The work update is
the direct-sum single-head SWAP sweep, with reversible endpoint bounce
collisions. Because this example has permutation gates, its exact zero
eigenspace can be counted by graph components and diagonal penalties;
no floating point diagonalization or preferred-history restriction is used.

This is a test of the interface and initialization mechanism, not the
universal preparation controller or a high-channel code.
"""

from itertools import product
import json
from pathlib import Path


CLOCK = ("L", "0", "1", "C0", "C1", "R0", "R1", "I", "Z", "A0", "A1", "D0", "D1", "E0", "E1")
COUNTER = ("C0", "C1", "R0", "R1")
SWEEPS = ("D0", "D1", "E0", "E1")
BITS = ("0", "1")
HEADS = set(COUNTER + SWEEPS + ("I",))
PAIRS = {("L", "L"), ("I", "Z"), ("Z", "Z")}
PAIRS |= {("L", h) for h in COUNTER + ("I",)}
PAIRS |= {(h, b) for h in COUNTER + SWEEPS for b in BITS}
PAIRS |= set(product(BITS, repeat=2))
PAIRS |= {("A" + a, "A" + b) for a in BITS for b in BITS}
PAIRS |= {("A" + a, h) for a in BITS for h in SWEEPS}
LEFT = {"L", "A0", "A1"} | HEADS
RIGHT = set(BITS) | HEADS | {"Z"}
RULES = [(("I", "Z"), ("L", "I"), "init")]
for b in BITS:
    RULES.extend([(("C1", b), ("L", "C" + b), "carry"), (("L", "R" + b), ("R0", b), "return")])
    for a in BITS:
        RULES.extend([(("D" + a, b), ("A" + a, "D" + b), "out"), (("A" + a, "E" + b), ("E" + a, b), "back")])


WORK = [(a,) for a in range(2)] + [(phase, direction, a) for phase in range(2) for direction in ("L", "R") for a in range(2)]


def g(pair):
    left, right = pair
    if len(left) == 3 and len(right) == 1:
        phase, direction, a = left
        b = right[0]
        if (phase, direction) == (0, "R"):
            return (b,), (1, "R", a)
        if (phase, direction) == (1, "L"):
            return (b,), (0, "L", a)
    if len(left) == 1 and len(right) == 3:
        a = left[0]
        phase, direction, b = right
        if (phase, direction) == (0, "L"):
            return (1, "L", b), (a,)
        if (phase, direction) == (1, "R"):
            return (0, "R", b), (a,)
    return pair


def f(symbol):
    if len(symbol) == 1:
        return symbol
    phase, direction, a = symbol
    return 1 - phase, ("R" if direction == "L" else "L"), a


def seed(symbol):
    if len(symbol) == 1:
        return (0, "L", symbol[0])
    if symbol[:2] == (0, "L"):
        return (symbol[2],)
    return symbol


def defects(clock):
    return (clock[0] not in LEFT, clock not in PAIRS, clock[1] not in RIGHT)


def clock_moves(clock):
    result = []
    for before, after, label in RULES:
        if clock == before:
            result.append((after, label, 0))
    for j, symbol in enumerate(clock):
        if symbol == "C0":
            after = list(clock)
            after[j] = "R1"
            result.append((tuple(after), "turn", j))
    if clock[0] in ("R0", "R1"):
        result.append((("D" + clock[0][1], clock[1]), "begin", 0))
    if clock[0] in ("E0", "E1"):
        result.append((("C" + clock[0][1], clock[1]), "finish", 0))
    if clock[1] == "I":
        result.append(((clock[0], "R0"), "activate", 1))
    if clock[1] in ("D0", "D1"):
        result.append(((clock[0], "E" + clock[1][1]), "bounce", 1))
    return result


def work_move(work, label):
    if label == "out":
        return g(work)
    if label == "back":
        return work[0], f(work[1])
    if label == "finish":
        return f(work[0]), work[1]
    if label == "activate":
        return work[0], seed(work[1])
    return work


def pin_penalty(clock, work):
    penalty = 0
    for j in range(2):
        if clock[j] == "I":
            penalty += len(work[j]) != 1
            # First input record is fixed to D_0. Second is a free bit.
            if j == 0:
                penalty += work[j][-1] != 0
    return penalty


def main():
    for value in WORK:
        assert f(f(value)) == seed(seed(value)) == value
        assert f(value)[-1] == seed(value)[-1] == value[-1]
    pairs = list(product(WORK, repeat=2))
    for pair in pairs:
        assert g(g(pair)) == pair
        assert sorted(x[-1] for x in g(pair)) == sorted(x[-1] for x in pair)
    assert len({g(pair) for pair in pairs}) == len(pairs)

    clocks = list(product(CLOCK, repeat=2))
    states = [(clock, work) for clock in clocks for work in pairs]
    index = {state: j for j, state in enumerate(states)}
    parent = list(range(len(states)))
    sizes = [1] * len(states)

    def root(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def unite(a, b):
        a, b = root(a), root(b)
        if a == b:
            return
        if sizes[a] < sizes[b]:
            a, b = b, a
        parent[b] = a
        sizes[a] += sizes[b]

    transitions = 0
    for j, (clock, work) in enumerate(states):
        for target_clock, label, site in clock_moves(clock):
            assert defects(target_clock) == defects(clock)
            target = (target_clock, work_move(work, label))
            unite(j, index[target])
            transitions += 1
    components = {}
    for j, (clock, work) in enumerate(states):
        component = components.setdefault(root(j), {"size": 0, "penalty": 0, "bad_clock": False})
        component["size"] += 1
        component["penalty"] += sum(defects(clock)) + pin_penalty(clock, work)
        component["bad_clock"] |= any(defects(clock))
    zero = [c for c in components.values() if not c["penalty"]]
    assert len(zero) == 2
    assert all(c["size"] == 30 and not c["bad_clock"] for c in zero)
    # All legal-clock histories exist before initialization; the pins leave two.
    legal_components = [c for c in components.values() if not c["bad_clock"]]
    assert len(legal_components) == 100
    assert all(c["size"] == 30 for c in legal_components)

    report = {
        "status": "passed",
        "scope": "all-sector two-site permutation-history toy core; not the universal controller",
        "physical_local_dimension": len(CLOCK) * len(WORK),
        "full_hilbert_space_dimension": len(states),
        "decorated_propagation_edges": transitions,
        "clock_work_components": len(components),
        "legal_clock_history_components_before_pins": len(legal_components),
        "unpinned_components_after_all_pins": len(zero),
        "exact_ground_space_dimension": len(zero),
        "vertices_in_each_surviving_history": 30,
        "all_zero_one_two_head_and_illegal_clock_states_included": True,
    }
    Path(__file__).with_name("tiny_history_results.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
