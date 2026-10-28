"""Exact finite rule checks for positive completion and readout absorption.

Standard library only. These are local grammar checks, not a CFT test.
"""

from itertools import product
from pathlib import Path
import json


ARROWS = {"#": (1, 0), "L": (0, 0), "X": (0, 1), "R": (1, 1)}


def bad(a, b):
    return int(ARROWS[a][1] != ARROWS[b][0])


def cells(word):
    starts = [j for j, a in enumerate(word) if a == "#"]
    assert starts[0] == 0
    return [word[j:k] for j, k in zip(starts, starts[1:]+[len(word)])]


def cell_bad(word):
    return sum(bad(a, b) for a, b in zip(word, word[1:]+("#",)))


def token(letter):
    return int(letter[1] is not None)


def incoming(letter):
    z, s, _ = letter
    if z == "B":
        return "R" if s is None else "L"
    return "R" if z == "R" else "L"


def outgoing(letter):
    z, s, _ = letter
    if z == "A":
        return "L" if s is None else "R"
    return "L" if z == "L" else "R"


def main():
    selector_cases = 0
    for before, after in [("LX", "XR"), ("#L", "R#")]:
        for a, b in product(ARROWS, repeat=2):
            w, z = a+before+b, a+after+b
            assert [bad(w[i], w[i+1]) for i in range(3)] == [
                bad(z[i], z[i+1]) for i in range(3)]
            selector_cases += 1

    # New pin support is exactly the two illegal selector pairs entering #.
    old = {("X", "#"), ("R", "#")}
    new = {(a, "#") for a in ARROWS}
    assert new-old == {("L", "#"), ("#", "#")}
    assert all(bad(*pair) for pair in new-old)

    fragment_cases = 0
    for length in range(1, 9):
        for tail in product("LXR", repeat=length-1):
            word = ("#",)+tail
            z = cell_bad(word)
            internal = [j for j in range(length-1) if bad(word[j], word[j+1])]
            cuts = [0]+[j+1 for j in internal]+[length]
            fragments = [word[i:j] for i, j in zip(cuts, cuts[1:])]
            assert len(fragments) <= z+1
            assert all(fragment.count("X") <= 1 for fragment in fragments)
            if z:
                # t=3,U=1 gives C=10 and K=21. The general algebra is exact.
                assert 21*z-10*len(fragments) >= 1
            fragment_cases += 1

    marker_cases = 0
    for l1, l2 in product(range(1, 5), repeat=2):
        for tail in product("LXR", repeat=l1+l2-2):
            word = ("#",)+tail[:l1-1]+("#",)+tail[l1-1:]
            before = [cell_bad(c) for c in cells(word)]
            length = len(word)
            for j in range(length):
                k = (j+1) % length
                pair = word[j]+word[k]
                if pair not in ("#L", "R#"):
                    continue
                out = list(word)
                out[j], out[k] = ("R", "#") if pair == "#L" else ("#", "L")
                root = 0
                if j == 0 and pair == "#L":
                    root = 1
                elif k == 0 and pair == "R#":
                    root = length-1
                out = tuple(out[root:]+out[:root])
                after = [cell_bad(c) for c in cells(out)]
                assert before == after
                marker_cases += 1

    species = ["v"]+[(r, s) for r in (1, 2) for s in (0, 1, -1)]
    excited = species[1:]
    rewrites = []
    for b, c in product((0, 1), repeat=2):
        for s in species:
            rewrites.append(((("d", s, b), ("R", None, c)),
                             (("L", None, b), ("d", s, c))))
        for s in excited:
            rewrites.append(((("x", s, b), ("R", None, c)),
                             (("L", None, b), ("x", s, c))))
    for b in (0, 1):
        for s in species:
            rewrites.append(((("A", s, None), ("R", None, b)),
                             (("A", None, None), ("d", s, b))))
        rewrites.append(((("d", "v", b), ("B", None, None)),
                         (("L", None, b), ("B", "v", None))))
        for s in excited:
            rewrites.append(((("x", s, b), ("B", None, None)),
                             (("L", None, b), ("B", s, None))))
            for am in [None]+species:
                rewrites.append(((("A", am, None), ("d", s, b)),
                                 (("A", am, None), ("x", s, b))))
    for initial, final in rewrites:
        assert sum(map(token, initial)) == sum(map(token, final))
        assert incoming(initial[0]) == incoming(final[0])
        assert outgoing(initial[1]) == outgoing(final[1])
    port_cases = 0
    for endpoint, bus, s in product(("A", "B"), (1, 2), (0, 1, -1)):
        vacuum, excited_letter = (endpoint, "v", None), (endpoint, (bus, s), None)
        assert token(vacuum) == token(excited_letter) == 1
        assert incoming(vacuum) == incoming(excited_letter)
        assert outgoing(vacuum) == outgoing(excited_letter)
        port_cases += 1

    result = dict(status="all exact local checks passed",
                  selector_outer_neighbor_cases=selector_cases,
                  selector_cell_fragment_cases=fragment_cases,
                  marker_reanchoring_cases=marker_cases,
                  core_pair_rewrites=len(rewrites), endpoint_port_cases=port_cases,
                  scope="Positive pin completion, selector fragmentation and token invariance")
    Path(__file__).with_name("illegal_checks.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
