"""Small controls of C58's new lemmas, not a Form W simulation."""
import itertools
import json
import math
from pathlib import Path
import numpy as np


def logz(h, beta):
    vals = -beta * np.linalg.eigvalsh(h)
    m = float(vals.max())
    return m + math.log(float(np.exp(vals - m).sum()))


def pin_checks():
    hop = 0.31 + 0.19j
    density, field = 0.23, -0.17
    bond = np.diag([0.0, field / 2, field / 2, field + density]).astype(complex)
    bond[1, 2], bond[2, 1] = hop, hop.conjugate()
    jnorm = float(np.linalg.norm(bond, 2))
    count, residual, slack = 0, 0.0, float("inf")
    for length in range(3, 8):
        for markers in range(1, length + 1):
            basis = [x for x in range(2**length) if x.bit_count() == markers]
            index = {x: j for j, x in enumerate(basis)}
            h = np.zeros((len(basis), len(basis)), complex)
            for col, x in enumerate(basis):
                for i in range(length):
                    k = (i + 1) % length
                    a, b = (x >> i) & 1, (x >> k) & 1
                    h[col, col] += density * a * b + field * (a + b) / 2
                    if a != b:
                        y = x ^ (1 << i) ^ (1 << k)
                        h[index[y], col] += hop if a == 1 else hop.conjugate()
            ps = [np.diag([(x >> i) & 1 for x in basis]) for i in range(length)]
            comp = sum(p @ h @ p for p in ps) / markers
            commsum = sum(p @ (p @ h - h @ p) - (p @ h - h @ p) @ p for p in ps)
            residual = max(residual, float(np.linalg.norm(comp - h + commsum / (2 * markers))))
            assert np.linalg.norm(comp - h, 2) <= 4 * jnorm * length / markers + 1e-12
            pinned = [i for i, x in enumerate(basis) if x & 1]
            hp = h[np.ix_(pinned, pinned)]
            for beta in (0.2, 1.0, 4.0):
                z, zp = logz(h, beta), logz(hp, beta)
                lower = z - zp - math.log(length / markers)
                upper = math.log(length) + 4 * beta * jnorm * length / markers + zp - z
                slack = min(slack, lower, upper)
                assert lower >= -1e-11 and upper >= -1e-11
                count += 1
    return {"trace_tests": count, "compression_residual": residual, "minimum_slack": slack}


ARROW = {"#": (1, 0), "L": (0, 0), "X": (0, 1), "R": (1, 1)}


def cell_counts(word, positions):
    length = len(word)
    counts = []
    for j, start in enumerate(positions):
        end = positions[(j + 1) % len(positions)]
        size = (end - start) % length or length
        counts.append(sum(ARROW[word[(start + k) % length]][1]
                          != ARROW[word[(start + k + 1) % length]][0]
                          for k in range(size)))
    return counts


def selector_checks():
    count = 0
    rewrites = {"LX": "XR", "XR": "LX", "R#": "#L", "#L": "R#"}
    for length in range(2, 7):
        for word in itertools.product("#LXR", repeat=length):
            positions = [i for i, c in enumerate(word) if c == "#"]
            if not positions:
                continue
            before = cell_counts(word, positions)
            for i in range(length):
                k = (i + 1) % length
                pair = word[i] + word[k]
                if pair not in rewrites:
                    continue
                new = list(word)
                new[i], new[k] = rewrites[pair]
                updated = positions.copy()
                if pair == "R#":
                    updated[positions.index(k)] = i
                elif pair == "#L":
                    updated[positions.index(i)] = k
                assert cell_counts(new, updated) == before
                count += 1
    return {"literal_selector_moves": count}


def kron_all(factors):
    out = np.ones((1, 1))
    for fac in factors:
        out = np.kron(out, fac)
    return out


def defect_checks():
    cells, b, penalty, slope, pin, coupling = 3, 2, 2.0, 0.4, 0.7, 0.005
    ident = np.eye(4)
    number = np.kron(np.diag([0.0, 1.0]), np.eye(2))
    helix_pin = np.kron(np.eye(2), np.diag([0.0, 1.0]))
    helix_x = np.kron(np.eye(2), np.array([[0.0, 1.0], [1.0, 0.0]]))
    lower = np.kron(np.array([[0.0, 1.0], [0.0, 0.0]]), np.eye(2))

    def at(op, site):
        return kron_all([op if j == site else ident for j in range(cells)])

    lengths = [b * np.eye(4**cells) + at(number, j) for j in range(cells)]
    pins = [at(helix_pin, j) for j in range(cells)]
    weak = sum(coupling * at(helix_x, j) for j in range(cells))
    for j in range(cells):
        edge = at(lower.T, j) @ at(lower, (j + 1) % cells)
        weak += coupling * (edge + edge.T)
    incidence = 3 * coupling
    epsilon = 0.02
    delta = penalty / 2 - epsilon / 2 - 4 * incidence
    amplitude = 1.0  # d_P=2, z_good^- <= 4 exp(beta epsilon/2)
    tests, slack = 0, float("inf")
    for tilt in (-epsilon, 0.0, epsilon):
        records = []
        for flags in itertools.product((0, 1), repeat=cells):
            ho = weak.copy() + tilt * sum(lengths)
            hc = weak.copy() + tilt * sum(lengths)
            for j in range(cells):
                if flags[j]:
                    bad = penalty * np.eye(4**cells) + slope * lengths[j]
                    ho += bad
                    hc += bad
                hc += pin * pins[j]
                if not flags[(j - 1) % cells]:
                    ho += pin * pins[j]
            records.append((flags, ho, hc))
        for beta in (2.0, 4.0, 8.0):
            legal = logz(records[0][1], beta)
            lo = np.array([logz(item[1], beta) for item in records])
            lc = np.array([logz(item[2], beta) for item in records])
            total_o = float(np.max(lo) + np.log(np.exp(lo - np.max(lo)).sum()))
            total_c = float(np.max(lc) + np.log(np.exp(lc - np.max(lc)).sum()))
            bound = cells * math.log1p(amplitude * math.exp(-beta * delta))
            assert legal <= total_c + 1e-11 <= total_o + 2e-11
            assert total_o - legal <= bound + 1e-11
            for value, (flags, _, _) in zip(lo, records):
                if sum(flags):
                    sector_bound = sum(flags) * (math.log(amplitude) - beta * delta)
                    assert value - legal <= sector_bound + 1e-11
            slack = min(slack, bound - (total_o - legal))
            tests += 1
    return {"noncommuting_comparisons": tests, "minimum_slack": slack,
            "matrix_dimension_per_flag_sector": 4**cells, "delta_bound": delta}


if __name__ == "__main__":
    results = {"pinning": pin_checks(), "selector": selector_checks(), "defects": defect_checks()}
    Path(__file__).with_name("checks.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))
