"""Recheck the published qubit fits and local rotation using saved data only."""
from pathlib import Path
import cmath
import hashlib
import json
import math
import re

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "d2_xy_yx_chirality_imps_20261118/tables/open_chain_ladders.txt"
source_bytes = SOURCE.read_bytes()
runs = {}
label = None
for line in source_bytes.decode().splitlines():
    if line.startswith("== "):
        label = line[3:].split()[0]
        runs[label] = []
    match = re.search(r"L=\s*(\d+).* Sbar=([-+0-9.e]+) dSbar\(chi\)=([-+0-9.e]+)", line)
    if match:
        length, entropy, change = int(match[1]), float(match[2]), float(match[3])
        coordinates = [math.log(2 * length / math.pi * math.sin(math.pi * cut / length)) / 6
                       for cut in range(1, length) if abs(cut - length / 2) <= length / 8]
        runs[label].append((length, sum(coordinates) / len(coordinates), entropy, change))

fits = []
for label in ["tg_d0.12500_h1.74050_rot", "tg_d0.18750_h1.78300_rot", "tg_d0.18750_h1.77200_rot"]:
    rows = [row for row in runs[label] if row[0] >= 96 and abs(row[3]) <= 0.02]
    x_mean = sum(row[1] for row in rows) / len(rows)
    y_mean = sum(row[2] for row in rows) / len(rows)
    xx = sum((row[1] - x_mean) ** 2 for row in rows)
    slope = sum((row[1] - x_mean) * (row[2] - y_mean) for row in rows) / xx
    intercept = y_mean - slope * x_mean
    residual = sum((row[2] - intercept - slope * row[1]) ** 2 for row in rows)
    error = math.sqrt(residual / (len(rows) - 2) / xx)
    fits.append({"run": label, "lengths": [row[0] for row in rows], "c_fit": slope,
                 "regression_error_only": error,
                 "largest_last_entropy_step": max(abs(row[3]) for row in rows)})


def kron(left, right):
    return [[left[i // len(right)][j // len(right[0])] * right[i % len(right)][j % len(right[0])]
             for j in range(len(left[0]) * len(right[0]))]
            for i in range(len(left) * len(right))]


def combine(*terms):
    return [[sum(coefficient * matrix[i][j] for coefficient, matrix in terms)
             for j in range(4)] for i in range(4)]


x, y, z = [[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]
xx, yy, xy, yx, zz = [kron(a, b) for a, b in [(x, x), (y, y), (x, y), (y, x), (z, z)]]
rotation = [cmath.exp(-1j * math.pi / 8), cmath.exp(1j * math.pi / 8)]
rotation_pair = [a * b for a in rotation for b in rotation]
delta = 3 / 16
original = combine(((1 + delta) / 2, xy), ((1 - delta) / 2, yx), (1, zz))
mapped = combine((0.5, xx), (-0.5, yy), (delta / 2, xy), (-delta / 2, yx), (1, zz))
rotation_error = max(abs(rotation_pair[i].conjugate() * original[i][j] * rotation_pair[j] - mapped[i][j])
                     for i in range(4) for j in range(4))
result = {
    "source": str(SOURCE.relative_to(HERE.parent.parent)),
    "source_sha256": hashlib.sha256(source_bytes).hexdigest(),
    "selection": {"minimum_length": 96, "maximum_absolute_last_entropy_step": 0.02},
    "fits": fits,
    "bond_rotation_max_entry_error": rotation_error,
    "displayed_full_field_coefficient_spreads": {"delta_3_16_h_1_783": 1.783 / (13 / 32),
                                                "delta_1_8_h_1_7405": 1.7405 / (7 / 16)},
    "scope": "Rounded-table regression and local matrix identity; no new many-body simulation.",
}
(HERE / "checks.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
