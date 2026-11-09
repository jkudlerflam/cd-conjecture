"""Pair energy shift from the complete leading facilitated extension.

This uses C175's existing, number-conserving three-body gap-coordinate
Hamiltonian, but derives the pair wavefunction and number-changing source
independently.  Output is a second-order perturbative coefficient in
alpha=lambda/Delta.  It is not a finite-alpha phase calculation.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import spsolve

from exact_pair import exact_pair

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "codex/work_C175"))
from nbody_box import build  # noqa: E402


def pair_source(cutoff: int, u: float = 1.0) -> tuple:
    info = exact_pair(u)
    z = complex(*info["z"])
    c = complex(*info["coefficient"])
    norm = info["unnormalized_norm_squared"]

    def f(r):
        # C175's original gauge and its leftmost-coordinate convention.
        return np.imag(c * z**r) * (1j) ** r / np.sqrt(norm)

    g, h = np.meshgrid(np.arange(2, cutoff + 1), np.arange(2, cutoff + 1), indexing="ij")
    source = 2 * (
        (g == 2) * (f(h) - f(g + h))
        + (h == 2) * (f(g + h) - f(g))
    )
    direct = 4 * abs(f(2)) ** 2
    return info, source.ravel(), float(direct)


def run() -> dict:
    output = {
        "u": 1.0,
        "K_original_gauge": float(np.pi),
        "perturbation_parameter": "alpha=lambda/Delta",
        "coefficient_definition": "E_pair(alpha)-E_pair(0)=coefficient*alpha^2+O(alpha^4)",
        "scope": "isolated two-body branch in a vacuum calculation; chemical potential held fixed",
        "rows": [],
    }
    for cutoff in (40, 80, 120, 160, 240):
        info, source, direct = pair_source(cutoff)
        h3 = build(3, np.pi, 1j, cutoff)
        assert np.linalg.norm(source.reshape(cutoff - 1, cutoff - 1) + source.reshape(cutoff - 1, cutoff - 1).T) < 1e-12
        for mu in (-2.4, -2.3, -2.27, -2.26):
            denominator = h3 + (-mu - info["energy"]) * eye(h3.shape[0])
            solution = spsolve(denominator.tocsc(), source)
            correction = np.vdot(source, solution)
            residual = np.linalg.norm(denominator @ solution - source)
            assert abs(correction.imag) < 1e-10
            assert residual < 1e-10
            row = {
                "cutoff_each_gap": cutoff,
                "dimension": h3.shape[0],
                "mu": mu,
                "direct_repulsion_coefficient": direct,
                "trimer_resolvent_coefficient": float(correction.real),
                "pair_shift_coefficient": float(direct - correction.real),
                "linear_solve_residual": float(residual),
                "source_norm_squared": float(np.vdot(source, source).real),
            }
            output["rows"].append(row)
            print(f"R={cutoff} mu={mu}: pair shift/alpha^2={row['pair_shift_coefficient']:.12f}", flush=True)
    return output


if __name__ == "__main__":
    result = run()
    target = Path(__file__).with_name("facilitated_pair_shift.json")
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(target)
