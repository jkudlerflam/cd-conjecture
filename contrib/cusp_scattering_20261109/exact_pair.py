"""Independent exact two-body solution of the symmetric Fibonacci cusp.

Requires numpy and scipy.  Run from the repository root with
  python contrib/cusp_scattering_20261109/exact_pair.py

The real gauge has nearest/next-nearest hops u,+1.  Its total K=0
is total K=pi in C173's original gauge.  No many-body phase is inferred.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.linalg import eig_banded
from scipy.optimize import brentq


def pair_polynomial(p: float, u: float) -> float:
    return (1 - p) * (3 - p) * (1 + p) ** 2 - u * u * p**3


def exact_pair(u: float) -> dict:
    """The unique subcontinuum K=0 bound state, for 0<u<4."""
    if not 0 < u < 4:
        raise ValueError("This certificate is stated only for 0 < u < 4.")
    p = brentq(pair_polynomial, 1e-14, 1, args=(u,), xtol=1e-14)
    s = -u * p / (1 + p)
    z = (s + 1j * np.sqrt(4 * p - s * s)) / 2
    c = np.conj(1 + u * z + z * z)
    norm2 = (
        abs(c) ** 2 * p**2 / (1 - p)
        - np.real(c * c * z**4 / (1 - z * z))
    ) / 2

    def shift_expectation(m: int) -> float:
        return float(
            (
                abs(c) ** 2 * np.real(z**m) * p**2 / (1 - p)
                - np.real(c * c * z ** (4 + m) / (1 - z * z))
            )
            / norm2
        )

    energy = -2 * (p + 1 / p + u * u * p / (1 + p) ** 2)
    threshold = -4 - u * u / 2
    c1, c2 = shift_expectation(1), shift_expectation(2)
    curvature = -u * c1 / 2 - 2 * c2
    r = np.arange(2, 514)
    wave = np.imag(c * z**r) / np.sqrt(norm2)
    ghost = lambda k: float(np.imag(c * z**k) / np.sqrt(norm2))
    residual = []
    for rr in range(2, 100):
        pot = -2 if rr == 2 else (-2 / 3 if rr == 3 else 0)
        actual = (pot - energy) * ghost(rr)
        for delta, hop in ((-2, 2), (-1, 2 * u), (1, 2 * u), (2, 2)):
            if rr + delta >= 2:
                actual += hop * ghost(rr + delta)
        residual.append(abs(actual))
    return {
        "u": u,
        "p": float(p),
        "z": [float(z.real), float(z.imag)],
        "coefficient": [float(c.real), float(c.imag)],
        "unnormalized_norm_squared": float(norm2),
        "energy": float(energy),
        "threshold": float(threshold),
        "binding": float(threshold - energy),
        "binding_positive_formula": float(
            (1 - p) ** 2 / (2 * p) * (4 - u * u * p / (1 + p) ** 2)
        ),
        "inverse_pair_mass": float(curvature),
        "free_pair_inverse_mass": float(4 - u * u / 4),
        "shift_one": c1,
        "shift_two": c2,
        "polynomial_residual": float(abs(pair_polynomial(p, u))),
        "recurrence_residual": float(max(residual)),
        "norm_through_r513": float(np.dot(wave, wave)),
    }


def box_energy(u: float, cutoff: int, total_k: float = 0) -> float:
    """Independent finite half-line diagonalization on 2<=r<=cutoff."""
    n = cutoff - 1
    bands = np.zeros((3, n))
    bands[0, 0] = -2
    bands[0, 1] = -2 / 3
    bands[1, : n - 1] = 2 * u * np.cos(total_k / 2)
    bands[2, : n - 2] = 2 * np.cos(total_k)
    return float(
        eig_banded(
            bands, lower=True, eigvals_only=True, select="i", select_range=(0, 0)
        )[0]
    )


def boundary_columns(qs, u: float, v2: float, v3: float) -> np.ndarray:
    z = np.exp(1j * np.asarray(qs))
    return np.array([2 * u * z + 2 - v2 * z**2, 2 * z - v3 * z**3])


def threshold_matrices(u: float, v2: float, v3: float) -> tuple:
    q = np.arccos(-u / 4)
    z = np.exp(1j * np.array([q, -q]))
    const = np.array([2 * u * z + 2 - v2 * z**2, 2 * z - v3 * z**3])
    linear = np.array([2 * u * z - 2 * v2 * z**2, 2 * z - 3 * v3 * z**3])
    raw = np.linalg.solve(const, linear)
    return raw, (raw + raw.conj().T) / 2


def flux_s_matrix(u: float, v2: float, v3: float, excess: float) -> np.ndarray:
    """S in incoming (q_low,-q_high), outgoing (q_high,-q_low) order.

    E(q)=4u*cos(q)+4*cos(2q).  Each incoming/outgoing state has unit
    relative-coordinate probability flux.  E=E_threshold+excess.
    """
    q_low = np.arccos((-u / 2 + np.sqrt(excess / 2)) / 2)
    q_high = np.arccos((-u / 2 - np.sqrt(excess / 2)) / 2)
    incoming = np.array([q_low, -q_high])
    outgoing = np.array([q_high, -q_low])
    velocity = lambda q: -4 * u * np.sin(q) - 8 * np.sin(2 * q)
    vin, vout = -velocity(incoming), velocity(outgoing)
    if min(vin.min(), vout.min()) <= 0:
        raise ValueError("Expected two incoming and two outgoing channels.")
    s = -np.linalg.solve(
        boundary_columns(outgoing, u, v2, v3),
        boundary_columns(incoming, u, v2, v3),
    )
    return np.sqrt(vout)[:, None] * s / np.sqrt(vin)[None, :]


def complex_list(a: np.ndarray) -> list:
    return [[[float(x.real), float(x.imag)] for x in row] for row in a]


def run_checks() -> dict:
    results = {
        "scope": "Two-body effective Hamiltonian only; no many-body c verdict.",
        "pair_rows": [],
        "flux_rows": [],
    }
    for u in (0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 3.9):
        row = exact_pair(u)
        cutoff = 2048 if u == 0.25 else 1024
        row["box_cutoff"] = cutoff
        row["box_energy"] = box_energy(u, cutoff)
        row["box_error"] = abs(row["box_energy"] - row["energy"])
        # A separately diagonalized finite K branch checks Hellmann-Feynman.
        dk = 1e-3
        row["box_curvature"] = (
            box_energy(u, cutoff, dk)
            + box_energy(u, cutoff, -dk)
            - 2 * row["box_energy"]
        ) / dk**2
        row["curvature_error"] = abs(
            row["box_curvature"] - row["inverse_pair_mass"]
        )
        assert row["binding"] > 0
        assert row["recurrence_residual"] < 1e-11
        assert row["box_error"] < 1e-10
        assert row["curvature_error"] < 1e-5
        results["pair_rows"].append(row)

    for u, v2, v3 in ((1.0, -2.0, -2 / 3), (1.0, 0.0, 0.0), (2.0, -2.0, -2 / 3)):
        raw, physical = threshold_matrices(u, v2, v3)
        row = {
            "u": u,
            "v2": v2,
            "v3": v3,
            "raw_matrix": complex_list(raw),
            "unit_flux_threshold_matrix": complex_list(physical),
            "eigenvalues": np.linalg.eigvalsh(physical).tolist(),
            "finite_energy": [],
        }
        for excess in (1e-4, 1e-6, 1e-8):
            s = flux_s_matrix(u, v2, v3, excess)
            delta = np.sqrt(excess / (8 - u * u / 2))
            finite_length = (s + np.eye(2)) / (2j * delta)
            # The O(delta) anti-Hermitian term is fixed by S unitarity;
            # compare the Hermitian part, whose error starts at O(delta^2).
            finite_hermitian = (finite_length + finite_length.conj().T) / 2
            error = float(np.linalg.norm(finite_hermitian - physical))
            unitarity = float(np.linalg.norm(s.conj().T @ s - np.eye(2)))
            assert unitarity < 1e-7
            row["finite_energy"].append(
                {"excess": excess, "threshold_error": error, "unitarity_error": unitarity}
            )
        assert row["finite_energy"][-1]["threshold_error"] < 3e-5
        results["flux_rows"].append(row)

    return results


if __name__ == "__main__":
    result = run_checks()
    destination = Path(__file__).with_name("results.json")
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Passed {len(result['pair_rows'])} exact-pair and {len(result['flux_rows'])} flux checks.")
    for row in result["pair_rows"]:
        print(
            f"u={row['u']:.2f} E_pair={row['energy']:.12f} "
            f"binding={row['binding']:.12f} inverse_mass={row['inverse_pair_mass']:.12f}"
        )
    print(destination)
