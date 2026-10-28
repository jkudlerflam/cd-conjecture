#!/usr/bin/env python3
"""Reproduce finite checks for SINGLE_CHANNEL_RG.md (requires NumPy).

These checks supplement the algebraic proofs; they do not test an infinite
renormalization-group flow or establish a central charge.
"""

import json
import math
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "contrib" / "routeA_demonstrator"))
from cell import cell_ops, kron, X, Z, I2  # noqa: E402


def symmetry_checks():
    h, _, code, charges = cell_ops(D0=1, D1=1, kap=0.4)
    charge = np.diag(charges)
    conjugation = kron(X, X, X, X, X)
    parity = kron(Z, Z, Z, Z, Z)
    reflection = np.zeros((32, 32))
    for i in range(32):
        bits = [(i >> (4 - j)) & 1 for j in range(5)]
        j = sum(b << (4 - k) for k, b in enumerate(reversed(bits)))
        reflection[j, i] = 1
    matrices = {
        "physical_C_commutes": h @ conjugation - conjugation @ h,
        "physical_Q_reverses": conjugation @ charge @ conjugation + charge,
        "physical_parity_reverses": conjugation @ parity @ conjugation + parity,
        "physical_code_intertwiner": conjugation @ code - code @ X,
        "physical_reflection_commutes": h @ reflection - reflection @ h,
        "physical_reflection_C_commutes": reflection @ conjugation - conjugation @ reflection,
    }

    # Five CAR modes represent one excited pair: f,a1,a2,b1,b2.
    annihilator = np.array([[0, 1], [0, 0]])
    ops = [kron(*([Z] * i + [annihilator] + [I2] * (4 - i))) for i in range(5)]
    numbers = [op.T @ op for op in ops]
    permutation = [0, 3, 4, 1, 2]
    swap = np.zeros((32, 32), dtype=int)
    for i in range(32):
        occupied = [j for j in range(5) if (i >> (4 - j)) & 1]
        image = [permutation[j] for j in occupied]
        inversions = sum(image[a] > image[b] for a in range(len(image)) for b in range(a + 1, len(image)))
        j = sum(1 << (4 - k) for k in image)
        swap[j, i] = (-1) ** inversions
    aux_parity = kron(I2, Z, Z, Z, Z)
    unitary = (ops[0] + ops[0].T) @ aux_parity @ swap
    vacuum = np.eye(32)[:, 0]
    embedding = np.column_stack([
        vacuum,
        ops[0].T @ vacuum,
        ops[1].T @ ops[2].T @ vacuum,
        ops[0].T @ ops[3].T @ ops[4].T @ vacuum,
    ])
    physical_C = np.kron(np.eye(2), X)
    matrices.update({
        "encoded_C_intertwiner": unitary @ embedding - embedding @ physical_C,
        "encoded_C_involution": unitary @ unitary - np.eye(32),
        "encoded_C_logical": unitary @ ops[0] @ unitary - ops[0].T,
        "encoded_C_a1": unitary @ ops[1] @ unitary - ops[3],
        "encoded_C_a2": unitary @ ops[2] @ unitary - ops[4],
    })
    for q in [-2, 0, 2]:
        qtilde = (0.5 * np.eye(32) - numbers[0] + numbers[1]
                  + (q - 1) * numbers[2] - numbers[3] + (1 - q) * numbers[4])
        target = np.diag([0.5, -0.5, 0.5 + q, -0.5 - q])
        matrices[f"charge_intertwiner_q{q}"] = qtilde @ embedding - embedding @ target
        matrices[f"charge_reversal_q{q}"] = unitary @ qtilde @ unitary + qtilde
    residuals = {name: float(np.max(np.abs(matrix))) for name, matrix in matrices.items()}
    assert max(residuals.values()) < 1e-12
    return residuals


def heat_kernel_checks():
    nodes, weights = np.polynomial.legendre.leggauss(500)
    scales, weights = (nodes + 1) / 2, weights / 2
    errors = []
    for energy in [0, 0.1, 0.7, 2]:
        for time in [-5, -2, -0.4, 0.4, 2, 5]:
            q = np.exp(-time * time / (4 * scales)) / np.sqrt(4 * math.pi * scales)
            integral = float(np.dot(weights, (time / (2 * scales) + energy) * q * np.exp(-scales * energy * energy)))
            if energy == 0:
                exact = 0.5 * math.copysign(1, time) * math.erfc(abs(time) / 2)
            else:
                exact = math.exp(-energy * time) * ((1 if time > 0 else 0) - 0.5 * math.erfc(energy - time / 2))
            errors.append(abs(integral - exact))
    assert max(errors) < 2e-13
    rows = {}
    for length in [4, 6, 10, 50, 101]:
        momenta = 2 * math.pi * np.arange(length) / length
        b = np.fft.ifft(np.exp(-(2 * np.cos(momenta)) ** 2))
        rows[length] = float(np.sum(np.abs(b)))
        assert rows[length] <= 1 + 2e-13
    return {"heat_kernel_cases": len(errors), "maximum_error": max(errors),
            "spatial_filter_absolute_sums": rows}


def matrix_uv_checks():
    h0 = 0.2 * np.array([[0.2, 0.25], [0.25, -0.3]])
    h1 = 0.2 * np.array([[1.0, 0.0], [0.17, 0.35]])
    row_bound = float(np.max(np.sum(np.abs(h0) + np.abs(h1) + np.abs(h1.T), axis=1)))
    filter_bound = math.exp(row_bound ** 2)
    nodes, weights = np.polynomial.legendre.leggauss(100)
    scales, weights = (nodes + 1) / 2, weights / 2
    errors = []
    for momentum in np.linspace(-math.pi, math.pi, 9):
        h = h0 + h1 * np.exp(1j * momentum) + h1.T * np.exp(-1j * momentum)
        values, vectors = np.linalg.eigh(h)
        for frequency in [0.2, 0.9, 1.7]:
            chi = (vectors * np.exp(-frequency ** 2 - values ** 2)) @ vectors.conj().T
            direct = (np.eye(2) - chi) @ np.linalg.inv(-1j * frequency * np.eye(2) + h)
            diagonal_integral = np.dot(weights, np.exp(-scales[:, None] * (frequency ** 2 + values[None, :] ** 2)))
            integral = (1j * frequency * np.eye(2) + h) @ ((vectors * diagonal_integral) @ vectors.conj().T)
            errors.append(float(np.max(np.abs(direct - integral))))
    assert max(errors) < 2e-13
    filter_rows = {}
    for length in [4, 7, 16, 51]:
        spectral_filter = []
        for momentum in 2 * math.pi * np.arange(length) / length:
            h = h0 + h1 * np.exp(1j * momentum) + h1.T * np.exp(-1j * momentum)
            values, vectors = np.linalg.eigh(h)
            spectral_filter.append((vectors * np.exp(-values ** 2)) @ vectors.conj().T)
        kernel = np.fft.ifft(np.array(spectral_filter), axis=0)
        actual = float(np.max(np.sum(np.abs(kernel), axis=(0, 2))))
        filter_rows[length] = actual
        assert actual <= filter_bound + 2e-13
    return {"matrix_identity_cases": len(errors), "maximum_error": max(errors),
            "A0": row_bound, "C_B_bound": filter_bound,
            "spatial_filter_absolute_row_sums": filter_rows}


def flow_sum_check():
    gamma, alpha, theta = 2, 0.2, 0.6
    finite = sum(gamma ** (-theta * n) * sum(gamma ** (alpha * r) for r in range(1, n + 1)) for n in range(1, 500))
    closed = 1 / ((1 - gamma ** (-theta)) * (gamma ** (theta - alpha) - 1))
    assert abs(finite - closed) < 1e-12
    return {"gamma": gamma, "alpha": alpha, "theta": theta,
            "finite_sum": finite, "closed_value": closed, "error": abs(finite - closed)}


if __name__ == "__main__":
    result = {"scope": "Finite algebra and covariance checks only; no infrared or central-charge claim.",
              "numpy_version": np.__version__, "symmetry": symmetry_checks(),
              "scalar_uv": heat_kernel_checks(), "matrix_uv": matrix_uv_checks(),
              "flow_sum": flow_sum_check()}
    output = HERE / "single_channel_checks.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"output": str(output), "passed": True,
                      "symmetry_checks": len(result["symmetry"]),
                      "scalar_cases": result["scalar_uv"]["heat_kernel_cases"],
                      "matrix_cases": result["matrix_uv"]["matrix_identity_cases"]}, indent=2))
