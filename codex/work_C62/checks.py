"""Finite controls for C62. No interacting XXZ RG calculation."""

from fractions import Fraction as F
from pathlib import Path
import json
import math
import numpy as np


def bilinear(state, create, annihilate):
    if not (state >> annihilate) & 1:
        return None
    sign = (-1) ** ((state & ((1 << annihilate) - 1)).bit_count())
    state ^= 1 << annihilate
    if (state >> create) & 1:
        return None
    sign *= (-1) ** ((state & ((1 << create) - 1)).bit_count())
    return state | (1 << create), sign


def current_columns(N, m):
    columns = []
    for state in range(1 << (2 * N)):
        column = {}
        for j in range(2 * N - abs(m)):
            create, annihilate = (j, j + m) if m > 0 else (j - m, j)
            out = bilinear(state, create, annihilate)
            if out:
                target, sign = out
                column[target] = column.get(target, 0) + sign
        columns.append({k: v for k, v in column.items() if v})
    return columns


def product_column(left, right, state):
    output = {}
    for intermediate, val in right[state].items():
        for target, coeff in left[intermediate].items():
            output[target] = output.get(target, 0) + val * coeff
    return output


def current_checks():
    cases = []
    for N in range(1, 5):
        for m in range(1, N + 1):
            plus, minus = current_columns(N, m), current_columns(N, -m)
            trace = 0
            max_defect = 0
            for state in range(1 << (2 * N)):
                pm = product_column(plus, minus, state)
                mp = product_column(minus, plus, state)
                comm = {k: pm.get(k, 0) - mp.get(k, 0) for k in pm.keys() | mp.keys()}
                comm = {k: v for k, v in comm.items() if v}
                lower = sum((state >> j) & 1 for j in range(m))
                upper = sum((state >> j) & 1 for j in range(2 * N - m, 2 * N))
                expected = lower - upper
                assert comm == ({state: expected} if expected else {})
                trace += expected
                max_defect = max(max_defect, abs(m - expected))
            sea = (1 << N) - 1
            sea_pm = product_column(plus, minus, sea).get(sea, 0)
            sea_mp = product_column(minus, plus, sea).get(sea, 0)
            assert sea_pm - sea_mp == m
            assert trace == 0 and max_defect == 2 * m
            cases.append({"N": N, "m": m, "dimension": 1 << (2 * N),
                          "commutator_trace": trace, "operator_norm_defect": max_defect})
    # Independent occupation evaluation of the free thermal boundary bound.
    max_ratio = 0.0
    for N in (2, 5, 12):
        for m in range(1, N + 1):
            for beta_spacing in (0.2, 0.7, 2.0):
                holes = sum(1 / (1 + math.exp(-beta_spacing * (j + 0.5)))
                            for j in range(-N, -N + m))
                particles = sum(1 / (1 + math.exp(beta_spacing * (j + 0.5)))
                                for j in range(N - m, N))
                bound = 2 * m * math.exp(-beta_spacing * (N - m + 0.5))
                ratio = (holes + particles) / bound
                assert ratio <= 1 + 1e-14
                max_ratio = max(max_ratio, ratio)
    return {"exact_matrix_cases": cases, "largest_free_thermal_bound_ratio": max_ratio}


def parameter_checks():
    cases = []
    v = F(7, 5)
    for tau in (F(-1, 5), F(-1, 17), F(0), F(1, 23), F(1, 4)):
        K = (1 - tau) / (1 + tau)
        A = v * (1 + tau * tau) / (1 - tau * tau)
        B = 2 * v * tau / (1 - tau * tau)
        assert A > abs(B)
        assert A * A - B * B == v * v
        assert (A - B) / (A + B) == K * K
        assert B / A == 2 * tau / (1 + tau * tau)
        assert B == (2 * v * tau) / (1 - tau * tau)
        cases.append({"tau": str(tau), "A": str(A), "B": str(B), "K": str(K)})
    return cases


def scale_checks():
    rows = []
    for sigma, theta, p in ((0.8, 0.2, 3), (0.4, 0.3, 3), (1.2, 0.6, 3)):
        bound = 1 / ((1 - 2 ** (-(sigma - theta))) * (1 - 2 ** (-p)))
        for depth in (0, 4, 16, 60, 120):
            normalized = math.fsum(
                2 ** (theta * r - p * s - sigma * max(r, s))
                for r in range(depth + 1) for s in range(depth + 1))
            assert normalized <= bound * (1 + 1e-14)
            rows.append({"sigma": sigma, "theta": theta, "p": p, "depth": depth,
                         "normalized_joint_sum": normalized, "bound": bound})
    # At theta=sigma, the s=0 subseries is depth+1.
    for depth in (4, 40, 120):
        assert sum(F(1) for _ in range(depth + 1)) == depth + 1
    return rows


def interpolation_checks():
    a, t = F(1, 20), F(2, 5)
    nu, r = -F(4, 7) * t * a * a, t * a
    rows = []
    for h in range(33):
        expected = -F(4, 7) * t * t * a * a * F(1, 4) ** h
        expected -= F(4, 7) * a * a * t * (1 - t) * 2 ** h
        assert nu == expected and r == t * a * F(1, 2) ** h
        if h in (0, 4, 16, 32):
            rows.append({"depth": h, "nu": str(nu), "r": str(r)})
        nu, r = 2 * nu + r * r, r / 2
    # Compare with the true critical endpoint under the same map.
    nu, r = -F(4, 7) * a * a, a
    for _ in range(33):
        assert nu == -F(4, 7) * r * r
        nu, r = 2 * nu + r * r, r / 2
    return rows


def reference_trace_checks():
    # L=2*pi, A=1.3, only transfer modes m=0,+1,-1 interact.
    # Assemble the finite fermion Hamiltonian, not truncated bosonic matrices.
    A, B = 1.3, 0.2
    u = math.sqrt(A * A - B * B)
    ground = u - A
    spectra = []
    for N in (1, 2):
        dim = 1 << (2 * N)
        current = np.zeros((dim, dim))
        for col, values in enumerate(current_columns(N, 1)):
            for row, coefficient in values.items():
                current[row, col] = coefficient
        charges = np.array([state.bit_count() - N for state in range(dim)])
        energies = np.array([sum(A * (j - N + 0.5)
                                 * (((state >> j) & 1) - (j < N))
                                 for j in range(2 * N)) for state in range(dim)])
        diagonal = (energies[:, None] + energies[None, :]
                    + B * charges[:, None] * charges[None, :]).reshape(-1)
        pair = np.kron(current, current)
        H = np.diag(diagonal) + B * (pair + pair.T)
        eig = np.linalg.eigvalsh(H)
        assert eig[0] >= ground - 1e-12
        spectra.append(eig)
    assert spectra[1][0] <= spectra[0][0] + 1e-12
    assert np.all(spectra[1][:len(spectra[0])] <= spectra[0] + 1e-12)
    rows = []
    for beta in (2.0, 4.0):
        # All omitted charge and oscillator terms here are below 1e-40.
        charge = math.fsum(math.exp(-beta * (A * (qp * qp + qm * qm) / 2 + B * qp * qm))
                           for qp in range(-20, 21) for qm in range(-20, 21))
        log_osc = -2 * math.log1p(-math.exp(-beta * u))
        log_osc -= 2 * math.fsum(math.log1p(-math.exp(-beta * A * m)) for m in range(2, 101))
        target = math.exp(-beta * ground + log_osc) * charge
        truncated = [float(np.exp(-beta * eig).sum()) for eig in spectra]
        assert truncated[0] <= truncated[1] <= target + 1e-11
        rows.append({"beta": beta, "fermion_traces_N1_N2": truncated,
                     "oscillator_trace": target,
                     "ground_energies_N1_N2": [float(eig[0]) for eig in spectra],
                     "oscillator_ground_energy": ground})
    return rows


if __name__ == "__main__":
    result = {"scope": "Finite algebra and positive scale sums only. TP remains unproved.",
              "current_algebra": current_checks(), "parameter_map": parameter_checks(),
              "joint_scales": scale_checks(), "critical_interpolation": interpolation_checks(),
              "reference_traces": reference_trace_checks()}
    Path(__file__).with_name("checks.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": "passed", "exact_current_cases": len(result["current_algebra"]["exact_matrix_cases"]),
                      "max_matrix_dimension": 256, "parameter_cases": len(result["parameter_map"]),
                      "joint_scale_cases": len(result["joint_scales"]),
                      "reference_trace_cases": len(result["reference_traces"])}, indent=2))
