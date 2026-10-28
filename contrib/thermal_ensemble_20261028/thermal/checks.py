"""Independent thermal controls; no c-d counterexample is asserted."""

from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np


def quadratic_control() -> dict:
    delta, eta, mu = 3.0, 0.8, 0.2
    kf = math.acos(-(mu + eta * eta / delta) / 2.0)
    z = 1.0 + (eta / delta) ** 2
    velocity = 2.0 * math.sin(kf) / z
    predicted = math.pi / (6.0 * velocity)
    missing_residue = math.pi / (12.0 * math.sin(kf))

    def energies(k: float) -> np.ndarray:
        return np.linalg.eigvalsh(
            [[-2.0 * math.cos(k) - mu, eta], [eta, delta]]
        )

    step = 1.0e-6
    derivative = (energies(kf + step)[0] - energies(kf - step)[0]) / (
        2.0 * step
    )
    assert abs(derivative - velocity) < 2.0e-9
    nodes, weights = np.polynomial.legendre.leggauss(128)

    def quadrature(function, left: float, right: float) -> float:
        halfwidth = (right - left) / 2.0
        midpoint = (right + left) / 2.0
        return halfwidth * sum(
            weight * function(midpoint + halfwidth * node)
            for node, weight in zip(nodes, weights)
        )

    rows = []
    for temperature in (0.08, 0.04, 0.02, 0.01, 0.005, 0.0025):
        width = 25.0 * temperature / velocity
        points = sorted(
            {
                0.0,
                math.pi,
                kf,
                max(0.0, kf - width),
                min(math.pi, kf + width),
            }
        )

        def integrand(k: float) -> float:
            return float(np.logaddexp(0.0, -np.abs(energies(k)) / temperature).sum())

        integral = sum(
            quadrature(integrand, left, right)
            for left, right in zip(points[:-1], points[1:])
        )
        coefficient = integral / (math.pi * temperature)
        aux_free_energy = -temperature * math.log1p(math.exp(-delta / temperature))
        aux_bound = temperature * math.exp(-delta / temperature)
        assert abs(aux_free_energy) <= aux_bound * (1.0 + 1.0e-14)
        rows.append(
            {
                "T": temperature,
                "minus_f_minus_e_over_T2": coefficient,
                "absolute_coefficient_error": abs(coefficient - predicted),
                "auxiliary_abs_free_energy": abs(aux_free_energy),
                "auxiliary_bound": aux_bound,
            }
        )
    assert rows[-1]["absolute_coefficient_error"] < 2.0e-6
    assert all(
        newer["absolute_coefficient_error"] < older["absolute_coefficient_error"]
        for older, newer in zip(rows[:-1], rows[1:])
    )
    assert abs(predicted - missing_residue) > 0.015
    return {
        "parameters": {"Delta": delta, "eta": eta, "mu": mu},
        "k_F": kf,
        "frequency_residue_Z": z,
        "velocity_formula": velocity,
        "velocity_full_matrix_difference": derivative,
        "coefficient_formula": predicted,
        "incorrect_coefficient_if_Z_dropped": missing_residue,
        "temperature_results": rows,
    }


def thermal_modulus_control() -> dict:
    """Check the exact geometric rate constant and a slow modulus example."""
    gamma, theta, exponent = 2.0, 0.25, 4.0
    prefactor = 1.0 / (1.0 - gamma ** (-(exponent - 2.0 - theta)))
    rows = []
    for cutoff in (4, 8, 16, 32, 64, 128):
        temperature = gamma ** (-cutoff)
        active_power = 0.0
        active_log = 0.0
        for h in range(-cutoff, 1):
            scale = gamma**h
            factor = (temperature / scale) ** (exponent - 2.0)
            active_power += factor * scale**theta
            active_log += factor / (1.0 + abs(math.log(scale)))
        bound = prefactor * temperature**theta
        assert active_power <= bound * (1.0 + 1.0e-14)
        rows.append(
            {
                "minus_log2_T": cutoff,
                "power_error_over_T2": active_power,
                "power_bound_over_T2": bound,
                "logarithmic_error_over_T2": active_log,
            }
        )
    assert all(
        new["logarithmic_error_over_T2"] < old["logarithmic_error_over_T2"]
        for old, new in zip(rows[:-1], rows[1:])
    )

    # Worst-case equality recursion with a genuinely relevant amplification.
    amplification, b, a = 2.05, 0.7, 0.1
    d_constant = max(b, b / (gamma ** (exponent - theta) - amplification))
    max_ratio = 0.0
    for cutoff in range(2, 60):
        temperature = gamma ** (-cutoff)
        difference = b * a * temperature**exponent
        for h in range(0, -cutoff - 1, -1):
            scale = gamma**h
            bound = d_constant * a * temperature**exponent * scale ** (theta - exponent)
            max_ratio = max(max_ratio, difference / bound)
            assert difference <= bound * (1.0 + 1.0e-13)
            difference = amplification * difference + b * a * temperature**exponent * scale ** (theta - exponent)
    return {
        "gamma": gamma,
        "theta": theta,
        "thermal_power_p": exponent,
        "modulus_results": rows,
        "causal_relevant_amplification": amplification,
        "max_recurrence_to_proved_bound_ratio": max_ratio,
    }


def localization_control() -> dict:
    """Exact rational checks of projection, remainder constants, and signs."""

    def evaluate(coefficients, value):
        return sum(c * value**degree for degree, c in enumerate(coefficients))

    cases = 0
    largest_a_ratio = Fraction(0)
    largest_b_ratio = Fraction(0)
    for denominator in (4, 10, 40, 100):
        a = Fraction(1, denominator)
        for trial in range(20):
            coefficients = [Fraction(((trial + 3) * (k + 5)) % 17 - 8, k + 1) for k in range(9)]
            values = {node: evaluate(coefficients, node * a) for node in (-3, -1, 1, 3)}
            a_value = (9 * (values[1] + values[-1]) - values[3] - values[-3]) / 16
            b_value = (27 * (values[1] - values[-1]) - values[3] + values[-3]) / (48 * a)
            m4 = sum(
                abs(coefficients[k]) * Fraction(math.factorial(k), math.factorial(k - 4)) * (3 * a) ** (k - 4)
                for k in range(4, len(coefficients))
            )
            a_bound = Fraction(15, 32) * m4 * a**4
            b_bound = Fraction(3, 16) * m4 * a**3
            a_error, b_error = abs(a_value - coefficients[0]), abs(b_value - coefficients[1])
            assert a_error <= a_bound and b_error <= b_bound
            largest_a_ratio = max(largest_a_ratio, a_error / a_bound)
            largest_b_ratio = max(largest_b_ratio, b_error / b_bound)
            cases += 1

        for degree in range(4):
            coefficients = [Fraction(int(k == degree)) for k in range(4)]
            values = {node: evaluate(coefficients, node * a) for node in (-3, -1, 1, 3)}
            a_value = (9 * (values[1] + values[-1]) - values[3] - values[-3]) / 16
            b_value = (27 * (values[1] - values[-1]) - values[3] + values[-3]) / (48 * a)
            assert a_value == int(degree == 0)
            assert b_value == int(degree == 1)
            cases += 1
    return {
        "exact_rational_cases": cases,
        "largest_A_error_to_bound_ratio": float(largest_a_ratio),
        "largest_B_error_to_bound_ratio": float(largest_b_ratio),
        "A_weight_absolute_sum": "5/4",
        "quadratic_thermal_comparison_power": 3,
        "quartic_thermal_comparison_power": 4,
    }


def main() -> None:
    data = {
        "scope": "Exact quadratic memory and abstract scalar RG checks only",
        "quadratic": quadratic_control(),
        "thermal_criterion": thermal_modulus_control(),
        "thermal_localization": localization_control(),
    }
    path = Path(__file__).with_name("checks.json")
    path.write_text(json.dumps(data, indent=2) + "\n")
    print("PASS: quadratic thermal coefficient, scalar bound, scale estimates, and exact localization checks")
    print(f"Saved {path}")


if __name__ == "__main__":
    main()
