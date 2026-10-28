"""Small controls of C59's partial lemmas, not an XXZ RG calculation."""
import json
import math
from pathlib import Path

from scipy.integrate import quad


def gaussian_symbol_check():
    # Phi(omega,k;a)=log(1+a exp(-omega^2-k^2)), q=1.
    # Rapid decay suffices for the integration-by-parts proof in Section 4.
    a = 0.23
    rows = []
    for tau in (0.5, 0.3, 0.2, 0.12):
        for derivative in range(3):
            def phi(omega, k):
                b = math.exp(-omega * omega - k * k)
                if derivative == 0:
                    return math.log1p(a * b)
                if derivative == 1:
                    return b / (1 + a * b)
                return -b * b / (1 + a * b) ** 2

            nmax = math.ceil(9 / (2 * math.pi * tau)) + 2
            thermal = tau * sum(
                quad(lambda k: phi((2 * n + 1) * math.pi * tau, k),
                     -9, 9, epsabs=2e-13, epsrel=2e-13)[0] / (2 * math.pi)
                for n in range(-nmax, nmax)
            )
            bulk_terms = []
            image_terms = []
            for r in range(max(1, derivative), 141):
                falling = math.prod(range(r - derivative + 1, r + 1))
                coeff = (-1) ** (r + 1) * falling * a ** (r - derivative) / r**2
                bulk_terms.append(coeff / (4 * math.pi))
                images = math.fsum(
                    2 * (-1) ** ell * math.exp(-ell * ell / (4 * r * tau**2))
                    for ell in range(1, 81)
                )
                image_terms.append(coeff * images / (4 * math.pi))
            bulk = math.fsum(bulk_terms)
            predicted = math.fsum(image_terms)
            residual = abs(thermal - bulk - predicted)
            assert residual < 2e-12, (tau, derivative, residual)
            rows.append({"T_over_q": tau, "derivative": derivative,
                         "poisson_residual": residual,
                         "thermal_difference": thermal - bulk})
    return {"maximum_residual": max(x["poisson_residual"] for x in rows),
            "rows": rows}


def history_check():
    p, theta, sigma, ell = 3.0, 0.25, 1.0, 2.0
    cepsilon, forcing = 0.08, 0.7
    r = 2 ** (-(p-theta)) * (
        ell + cepsilon / (1 - 2 ** (-(sigma+p-theta))))
    bound = forcing * 2 ** (-(p-theta)) / (1-r)
    assert r < 1
    maxima = []
    for depth in (12, 40, 100):
        T = 2.0 ** (-depth)
        values = [0.0]
        for h in range(depth):
            q = 2.0 ** (-h)
            weight = q**theta * (T/q)**p
            values.append(ell * values[h] + cepsilon * math.fsum(
                2 ** (-sigma*(h-k)) * values[k] for k in range(h+1))
                + forcing * weight)
        ratios = [value / ((2.0**(-h))**theta * (T/(2.0**(-h)))**p)
                  for h, value in enumerate(values)]
        maximum = max(ratios)
        assert maximum <= bound * (1 + 2e-14)
        maxima.append({"depth": depth, "maximum_normalized_difference": maximum})
    return {"contraction_factor": r, "proved_bound": bound, "runs": maxima}


def normalization_check():
    theta, B, eta = 0.4, 0.03, 0.06
    log_c = B / (1-2**(-theta))
    max_residual = 0.0
    for h in range(80):
        # Work with logarithms, since individual Z are allowed to diverge.
        log_x = eta*h + B * (1-2**(-theta*h)) / (1-2**(-theta))
        log_y = eta*h
        predicted = -B * 2**(-theta*h) / (1-2**(-theta))
        residual = abs((log_x-log_y-log_c)-predicted)
        max_residual = max(max_residual, residual)
    assert max_residual < 2e-15
    return {"maximum_log_ratio_residual": max_residual}


if __name__ == "__main__":
    result = {"scope": "Analytic lemma controls only. No interacting XXZ computation.",
              "gaussian_symbol": gaussian_symbol_check(),
              "history": history_check(),
              "normalization": normalization_check()}
    Path(__file__).with_name("checks.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": "passed",
                      "poisson_max_residual": result["gaussian_symbol"]["maximum_residual"],
                      "history_contraction": result["history"]["contraction_factor"],
                      "normalization_max_residual": result["normalization"]["maximum_log_ratio_residual"]},
                     indent=2))
