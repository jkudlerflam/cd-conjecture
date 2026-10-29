"""[Numerical] Small controls for C90, not a microscopic infrared test."""

from pathlib import Path
import json
import math

import numpy as np
from scipy.integrate import quad


ROOT = Path(__file__).resolve().parent
assert Path.cwd().resolve() == ROOT, "Run only inside work_C90."


def reference_checks():
    t = 1.3
    rows = []
    for rho in (0.25, 0.4, 0.5, 0.65, 0.8):
        r = 1.0 / rho - 1.0

        def e_gap(x):
            return -2 * t / math.pi * (1 + x) * math.sin(math.pi / (1 + x))

        step = 1.0e-3
        curvature_fd = (
            -e_gap(r + 2 * step) + 16 * e_gap(r + step)
            - 30 * e_gap(r) + 16 * e_gap(r - step) - e_gap(r - 2 * step)
        ) / (12 * step**2)
        curvature = 2 * math.pi * t * rho**3 * math.sin(math.pi * rho)
        v_gap = 2 * t * rho * math.sin(math.pi * rho)
        k_from_fd = math.pi * v_gap / curvature_fd
        assert abs(curvature_fd / curvature - 1) < 2e-7
        assert abs(k_from_fd * rho**2 - 1) < 2e-7
        for b in (14, 15, 24):
            u = (b / rho) * (v_gap / b)
            assert abs(u - 2 * t * math.sin(math.pi * rho)) < 1e-13
        rows.append({"rho": rho, "K_from_curvature": k_from_fd,
                     "K_exact": rho**-2, "curvature_relative_error":
                     abs(curvature_fd / curvature - 1)})

    # Independent physical free-fermion pressure at mu=0, hence rho=1/2.
    e0 = -2 * t / math.pi
    predicted = math.pi / (12 * t)
    thermal = []
    for temperature in (0.08, 0.04, 0.02, 0.01):
        integral = quad(
            lambda k: float(np.logaddexp(0.0, 2*t*math.cos(k)/temperature)),
            0.0, math.pi, points=[math.pi/2], epsabs=2e-11, epsrel=2e-12,
        )[0]
        f = -temperature * integral / math.pi
        measured = (e0 - f) / temperature**2
        thermal.append({"T": temperature, "coefficient": measured,
                        "relative_error": abs(measured / predicted - 1)})
    assert thermal[-1]["relative_error"] < 1e-4
    assert all(thermal[j+1]["relative_error"] < thermal[j]["relative_error"]
               for j in range(len(thermal)-1))
    return {"curvatures": rows, "thermal": thermal,
            "thermal_predicted": predicted,
            "half_density_contact": 0.5 - 2 / math.pi**2}


def pinning_checks():
    # Number-conserving spin chain with complex hopping, density interaction,
    # and a chemical potential. Test every nonempty number sector.
    sigma = np.array([[0, 1], [0, 0]], dtype=complex)
    number = np.diag([0.0, 1.0])
    ident = np.eye(2)
    phase = complex(math.cos(0.37), math.sin(0.37))
    h = (-0.7 * (phase * np.kron(sigma.conj().T, sigma)
                 + phase.conjugate() * np.kron(sigma, sigma.conj().T))
         + 0.23 * np.kron(number, number)
         - 0.11 * (np.kron(number, ident) + np.kron(ident, number)))
    bond_norm = np.max(np.abs(np.linalg.eigvalsh(h)))
    residual = 0.0
    count = 0

    def logtrace(matrix, beta):
        ev = np.linalg.eigvalsh(matrix)
        z = -beta * ev
        return float(np.logaddexp.reduce(z))

    for length in (3, 4, 5, 6):
        full = np.zeros((2**length, 2**length), dtype=complex)
        for word in range(2**length):
            for i in range(length):
                j = (i + 1) % length
                old = 2 * ((word >> i) & 1) + ((word >> j) & 1)
                stripped = word & ~(1 << i) & ~(1 << j)
                for new in range(4):
                    out = stripped | ((new >> 1) << i) | ((new & 1) << j)
                    full[out, word] += h[new, old]
        for occupied in range(1, length+1):
            basis = [w for w in range(2**length) if w.bit_count() == occupied]
            mat = full[np.ix_(basis, basis)]
            projectors = [np.diag([(w >> i) & 1 for w in basis]).astype(complex)
                          for i in range(length)]
            compressed = sum(p @ mat @ p for p in projectors) / occupied
            commsum = sum(p @ (p @ mat - mat @ p) - (p @ mat - mat @ p) @ p
                          for p in projectors)
            residual = max(residual, float(np.max(np.abs(
                compressed - mat + commsum/(2*occupied)))))
            pin = [a for a, w in enumerate(basis) if w & 1]
            pinned = mat[np.ix_(pin, pin)]
            for beta in (0.5, 2.0, 7.0):
                logz = logtrace(mat, beta)
                logpin = logtrace(pinned, beta)
                assert logz + 1e-10 >= math.log(length/occupied) + logpin
                assert logz <= math.log(length) + 4*beta*bond_norm*length/occupied + logpin + 1e-10
                count += 1
    assert residual < 1e-12
    return {"number_of_trace_tests": count, "largest_matrix_dimension": 64,
            "compression_identity_residual": residual}


def penalty_and_harmonic_checks():
    # Synthetic fixed constants, not a value of the unknown C88 C3.
    t, tau, c3, lam = 1.0, 1e-7, 5.0, 0.01
    budget_rows = []
    for coupling in (0.1, 0.01, 0.001):
        nu = 4*t + 8192*coupling
        a2 = 2*t*t + 3*coupling**2 + 6*coupling**2*lam**2
        k = 4 * max(2*a2/(coupling**2*tau),
                    math.sqrt(2*c3*nu**3/(coupling**3*tau)),
                    100*nu/coupling)
        explicit = a2/(coupling*k)
        remainder = c3*nu**3/(coupling*k)**2
        assert explicit < coupling*tau/2
        assert remainder < coupling*tau/2
        budget_rows.append({"J": coupling, "synthetic_K": k,
                            "explicit_over_J_tau": explicit/(coupling*tau),
                            "remainder_over_J_tau": remainder/(coupling*tau)})
    fixed_k = 1e12
    ratios = [t*t/(j*j*fixed_k) for j in (0.1, 0.01, 0.001)]
    assert abs(ratios[1]/ratios[0] - 100) < 1e-10
    assert abs(ratios[2]/ratios[1] - 100) < 1e-10
    margin = 3.2513114e-6 / 30
    assert margin > 1.0837704e-7
    assert 3*math.sqrt(3/5) > 2
    assert np.allclose(np.exp(1j*math.pi*np.array([-1, 1])/2), [-1j, 1j])
    return {"synthetic_budget": budget_rows,
            "fixed_K_color_over_J": ratios, "C88_physical_margin_lower_bound": margin,
            "dimension_three_lower_bound": 3*math.sqrt(3/5)}


result = {"status": "numerical controls passed; no microscopic infrared test",
          "reference": reference_checks(), "pinning": pinning_checks(),
          "penalty_and_harmonics": penalty_and_harmonic_checks()}
(ROOT / "checks.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
