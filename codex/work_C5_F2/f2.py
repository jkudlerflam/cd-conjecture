#!/usr/bin/env python3
"""Closed-book checks for F2. No input files, network, or external datasets.

Run: python3 f2.py
Writes f2_checks.json beside this script. These are finite-size/algebraic checks,
not a numerical proof of the c-d conjecture. Only numpy/scipy are required.
"""

import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import spence


def spin_hamiltonian(length, jxy, jz, vacancy=None):
    """Spin-1/2 ring, with an optional deleted site and its two bonds."""
    active = [j for j in range(length) if j != vacancy]
    position = {j: p for p, j in enumerate(active)}
    dim = 1 << len(active)
    h = np.zeros((dim, dim), dtype=float)
    for site in range(length):
        other = (site + 1) % length
        if site not in position or other not in position:
            continue
        p, q = position[site], position[other]
        for state in range(dim):
            bp, bq = (state >> p) & 1, (state >> q) & 1
            h[state, state] += jz * (bp - 0.5) * (bq - 0.5)
            if bp != bq:
                flipped = state ^ (1 << p) ^ (1 << q)
                h[flipped, state] += 0.5 * jxy
    return h


def contract_site(psi, length, site, bit):
    """Contract a local spin. Fermion string at fixed single occupancy is a sign."""
    result = np.zeros(1 << (length - 1), dtype=complex)
    for state, amplitude in enumerate(psi):
        if ((state >> site) & 1) != bit:
            continue
        low = state & ((1 << site) - 1)
        high = (state >> (site + 1)) << site
        result[low | high] = ((-1) ** site) * amplitude
    return result


def electron_green_check(length=6, u=4.0):
    output = []
    for name, jxy, jz in [("critical_XX", 1.0, 0.0),
                           ("gapped_Ising", 0.0, 1.0)]:
        h = spin_hamiltonian(length, jxy, jz)
        if jxy:
            energies, vectors = np.linalg.eigh(h)
            e0, psi = energies[0], vectors[:, 0]
        else:
            # A spin-flip invariant ground state. Its degeneracy is finite,
            # not an extensive set of zero-energy degrees of freedom.
            state = sum(1 << j for j in range(0, length, 2))
            complement = ((1 << length) - 1) ^ state
            psi = np.zeros(1 << length)
            psi[state] = psi[complement] = 1 / math.sqrt(2)
            e0 = float(psi @ h @ psi)
        vacancy = 0
        ev, vv = np.linalg.eigh(spin_hamiltonian(length, jxy, jz, vacancy))
        addition = contract_site(psi, length, vacancy, 0)  # down -> double
        removal = contract_site(psi, length, vacancy, 1)   # up -> empty
        wp = np.abs(vv.T.conj() @ addition) ** 2
        wm = np.abs(vv.T.conj() @ removal) ** 2
        gaps = u + ev - e0

        def green(z):
            return np.sum(wp / (z - gaps) + wm / (z + gaps))

        frequencies = [0.0, 0.1, 0.5, 2.0]
        values = [green(1j * omega) for omega in frequencies]
        # Local number conservation proves spatial off-diagonal entries vanish.
        # The numerical calculation independently checks the onsite cancellation.
        for value in values[1:]:
            assert abs(value.real) < 1e-11
            assert value.imag < 0
        assert abs(values[0]) < 1e-11
        assert abs(np.sum(wp + wm) - 1) < 1e-11
        odd_error = max(abs(green(1j*w) + green(-1j*w))
                        for w in frequencies[1:])
        assert odd_error < 1e-11
        bond = spin_hamiltonian(2, jxy, jz) / 2  # ring counts this bond twice
        bond_norm = float(np.max(np.abs(np.linalg.eigvalsh(bond))))
        output.append({
            "model": name, "length": length, "U": u,
            "spin_ground_energy": float(e0),
            "ground_residual": float(np.linalg.norm(h @ psi - e0 * psi)),
            "rigorous_charge_gap_lower_bound": u - 2 * bond_norm,
            "minimum_one_vacancy_energy": float(np.min(gaps)),
            "spectral_weight": float(np.sum(wp + wm)),
            "G_zero_abs": float(abs(values[0])),
            "frequency_oddness_error": float(odd_error),
            "G_iomega": [{"omega": w, "real": float(g.real),
                           "imag": float(g.imag)}
                          for w, g in zip(frequencies, values)],
            "central_charge_from_analytic_argument": 1 if jxy else 0,
        })
    return output


def dressed_sea_check(kappa=0.35, mu=-0.2):
    """Nonzero-kernel scalar TBA. This is a consistency model, not a claimed
    Bethe solution of a microscopic NN chain. K=kappa/(2*pi), a=1/(2*pi).
    """
    def a_zero(m):
        q = math.acos(-m)
        return (math.sin(q) + m*q) / math.pi

    m0 = brentq(lambda m: m-mu-kappa*a_zero(m), -0.9, 0.9, xtol=1e-14)
    q = math.acos(-m0)
    density = 1 / (2*math.pi*(1-kappa*q/math.pi))
    velocity = math.sin(q) / (2*math.pi*density)
    f0 = -a_zero(m0)
    expected = math.pi / (6*velocity)
    results = []
    for temp in [0.04, 0.02, 0.01, 0.005]:
        def a_temp(m):
            integrand = lambda x: temp*np.logaddexp(0.0, (math.cos(x)+m)/temp)
            return quad(integrand, 0, math.pi, epsabs=1e-13,
                        epsrel=1e-12, points=[q])[0] / math.pi

        mt = brentq(lambda m: m-mu-kappa*a_temp(m), -0.9, 0.9, xtol=1e-14)
        free_energy = -a_temp(mt)
        coefficient = -(free_energy-f0)/temp**2
        results.append({"T": temp, "coefficient": coefficient,
                        "c_inferred_with_dressed_velocity": coefficient/expected})
    assert abs(results[-1]["c_inferred_with_dressed_velocity"]-1) < 2e-4
    return {"kernel": "kappa/(2*pi)", "kappa": kappa, "mu": mu,
            "Q": q, "rho_total": density, "dressed_velocity": velocity,
            "predicted_coefficient": expected, "results": results}


def rogers_l(x):
    return float(spence(1-x) + 0.5*math.log(x)*math.log1p(-x))


def parafermion_check():
    """Check, but do not prove, the stated finite A-type dilogarithm identity."""
    output = []
    for k in range(2, 13):
        theta = math.pi/(k+2)
        xs = [math.sin(theta)**2 / math.sin((a+1)*theta)**2
              for a in range(1, k)]
        ys = [0.0] + [1/x-1 for x in xs] + [0.0]
        recurrence_error = max(abs(ys[a]**2-(1+ys[a-1])*(1+ys[a+1]))
                               for a in range(1, k))
        central = 6/math.pi**2 * sum(rogers_l(x) for x in xs)
        expected = 2*(k-1)/(k+2)
        assert abs(central-expected) < 1e-12
        assert recurrence_error < 1e-9
        output.append({"k": k, "auxiliary_nodes": k-1,
                       "c_parafermion_dilog": central,
                       "c_parafermion_formula": expected,
                       "c_boson_plus_parafermion": 1+central,
                       "c_SU2k": 3*k/(k+2),
                       "Y_system_max_error": recurrence_error})
    return output


def permutation_schmidt_check():
    results = []
    for d in [2, 3, 4, 5]:
        p = np.zeros((d, d, d, d))
        for i in range(d):
            for j in range(d):
                p[i, j, j, i] = 1
        identity = np.eye(d*d).reshape(d, d, d, d)
        traceless = p - identity/d
        # Reshuffle (output A, output B, input A, input B) to operator bipartition.
        r = p.transpose(0, 2, 1, 3).reshape(d*d, d*d)
        rt = traceless.transpose(0, 2, 1, 3).reshape(d*d, d*d)
        rank = int(np.linalg.matrix_rank(r))
        rank_t = int(np.linalg.matrix_rank(rt))
        assert rank == d*d and rank_t == d*d-1
        results.append({"d": d, "swap_operator_Schmidt_rank": rank,
                        "traceless_operator_Schmidt_rank": rank_t,
                        "SUd_Lie_rank": d-1, "ULS_c_given": d-1})
    return results


def neutral_susceptibility_check(j=1.0):
    """Thermodynamic XX energy in a field h coupled to S^z, within the
    exactly singly occupied sector. This is not a finite-size derivative.
    """
    def energy(h):
        return -(math.sqrt(j*j-h*h)+h*math.asin(h/j))/math.pi

    results = []
    for field in [0.1, 0.03, 0.01, 0.003]:
        chi = -(energy(field)+energy(-field)-2*energy(0))/field**2
        results.append({"field": field, "spin_susceptibility": chi})
    exact = 1/(math.pi*j)
    assert abs(results[-1]["spin_susceptibility"]/exact-1) < 1e-6
    # Q_up=N/2+S^z, Q_down=N/2-S^z. Total physical charge is in a plateau.
    matrix = exact*np.array([[1.0, -1.0], [-1.0, 1.0]])
    assert np.linalg.matrix_rank(matrix) == 1
    return {"thermodynamic_spin_susceptibility": exact,
            "physical_total_charge_susceptibility": 0.0,
            "spin_resolved_charge_matrix": matrix.tolist(),
            "active_rank": 1, "results": results}


def main():
    started = time.monotonic()
    result = {
        "scope": "Finite checks only; no conjecture proof or CFT extraction by ED.",
        "electron_green": electron_green_check(),
        "regular_TBA": dressed_sea_check(),
        "parafermion_dilog": parafermion_check(),
        "bond_ranks": permutation_schmidt_check(),
        "neutral_susceptibility": neutral_susceptibility_check(),
    }
    result["runtime_seconds"] = time.monotonic()-started
    destination = Path(__file__).with_name("f2_checks.json")
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "status": "all assertions passed",
        "runtime_seconds": result["runtime_seconds"],
        "max_abs_G_zero": max(x["G_zero_abs"] for x in result["electron_green"]),
        "TBA_lowest_T_c": result["regular_TBA"]["results"][-1]["c_inferred_with_dressed_velocity"],
        "output": str(destination),
    }, indent=2))


if __name__ == "__main__":
    main()
