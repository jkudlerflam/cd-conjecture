"""Finite-matrix sanity check, not a locality proof or vacancy simulation."""
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import expm
from scipy.special import erfc


def main():
    n = 4
    eye = np.eye(2)
    pauli_x = np.array([[0, 1], [1, 0]], dtype=complex)
    pauli_z = np.diag([1, -1]).astype(complex)

    def embed(ops):
        ans = np.array([[1.0]], dtype=complex)
        for j in range(n):
            ans = np.kron(ans, ops.get(j, eye))
        return ans

    terms = []
    supports = []
    for j in range(n):
        k = (j + 1) % n
        terms.append(0.4 * embed({j: pauli_x, k: pauli_x})
                     + 0.13 * embed({j: pauli_z, k: pauli_z}))
        supports.append({j, k})
    h0 = sum(terms)
    vv = [0.7 * embed({j: pauli_z}) + 0.11 * embed({j: pauli_x})
          for j in range(n)]
    v = sum(vv)
    delta = 0.8
    cutoff_time = delta ** -2
    # Radius one is chosen to retain nonzero boundary terms in this algebra
    # check. It is not the large radius prescribed for the theorem's bound.
    nodes, weights = leggauss(96)
    times = np.concatenate((cutoff_time * (nodes - 1) / 2,
                            cutoff_time * (nodes + 1) / 2))
    quad_weights = np.tile(cutoff_time * weights / 2, 2)
    ss, ww, ee = [], [], []
    for j in range(n):
        region = {(j - 1) % n, j, (j + 1) % n}
        inside = sum(a for a, support in zip(terms, supports)
                     if support <= region)
        boundary = sum(a for a, support in zip(terms, supports)
                       if support & region and not support <= region)
        eig, basis = np.linalg.eigh(inside)
        mat = basis.conj().T @ vv[j] @ basis
        differences = eig[:, None] - eig[None, :]

        def tau(time):
            return basis @ (np.exp(1j * time * differences) * mat) @ basis.conj().T

        s_local = np.zeros_like(h0)
        w_local = np.zeros_like(h0)
        e_local = np.zeros_like(h0)
        for time, weight in zip(times, quad_weights):
            evolved = tau(time)
            ff = -0.5j * np.sign(time) * erfc(delta * abs(time) / 2)
            gg = delta / (2 * np.sqrt(np.pi)) * np.exp(-delta**2 * time**2 / 4)
            s_local += weight * ff * evolved
            w_local += weight * gg * evolved
            e_local += weight * ff * (evolved @ boundary - boundary @ evolved)
        ff_pos = -0.5j * erfc(delta * cutoff_time / 2)
        ff_neg = -ff_pos
        e_local += 1j * ff_pos * tau(cutoff_time) - 1j * ff_neg * tau(-cutoff_time)
        ss.append(s_local)
        ww.append(w_local)
        ee.append(e_local)

    generator = sum(ss)
    derivative = v + generator @ h0 - h0 @ generator
    first_order = sum(ww) + sum(ee)
    identity_error = float(np.linalg.norm(derivative - first_order, 2))
    antihermitian_error = float(max(np.linalg.norm(a + a.conj().T, 2) for a in ss))
    trials = []
    for coupling in (0.02, 0.01, 0.005, 0.0025):
        unitary = np.eye(2**n, dtype=complex)
        for generator_local in ss:
            unitary = expm(coupling * generator_local) @ unitary
        transformed = unitary @ (h0 + coupling * v) @ unitary.conj().T
        remainder = transformed - h0 - coupling * first_order
        norm = float(np.linalg.norm(remainder, 2))
        trials.append({"s": coupling, "remainder_norm": norm,
                       "remainder_over_s_squared": norm / coupling**2})
    assert identity_error < 1e-11
    assert antihermitian_error < 1e-12
    ratios = [row["remainder_over_s_squared"] for row in trials]
    assert max(ratios) / min(ratios) < 1.1
    result = {
        "scope": "Finite matrix algebra sanity check only; no vacancy or locality conclusion",
        "sites": n, "delta": delta, "truncation_time": cutoff_time,
        "generator_antihermitian_error": antihermitian_error,
        "first_derivative_identity_error": identity_error,
        "quadratic_remainder_trials": trials,
    }
    output = Path(__file__).with_name("checks.json")
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
