"""Finite algebra controls for SEPARATE_JW.md. NumPy only."""
from itertools import product
import json
from pathlib import Path

import numpy as np


def kron(*operators):
    result = np.ones((1, 1))
    for operator in operators:
        result = np.kron(result, operator)
    return result


def error(matrix):
    return float(np.max(np.abs(matrix)))


def main():
    i2 = np.eye(2)
    i4 = np.eye(4)
    i16 = np.eye(16)
    a = np.array([[0., 1.], [0., 0.]])
    z = np.diag([1., -1.])
    x = a + a.T
    # On-cell order: payload qubit, length qubit, auxiliary payload and
    # length charge qubits. The auxiliary vacuum is |00>.
    local_c = kron(a, i2, i4)
    local_b = kron(i2, a, i4)
    local_p = kron(z, i2, i4)
    local_q = kron(i2, z, i4)
    local_pi = kron(i2, i2, z, i2)
    local_lam = kron(i2, i2, i2, z)
    full_i = np.eye(256)

    def at(operator, j):
        return kron(operator, i16) if j == 0 else kron(i16, operator)

    c = [at(local_c, j) for j in range(2)]
    b = [at(local_b, j) for j in range(2)]
    p = [at(local_p, j) for j in range(2)]
    q = [at(local_q, j) for j in range(2)]
    pi = [at(local_pi, j) for j in range(2)]
    lam = [at(local_lam, j) for j in range(2)]
    ptot = [p[j] @ pi[j] for j in range(2)]
    qtot = [q[j] @ lam[j] for j in range(2)]
    prefix_p = [full_i, ptot[0]]
    prefix_q = [full_i, qtot[0]]
    klein = ptot[0] @ ptot[1]
    psi = [prefix_p[j] @ c[j] for j in range(2)]
    chi = [klein @ prefix_q[j] @ b[j] for j in range(2)]
    logical = psi + chi

    def xi(operator, j, r, s):
        result = at(kron(i2, i2, operator), j)
        if r:
            result = prefix_p[j] @ p[j] @ result
        if s:
            result = klein @ pi[j] @ prefix_q[j] @ q[j] @ result
        return result

    max_car = 0.
    for j, left in enumerate(logical):
        for k, right in enumerate(logical):
            max_car = max(max_car, error(left @ right + right @ left),
                          error(left @ right.T + right.T @ left
                                - (full_i if j == k else 0)))
    assert max_car < 1e-12

    # Auxiliary generators with all homogeneous parity possibilities.
    aux_local = [(kron(x, i2), 1, 0), (kron(i2, x), 0, 1),
                 (kron(x, x), 1, 1), (kron(z, i2), 0, 0)]
    aux_images = [[xi(op, j, r, s) for op, r, s in aux_local]
                  for j in range(2)]
    max_graded = 0.
    for j in range(2):
        for auxiliary, (_, r, s) in zip(aux_images[j], aux_local):
            sign = (-1) ** (r+s)
            for field in logical:
                max_graded = max(max_graded,
                                 error(field @ auxiliary-sign*auxiliary @ field))
    for left, (_, r, s) in zip(aux_images[0], aux_local):
        for right, (_, u, v) in zip(aux_images[1], aux_local):
            sign = (-1) ** ((r+s)*(u+v))
            max_graded = max(max_graded, error(left @ right-sign*right @ left))
    assert max_graded < 1e-12

    max_aux_product = 0.
    for j in range(2):
        for op, r, s in aux_local:
            for op2, u, v in aux_local:
                max_aux_product = max(max_aux_product,
                    error(xi(op, j, r, s) @ xi(op2, j, u, v)
                          - xi(op @ op2, j, r ^ u, s ^ v)))
    assert max_aux_product < 1e-12

    # Full vertices can have odd logical degree, compensated by an
    # auxiliary charge lowering operator. They remain local and commute
    # with every outside field/auxiliary generator after rewriting.
    pay_exchange = at(kron(a.T, i2, a, i2), 1)
    len_exchange = at(kron(i2, a.T, i2, a), 1)
    max_locality = 0.
    for exchange in (pay_exchange, len_exchange):
        vertex = exchange + exchange.T
        for outside in [psi[0], chi[0]] + aux_images[0]:
            max_locality = max(max_locality,
                              error(vertex @ outside-outside @ vertex))
    assert max_locality < 1e-12

    # The Gaussian reference has the same spectrum in every auxiliary
    # parity sector, with precisely the internal-prefix gauge of (10).
    hfull = .7*(psi[0].T @ psi[1] + psi[1].T @ psi[0])
    hfull += 1.2*(chi[0].T @ chi[1] + chi[1].T @ chi[0])
    hfull += .13*(psi[0].T @ psi[0])-.09*(chi[1].T @ chi[1])
    logical_words = list(product(range(2), repeat=4))

    def indices(aux0, aux1):
        return [(pay0*8+len0*4+aux0)*16 + pay1*8+len1*4+aux1
                for pay0, len0, pay1, len1 in logical_words]

    idx0 = indices(0, 0)
    hcode = hfull[np.ix_(idx0, idx0)]
    max_gauge = 0.
    norm_errors = []
    reference_norm = np.max(np.abs(np.linalg.eigvalsh(hcode)))
    for aux0, aux1 in product(range(4), repeat=2):
        alpha0 = (-1) ** (aux0 // 2)
        beta0 = (-1) ** (aux0 % 2)
        gauge = np.array([alpha0 ** pay1 * beta0 ** len1
                          for pay0, len0, pay1, len1 in logical_words])
        idx = indices(aux0, aux1)
        block = hfull[np.ix_(idx, idx)]
        max_gauge = max(max_gauge, error(block-gauge[:, None]*hcode*gauge))
        norm_errors.append(abs(np.max(np.abs(np.linalg.eigvalsh(block)))-reference_norm))
    assert max_gauge < 1e-12 and max(norm_errors) < 1e-12

    output = {"status": "all algebra checks passed", "full_dimension": 256,
              "logical_modes": 4, "auxiliary_parity_sectors": 16,
              "CAR_error": max_car, "graded_relation_error": max_graded,
              "auxiliary_product_error": max_aux_product,
              "charge_exchange_locality_error": max_locality,
              "reference_gauge_error": max_gauge,
              "reference_norm_error": float(max(norm_errors)),
              "scope": "Separate JW with two auxiliary charges; no infrared inference"}
    Path(__file__).with_name("jw_checks.json").write_text(json.dumps(output, indent=2)+"\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
