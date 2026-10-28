"""Independent finite-matrix controls of TILTED_EMBEDDING (numpy only)."""

from pathlib import Path
import json
import math

import numpy as np


def logtrace(energies, beta):
    z = -beta * energies
    m = float(z.max())
    return m + math.log(float(np.exp(z - m).sum()))


def one_site(op, j, count):
    out = np.ones((1, 1))
    for k in range(count):
        out = np.kron(out, op if k == j else np.eye(len(op)))
    return out


def chain(onsite, bond, count):
    d = len(onsite)
    result = sum(one_site(onsite, j, count) for j in range(count))
    for j in range(count - 1):
        result = result + np.kron(
            np.eye(d**j), np.kron(bond, np.eye(d ** (count - j - 2)))
        )
    return result


def main():
    rng = np.random.default_rng(20261028)
    b, delta, gamma, d_p = 2, 2.0, 0.5, 2
    lengths = np.array([2, 3, 1, 4])
    energies = np.array([0, 0, 2, 3.5])
    pair_length = (lengths[:, None] + lengths[None, :]).ravel()
    raw = rng.normal(size=(16, 16)) + 1j * rng.normal(size=(16, 16))
    bond = (raw + raw.conj().T) * (pair_length[:, None] == pair_length[None, :])
    bond *= 0.04 / np.linalg.norm(bond, 2)

    # Basis: logical bit first, then auxiliary vacuum or either excluded state.
    embedding = np.zeros((6, 4))
    embedding[[0, 3, 1, 2], np.arange(4)] = 1
    pair_embedding = np.kron(embedding, embedding)
    enlarged_bond = pair_embedding @ bond @ pair_embedding.T
    aux_a = np.array([0, 2, 3.5, 0, 2, 3.5])
    aux_q = np.array([0, -1, 2, 0, -1, 2])
    len_bit = np.array([0, 0, 0, 1, 1, 1])
    enlarged_length = b + len_bit + aux_q
    assert np.max(np.abs(np.diag(enlarged_length) @ embedding
                         - embedding @ np.diag(lengths))) == 0
    s0 = min(gamma / 4, delta / (4 * b))
    rows = []
    worst = 0.0
    for count in (2, 3):
        # Open chain. At three cells the middle cell touches both bonds.
        w = (1 if count == 2 else 2) * np.linalg.norm(bond, 2)
        for s in (-s0, 0.0, s0):
            h = chain(np.diag(energies + s * lengths), bond, count)
            he = chain(np.diag(aux_a + s * enlarged_length), enlarged_bond, count)
            eig = np.linalg.eigvalsh(h)
            eige = np.linalg.eigvalsh(he)
            for beta in (1.0, 3.0, 7.0):
                logratio = logtrace(eige, beta) - logtrace(eig, beta)
                zhat = sum(math.exp(-beta * (aux_a[j]
                           + s * (enlarged_length[j] - b - 0.5)))
                           for j in (4, 5))
                upper = count * math.log1p(math.exp(beta * w) * zhat / d_p)
                violation = max(-logratio, logratio - upper, 0.0)
                worst = max(worst, violation)
                assert violation < 1e-11
                rows.append(dict(cells=count, tilt=s, beta=beta,
                                 log_ratio=logratio, upper_bound=upper))

    derivative_ratios = []
    a_q = energies[2:]
    q = lengths[2:] - b
    for s in (-s0, 0.0, s0):
        for time in (1.0, 4.0):
            for order in range(5):
                left = float(np.sum(np.abs(time * q) ** order
                                    * np.exp(-time * (a_q + s * q))))
                right = math.factorial(order) * (time * b + 4 / gamma) ** order
                right *= float(np.exp(-time * a_q / 4).sum())
                derivative_ratios.append(left / right)
                assert left <= right + 1e-12

    result = dict(status="passed", description="tilted invalid-state and heat-trace controls",
                  trace_cases=len(rows), worst_trace_violation=worst,
                  derivative_cases=len(derivative_ratios),
                  max_derivative_bound_ratio=max(derivative_ratios),
                  rows=rows)
    Path(__file__).with_name("tilted_embedding_checks.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps({k: v for k, v in result.items() if k != "rows"}, indent=2))


if __name__ == "__main__":
    main()
