"""Finite-matrix corroboration of the entropy step in CONTACT_BOUND.md.

Uses NumPy only. The finite truncated lowering channel is trace preserving,
so no Hamiltonian truncation or thermodynamic extrapolation enters this check.
The analytic proof, not the random cases, establishes the general inequality.
"""

from pathlib import Path
import json
import numpy as np


def entropy(rho):
    eig = np.linalg.eigvalsh((rho + rho.conj().T) / 2)
    assert eig.min() > -2e-12
    eig = eig[eig > 1e-15]
    return float(-np.dot(eig, np.log(eig)))


def check(rho, gap_dim, aux_dim):
    p = np.zeros((gap_dim, gap_dim), complex)
    p[0, 0] = 1
    lower = np.diag(np.ones(gap_dim - 1), 1).astype(complex)
    kraus = [np.kron(p, np.eye(aux_dim)), np.kron(lower, np.eye(aux_dim))]
    image = sum(k @ rho @ k.conj().T for k in kraus)
    env = np.array([
        [np.trace(a @ rho @ b.conj().T) for b in kraus]
        for a in kraus
    ])
    a = float(np.trace(kraus[0] @ rho).real)
    h = 0.0 if a <= 0 or a >= 1 else -a * np.log(a) - (1 - a) * np.log1p(-a)
    loss = entropy(rho) - entropy(image)
    env_entropy = entropy(env)
    assert abs(np.trace(image) - 1) < 2e-12
    assert abs(env[0, 0] - a) < 2e-12
    assert abs(env[1, 1] - (1 - a)) < 2e-12
    assert loss <= env_entropy + 2e-12
    assert env_entropy <= h + 2e-12
    return dict(
        a=a, entropy_loss=loss, environment_entropy=env_entropy,
        binary_entropy=float(h), subadditivity_slack=env_entropy - loss,
        dephasing_slack=float(h - env_entropy),
    )


def main():
    rng = np.random.default_rng(20261109)
    records = []
    for gap_dim in (2, 3, 5, 8):
        for aux_dim in (1, 2, 3):
            dim = gap_dim * aux_dim
            for rank in (1, dim):
                for bias in (0.01, 1.0, 10.0):
                    z = rng.normal(size=(dim, rank)) + 1j * rng.normal(size=(dim, rank))
                    z[:aux_dim] *= bias
                    rho = z @ z.conj().T
                    rho /= np.trace(rho)
                    r = check(rho, gap_dim, aux_dim)
                    r.update(gap_dim=gap_dim, aux_dim=aux_dim, rank=rank, bias=bias)
                    records.append(r)

    # Both input levels map to |0>: entropy loss is exactly h_2(a).
    saturated = []
    for a in (0.0, 1 / 4096, 1 / 264, 0.1, 0.5, 0.9, 1.0):
        r = check(np.diag([a, 1 - a]).astype(complex), 2, 1)
        assert abs(r["entropy_loss"] - r["binary_entropy"]) < 2e-12
        saturated.append(r)

    result = dict(
        status="all checks passed",
        random_cases=len(records),
        sharp_classical_cases=len(saturated),
        minimum_subadditivity_slack=min(r["subadditivity_slack"] for r in records),
        minimum_dephasing_slack=min(r["dephasing_slack"] for r in records),
        sharp_classical_examples=saturated,
        scope="entropy-channel inequality only; no Hamiltonian or CFT test",
    )
    text = json.dumps(result, indent=2) + "\n"
    (Path(__file__).parent / "entropy_checks.json").write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
