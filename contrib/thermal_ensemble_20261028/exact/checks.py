"""Independent finite controls for LENGTH_FUGACITY.md; numpy only."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from numpy.linalg import eigvalsh


def logsumexp(values):
    largest = np.max(values)
    return largest + np.log(np.exp(values-largest).sum())


def bisect_zero(function, left, right):
    """Zero of a continuous decreasing function with a supplied bracket."""
    assert function(left) * function(right) <= 0
    for _ in range(100):
        middle = (left+right)/2
        if function(middle) > 0:
            left = middle
        else:
            right = middle
    return (left+right)/2


def kron_all(items):
    out = np.ones((1, 1))
    for item in items:
        out = np.kron(out, item)
    return out


def logtrace(eigenvalues, beta):
    return float(logsumexp(-beta * eigenvalues))


def thermal_controls():
    rng = np.random.default_rng(20261028)
    # Two code states, followed by five excited states. Length is conserved,
    # whereas the number of excited cells need not be conserved.
    lengths = np.array([2, 3, 1, 2, 3, 4, 5])
    energies = np.array([0, 0, .8, 1, 1.2, 2, 3])
    d, dp = len(lengths), 2
    pair_length = (lengths[:, None] + lengths).ravel()
    is_code = np.zeros(d, dtype=bool)
    is_code[:dp] = True
    pair_code = np.logical_and(is_code[:, None], is_code).ravel()
    bond = np.zeros((d * d, d * d))
    for ell in np.unique(pair_length):
        for code in (False, True):
            idx = np.flatnonzero((pair_length == ell) & (pair_code == code))
            block = rng.normal(size=(len(idx), len(idx)))
            block = (block + block.T) / 2
            bond[np.ix_(idx, idx)] = block
    bond *= .03 / max(abs(eigvalsh(bond)))
    identity = np.eye(d)
    all_cases = []
    nonconservation = 0.
    for m in (2, 3):
        configurations = np.array(list(np.ndindex(*(d,) * m)))
        total_length = lengths[configurations].sum(axis=1)
        total_q = (~is_code[configurations]).sum(axis=1)
        ham = np.diag(energies[configurations].sum(axis=1))
        for j in range(m - 1):
            ham += kron_all([identity] * j + [bond] + [identity] * (m - j - 2))
        nonconservation = max(nonconservation,
                             float(np.linalg.norm((total_q[:, None] - total_q) * ham)))
        assert np.linalg.norm((total_length[:, None] - total_length) * ham) < 1e-12
        code_idx = np.flatnonzero(np.all(is_code[configurations], axis=1))
        # A uniform upper bound on the incident interaction norm, including
        # boundary cells. The proof permits this conservative choice.
        j_bound = .06
        for s in (-.02, 0., .02):
            shifted = ham + np.diag(s * total_length)
            ev = eigvalsh(shifted)
            pev = eigvalsh(shifted[np.ix_(code_idx, code_idx)])
            for beta in (1., 4., 10.):
                delta = logtrace(ev, beta) - logtrace(pev, beta)
                t = np.exp(-beta * (energies[dp:] + s * (lengths[dp:] - 2.5)
                                    - 2 * j_bound)).sum()
                upper = m * np.log1p(t / dp)
                assert delta >= -1e-10 and delta <= upper + 1e-10
                all_cases.append({"M": m, "s": s, "beta": beta,
                                  "log_ratio": delta, "upper": float(upper)})
    assert nonconservation > 1e-4
    return {"cases": all_cases, "defect_count_commutator_frobenius": nonconservation,
            "passed": len(all_cases)}


def renewal_controls():
    # Exact all-length independent cells. Code lengths 2 and 3 each have
    # one state at zero energy. Wrong lengths have exponentially small
    # weights but are not truncated in the pressure equation.
    # Excited labels: ell=1 energy 2, ell=2,3 energy 2;
    # ell>=4 energy ell, one state per length.
    def q_selected(s, beta):
        return np.exp(-beta * 2 * s) + np.exp(-beta * 3 * s)

    def q_excited(s, beta):
        r = np.exp(-beta * (1 + s))
        return (np.exp(-beta * (2 + s)) + np.exp(-beta * (2 + 2 * s))
                + np.exp(-beta * (2 + 3 * s)) + r ** 4 / (1 - r))

    output = []
    for beta in (2., 4., 6., 10.):
        selected = bisect_zero(lambda s: q_selected(s, beta) - 1, 0., 1.)
        full = bisect_zero(lambda s: q_selected(s, beta) + q_excited(s, beta) - 1, 0., 1.)
        # Bound (17) with exact on-cell traces, using J=0 and a uniform
        # bound over the tiny interval between the two pressure roots.
        # t_s=sum_exc exp[-beta(A+s(ell-2.5))]. Convexity means its maximum
        # on this compact interval is at an endpoint.
        def t(s):
            return np.exp(beta * 2.5 * s) * q_excited(s, beta)
        delta = np.log1p(max(t(selected), t(full)) / 2) / beta
        assert 0 <= full - selected <= delta / 2 + 1e-12

        # Direct fixed-L partition recursion, rescaled by the predicted
        # pressure to avoid overflow. Every cell length through L is kept.
        limit = 600
        weights = np.zeros(limit + 1)
        weights[1] = np.exp(-beta * (2 + full))
        weights[2] = np.exp(-beta * 2 * full) * (1 + np.exp(-2 * beta))
        weights[3] = np.exp(-beta * 3 * full) * (1 + np.exp(-2 * beta))
        lengths = np.arange(4, limit + 1)
        weights[4:] = np.exp(-beta * (1 + full) * lengths)
        z = np.zeros(limit + 1)
        z[0] = 1.
        for ell in range(1, limit + 1):
            z[ell] = np.dot(weights[1:ell + 1], z[ell - 1::-1])
        finite_pressure = full + np.log(z[limit]) / (beta * limit)
        assert abs(finite_pressure - full) < 1 / limit
        output.append({"beta": beta, "selected_pressure": float(selected),
                       "full_pressure": float(full), "difference": float(full-selected),
                       "bound": float(delta/2), "L": limit,
                       "direct_fixed_L_pressure": float(finite_pressure)})
    return output


if __name__ == "__main__":
    result = {"thermal": thermal_controls(), "renewal": renewal_controls()}
    path = Path(__file__).with_name("checks.json")
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"thermal_cases_passed": result["thermal"]["passed"],
                      "defect_count_commutator_frobenius":
                      result["thermal"]["defect_count_commutator_frobenius"],
                      "renewal_controls": len(result["renewal"]),
                      "output": str(path)}, indent=2))
